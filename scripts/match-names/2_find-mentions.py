"""Find mentions of molecules in a set of text files

Assumes text files are each just one sentence per line.
"""

import sys
from collections import defaultdict
from dataclasses import dataclass, field
from hashlib import sha256
from queue import Queue
from threading import Event
from typing import List, Tuple, Dict
from pathlib import Path
from datetime import datetime
from glob import glob
import argparse
import logging
import gzip
import json

from colmena.queue import PipeQueues
from colmena.task_server import ParslTaskServer
from colmena.thinker import BaseThinker, agent, ResourceCounter
from more_itertools import batched
from parsl import Config, HighThroughputExecutor
from parsl.addresses import address_by_hostname
from parsl.launchers import AprunLauncher
from parsl.providers import CobaltProvider
from proxystore.store.file import FileStore, FileStoreKey
from proxystore.store.utils import get_key
from pymongo.database import Database
from pymongo import MongoClient
import proxystore as ps
import pandas as pd
import numpy as np


def parsl_config(name: str, run_path: Path) -> Tuple[Config, int]:
    """Make compute resource configuration

    Args:
        name: Name of the desired configuration
        run_path: Output path for th erun
    Returns:
        - Parsl compute configuration
        - Number of compute slots: Includes execution slots and pre-fetch buffers
    """

    run_path = str(run_path / 'runinfo')
    if name == 'local':
        return Config(
            executors=[
                HighThroughputExecutor(max_workers=1, prefetch_capacity=1)
            ],
            run_dir=run_path
        ), 1
    elif name == 'theta-debug':
        return Config(
            retries=2,
            executors=[HighThroughputExecutor(
                address=address_by_hostname(),
                label="debug",
                max_workers=64,
                prefetch_capacity=64,
                cpu_affinity='block',
                provider=CobaltProvider(
                    account='CSC249ADCD08',
                    queue='debug-flat-quad',
                    nodes_per_block=8,
                    scheduler_options='#COBALT --attrs enable_ssh=1',
                    walltime='00:60:00',
                    init_blocks=0,
                    max_blocks=1,
                    cmd_timeout=360,
                    launcher=AprunLauncher(overrides='-d 64 --cc depth -j 1'),
                    worker_init='''
module load miniconda-3
conda activate /lus/eagle/projects/ExaLearn/carbon-free-ldrd/env-cpu''',
                ),
            )],
            run_dir=run_path
        ), 64 * 8 * 3
    elif name == 'theta-full':
        return Config(
            retries=2,
            executors=[HighThroughputExecutor(
                address=address_by_hostname(),
                label="full",
                max_workers=64,
                prefetch_capacity=64,
                cpu_affinity='block',
                provider=CobaltProvider(
                    account='SuperBERT',
                    nodes_per_block=900,
                    scheduler_options='#COBALT --attrs enable_ssh=1',
                    walltime='20:00:00',
                    init_blocks=0,
                    max_blocks=1,
                    cmd_timeout=360,
                    launcher=AprunLauncher(overrides='-d 64 --cc depth -j 2'),
                    worker_init=f'''
module load miniconda-3
conda activate /lus/eagle/projects/ExaLearn/carbon-free-ldrd/env
export PYTHONPATH="$PYTHONPATH:{Path(__file__).parent}"
pwd 
which python
''',
                ),
            )],
            run_dir=run_path
        ), 900 * 64 * 4
    else:
        raise ValueError(f'Configuration not defined: {name}')


@dataclass(order=True)
class _PrioritizedItem:
    priority: float
    item: str = field(compare=False)


class ScreenEngine(BaseThinker):
    """Screening engine that screens molecules in parallel

    Reads from a list of molecules too large to fit in memory and sends it out gradually in small chunks to be evaluated

    Parameters:
        queues: List of queues
        text_paths: Path to text mentioning molecules
        store: Store used to pass inputs from thinker to worker
        subsets: List of molecule sys
        output_dir: Path to the output directory
        slot_count: Number of execution slots
        molecules_per_task: Number of molecules per chunk
        lines_per_task: Number of lines of text per chunk
    """

    def __init__(self,
                 queues: PipeQueues,
                 store: FileStore,
                 database: Database,
                 text_paths: List[Path],
                 subsets: List[str],
                 output_dir: Path,
                 slot_count: int,
                 molecules_per_task: int,
                 lines_per_task: int,
                 ):
        # Init the base class
        super().__init__(queues, ResourceCounter(slot_count, ['screen']))
        self.rec.reallocate(None, 'screen', 'all')

        # Store the input and output information
        self.screen_paths = sorted(text_paths)
        self.output_dir = output_dir
        self.molecules_per_task = molecules_per_task
        self.lines_per_task = lines_per_task
        self.molecules = database['molecule_record']
        self.mentions = database['mention']
        self.store = store
        self.subsets = list(subsets)

        # Queue to store ready-to-compute chunks
        self.queue_depth = slot_count * 3  # One in queue for each running already
        self.screen_queue = Queue(self.queue_depth)

        # Queue to store results ready to add to the list
        self.result_queue = Queue()

        # Tracking when we are done with a specific proxy
        self.proxy_key: Dict[str, FileStoreKey] = {}  # Key belonging to a proxy. Map of (filename, start line) to whatever ProxyStore uses
        self.all_submitted: Dict[str, Event] = {}  # Whether all molecules for this have been submitted
        self.tasks_remaining: Dict[str, int] = {}  # How many tasks for this file are left. Key: filename, start line
        self.stores_remaining: Dict[str, int] = {}  # How many tasks for this file are left. Key: filename, start line

        # Global tracking if we're done with a full file
        self.completed_path = Path("files.complete")

        # Things to know if we are done
        self.all_read = Event()
        self.queue_filled = Event()
        self.total_molecule_names = 0
        self.total_chunks = 0
        self.total_matches = 0
        self.total_pairs = 0

        # Recording for performance information
        self.start_time = np.inf  # Use the time the first compute starts

    @agent(startup=True)
    def prepare_inputs(self):
        """Read chunks of molecules and text files to be screened"""

        # Load in the list of file chunks which have been processed fully
        if not self.completed_path.exists():
            self.completed_path.write_text("filename,lines_per_task,start_line\n")
            already_ran = []
        else:
            already_ran = pd.read_csv(self.completed_path).query(f'lines_per_task=={self.lines_per_task}').apply(tuple, axis=1)

        # Get a list of all the molecules we could run
        cursor = self.molecules.find({'subsets': {'$in': self.subsets}}, projection=['names'])

        def yield_names(cursor):
            for record in cursor:
                for name in record['names']:
                    self.total_molecule_names += 1
                    yield record['_id'], name

        mol_batch_sizes = []
        to_match = []
        for batch_id, molecule_batch in enumerate(batched(yield_names(cursor), self.molecules_per_task)):
            mol_batch_sizes.append(len(molecule_batch))
            to_match.append(self.store.proxy(molecule_batch))

        self.logger.info(f'Prepared a list of {self.total_molecule_names} names to match to the text. '
                         f'Will run them in {len(to_match)} batches')

        # Loop over files, as we want to minimize the time a proxy exists
        for file_id, filename in enumerate(self.screen_paths):
            self.logger.info(f'Screening file {file_id + 1}/{len(self.screen_paths)}.')

            with open(filename) as fp:
                submitted_chunks = 0
                for lines_id, lines in enumerate(batched(fp, self.lines_per_task)):
                    # Store which line this chunk starts at
                    start_line = lines_id * self.lines_per_task
                    key = f'{filename.absolute()}-{start_line}'

                    # Check if we've done everything already
                    if (str(filename.absolute()), self.lines_per_task, start_line) in already_ran:
                        self.logger.info(f'Already ran {filename} starting at {start_line}. Skipping')

                    # Proxy the section of file
                    lines_proxy = self.store.proxy(lines)

                    # Prepare to track when we are done with it
                    self.proxy_key[key] = get_key(lines_proxy)
                    self.all_submitted[key] = Event()
                    self.tasks_remaining[key] = 0
                    self.stores_remaining[key] = 0

                    # Submit the molecules for this task
                    for mol_id, molecules in enumerate(to_match):
                        submitted_chunks += 1
                        n_pairs = len(lines) * mol_batch_sizes[mol_id]
                        self.total_pairs += n_pairs
                        self.screen_queue.put((key, str(filename.absolute()), start_line, n_pairs, lines_proxy, molecules))
                        self.total_chunks += 1
                        self.tasks_remaining[key] += 1
                        self.stores_remaining[key] += 1

                        # Check if we have submitted enough task to start
                        if self.total_chunks == self.queue_depth:
                            self.logger.info('Queue is filled. Submit tasks now!')
                            self.queue_filled.set()

                    # Mark that we have finished submitting all tasks for this file
                    self.all_submitted[key].set()

                # Message that we're done with this file
                self.logger.info(f'Submitted {submitted_chunks} chunks for {filename}')

        # Put a None at the end to signal we are done
        self.queue_filled.set()  # Make sure the task submitter starts
        self.screen_queue.put(None)

        # Mark that we are done reading
        self.logger.info(f'Finished submitting screening tasks and submitting {self.total_chunks} blocks')
        self.all_read.set()

    @agent(startup=True)
    def submit_task(self):
        """Submit chunks of molecules to be screened"""

        # Wait until the queue is filled
        self.queue_filled.wait()
        self.logger.info('Starting to submit tasks.')

        while True:
            # Get the next chunk and, if None, break
            msg = self.screen_queue.get()
            if msg is None:
                break

            # Create a proxy for the chunk data
            key, filename, start_line, total_matches, lines, names = msg

            # Submit once we have resources available
            self.rec.acquire("screen", 1)
            self.queues.send_inputs(lines, names, method='find_matches',
                                    task_info={'key': key, 'filename': filename, 'start_line': start_line,
                                               'total_matches': total_matches})

    @agent(startup=True)
    def receive_results(self):
        """Mark results complete"""

        # Open the output file
        num_recorded = 0
        with gzip.open(self.output_dir / 'inference-results.json.gz', 'wt') as fq:
            while not self.done.is_set():
                # Stop when all have been recorded
                if self.all_read.is_set() and num_recorded >= self.total_chunks:
                    break

                # Wait for a chunk to be received
                result = self.queues.get_result()
                self.rec.release("screen", 1)
                num_recorded += 1

                # Save the JSON result
                print(result.json(exclude={'inputs', 'value'}), file=fq)

                # Update the start time
                self.start_time = min(self.start_time, result.time_compute_started)

                # Remove the chunk from the proxystore if we can
                chunk_key = result.task_info['key']
                self.tasks_remaining[chunk_key] -= 1
                if self.tasks_remaining[chunk_key] == 0 and self.all_submitted[chunk_key].is_set():
                    self.logger.info(f'Done with all work on {chunk_key}. Evicting from store')
                    proxy_key = self.proxy_key.pop(chunk_key)
                    self.store.evict(proxy_key)

                # If unsuccessful, just skip this one
                if not result.success:
                    self.logger.error(f'Task failed. Traceback: {result.failure_info.traceback}')
                else:
                    # Push the result to be processed
                    self.result_queue.put((chunk_key, result.task_info['filename'], result.task_info['start_line'], result.value))

                # Print a status message
                if self.all_read.is_set():
                    self.logger.info(f'Received task {num_recorded}/{self.total_chunks}. Processing backlog: {self.result_queue.qsize()}')
                else:
                    self.logger.info(f'Received task {num_recorded}/???. Processing backlog: {self.result_queue.qsize()}')

        # Send a "done" message
        self.result_queue.put(None)

        # Print the final status
        run_time = datetime.now().timestamp() - self.start_time
        self.logger.info(f'Completed storing all results. {run_time:.2f} s.')

    @agent()
    def store_results(self):
        """Store the matches to a database"""

        # Loop while the simulation hasn't been killed
        count = 0
        while not self.done.is_set():
            count += 1
            # Get the next chunk ready for processing
            msg: Tuple[str, str, int, Dict[str, Tuple[str, int, str]]] = self.result_queue.get()

            # If it is `None` we are done
            if msg is None:
                break

            # Otherwise, store results
            chunk_key, filename, start_line, match_result = msg
            for key, matches in match_result.items():  # match_result is (inchi_key, match infomration)
                self.total_matches += len(matches)

                # Format the match as a dictionary
                matches_by_line = defaultdict(list)  # Map of line-id -> molecule names
                texts = {}
                for name, line_id, text in matches:  # matches is a list of (name of molecule, line_id of match, text of match)
                    matches_by_line[line_id + start_line].append({
                        'key': key,
                        'name': name,
                    })
                    texts[line_id + start_line] = text

                for line_id, hits in matches_by_line.items():
                    self.mentions.update_one({'filename': filename, 'line': line_id},
                                             {'$addToSet': {'matches': {'$each': hits}},
                                              '$setOnInsert': {'text': texts[line_id]}}, upsert=True)

            # Mark that we're done processing this result, and the whole file if we're done
            self.stores_remaining[chunk_key] -= 1
            if self.stores_remaining[chunk_key] == 0:
                self.logger.info(f'Finished storing all results from {self.stores_remaining}. Writing it to disk')
                with self.completed_path.open("a") as fp:
                    print(f'{filename},{self.lines_per_task},{start_line}', file=fp)

            # Print a status message
            status = f'Found matches for {len(match_result)} molecules in {filename}. Total found: {self.total_matches}'
            if self.all_read.is_set():
                self.logger.info(f'Processed task {count}/{self.total_chunks}. {status}')
            else:
                self.logger.info(f'Processed task {count}/???. {status}')

        # Write out the total run time.
        run_time = datetime.now().timestamp() - self.start_time
        self.logger.info(f'Runtime {run_time:.2f} s. Total matches made: {self.total_matches}. Total screened: {self.total_pairs}.'
                         f' Overall evaluation rate: {self.total_pairs / run_time:.3e} pairs/s')


def find_matches(lines: List[str], names: List[Tuple[str, str]]) -> Dict[str, List[Tuple[str, int, str]]]:
    """Find occurrences of different names in a file

    Args:
        lines: Lines from a file to be matched
        names: Pairs of molecule key and name
    Returns:
        - List of matches: molecule key -> [molecule name, line number, line content)
    """
    from collections import defaultdict
    from itertools import product
    output = defaultdict(list)

    # Beginning and ends of strings
    start_chars = ' ('
    end_chars = ' ).;,:'
    caps = list(product(start_chars, end_chars))

    # Lower all names
    names_lower = [n.lower() for _, n in names]

    # Match the name as a full token. It must be preceded by a parenthesis, comma, space,
    #  and followed by a comma, period, parenthesis; unless it as the beginning or end of line
    for line_id, line in enumerate(lines):
        test_line = ' ' + line.lower() + ' '  # Pad with spaces, to simplify matching logic
        for (key, name), name_lower in zip(names, names_lower):
            if name_lower in test_line and any((s + name_lower + e) in test_line for s, e in caps):
                output[key].append((name, line_id, line))
    return output


if __name__ == '__main__':
    # User inputs
    parser = argparse.ArgumentParser()
    parser.add_argument('--overwrite', action='store_true',
                        help='Whether to overwrite a previous run')

    group = parser.add_argument_group(title='Compute Configuration',
                                      description='Settings related to how compute is executed')
    group.add_argument("--molecules-per-chunk", default=2000, type=int, help="Number molecules per screening task")
    group.add_argument("--lines-per-chunk", default=100000, type=int, help="Number lines of files to run per task")
    group.add_argument('--compute', default='local', help='Resource configuration')
    group.add_argument('--proxy-store', default='file', help='What kind of proxy store to use')

    group = parser.add_argument_group(title='Search Configuration', description='How we screen for molecules')
    group.add_argument('--search-space', nargs='+', help='Path to molecules to be screened. If only one path provided, we assume it is a glob string')
    group.add_argument('--subsets', default=['known'], nargs='+', help='Which subsets of molecules to run')

    group = parser.add_argument_group(title='Database Configuration', description='Things associated with data storage')
    group.add_argument('--mongo-port', default=27856, type=int, help='Port on which to connect to MongoDB')

    # Parse the arguments
    args = parser.parse_args()
    run_params = args.__dict__.copy()
    run_params.pop('overwrite')

    # Connect to Mongo
    client = MongoClient(port=args.mongo_port)
    db = client['cfree']
    coll = db['molecule_record']

    # Make sure at least one molecule is available for matching from each subset
    for subset in args.subsets:
        mol_count = coll.count_documents({'subsets': subset}, limit=1)
        assert mol_count > 0, f'Did not match any molecules for {subset}'

    # Check that the search path exists
    if len(args.search_space) == 1:
        search_paths = [Path(x) for x in glob(args.search_space[0], recursive=True)]
    else:
        search_paths = [Path(x) for x in args.search_space]
    assert all(x.is_file() for x in search_paths)

    # Create an output directory with the name of the directory
    run_hash = sha256()
    run_hash.update(str(run_params).encode())
    subset_names = "_".join(args.subsets)
    out_path = Path().joinpath('runs', f'matching-{subset_names}-mols-{run_hash.hexdigest()[:6]}')
    out_path.mkdir(parents=True, exist_ok=args.overwrite)

    # Store the run parameters
    with open(out_path / 'run-config.json', 'w') as fp:
        json.dump(run_params, fp, indent=2)

    # Set up the logging
    handlers = [logging.FileHandler(out_path / 'runtime.log', mode='w'), logging.StreamHandler(sys.stdout)]


    class ParslFilter(logging.Filter):
        """Filter out Parsl debug logs"""

        def filter(self, record):
            return not (record.levelno == logging.DEBUG and '/parsl/' in record.pathname)


    for h in handlers:
        h.addFilter(ParslFilter())

    logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                        level=logging.INFO, handlers=handlers)

    # Make Parsl engine
    config, n_slots = parsl_config(args.compute, out_path)

    # Configure the file
    if args.proxy_store == 'file':
        ps_file_dir = out_path / 'file-store'
        ps_file_dir.mkdir(exist_ok=True)
        store = FileStore(name='file', store_dir=str(ps_file_dir))
    else:
        raise ValueError('ProxyStore config not recognized: {}')

    # Make the task queues and task server
    queues = PipeQueues(proxystore_name=args.proxy_store, proxystore_threshold=1000)
    task_server = ParslTaskServer([find_matches], queues, config)

    # Make the thinker
    thinker = ScreenEngine(queues, store, db, search_paths, args.subsets,
                           out_path, n_slots, args.molecules_per_chunk, args.lines_per_chunk)

    # Run the program
    try:
        task_server.start()
        thinker.run()
    finally:
        queues.send_kill_signal()

    task_server.join()
    store.close()
