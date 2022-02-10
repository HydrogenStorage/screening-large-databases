import sys
from dataclasses import dataclass, field
from functools import partial, update_wrapper
from hashlib import sha256
from queue import Queue
from threading import Event
from typing import List, Tuple
from pathlib import Path
from csv import DictWriter
from datetime import datetime
import argparse
import logging
import heapq
import json

from colmena.redis.queue import make_queue_pairs, ClientQueues
from colmena.task_server import ParslTaskServer
from colmena.thinker import BaseThinker, agent, ResourceCounter
from parsl import Config, HighThroughputExecutor
from parsl.addresses import address_by_hostname
from parsl.launchers import AprunLauncher
from parsl.providers import CobaltProvider
import proxystore as ps
import numpy as np


def parsl_config(name: str) -> Tuple[Config, int]:
    """Make the compute resource configuration

    Args:
        name: Name of the diesred configuration
    Returns:
        - Parsl compute configuration
        - Number of compute slots: Includes execution slots and pre-fetch buffers
    """

    if name == 'local':
        return Config(
            executors=[
                HighThroughputExecutor(max_workers=16, prefetch_capacity=1)
            ]
        ), 64
    elif name == 'theta-debug':
        return Config(
            retries=16,
            executors=[HighThroughputExecutor(
                    address=address_by_hostname(),
                    label="debug",
                    max_workers=64,
                    prefetch_capacity=2,
                    cpu_affinity='block',
                    provider=CobaltProvider(
                        account='redox_adsp',
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
conda activate /lus/theta-fs0/projects/CSC249ADCD08/edw/env''',
                    ),
                )]
            ), 64 * 8 * 4
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
        screen_path: Path to molecules to be screened
        store: Store used to pass inputs from thinker to worker
        output_dir: Path to the output directory
        slot_count: Number of execution slots
        chunk_size: Number of molecules per chunk
    """

    def __init__(self,
                 queues: ClientQueues,
                 store: ps.store.remote.RemoteStore,
                 screen_path: Path,
                 target_count: int,
                 output_dir: Path,
                 slot_count: int,
                 chunk_size: int
                 ):
        # Init the base class
        super().__init__(queues, ResourceCounter(slot_count, ['screen']))
        self.rec.reallocate(None, 'screen', 'all')

        # Store the input and output information
        self.screen_path = screen_path
        self.output_dir = output_dir
        self.chunk_size = chunk_size
        self.target_count = target_count
        self.store = store

        # Queue to store ready-to-compute chunks
        self.screen_queue = Queue(slot_count * 2)

        # Queue to store results ready to add to the list
        self.result_queue = Queue()

        # Variables to hold the outputs
        self.best_mols: List[_PrioritizedItem] = []  # Will be used as a heap
        self.current_threshold = 0

        # Things to know if we are done
        self.all_read = Event()
        self.total_chunks = 0
        self.total_molecules = 0

    @agent(startup=True)
    def read_chunks(self):
        """Read chunks to be screened"""
        
        # Make a function to push a chunk of molecules to the execution queue
        def _proxy_and_push(chunk):
            # Submit the object to the proxy store
            key = f'chunk-{self.total_chunks}'
            chunk_proxy = self.store.proxy(chunk, key=key)
            
            # Increment the count
            self.total_chunks += 1
            
            # Put the proxy and the key in the queue
            self.screen_queue.put((chunk_proxy, key))

        with self.screen_path.open() as fp:
            chunk = []
            for line in fp:
                try:
                    # The line is comma-separated with the last entry as the string
                    _, smiles = line.strip().rsplit(",", 1)
                except:
                    continue
                self.total_molecules += 1

                # Add to the chunk and submit if we hit the target size
                chunk.append(smiles)
                if len(chunk) >= self.chunk_size:
                    _proxy_and_push(chunk)
                    chunk = []
                    
        # Submit whatever remains
        _proxy_and_push(chunk)

        # Put a None at the end to signal we are done
        self.screen_queue.put(None)

        # Mark that we are done reading
        self.logger.info(f'Finished reading {self.total_molecules} molecules and submitting {self.total_chunks} blocks')
        self.all_read.set()

    @agent(startup=True)
    def submit_task(self):
        """Submit chunks of molecules to be screened"""
        while True:
            # Get the next chunk and, if None, break
            msg = self.screen_queue.get()
            if msg is None:
                break
                
            # Create a proxy for the chunk data
            chunk, key = msg

            # Submit once we have resources available
            self.rec.acquire("screen", 1)
            self.queues.send_inputs(chunk, self.current_threshold, method='screen_molecules', task_info={'key': key})

    @agent(startup=True)
    def receive_results(self):
        """Mark results complete"""

        # Open the output file
        start_time = np.inf  # Use the time the first compute starts
        num_recorded = 0
        with open(self.output_dir / 'inference-results.json', 'w') as fq:
            while True:
                # Stop when all have been recorded
                if self.all_read.is_set() and num_recorded >= self.total_chunks:
                    break

                # Wait for a chunk to be received
                result = self.queues.get_result()
                self.rec.release("screen", 1)
                if not result.success:
                    self.logger.error(f'Task failed. Traceback: {result.failure_info.traceback}')
                    raise ValueError('Failed task')

                # Push the result to be processed
                self.result_queue.put(result.value)
                num_recorded += 1
                    
                # Update the start time
                start_time = min(start_time, result.time_compute_started)
                
                # Remove the chunk from the proxystore
                self.store.evict(result.task_info['key'])
                
                # Save the JSON result
                print(result.json(exclude={'inputs', 'value'}), file=fq)

                # Print a status message
                if self.all_read.is_set():
                    self.logger.info(f'Recorded task {num_recorded}/{self.total_chunks}. Processing backlog: {self.result_queue.qsize()}')
                else:
                    self.logger.info(f'Recorded task {num_recorded}/???. Processing backlog: {self.result_queue.qsize()}')

        # Send a "done" message
        self.result_queue.put(None)

        # Print the final status
        run_time = datetime.now().timestamp() - start_time
        self.logger.info(f'Runtime {run_time:.2f} s. Evaluation rate: {self.total_molecules / run_time:.3e} mol/s')

    @agent()
    def store_results(self):
        """Keep a running tally of the best molecules"""

        # Keep track of whether we've reached the target size or not
        reached_target_size = False

        # Loop while the simulation hasn't been killed
        count = 0
        while not self.done.is_set():
            # Get the next chunk ready for processing
            result: List[Tuple[float, str]] = self.result_queue.get()
            count += 1

            # If it is `None` we are done
            if result is None:
                break

            # Add objects to the priority queue heap
            for score, entry in result:
                item = _PrioritizedItem(score, entry)
                if reached_target_size:
                    # Maintain the heap size
                    removed = heapq.heappushpop(self.best_mols, item)
                    assert removed.priority <= item.priority, (removed.priority, item.priority)
                else:
                    # Add to the queue without removing
                    heapq.heappush(self.best_mols, item)
                    reached_target_size = len(self.best_mols) >= self.target_count
                    if reached_target_size:
                        self.logger.info(f'We have filled the queue of {len(self.best_mols)} best molecules')

            # Update the minimum value required to be on
            if reached_target_size:
                self.current_threshold = self.best_mols[0].priority

            # Print a status message
            if self.all_read.is_set():
                self.logger.info(f'Recorded task {count}/{self.total_chunks}. Current threshold: {self.current_threshold:.4f}')
            else:
                self.logger.info(f'Recorded task {count}/???. Current threshold: {self.current_threshold:.4f}')

        # Write the data out to disk
        self.logger.info(f'Completed processing all results. Output list size: {len(self.best_mols)}')

        with open(self.output_dir / 'best_molecules.csv', 'w') as fp:
            writer = DictWriter(fp, ['smiles', 'inchi', 'score'])
            writer.writeheader()
            for record in sorted(self.best_mols, reverse=True):
                entry = json.loads(record.item)
                entry['score'] = record.priority
                writer.writerow(entry)

        self.logger.info('Finished writing results to disk')


def screen_molecules(to_screen: List[str], min_similarity: float, to_compare: Tuple[str], radius: int = 4) -> List[Tuple[float, str]]:
    """Compute the difference betwee

    Args:
        to_screen: List of SMILES strings to string
        to_compare: List of SMILES strings to compare against
        min_similarity: Minimum similarity to be returned by this function
        radius: Radius of the fingerprint
    Returns:
        List of tuples of the similarity score and molecule IDs for each molecule above min_similarity
    """
    import json
    from rdkit import Chem
    from rdkit.Chem import DataStructs, AllChem

    # Check if we have already parsed the comparison set
    my_globals = globals()
    compare_mols = my_globals.get('compare_mols', None)
    compare_mols_hash = my_globals.get('compare_mols_hash', None)
    if (hash(to_compare) + radius) != compare_mols_hash:
        compare_mols = None

    # Parse the molecules, if needed
    if compare_mols is None:
        # Parse them and compute fingerprints
        compare_mols = []
        for mol_str in to_compare:
            mol = Chem.MolFromSmiles(mol_str)
            if mol is None:
                raise ValueError(f'Failed to parse: {mol_str}')
            compare_mols.append(AllChem.GetMorganFingerprint(mol, radius))

        # Store them for later use
        my_globals['compare_mols'] = compare_mols
        my_globals['compare_mols_hash'] = hash(to_compare) + radius

    # Compute max similarity to any of the comparison set, return only those above threshold
    passed = []
    for smiles in to_screen:
        mol = Chem.MolFromSmiles(smiles)

        # Skip if molecule does not parse
        if mol is None:
            continue

        # Compute the fingerprints
        fp = AllChem.GetMorganFingerprint(mol, radius)

        # Compute the maximum similarity to the comparison set
        mol_similarity = max(DataStructs.BulkTanimotoSimilarity(fp, compare_mols))

        # If it doesn't meet the threshold, don't return it
        if mol_similarity < min_similarity:
            continue

        # If not, add it to the output dictionary
        passed.append((mol_similarity, json.dumps({'smiles': smiles, 'inchi': Chem.MolToInchi(mol)})))

    return passed


if __name__ == '__main__':
    # User inputs
    parser = argparse.ArgumentParser()

    group = parser.add_argument_group(title='Compute Configuration', description='Settings related to how compute is executed')

    group.add_argument("--redishost", default="127.0.0.1", help="Address at which the redis server can be reached")
    group.add_argument("--redisport", default="6379", help="Port on which redis is available")
    group.add_argument("--molecules-per-chunk", default=50000, type=int, help="Number molecules per screening task")
    group.add_argument('--compute', default='local', help='Resource configuration')
    group.add_argument('--proxy-store', default='file', help='What kind of proxy store to use')
    group.add_argument('--overwrite', action='store_true', help='Whether to overwrite a previous run')

    group = parser.add_argument_group(title='Search Coonfiguration', description='How we screen for the top molecules')
    group.add_argument('--search-space', help='Path to molecules to be screened', required=True)
    group.add_argument('--num-top', default=1000000, type=int, help='Number of most-similar molecules to find')
    group.add_argument('--comparison-molecules', required=True, help='Path to the molecule(s) to compare against. Must be line-delimited SMILES')
    group.add_argument('--radius', default=4, type=int, help='Radius of the fingerprint')

    # Parse the arguments
    args = parser.parse_args()
    run_params = args.__dict__

    # Check that the search path exists
    search_path = Path(args.search_space)
    assert search_path.is_file()

    # Load in the molecules to be screened
    with open(args.comparison_molecules) as fp:
        to_compare = [x.strip() for x in fp]
    
    # Create an output directory with the name of the directory
    run_hash = sha256()
    run_hash.update(str(run_params).encode())
    run_hash.update(str(to_compare).encode())
    out_path = Path().joinpath('runs', f'{search_path.name[:-4]}-top{args.num_top}-{run_hash.hexdigest()[:6]}')
    out_path.mkdir(parents=True, exist_ok=args.overwrite)
    
    # Store the run parameters
    with open(out_path / 'run-config.json', 'w') as fp:
        json.dumps(run_params, indent=2)
    
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

    # Log the run information
    logging.info(f'Finding the {args.num_top} molecules from {search_path} closest to {len(to_compare)} molecules. Saving to {out_path}')

    # Prepare the screening function
    screen_fun = partial(screen_molecules, to_compare=tuple(to_compare), radius=args.radius)
    update_wrapper(screen_fun, screen_molecules)

    # Make Parsl engine
    config, n_slots = parsl_config(args.compute)

    # Configure the file 
    if args.proxy_store == 'file':
        ps_file_dir = out_path / 'file-store'
        ps_file_dir.mkdir(exist_ok=True)
        store = ps.store.init_store(ps.store.STORES.FILE, name='file', store_dir=str(ps_file_dir))
    elif args.proxy_store == 'redis':
        store = ps.store.init_store(ps.store.STORES.REDIS, name='redis', hostname=args.redishost, port=args.redisport)
    else:
        raise ValueError('ProxyStore config not recognized: {}')

    # Make the task queues and task server
    client_q, server_q = make_queue_pairs(args.redishost, args.redisport, keep_inputs=False, serialization_method='pickle',
                                          proxystore_threshold=1000, proxystore_name=store.name)
    task_server = ParslTaskServer([screen_fun], server_q, config)

    # Make the thinker
    thinker = ScreenEngine(client_q, store, search_path, args.num_top, out_path, n_slots, args.molecules_per_chunk)

    # Run the program
    try:
        task_server.start()
        thinker.run()
    finally:
        client_q.send_kill_signal()

    task_server.join()
    store.cleanup()
