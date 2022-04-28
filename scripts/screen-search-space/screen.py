import sys
from dataclasses import dataclass, field
from functools import partial, update_wrapper
from hashlib import sha256
from queue import Queue
from threading import Event
from typing import List, Tuple, Union
from pathlib import Path
from csv import DictWriter
from datetime import datetime
from glob import glob
import argparse
import logging
import heapq
import gzip
import json

from colmena.redis.queue import make_queue_pairs, ClientQueues
from colmena.task_server import ParslTaskServer
from colmena.thinker import BaseThinker, agent, ResourceCounter
from parsl import Config, HighThroughputExecutor
from parsl.addresses import address_by_hostname
from parsl.launchers import AprunLauncher
from parsl.providers import CobaltProvider
from tqdm import tqdm
import proxystore as ps
import numpy as np


def parsl_config(name: str) -> Tuple[Config, int]:
    """Make the compute resource configuration

    Args:
        name: Name of the desired configuration
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
            retries=2,
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
                        launcher=AprunLauncher(overrides='-d 64 --cc depth -j 2'),
                        worker_init='''
module load miniconda-3
conda activate /lus/theta-fs0/projects/CSC249ADCD08/carbon-free-ldrd/env''',
                    ),
                )]
            ), 64 * 8 * 2
    elif name == 'theta-full':
        return Config(
            retries=2,
            executors=[HighThroughputExecutor(
                    address=address_by_hostname(),
                    label="full",
                    max_workers=64,
                    prefetch_capacity=192,
                    cpu_affinity='block',
                    provider=CobaltProvider(
                        account='CSC249ADCD08',
                        nodes_per_block=512,
                        scheduler_options='#COBALT --attrs enable_ssh=1',
                        walltime='09:00:00',
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
                )]
            ), 512 * 64 * 5
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
        screen_paths: Path to molecules to be screened
        store: Store used to pass inputs from thinker to worker
        output_dir: Path to the output directory
        slot_count: Number of execution slots
        chunk_size: Number of molecules per chunk
        read_on_nodes: Whether to read files on computes nodes rather than main
    """

    def __init__(self,
                 queues: ClientQueues,
                 store: ps.store.remote.RemoteStore,
                 screen_paths: List[Path],
                 target_count: int,
                 output_dir: Path,
                 slot_count: int,
                 chunk_size: int,
                 read_on_nodes: bool
                 ):
        # Init the base class
        super().__init__(queues, ResourceCounter(slot_count, ['screen']))
        self.rec.reallocate(None, 'screen', 'all')

        # Store the input and output information
        self.read_on_nodes = read_on_nodes
        self.screen_paths = screen_paths
        self.output_dir = output_dir
        self.chunk_size = chunk_size
        self.target_count = target_count
        self.store = store

        # Queue to store ready-to-compute chunks
        self.queue_depth = slot_count * 8
        self.screen_queue = Queue(self.queue_depth)

        # Queue to store results ready to add to the list
        self.result_queue = Queue()

        # Variables to hold the outputs
        self.best_mols: List[_PrioritizedItem] = []  # Will be used as a heap
        self.current_threshold = 0

        # Things to know if we are done
        self.all_read = Event()
        self.queue_filled = Event()
        self.total_chunks = 0
        self.total_molecules = 0
        
        # Recording for performance information
        self.start_time = np.inf  # Use the time the first compute starts

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
            
            # Mark that the queue is completely filled
            if self.total_chunks == self.queue_depth:
                self.logger.info('Queue is filled. Submit tasks now!')
                self.queue_filled.set()
            
        # Loop over all files
        chunk = []
        for i, path in enumerate(self.screen_paths): 
            self.logger.info(f'Reading from file {i+1}/{len(self.screen_paths)}: {path}')

            if self.read_on_nodes:
                # If read on nodes, count the number of molecules so we have some performance tracking
                with path.open() as fp:
                    for _ in fp:
                        self.total_molecules += 1
                self.total_chunks += 1

                # Submit only the path the queue
                self.screen_queue.put((path.absolute(), f'file-{self.total_chunks}'))

                # If needed Mark that the queue is completely filled
                if self.total_chunks == self.queue_depth:
                    self.logger.info('Queue is filled. Submit tasks now!')
                    self.queue_filled.set()
            else:
                # If not, read the file and create chunks
                with path.open() as fp:
                    for line in fp:
                        try:
                            if path.suffix == '.csv':
                                # The line is comma-separated with the last entry as the string
                                _, smiles = line.strip().rsplit(",", 1)
                            elif path.suffix == '.smi':
                                # The line is just a SMILES string
                                smiles = line.strip()
                            else:
                                raise ValueError(f'Extension "{path.suffix}" not recognized for "')
                        except Exception:
                            continue
                        self.total_molecules += 1

                        # Add to the chunk and submit if we hit the target size
                        chunk.append(smiles)
                        if len(chunk) >= self.chunk_size:
                            _proxy_and_push(chunk)
                            chunk = []

        if not self.read_on_nodes:
            # Submit whatever remains
            _proxy_and_push(chunk)

        # Put a None at the end to signal we are done
        self.queue_filled.set()  # Make sure the task submitter starts
        self.screen_queue.put(None)

        # Mark that we are done reading
        self.logger.info(f'Finished reading {self.total_molecules} molecules and submitting {self.total_chunks} blocks')
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
            chunk, key = msg

            # Submit once we have resources available
            self.rec.acquire("screen", 1)
            self.queues.send_inputs(chunk, self.current_threshold, method='screen_molecules',
                                    task_info={'key': key, 'threshold': self.current_threshold})

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
                if not result.success:
                    self.logger.error(f'Task failed. Traceback: {result.failure_info.traceback}')
                    raise ValueError('Failed task')

                # Push the result to be processed
                self.result_queue.put(result.value)
                num_recorded += 1
                    
                # Update the start time
                self.start_time = min(self.start_time, result.time_compute_started)
                
                # Remove the chunk from the proxystore
                self.store.evict(result.task_info['key'])
                
                # Save the JSON result
                print(result.json(exclude={'inputs', 'value'}), file=fq)

                # Print a status message
                if self.all_read.is_set():
                    self.logger.info(f'Received task {num_recorded}/{self.total_chunks}. Processing backlog: {self.result_queue.qsize()}')
                else:
                    self.logger.info(f'Received task {num_recorded}/???. Processing backlog: {self.result_queue.qsize()}')

        # Send a "done" message
        self.result_queue.put(None)

        # Print the final status
        run_time = datetime.now().timestamp() - self.start_time
        self.logger.info(f'Completed storing all results. {run_time:.2f} s. Overall evaluation rate: {self.total_molecules / run_time:.3e} mol/s')

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
            status = f'Number received: {len(result)}. New threshold: {self.current_threshold:.4f}'
            if self.all_read.is_set():
                self.logger.info(f'Processed task {count}/{self.total_chunks}. {status}')
            else:
                self.logger.info(f'Processed task {count}/???. {status}')

        # Write the data out to disk
        self.logger.info(f'Completed processing all results. Output list size: {len(self.best_mols)}. Sorting results now')
        
        self.best_mols.sort(reverse=True)
        self.logger.info('Finished sorting. Writing to disk')

        with gzip.open(self.output_dir / 'best_molecules.csv.gz', 'wt') as fp:
            writer = DictWriter(fp, ['smiles', 'inchi', 'score', 'similarities'])
            writer.writeheader()
            for record in self.best_mols:
                entry = json.loads(record.item)
                entry['score'] = record.priority
                writer.writerow(entry)
        self.logger.info('Finished writing results to disk')
        
        # Write out the total run time.
        run_time = datetime.now().timestamp() - self.start_time
        self.logger.info(f'Runtime {run_time:.2f} s. Overall evaluation rate: {self.total_molecules / run_time:.3e} mol/s')


def screen_molecules(to_screen: Union[List[str], Path], min_similarity: float, to_compare: Tuple[str], radius: int = 4) -> List[Tuple[float, str]]:
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
    from pathlib import Path
    from rdkit import Chem
    from rdkit.Chem import DataStructs, AllChem
    
    # Turn off logging
    import rdkit.rdBase as rkrb
    import rdkit.RDLogger as rkl
    logger = rkl.logger()
    logger.setLevel(rkl.ERROR)
    rkrb.DisableLog('rdApp.error')

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

    # If `to_screen` is a file, read the molecules from it
    if isinstance(to_screen, Path):
        assert to_screen.suffix == '.smi', 'File reading only supported from .smi files'
        with to_screen.open() as fp:
            to_screen = [x.strip() for x in fp.readlines()]

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
        mol_similarity = DataStructs.BulkTanimotoSimilarity(fp, compare_mols)
        max_similarity = max(mol_similarity)
        
        # If it is an exact match, skip it
        if max_similarity == 1.0:
            continue

        # If it doesn't meet the threshold, don't return it
        if max_similarity < min_similarity:
            continue

        # If not, add it to the output dictionary
        passed.append((max_similarity, json.dumps({
            'smiles': smiles, 'inchi': Chem.MolToInchi(mol),
            'similarities': json.dumps(mol_similarity)
        })))

    return passed


if __name__ == '__main__':
    # User inputs
    parser = argparse.ArgumentParser()
    parser.add_argument('--overwrite', action='store_true', help='Whether to overwrite a previous run')

    group = parser.add_argument_group(title='Compute Configuration', description='Settings related to how compute is executed')
    group.add_argument("--redishost", default="127.0.0.1", help="Address at which the redis server can be reached")
    group.add_argument("--redisport", default="6379", help="Port on which redis is available")
    group.add_argument("--molecules-per-chunk", default=50000, type=int, help="Number molecules per screening task")
    group.add_argument('--compute', default='local', help='Resource configuration')
    group.add_argument('--proxy-store', default='file', help='What kind of proxy store to use')
    group.add_argument('--read-on-nodes', action='store_true', help='Read the molecule files on nodes instead of master')

    group = parser.add_argument_group(title='Search Configuration', description='How we screen for the top molecules')
    group.add_argument('--search-space', nargs='+', help='Path to molecules to be screened. If only one path provided, we assume it is a glob string')
    group.add_argument('--name', help='Name for the screening selection')
    group.add_argument('--num-top', default=1000000, type=int, help='Number of most-similar molecules to find')
    group.add_argument('--comparison-molecules', required=True, help='Path to the molecule(s) to compare against. Must be line-delimited SMILES')
    group.add_argument('--radius', default=4, type=int, help='Radius of the fingerprint')

    # Parse the arguments
    args = parser.parse_args()
    run_params = args.__dict__

    # Check that the search path exists
    if len(args.search_space) == 1:
        search_paths = [Path(x) for x in glob(args.search_space[0], recursive=True)]
    else:
        search_paths = [Path(x) for x in args.search_space]

    # Load in the molecules to be screened
    with open(args.comparison_molecules) as fp:
        to_compare = [x.strip() for x in fp]
    
    # Create an output directory with the name of the directory
    run_hash = sha256()
    run_hash.update(str(run_params).encode())
    run_hash.update(str(to_compare).encode())
    name = args.name if args.name is not None else search_paths[0].name[:-4]
    out_path = Path().joinpath('runs', f'{name}-top{args.num_top}-{run_hash.hexdigest()[:6]}')
    out_path.mkdir(parents=True, exist_ok=args.overwrite)
    
    # Store the run parameters
    with open(out_path / 'run-config.json', 'w') as fp:
        json.dump(run_params, fp, indent=2)
    with open(out_path / 'to-compare.smi', 'w') as fp:
        for s in to_compare:
            print(s, file=fp)
    
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
    logging.info(f'Finding the {args.num_top} molecules that are closest to {len(to_compare)} examples out of {len(search_paths)} files. Saving to {out_path}')

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
    thinker = ScreenEngine(client_q, store, search_paths, args.num_top, out_path, n_slots, args.molecules_per_chunk, args.read_on_nodes)

    # Run the program
    try:
        task_server.start()
        thinker.run()
    finally:
        client_q.send_kill_signal()

    task_server.join()
    store.cleanup()
