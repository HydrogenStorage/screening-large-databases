# Screen Search Space

We plan to use large resources of molecules as the basis for a search and want to downselect the list of molecules to only those that are most interesting for later analysis.
Our approach is to find the molecules that are most similar to known molecules using Tanimoto similarity with Morgan fingerprints.

Running the screens on large search spaces is easy with chemoinformatics toolkits, like RDKit, but running them in parallel has a few challenges.
The largest is that storing the entire search space in memory is problematic.
So, we implement an out-of-core strategy where we gradually read molecules from the search space and send them out to workers as tasks complete.
That way, we do not have the whole dataset in memory at any one time.

## Running the Program

The `screen.py` program implements the parallel screening algorithm with [Colmena](https://colmena.rtfd.org). 
It requires a path to the search space of molecules and a description of the molecules we want to compare to as a `.smi` file.

Call `python screen.py --help` for the full options. 
