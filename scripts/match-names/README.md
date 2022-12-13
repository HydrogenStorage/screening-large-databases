# Finding Mentions of Molecules in Papers

These notebooks show matching molecules to their mentions in the literature. 
There are several steps:

1. Getting the names of _all_ molecules by downloading them from PubChem
2. Marking which molecules of those 110M are ones we care about
3. Finding mentions of those molecules in a pile of papers
4. Saving the databases

We store the results of each of these steps in a MongoDB.
