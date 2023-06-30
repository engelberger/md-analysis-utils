# Documentation for Molecular Dynamics Analysis Code

This code is designed to analyze molecular dynamics (MD) simulations, specifically focusing on distance calculations between specified atoms in the system. It uses the MDAnalysis library to handle trajectory and topology data from MD simulations. 

The code is organized into several functions, each with a specific role in the analysis process. Here is a brief overview of each function:

1. `createSystemsDict(path, tag)`: This function creates a dictionary of systems given a path and a tag. The dictionary keys are the system names and the values are the system paths.

2. `createUniverseDict(systemsDict, systemName)`: This function creates a dictionary of MDAnalysis Universe objects for each replica of a given system.

3. `createSelection(replicaUniverse, residueNumber1, residueNumber2, atomList1, atomList2)`: Given a Universe, two residue numbers and two atom lists, this function creates a MDAnalysis selection.

4. `calculateDistances(selectionDict)`: This function calculates the distances between the first atom in the list and the rest of the atoms in the list, using the selections created by `createSelection()`.

5. `plotDistances(distancesDict, systemName, replica)`: This function plots the calculated distances for each replica of a system.

6. `calculateDistancesParallelFunction(systemName, replicaNumber, serineResidueNumber, petResidueNumber, serine_atom_list, pet_atoms_list, systemsFolder, figuresFolder)`: This function is a parallelized version of the distance calculation and plotting process. It uses multiprocessing to perform these operations in parallel for each replica of a system.

7. `loadDistances(systemsFolder,systemName, replicaNumber)`: This function loads the calculated distances from a pickle file and returns a dictionary of distances.

8. `createDataframeFromDistancesDict(distancesDict)`: This function creates a pandas DataFrame from the distances dictionary.

9. `plotDataFrameDistances(df, systemName, replicaNumber)`: This function plots the distances with rolling average and standard deviation for a system and a replica.

10. `loadAndPlotReplicaSystemFromPickle(total_time, systemName, dfList, replicaNumber)`: This function loads the distances from a pickle file, creates a DataFrame, calculates rolling averages and standard deviations, and plots the distances for each system and replica.

The script also includes a main section that uses multiprocessing to perform the distance calculations and plotting for each system and replica in parallel.

The output of this script includes plots of distances for each system and replica, saved in the specified figures folder, and pickle files containing the calculated distances, saved in the specified pickle folder. 

To use this script, you will need to have MDAnalysis installed and your MD simulation data in the correct format and location. You will also need to specify the correct atom and residue numbers for your system.



## Requirements

```bash
pip install -r requirements.txt
```

Please note that `re`, `os`, and `math` are part of the Python Standard Library, so they don't need to be included in the requirements file.

Also, `%matplotlib inline` is a magic command for Jupyter notebooks and not a package, so it doesn't need to be included in the requirements file.

To install these packages using the `requirements.txt` file, you can use the following command: