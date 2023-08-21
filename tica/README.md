# Protein Dynamics Analysis Pipeline

This Python script provides a pipeline for analysis of protein dynamics. It uses Principal Component Analysis (PCA) and Time-lagged Independent Component Analysis (TICA) to analyze the trajectory data of a protein. The data is then clustered using the K-means algorithm, and the implied timescales are computed.

## Features

- Load trajectory data for a given protein structure and feature selection
- Perform Principal Component Analysis (PCA) on the trajectory data
- Perform Time-lagged Independent Component Analysis (TICA) on the trajectory data
- Cluster the trajectory data using the K-means algorithm
- Compute the implied timescales for a given list of lag times
- Save and load data from pickle files
- Save the first two TICA components to text files
- Parse command-line arguments and load configuration from a yaml file

## Usage

1. Prepare a configuration file in yaml format. Here is an example:

    ```yaml
    pdb_file: /path/to/pdb_file
    xtc_files: 
      - /path/to/xtc_file1
      - /path/to/xtc_file2
    sysName: 140s_LightChain_amd
    lag: 100
    n_components: 2
    n_clusters: 200
    lags: [100, 200, 300, 400, 500, 1000]
    selection: "name CA and resid 140 to 156 and not name H"
    ```

    The `selection` field is a string representing the protein residues/features to be considered.

2. Run the script with the configuration file as a command-line argument:

    ```bash
    python protein_dynamics.py config.yml
    ```

3. The results will be saved in pickle files and text files. The pickle files contain the PCA and TICA transformed trajectory data, and the text files contain the first two TICA components.

## Requirements

- Python 3
- tqdm
- pyemma
- deeptime
- sklearn
- pickle
- yaml
- argparse
- numpy
- matplotlib

## Author

Felipe Engelberger-Aliaga