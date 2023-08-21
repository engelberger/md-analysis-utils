import logging
from tqdm.auto import tqdm
from typing import Tuple, List, Any
import matplotlib.pyplot as plt
import numpy as np
import pickle
import pyemma
from deeptime.clustering import KMeans
from deeptime.decomposition import TICA
from deeptime.markov import TransitionCountEstimator, MaximumLikelihoodMSM
from deeptime.plots import plot_implied_timescales
from deeptime.util.validation import implied_timescales
from pyemma.util.contexts import settings
from sklearn.decomposition import PCA
import argparse
import yaml

def load_data(pdb_file: str, xtc_files: List[str], selection: str) -> Any:
    """
    Load trajectory data for a given protein structure and feature selection.

    Parameters
    ----------
    pdb_file : str
        Path to the pdb file.
    xtc_files : List[str]
        List of paths to the xtc files.
    selection : str
        Selection string for the protein residues/features to be considered.

    Returns
    -------
    Any
        Loaded trajectory data.
    """
    feat = pyemma.coordinates.featurizer(pdb_file)
    selection = feat.select(selection)
    feat.add_distances(selection, periodic=False)
    data = pyemma.coordinates.load(xtc_files, features=feat)
    return data

def perform_pca(data: Any, n_components: int = 2) -> Tuple[PCA, List[np.ndarray]]:
    """
    Perform Principal Component Analysis (PCA) on the trajectory data.

    Parameters
    ----------
    data : Any
        Trajectory data.
    n_components : int, optional
        Number of principal components to be computed, by default 2.

    Returns
    -------
    Tuple[PCA, List[np.ndarray]]
        Fitted PCA model and transformed trajectory data.
    """
    pca = PCA(n_components=n_components)
    pca_output = [pca.fit_transform(traj) for traj in data]
    return pca, pca_output

def perform_tica(data: Any, lag: int, n_components: int = 2) -> Tuple[TICA, List[np.ndarray]]:
    """
    Perform Time-lagged Independent Component Analysis (TICA) on the trajectory data.

    Parameters
    ----------
    data : Any
        Trajectory data.
    lag : int
        Time lag to be used for the TICA computation.
    n_components : int, optional
        Number of independent components to be computed, by default 2.

    Returns
    -------
    Tuple[TICA, List[np.ndarray]]
        Fitted TICA model and transformed trajectory data.
    """
    tica_estimator = TICA(lagtime=lag, dim=n_components)
    tica = tica_estimator.fit_fetch(data)
    tica_output = [tica.transform(traj) for traj in data]
    return tica, tica_output

def cluster_data(data: List[np.ndarray], n_clusters: int, max_iter: int = 100, n_jobs: int = 2) -> Tuple[KMeans, List[int]]:
    """
    Cluster the trajectory data using K-means algorithm.

    Parameters
    ----------
    data : List[np.ndarray]
        Trajectory data.
    n_clusters : int
        Number of clusters to be formed.
    max_iter : int, optional
        Maximum number of iterations for the K-means algorithm, by default 100.
    n_jobs : int, optional
        Number of jobs to be used for the computation, by default 2.

    Returns
    -------
    Tuple[KMeans, List[int]]
        Fitted KMeans model and list of cluster indices for each data point.
    """
    cls = KMeans(n_clusters, max_iter=max_iter, n_jobs=n_jobs).fit(np.concatenate(data)[::10]).fetch_model()
    dtrajs = [cls.transform(traj) for traj in data]
    return cls, dtrajs

def compute_implied_timescales(dtrajs: List[int], lags: List[int]) -> np.ndarray:
    """
    Compute the implied timescales for a given list of lag times.

    Parameters
    ----------
    dtrajs : List[int]
        List of cluster indices for each data point.
    lags : List[int]
        List of lag times to be used for the computation.

    Returns
    -------
    np.ndarray
        Computed implied timescales.
    """
    return implied_timescales([MaximumLikelihoodMSM(lagtime=lag).fit_fetch(dtrajs) for lag in lags])

def save_to_pickle(data: Any, filename: str) -> None:
    """
    Save data to a pickle file.

    Parameters
    ----------
    data : Any
        Data to be saved.
    filename : str
        Path to the pickle file.
    """
    with open(filename, 'wb') as f:
        pickle.dump(data, f)

def load_from_pickle(filename: str) -> Any:
    """
    Load data from a pickle file.

    Parameters
    ----------
    filename : str
        Path to the pickle file.

    Returns
    -------
    Any
        Loaded data.
    """
    with open(filename, 'rb') as f:
        return pickle.load(f)

def save_tics_to_txt(tica_output: List[np.ndarray], sysName: str) -> None:
    """
    Save the first two TICA components to text files.

    Parameters
    ----------
    tica_output : List[np.ndarray]
        Transformed trajectory data after TICA.
    sysName : str
        System name to be used for the output file names.
    """
    # check if the tics folder exists 
    # ODO
    np.savetxt(f'./tics/{sysName}_tica_tic1.txt',np.concatenate(tica_output)[:, 0], delimiter='\n', fmt='%.10f')
    np.savetxt(f'./tics/{sysName}_tica_tic2.txt',np.concatenate(tica_output)[:, 1], delimiter='\n', fmt='%.10f')

def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("config", help="Path to the configuration file")
    return parser.parse_args()

def load_config(config_path: str) -> dict:
    """
    Load configuration from a yaml file.

    Parameters
    ----------
    config_path : str
        Path to the yaml configuration file.

    Returns
    -------
    dict
        Loaded configuration.
    """
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

if __name__ == "__main__":
    args = parse_args()
    config = load_config(args.config)

    pdb_file = config['pdb_file']
    xtc_files = config['xtc_files']
    sysName = config['sysName']
    lag = config['lag']
    n_components = config['n_components']
    n_clusters = config['n_clusters']
    lags = config['lags']
    selection = config['selection']

    logging.info("Loading data...")
    data = load_data(pdb_file, xtc_files, selection)
    save_to_pickle(data, f'{sysName}.pickle')

    logging.info("Performing PCA...")
    pca, pca_output = perform_pca(data, n_components)
    logging.info("Performing TICA...")
    tica, tica_output = perform_tica(data, lag, n_components)

    logging.info("Clustering data...")
    cls_pca, dtrajs_pca = cluster_data(pca_output, n_clusters)
    cls_tica, dtrajs_tica = cluster_data(tica_output, n_clusters)

    logging.info("Computing implied timescales...")
    its_pca = compute_implied_timescales(dtrajs_pca, lags)
    its_tica = compute_implied_timescales(dtrajs_tica, lags)

    save_to_pickle(tica_output, f'{sysName}_tica_output.pickle')
    save_to_pickle(pca_output, f'{sysName}_pca_output.pickle')

    save_tics_to_txt(tica_output, sysName)
