import pandas as pd
import numpy as np
import json
import os
from os.path import exists, join
from sklearn.decomposition import NMF, non_negative_factorization
from tqdm import tqdm
from pathlib import Path
import warnings
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.ndimage.filters import maximum_filter1d

from decoden.utils import get_blacklisted_regions_mask, load_files, adjust_matrices, get_genome_size
from decoden.constants import *

def relax_mixing_matrix(data, mmat, smat, alpha_H=0.01, alpha_W=0.001):
    n_comps = len(mmat)
    mmat_1, smat_1, n_iter1 = non_negative_factorization(data.T,
                                                        n_components=n_comps, init="custom",
                                                        W=mmat.T, H=smat.T, update_H=False, alpha_W=alpha_H,
                                                        beta_loss=1, solver="mu")
    smat_2, mmat_2, n_iter2 = non_negative_factorization(data,
                                                        n_components=n_comps, init="custom",
                                                        W=smat_1.T, H=mmat_1.T, update_H=False, alpha_W=alpha_W,
                                                        beta_loss=1, solver="mu")
    return mmat_2, smat_2



def extract_mixing_matrix(data_df, conditions_list, conditions_counts_ref, alpha_W=0.01, alpha_H=0.001,
                          control_cov_threshold=1.0, n_train_bins=300000, seed=42):
    """Extract mixing matrix in the NMF step of DecoDen

    Args:
        data_df (np.array): Data matrix of all samples and conditions
        conditions_list (list): List of different experimental conditions
        conditions_counts_ref (list): Counts of different conditions
        alpha_W (float, optional): Regularisation for the signal matrix. Defaults to 0.01.
        alpha_H (float, optional): Regularisation for the mixing matrix. Defaults to 0.001.
        control_cov_threshold (float, optional): Minimum coverage for the training data for the NMF. Defaults to 1.0.
        n_train_bins (int, optional): Number of training bins for the extraction of the mixing matrix. Defaults to 300000.
        seed (int, optional): Random state for reproductibility. Defaults to 42.

    Returns:
        numpy.array: Mixing matrix from NMF
    """
    # Filter data to have sufficient control coverage
    dsel = data_df[(data_df[[c for c in data_df.columns if c.startswith(
        conditions_list[0])]] > control_cov_threshold).any(axis=1)]
    train_data = dsel.sample(n_train_bins, random_state=seed)
    

    # Extract unspecific signal from control samples
    model = NMF(n_components=1, init='random', random_state=0,
                beta_loss=1, solver="mu", alpha_W=alpha_W, alpha_H=alpha_H)
    control_cols = [
        c for c in data_df.columns if c.startswith(conditions_list[0])]
    treatment_cols = [c for c in data_df.columns if c not in control_cols]
    W_unspec = model.fit_transform(train_data.loc[:, control_cols])
    H_unspec = model.components_

    # Calculate unspecific signal coefficients for treatment
    # I swapped and transposed the matrices to make use of the update_H parameter
    W_treat_uns, H_treat_uns, n_iter = non_negative_factorization(train_data.loc[:, treatment_cols].T,
                                                                  n_components=1, init="custom",
                                                                  H=W_unspec.reshape(1, -1), update_H=False, alpha_W=alpha_W,
                                                                  beta_loss=1, solver="mu")
    H_unspec_coefs = np.concatenate(
        (H_unspec.flatten(), W_treat_uns.T.flatten()))

    # Subtract the unspecific contribution from the data
    # Cap the minimum value to 0 to account for predictions higher than the signal
    specific_data_mat = np.maximum(
        train_data.values - W_unspec.dot(H_unspec_coefs.reshape(1, -1)), 0)

    # For each modification, extract the specific signal components
    treatment_conditions = conditions_list[1:]
    n_replicates = [conditions_counts_ref[c] for c in conditions_list]
    n_control_replicates = n_replicates[0]
    histone_n_replicates = n_replicates[1:]
    signal_matrix = [W_unspec]
    mixing_matrix = np.zeros((len(conditions_list), np.sum(n_replicates)))
    mixing_matrix[0, :] = H_unspec_coefs

    ix = n_control_replicates

    cross_data_mat = []
    for i, modif in tqdm(enumerate(treatment_conditions)):

        c = histone_n_replicates[i]
        W_spec, H_spec, n_iter = non_negative_factorization(specific_data_mat[:, ix:ix+c],
                                                            n_components=1, beta_loss=1, solver="mu", alpha_W=alpha_W, alpha_H=alpha_H)

        mixing_matrix[i+1, ix:ix+c] = H_spec

        cross_data_mat.append(
            np.clip(specific_data_mat[:, ix:ix+c] - W_spec.dot(H_spec), 0, None)
        )
        signal_matrix.append(W_spec)
        ix += c
    
    signal_matrix = np.hstack(signal_matrix)
    cross_data_mat = np.hstack(cross_data_mat)
    ix = n_control_replicates

    if len(treatment_conditions) > 1:
        for i, modif in tqdm(enumerate(treatment_conditions)):
            c = histone_n_replicates[i]
            excl_idxs = list(range(n_control_replicates)) + list(range(ix, ix+c))
            cross_samples_idxs = np.array([i for i in range(mixing_matrix.shape[1]) if i not in excl_idxs]).astype(int)
            

            cross_cond_idxs = np.array([j for j in range(len(treatment_conditions)+1) if j>0 and j!=i+1]).astype(int)

            cond_ix = ix - n_control_replicates
            cross_mixing_coefs, cross_signal_mat, n_iter = non_negative_factorization(cross_data_mat.T[cond_ix:cond_ix+c,:],
                                                                        n_components=len(cross_cond_idxs), init="custom",
                                                                        H=signal_matrix[:,cross_cond_idxs].T , update_H=False, alpha_W=alpha_W,
                                                                        beta_loss=1, solver="mu")
            mixing_matrix[cross_cond_idxs, ix:ix+c] = cross_mixing_coefs.T
            ix += c
    
    # Add relaxation step, where we allow the matrices to vary jointly across modifications
    mixing_matrix, signal_matrix = relax_mixing_matrix(train_data, mixing_matrix, signal_matrix, alpha_H=alpha_H, alpha_W=alpha_W)
    
    mm = pd.DataFrame(mixing_matrix, index=[
                      UNSPECIFIC_SIGNAL_LABEL]+treatment_conditions, columns=train_data.columns)
    
    return mm


def extract_signal(data_df, mmatrix, conditions_list, chunk_size=100000, alpha_W=0.01, seed=42):
    """Use the mixing matrix to extract the signals. Can be done in chunks to fit in memory

    Args:
        data_df (np.array): Data matrix of all samples and conditions
        mmatrix (np.array): Mixing matrix
        conditions_list (list): List of different experimental conditions
        chunk_size (int, optional): Number of genomic bins to process in one chunk. Defaults to 100000.
        alpha_W (float, optional): Regularisation for the signal matrix. Defaults to 0.01.
        seed (int, optional): Random state for reproductibility. Defaults to 42.

    Returns:
        np.array: signal matrix from NMF
    """
    # Use the mixing matrix to extract the signals. Can be done in chunks to fit in memory
    processed_W = []
    for ck in tqdm(np.split(data_df.values, range(chunk_size, len(data_df)+chunk_size, chunk_size))):
        if len(ck) < 1:
            # print(ck)
            break
        W = np.random.uniform(0, 0.1, size=(len(ck), len(mmatrix)))
        ck_W, H, n_iter = non_negative_factorization(ck, W, mmatrix.values,
                                                     n_components=len(mmatrix), init='custom', random_state=seed, update_H=False,
                                                     beta_loss="kullback-leibler", solver="mu", alpha_W=alpha_W, max_iter=500
                                                     )
        processed_W.append(ck_W)
    processed_W = np.vstack(processed_W)
    processed_W = pd.DataFrame(
        processed_W, index=data_df.index, columns=[UNSPECIFIC_SIGNAL_LABEL]+conditions_list[1:])
    return processed_W


def run_NMF(files_reference,
        bin_size,
        genome_size, 
        conditions, 
        output_folder, 
        blacklist_file=None,
        alpha_W=0.01, 
        alpha_H=0.001, 
        control_cov_threshold=0.1, 
        n_train_bins=50000, 
        chunk_size=50000, 
        seed=0,
        plotting=True):
    
    
    """`main` function that runs the internal pipeline for DecoDen

    Args:
        files_reference: Path to JSON file with experiment conditions. 
                        If you used DecoDen for pre-processing, use the `experiment_conditions.json` file
        bin_size: bin size
        genome_size: genome size to calculate background
        conditions: list of experimental conditions. First condition MUST correspond to the control/input samples.
        output_folder: Path to output directory
        blacklist_file: Path to blacklist file. Make sure to use the blacklist that is appropriate for the genome assembly/organism.
        control_cov_threshold : Threshold for coverage in control samples. Only genomic bins above this threshold will be used. 
                                It is recommended to choose a value larger than 1/bin_size.
        n_train_bins: Number of genomic bins to be used for training
        chunk_size: Chunk size for processing the signal matrix. Should be smaller than `n_train_bins`
        alpha_W: Regularisation for the signal matrix
        alpha_H: Regularisation for the mixing matrix
    """

    # validate arguments
    assert control_cov_threshold > 0, 'Control coverage threshold must be greater than 0'
    assert n_train_bins > 0, 'Number of training bins must be greater than 0'
    assert chunk_size > 0, 'Chunk size must be greater than 0'
    assert alpha_W > 0, '`alpha_W` must be greater than 0'
    assert alpha_H > 0, '`alpha_H` must be greater than 0'
    assert genome_size is not None, '`genome_size` cannnot be None'

    with open(files_reference, "r") as f:
        files = json.load(f)
        
    data_folder = os.path.dirname(files_reference) 
    assert set(conditions) == set([v['condition'] for v in files.values()]), 'Conditions do not match conditions in reference file. Perhaps there is an error in `--conditions` argument?'
    if not exists(output_folder):
        os.makedirs(output_folder)

    # Load data
    data, conditions_counts = load_files(files, data_folder, conditions)

    # Get bin_size
    a_key = list(files.keys())[0]
    bin_size = files[a_key]['bin_size']

    # Filter BL regions
    if blacklist_file is not None:
        bl_regions = pd.read_csv(blacklist_file,
                            sep="\t", header=None, names=["seqnames", "start", "end", "type"])
        mask = get_blacklisted_regions_mask(data, bl_regions)
    else:
        warnings.warn('No blacklist supplied, using all data for next steps')
        mask = np.ones(len(data)).astype(bool)
    data_noBL = data[mask]
    nmf_folder = join(output_folder, "NMF")

    # skip NMF step if already computed
    if Path(join(nmf_folder, "mixing_matrix.csv")).exists()  and Path(join(nmf_folder, "signal_matrix.ftr")).exists():
        print('`NMF` directory found. Using existing NMF results')
        mmatrix = pd.read_csv(join(nmf_folder, "mixing_matrix.csv"), index_col=0)
        wmatrix = pd.read_feather(join(nmf_folder, "signal_matrix.ftr"))
        wmatrix.set_index(["seqnames", "start", "end"], inplace=True)

    else:
        os.makedirs(nmf_folder, exist_ok=True)

        # Extract mixing matrix
        mmatrix = extract_mixing_matrix(data_noBL, conditions, conditions_counts, alpha_W=alpha_W, 
                                    alpha_H=alpha_H, control_cov_threshold=control_cov_threshold, 
                                    n_train_bins=n_train_bins, seed=seed)
        
        # Extract signal matrix
        wmatrix = extract_signal(data, mmatrix, conditions, chunk_size=chunk_size, alpha_W=alpha_W, seed=seed)
        
        # Rescale matrixes to have comparable signals
        # mmatrix, wmatrix = adjust_matrices(mmatrix, wmatrix, q=0.98)

        # Smooth the chromatin bias (unspecific signal)
        chrom_bias = wmatrix.iloc[ : , 0]
        slocal_background = maximum_filter1d(chrom_bias, size=int(1000/bin_size), mode='nearest') * (bin_size / 1000)
        llocal_background = maximum_filter1d(chrom_bias, size=int(10000/bin_size), mode='nearest') * (bin_size / 10000)
        genome_background = 1/get_genome_size(genome_size)
        chrom_bias = np.maximum(chrom_bias, slocal_background)
        chrom_bias = np.maximum(chrom_bias, llocal_background)
        chrom_bias = np.maximum(chrom_bias, genome_background)
        wmatrix.iloc[ : , 0] = chrom_bias

        # Save results
        mmatrix.to_csv(join(nmf_folder, "mixing_matrix.csv"))
        wmatrix.reset_index().to_feather(join(nmf_folder, "signal_matrix.ftr"))
        
        
        if plotting:
            # visualise the mixing matrix
            plt.ioff()
            fig, ax = plt.subplots(figsize=(len(files),len(conditions)))
            sns.heatmap(mmatrix, annot=True, fmt=".2f", ax=ax)
            fig.savefig(join(nmf_folder, "mixing_matrix.pdf"), bbox_inches="tight")
            plt.close(fig)

            # Sanity check plots
            plt.ioff()
            fig, axs = plt.subplots(1, 3, figsize=(10,20))
            np.random.seed(seed)
            ixs = np.sort(np.random.choice(range(len(data)), size=500, replace=False))
            sns.heatmap(data.values[ixs], ax=axs[0])
            axs[0].set_title("Data Matrix")
            sns.heatmap(wmatrix.values.dot(mmatrix.values)[ixs], ax=axs[1])
            axs[1].set_title("Reconstruction")
            sns.heatmap(wmatrix.values[ixs], ax=axs[2])
            axs[2].set_title("Signal Matrix")
            axs[0].set_yticklabels([])
            axs[0].set_yticks([])
            axs[1].set_yticklabels([])
            axs[1].set_yticks([])
            axs[2].set_yticklabels([])
            axs[2].set_yticks([])
            fig.savefig(join(nmf_folder, "signal_matrix_sample.pdf"), bbox_inches="tight")
            plt.close(fig)

    return wmatrix, mmatrix, data_noBL, mask, conditions_counts