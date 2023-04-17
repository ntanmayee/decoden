#!/usr/bin/env python3

"""This module runs the DecoDen pipeline for peak consolidation
"""

import pandas as pd
import numpy as np
import json
from argparse import ArgumentParser
import os
from os.path import join, exists
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
from decoden.utils import get_blacklisted_regions_mask, load_files, print_message
from decoden.denoising.nmf import extract_mixing_matrix, extract_signal
from decoden.denoising.hsr import run_HSR


def run_consolidate(data_folder, 
                    files_reference, 
                    conditions, 
                    output_folder, 
                    blacklist_file,
                    alpha_W=0.01, 
                    alpha_H=0.001, 
                    control_cov_threshold=0.1, 
                    n_train_bins=50000, 
                    chunk_size=50000, 
                    seed=0):
    
    
    """`main` function that runs the internal pipeline for DecoDen

    Args:
        data_folder: Path to preprocessed data files in BED format
        files_reference: Path to JSON file with experiment conditions. 
                        If you used DecoDen for pre-processing, use the `experiment_conditions.json` file
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

    with open(files_reference, "r") as f:
        files = json.load(f)
    assert set(conditions) == set([files[i] for i in files]), 'Conditions do not match conditions in reference file. Perhaps there is an error in `--conditions` argument?'
    if not exists(output_folder):
        os.makedirs(output_folder)

    # with open(join(output_folder, "config.json"), "w") as f:
    #     json.dump(vars(args), f, indent=2)


    # Load data
    data, conditions_counts = load_files(files, data_folder, conditions)

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
    if Path(join(nmf_folder, "mixing_matrix.csv")).exists() and Path(join(nmf_folder, "mixing_matrix.pdf")).exists():
        print('`NMF` directory found. Using existing NMF results')
        mmatrix = pd.read_csv(join(nmf_folder, "mixing_matrix.csv"), index_col=0)
        wmatrix = pd.read_csv(join(nmf_folder, "signal_matrix.csv"))
        wmatrix.set_index(["seqnames", "start", "end"], inplace=True)

    else:
        if not exists(nmf_folder):
            os.makedirs(nmf_folder)

        # Extract mixing matrix
        mmatrix = extract_mixing_matrix(data_noBL, conditions, conditions_counts, alpha_W=alpha_W, 
                                    alpha_H=alpha_H, control_cov_threshold=control_cov_threshold, 
                                    n_train_bins=n_train_bins, seed=seed)
        mmatrix.to_csv(join(nmf_folder, "mixing_matrix.csv"))

        # visualise the mixing matrix
        plt.ioff()
        fig, ax = plt.subplots(figsize=(len(files),len(conditions)))
        sns.heatmap(mmatrix, annot=True, fmt=".2f", ax=ax)
        fig.savefig(join(nmf_folder, "mixing_matrix.pdf"), bbox_inches="tight")
        plt.close(fig)

        # Extract signal matrix
        wmatrix = extract_signal(data, mmatrix, conditions, chunk_size=chunk_size, alpha_W=alpha_W, seed=seed)
        wmatrix.reset_index().to_feather(join(nmf_folder, "signal_matrix.ftr"))

        # Sanity check plots
        plt.ioff()
        fig, axs = plt.subplots(1, 3, figsize=(10,20))
        np.random.seed(42)
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

        del data, data_noBL

    # Perform HSR to remove multiplicative noise
    hsr_df = run_HSR(wmatrix, mask, conditions)
    hsr_df.reset_index().to_feather(join(output_folder, "HSR_results.ftr"))
    
    print("\nDecoDen complete!")


if __name__=="__main__":
    print_message()

    parser = ArgumentParser()

    parser.add_argument("--data_folder", type=str, help="path to preprocessed data files in BED format", required=True)
    parser.add_argument("--output_folder", type=str, help='path to output directory', required=True)
    parser.add_argument("--files_reference", type=str, help='path to JSON file with experiment conditions. If you used DecoDen for pre-processing, use the `experiment_conditions.json` file', required=True)

    # First experimental condition should always correspond to the control/input samples
    parser.add_argument("--conditions", type=str, nargs="+", help='list of experimental conditions. First condition MUST correspond to the control/input samples.', required=True)

    parser.add_argument("--blacklist_file", type=str, help='path to blacklist file. Make sure to use the blacklist that is appropriate for the genome assembly/organism.')
    parser.add_argument("--control_cov_threshold", type=float, default=1.0, help='Threshold for coverage in control samples. Only genomic bins above this threshold will be used. It is recommended to choose a value larger than 1/bin_size.')
    parser.add_argument("--n_train_bins", type=int, default=300000, help='Number of genomic bins to be used for training')
    parser.add_argument("--chunk_size", type=int, default=100000, help='Chunk size for processing the signal matrix. Should be smaller than `n_train_bins`')
    parser.add_argument("--seed", type=int, default=42, help='Random state for reproducibility')

    parser.add_argument("--alpha_W", type=float, default=0.01, help='Regularisation for the signal matrix')
    parser.add_argument("--alpha_H", type=float, default=0.001, help='Regularisation for the mixing matrix')

    args = parser.parse_args()
    run_consolidate(args)
