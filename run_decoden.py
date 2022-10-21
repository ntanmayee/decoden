"""This module runs the DecoDen pipeline
"""

import pandas as pd
import numpy as np
import json
from argparse import ArgumentParser
import os
from os.path import join, exists
from tqdm import tqdm
from decoden.utils import get_blacklisted_regions_mask, load_files, print_message
from decoden.functions import extract_mixing_matrix, extract_signal, run_HSR

import matplotlib.pyplot as plt
import seaborn as sns

def main(args):
    """`main` function that runs the internal pipeline for DecoDen

    Args:
        args : arguments from ArgumentParser
    """

    with open(args.files_reference, "r") as f:
        files = json.load(f)
    if not exists(args.output_folder):
        os.makedirs(args.output_folder)

    with open(join(args.output_folder, "config.json"), "w") as f:
        json.dump(vars(args), f, indent=2)


    # Load data
    conditions = args.conditions
    data, conditions_counts = load_files(files, args.data_folder, conditions)

    # Filter BL regions
    if args.blacklist_file is not None:
        bl_regions = pd.read_csv(args.blacklist_file,
                            sep="\t", header=None, names=["seqnames", "start", "end", "type"])
        mask = get_blacklisted_regions_mask(data, bl_regions)
    else:
        mask = np.ones(len(data)).astype(bool)
    data_noBL = data[mask]


    nmf_folder = join(args.output_folder, "NMF")
    if not exists(nmf_folder):
        os.makedirs(nmf_folder)


    # Extract mixing matrix
    mmatrix = extract_mixing_matrix(data_noBL, conditions, conditions_counts, alpha_W=args.alpha_W, 
                                alpha_H=args.alpha_H, control_cov_threshold=args.control_cov_threshold, 
                                n_train_bins=args.n_train_bins, seed=args.seed)
    mmatrix.to_csv(join(nmf_folder, "mixing_matrix.csv"))
    # mmatrix = pd.read_csv(join(nmf_folder, "mixing_matrix.csv"), index_col=0)

    plt.ioff()
    fig, ax = plt.subplots(figsize=(len(files),len(conditions)))
    sns.heatmap(mmatrix, annot=True, fmt=".2f", ax=ax)
    fig.savefig(join(nmf_folder, "mixing_matrix.pdf"), bbox_inches="tight")
    plt.close(fig)


    # Extract signal matrix
    wmatrix = extract_signal(data, mmatrix, conditions, chunk_size=args.chunk_size, alpha_W=args.alpha_W, seed=args.seed)
    wmatrix.reset_index().to_feather(join(nmf_folder, "signal_matrix.ftr"))

    # wmatrix = pd.read_csv(join(nmf_folder, "signal_matrix.csv"))
    # wmatrix.set_index(["seqnames", "start", "end"], inplace=True)

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
    hsr_df.reset_index().to_feather(join(args.output_folder, "HSR_results.ftr"))
    
    print("DecoDen complete!")





if __name__=="__main__":
    print_message()

    parser = ArgumentParser()

    parser.add_argument("--data_folder", type=str, default="../DecoDen_GV/data/shallow_e114_200bp_bedGraph_files/", help="path to preprocessed data files in BED format")
    parser.add_argument("--output_folder", type=str, default="../DecoDen_GV/outputs/shallow_e114_200bp_results_centering", help='path to output directory')
    parser.add_argument("--files_reference", type=str, default="../DecoDen_GV/data/shallow_e114_200bp_bedGraph_files/sample_files.json", help='path to JSON file with experiment conditions. If you used DecoDen for pre-processing, use the `experiment_conditions.json` file')
    parser.add_argument("--blacklist_file", type=str, default="../DecoDen_GV/data/annotations/hg19-blacklist.v2.bed", help='path to blacklist file. Make sure to use the blacklist that is appropriate for the genome assembly/organism.')

    # First experimental condition should always correspond to the control/input samples
    parser.add_argument("--conditions", type=str, nargs="+", default=["control", "H3K27me3", "H3K4me3"], help='list of experimental conditions. First condition MUST correspond to the control/input samples, otherwise DecoDen will not work correctly.')

    parser.add_argument("--control_cov_threshold", type=float, default=1.0, help='Threshold for coverage in control samples. Only genomic bins above this threshold will be used. It is recommended to choose a value larger than 1/bin_size.')
    parser.add_argument("--n_train_bins", type=int, default=300000, help='Number of genomic bins to be used for training')
    parser.add_argument("--chunk_size", type=int, default=100000, help='Chunk size for processing the signal matrix. Should be smaller than `n_train_bins`')
    parser.add_argument("--seed", type=int, default=42, help='Random state for reproducability')

    parser.add_argument("--alpha_W", type=float, default=0.01, help='Regularisation for the signal matrix')
    parser.add_argument("--alpha_H", type=float, default=0.001, help='Regularisation for the mixing matrix')

    args = parser.parse_args()
    main(args)