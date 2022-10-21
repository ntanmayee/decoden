import argparse
from argparse import Namespace
from pathlib import Path
import json
import os
from preprocess.logger import logger
from decoden.utils import print_message
from run_preprocess import run
from run_decoden import main


def extract_conditions(json_file):
    """Extract list of different experimental conditions, to pass to run_decoden.py. Assumes that `control` is the first condition.

    Args:
        json_file (string): path to `experiment_conditions.json` generated from run_preprocess.py

    Returns:
        list: list of different experimental conditions
    """
    json_object = json.load(open(json_file))
    
    conditions = []
    for key in json_object:
        value = json_object[key]
        if value not in conditions:
            conditions.append(value)
    logger.info(conditions)
    return conditions 


if __name__ == '__main__':
    print_message()
    
    parser = argparse.ArgumentParser()

    # arguments for preprocessing
    parser.add_argument('-i', "--input_csv", required=True, help='path to CSV file with information about experimental conditions. Must contain `filepath`, `exp_name` and `is_control` columns. Control/input should be the first condition. Input files can be in BED/BAM format.')
    parser.add_argument('-bs', "--bin_size", default=200, type=int, help='size of genomic bin for tiling. Recommended value is 10-200. Smaller bin size increases space and runtime, larger binsizes may occlude small variations. Default: 200')
    parser.add_argument('-n', "--num_jobs", default=1, type=int, help='Number of parallel jobs for preprocessing. Default: 1')
    parser.add_argument('-o', "--out_dir", required=True, help='path to directory where all output files will be written')

    # arguments from decoden
    parser.add_argument('-bl', "--blacklist_file", help='path to blacklist file')
    parser.add_argument("--control_cov_threshold", type=float, default=1.0, help='Threshold for coverage in control samples. Only genomic bins above this threshold will be used. It is recommended to choose a value larger than 1/bin_size.')
    parser.add_argument("--n_train_bins", type=int, default=300000, help='Number of genomic bins to be used for training')
    parser.add_argument("--chunk_size", type=int, default=100000, help='Chunk size for processing the signal matrix. Should be smaller than `n_train_bins`')
    parser.add_argument("--seed", type=int, default=42, help='Random state for reproducability')

    parser.add_argument("--alpha_W", type=float, default=0.01, help='Regularisation for the signal matrix')
    parser.add_argument("--alpha_H", type=float, default=0.001, help='Regularisation for the mixing matrix')
    
    args = parser.parse_args()
    
    logger.info('Checking for output directory...')
    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    
    tiled_files = run(args.input_csv, args.bin_size, args.num_jobs, args.out_dir)
    
    logger.info("Running decoden...")
    decoden_args = Namespace(
        data_folder=os.path.join(args.out_dir),
        output_folder=args.out_dir,
        files_reference=os.path.join(args.out_dir, 'experiment_conditions.json'),
        blacklist_file=args.blacklist_file,
        conditions=extract_conditions(os.path.join(args.out_dir, 'experiment_conditions.json')),
        control_cov_threshold=args.control_cov_threshold,
        n_train_bins=args.n_train_bins,
        chunk_size=args.chunk_size,
        seed=args.seed,
        alpha_W=args.alpha_W,
        alpha_H=args.alpha_H
    )

    main(decoden_args)
