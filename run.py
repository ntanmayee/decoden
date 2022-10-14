import argparse
from argparse import Namespace
from pathlib import Path
import json
import os
from preprocess.logger import logger
from run_preprocess import run
from run_decoden import main

def print_message():
    with open('utils/message.txt') as fp:
        print(fp.read())

def extract_conditions(json_file):
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
    parser.add_argument('-i', "--input_csv", required=True)
    parser.add_argument('-bs', "--bin_size", default=200, type=int)
    parser.add_argument('-n', "--num_jobs", default=1, type=int)
    parser.add_argument('-o', "--out_dir", required=True)

    # arguments from decoden
    parser.add_argument('-bl', "--blacklist_file")
    parser.add_argument("--control_cov_threshold", type=float, default=1.0)
    parser.add_argument("--n_train_bins", type=int, default=300000)
    parser.add_argument("--chunk_size", type=int, default=100000)
    parser.add_argument("--seed", type=int, default=42)

    parser.add_argument("--alpha_W", type=float, default=0.01)
    parser.add_argument("--alpha_H", type=float, default=0.001)
    
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
