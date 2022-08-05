import argparse
import pandas as pd
from joblib import Parallel, delayed
from pathlib import Path
import json
import os
from preprocess.macs_steps import run_pipeline
from preprocess.logger import logger

def read_csv(input_csv_filepath):
    logger.info(f'Reading CSV file {input_csv_filepath}')
    input_csv = pd.read_csv(input_csv_filepath)
    
    if 'filepath' not in input_csv.columns:
        logger.error('`filepath` column not found. Check input CSV.')
        raise ValueError("`filepath` column not found. Check input CSV.")
    if 'exp_name' not in input_csv.columns:
        logger.error('`exp_name` column not found. Check input CSV.')
        raise ValueError("`exp_name` column not found. Check input CSV.")
    if 'is_control' not in input_csv.columns:
        logger.error('`is_control` column not found. Check input CSV.')
        raise ValueError("`is_control` column not found. Check input CSV.")
    
    return input_csv

def make_args(input_csv, out_dir, bin_size):
    logger.info(f'Making arguments for parallelization...')
    arg_list = []
    for _, row in input_csv.iterrows():
        arg_list.append(
            (
                row.filepath,
                row.exp_name,
                out_dir,
                row.is_control,
                bin_size
            )
        )
        
    logger.debug(arg_list)
    return arg_list

def run_single(args):
    logger.info(f'Running pipeline for {args}')
    input_filepath, name, out_dir, is_control, bin_size = args
    tiled_filepath = run_pipeline(input_filepath, name, out_dir, is_control, bin_size)

    return (tiled_filepath, name)

def write_json(tiled_files, out_dir):
    out_filename = os.path.join(out_dir, 'experiment_conditions.json')
    logger.info(f'Writing json file in {out_filename}')
    json_obj = {a[0]: a[1] for a in tiled_files}
    json.dump(json_obj, out_filename, indent=1)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', "--input_csv", required=True)
    parser.add_argument('-bs', "--bin_size", default=200, type=int)
    parser.add_argument('-n', "--num_jobs", default=2, type=int)
    parser.add_argument('-o', "--out_dir", required=True)
    
    logger.info('Parsing arguments...')
    args = parser.parse_args()
    
    logger.info('Checking for output directory...')
    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    
    input_csv = read_csv(args.input_csv)
    arg_list = make_args(input_csv, args.out_dir, args.bin_size)
    
    logger.info(f'Running {args.num_jobs} in parallel...')
    tiled_files = Parallel(n_jobs=args.num_jobs)(
        delayed(run_single)(args) for args in arg_list 
    )
    logger.info(f'Parallel jobs completed.')
    
    write_json(tiled_files, args.out_dir)
