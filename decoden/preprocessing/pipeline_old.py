"""This module runs the preprocessing pipeline. 
"""

import argparse
import pandas as pd
from joblib import Parallel, delayed
from pathlib import Path
import json
import os
from os.path import isabs, join
from refact.decoden.preprocessing.macs_steps_old import run_pipeline
from decoden.preprocessing.logger import logger
from decoden.utils import print_message

def read_csv(input_csv_filepath):
    """Read in CSV of different samples. CSV should contain `filepath`, `exp_name` and `is_control` columns. Control should be the first condition.

    Args:
        input_csv_filepath (string): path to CSV file

    Returns:
        pandas.DataFrame: DataFrame of samples
    """
    logger.info(f'Reading CSV file {input_csv_filepath}')
    
    
    assert os.path.exists(input_csv_filepath), f"Annotation file {input_csv_filepath} not found"
    
    input_csv = pd.read_csv(input_csv_filepath)

    assert len(input_csv)>0, "Annotation .csv empty"
    
    n_cols = len(input_csv.columns)
    for col in ["filepath", "exp_name", "is_control", "replicate", "cell_type"]:
        assert col in input_csv.columns, f"Annotation .csv requires `{col}` column"
    if n_cols==6:
        assert "sample_label" in input_csv.columns, "Optional labelling column must be `sample_label`"
        
    assert len(input_csv[input_csv["is_control"]==1]["exp_name"].unique())>0, "Specify at least one control condition"
    assert len(input_csv[input_csv["is_control"]==1]["exp_name"].unique())<2, "Multiple control labels not allowed"
    
    return input_csv

def make_args(input_csv, out_dir, bin_size):
    """Make arguments for parallelization. 

    Args:
        input_csv (string, Path): path to CSV with details about experiment conditions and files.
        out_dir (string): path to write processed files
        bin_size (int): bin size for tiling

    Returns:
        list: list of args for each sample. The parameters are in this order - filepath, exp_name, out_dir, is_control, bin_size
    """
    input_data = read_csv(input_csv)
    input_dirname = os.path.dirname(input_csv)

    logger.info(f'Making arguments for parallelization...')
    arg_list = []
    
    for _, row in input_data.iterrows():
        filepath = row.filepath 
        if not isabs(filepath):
            filepath = join(input_dirname, filepath)
        
        if "sample_label" in input_data.columns:
            sample_name = row.sample_label + "_rep" + str(row.replicate)
        else:
            sample_name = row.exp_name + "_"+ row.cell_type + "_rep" + str(row.replicate)
        sample_name = sample_name.replace(" ", "_")
        arg_list.append(
            {
                "filepath": filepath,
                "exp_name": row.exp_name,
                "sample_name": sample_name,
                "out_dir": out_dir,
                "is_control": row.is_control,
                "bin_size": bin_size,
                "replicate": row.replicate
            }
        )

    logger.debug(arg_list)
    return arg_list



def write_json(tiled_files, out_dir):
    """write `experiment_conditions.json` to output directory. This file is required to run DecoDen later.

    Args:
        tiled_files (list): list of tuples from run_single() function. Each tuple should contain `(tiled_filepath, name)`
        out_dir (string): path to output directory
    """
    out_filename = os.path.join(out_dir, 'experiment_conditions.json')
    logger.info(f'Writing json file in {out_filename}')
    json_obj = {os.path.join('data', os.path.basename(a[0])): a[1] for a in tiled_files}
    json.dump(json_obj, open(out_filename, 'w'), indent=1)


def run_single(args):
    """Run preprocessing pipeline for single sample

    Args:
        args (tuple): list of args for each sample. The parameters are in this order - filepath, exp_name, out_dir, is_control, bin_size

    Returns:
        tuple: (tiled_filepath, name) - filepath to the processed file and the unique name of the sample.
    """
    logger.info(f'Running pipeline for {args}')
    input_filepath = args["filepath"]
    exp_name = args["exp_name"]
    sample_name = args["sample_name"]
    out_dir = args["out_dir"]
    is_control = args["is_control"]
    bin_size = args["bin_size"]
    replicate = args["replicate"]
    
    files = os.listdir(os.path.join(out_dir, 'data'))
    tiled_filepath = sample_name + "_" + exp_name + "_filterdup_pileup_tiled.bed"
    
    if tiled_filepath not in files: 
        assert os.path.exists(input_filepath), f"File {input_filepath} not found"
        tiled_filepath = run_pipeline(input_filepath, sample_name + "_" + exp_name, out_dir, is_control, bin_size)
    else:
        logger.info(f'Skipping {input_filepath} as output already exists.')

    return (tiled_filepath, [exp_name, replicate, sample_name])


def run_preprocessing(input_csv, bin_size, num_jobs, out_dir):
    """Run DecoDen for all samples. This is supposed to be parallel, but currently it is not working. 

    Args:
        input_csv (string): path to CSV with details about experiment conditions and files. CSV should contain `filepath`, `exp_name` and `is_control` columns. Control should be the first condition.
        bin_size (int): width of bin for tiling. Recommended to choose a bin width from 10 - 200 bp. Smaller bin  width increases run time.
        num_jobs (int): Number of parallel jobs
        out_dir (string): path to output directory 

    Returns:
        list: list of tuples (tiled_filepath, name). `tiled_filepath` is the path to the processed file.
    """  
    
    Path(out_dir, 'data').mkdir(parents=True, exist_ok=True)
    
    arg_list = make_args(input_csv, out_dir, bin_size)

    
    logger.info(f'Running {num_jobs} in parallel...')
    tiled_files = Parallel(n_jobs=num_jobs)(
        delayed(run_single)(args) for args in arg_list 
    )
    logger.info(f'Parallel jobs completed.')
    
    write_json(tiled_files, out_dir)

    return tiled_files

if __name__ == '__main__':
    print_message()

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', "--input_csv", required=True, help='path to CSV file with information about experimental conditions. Must contain `filepath`, `exp_name` and `is_control` columns. Control/input should be the first condition. Input files can be in BED/BAM format.')
    parser.add_argument('-bs', "--bin_size", default=200, type=int, help='size of genomic bin for tiling. Recommended value is 10-200. Smaller bin size increases space and runtime, larger binsizes may occlude small variations. Default: 200')
    parser.add_argument('-n', "--num_jobs", default=1, type=int, help='Number of parallel jobs for preprocessing. Default: 1')
    parser.add_argument('-o', "--out_dir", required=True, help='path to directory where all output files will be written')
    
    logger.info('Parsing arguments...')
    args = parser.parse_args()
    
    logger.info('Checking for output directory...')
    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    
    _ = run_preprocessing(args.input_csv, args.bin_size, args.num_jobs, args.out_dir)
