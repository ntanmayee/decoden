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
import subprocess
import deeptools.countReadsPerBin as crpb
import pysam


def get_fragment_length(list_of_filepaths, output_dir):
    """Estimate fragment length. Internally uses `macs2 predictd`.

    Args:
        list_of_filepaths (list): list of strings to path to file with reads

    Raises:
        Exception: Unable to compute fragment length

    Returns:
        int: estimated fragment length
    """
    logger.info(f'Getting fragment length for {list_of_filepaths}')

    for filepath in list_of_filepaths:    
        assert os.path.exists(filepath), f"File {filepath} not found"

    result = subprocess.run(f'macs2 predictd -i {" ".join(list_of_filepaths)} -g hs -m 5 50 --outdir {output_dir}', capture_output=True, text=True, shell=True)
    try:
        fragment_length = int([s for s in result.stderr.split('\n') if 'tag size is' in s][0].split()[-2])
    except:
        logger.error(f'Unable to compute fragment length for {list_of_filepaths}')
        raise Exception(f'Unable to compute fragment length for {list_of_filepaths}')
    return fragment_length


class Preprocessor(object):
    def __init__(self, input_csv_filepath, bin_size, num_jobs, out_dir, organism_name) -> None:
        self.input_csv_filepath = input_csv_filepath
        self.bin_size = bin_size
        self.num_jobs = num_jobs
        self.out_dir = out_dir
        self.organism_name = organism_name
        self.fragment_lengths = {}

        self.read_csv()
        self.init_chrom_sizes()

    def init_chrom_sizes(self):
        if os.path.exists(os.path.dirname(self.organism_name)):
            try:
                self.chrom_sizes = pd.read_csv(self.organism_name, sep='\t', names=['chr_name', 'length'])
            except:
                raise ValueError(f'Unable to read {self.organism_name}. Tab separated file must contain 2 columns with chromosome name and length.')
        else:
            try:
                url = f'https://hgdownload-test.gi.ucsc.edu/goldenPath/{self.organism_name}/bigZips/{self.organism_name}.chrom.sizes'
                self.chrom_sizes = pd.read_csv(url, sep='\t', names=['chr_name', 'length'])
            except: 
                raise ValueError(f'Unknown assembly {self.organism_name}')
        
        # create a bed file of chrom sizes to pass to deeptools later
        self.chrom_sizes_path = join(self.out_dir, 'chrom_sizes.bed')
        self.chrom_sizes['start'] = 0
        self.chrom_sizes[['chr_name', 'start', 'length']].to_csv(self.chrom_sizes_path, sep='\t', header=False, index=False)

    def read_csv(self):
        """Read in CSV of different samples. CSV should contain `filepath`, `exp_name` and `is_control` columns. Control should be the first condition.

        Args:
            input_csv_filepath (string): path to CSV file

        Returns:
            pandas.DataFrame: DataFrame of samples
        """
        logger.info(f'Reading CSV file {self.input_csv_filepath}')
        assert os.path.exists(self.input_csv_filepath), f"Annotation file {self.input_csv_filepath} not found"    
        input_csv = pd.read_csv(self.input_csv_filepath)
        assert len(input_csv)>0, "Annotation .csv empty"
        
        n_cols = len(input_csv.columns)
        for col in ["filepath", "exp_name", "is_control", "replicate", "cell_type"]:
            assert col in input_csv.columns, f"Annotation .csv requires `{col}` column"
        if n_cols==6:
            assert "sample_label" in input_csv.columns, "Optional labelling column must be `sample_label`"
            
        assert len(input_csv[input_csv["is_control"]==1]["exp_name"].unique())>0, "Specify at least one control condition"
        assert len(input_csv[input_csv["is_control"]==1]["exp_name"].unique())<2, "Multiple control labels not allowed"
        
        self.input_csv = input_csv

    def preprocess_single(self, condition, group):
        # estimate fragment length for condition
        list_of_filepaths = list(group['filepath'].unique())
        fragment_length = get_fragment_length(list_of_filepaths, self.out_dir)
        self.fragment_lengths[condition] = fragment_length

        assert len(group.is_control.unique()) == 1, 'BAM files from same condition cannot have mixed `is_control` column'
        is_control = True if 1 in group.is_control.unique() else False

        # count reads
        readcount_object = crpb.CountReadsPerBin(list_of_filepaths, binLength=self.bin_size, bedFile=self.chrom_sizes_path, 
                                                 stepSize=self.bin_size, bed_and_bin=True, ignoreDuplicates=True, 
                                                 extend_reads=fragment_length, verbose=True, center_read=is_control)
        processed_reads = readcount_object.run()
        return processed_reads

    def run(self):
        """Run DecoDen for all samples 

        Args:
            input_csv (string): path to CSV with details about experiment conditions and files. CSV should contain `filepath`, `exp_name` and `is_control` columns. Control should be the first condition.
            bin_size (int): width of bin for tiling. Recommended to choose a bin width from 10 - 200 bp. Smaller bin width increases run time.
            num_jobs (int): Number of parallel jobs
            out_dir (string): path to output directory 

        Returns:
            list: list of tuples (tiled_filepath, name). `tiled_filepath` is the path to the processed file.
        """
        for condition, group in self.input_data.groupby('exp_name'):
            processed_reads = self.preprocess_single(condition, group)
            processed_reads.save(join(f'{condition}_reads.npy'))

def run_preprocessing(input_csv, bin_size, num_jobs, out_dir, organism_name):
    preprocess_object = Preprocessor(input_csv, bin_size, num_jobs, out_dir, organism_name)
    preprocess_object.run()
