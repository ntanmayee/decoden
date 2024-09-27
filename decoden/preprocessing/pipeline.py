"""This module runs the preprocessing pipeline. 
"""

import argparse
import pandas as pd
from pathlib import Path
import json
import os
from os.path import isabs, join
from decoden.preprocessing.logger import logger
from decoden.utils import print_message
import subprocess
import deeptools.countReadsPerBin as crpb
from deeptools import bamHandler
from deeptools.utilities import getCommonChrNames
import pysam
import numpy as np
from tqdm import tqdm
from functools import reduce


def get_fragment_length(list_of_filepaths, output_dir, genome_size):
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

    result = subprocess.run(f'macs2 predictd -i {" ".join(list_of_filepaths)} -g {genome_size} -m 5 50 --outdir {output_dir}', capture_output=True, text=True, shell=True)
    try:
        fragment_length = int([s for s in result.stderr.split('\n') if 'tag size is' in s][0].split()[-2])
    except:
        logger.error(f'Unable to compute fragment length for {list_of_filepaths}')
        raise Exception(f'Unable to compute fragment length for {list_of_filepaths}')
    return fragment_length


class Preprocessor(object):
    def __init__(self, input_csv_filepath, bin_size, num_jobs, out_dir, genome_size) -> None:
        self.input_csv_filepath = input_csv_filepath
        self.bin_size = bin_size
        self.num_jobs = num_jobs
        self.out_dir = out_dir
        self.genome_size = genome_size
        self.fragment_lengths = {}
        self.experiment_conditions = {}

        Path(self.out_dir).mkdir(parents=True, exist_ok=True)
        self.read_csv()

    def get_genome_size(self):
        if self.genome_size == 'hs':
            return 3137161264
        if self.genome_size == 'mm':
            return 2725765481
        try:
            genome_size = int(self.genome_size)
            assert genome_size > 0
            return genome_size
        except:
            raise NotImplementedError(f'{self.genome_size} is unknown. Genome size should be an integer or `hs` (Homo sapien) or `mm` (Mus musculus). Please raise an issue if this is a mistake.')

    def init_chrom_sizes(self):
        logger.info('Initialising chrom.sizes...')
   
        bam_handles = []
        logger.info('Indexing bam files..')
        for filename in tqdm(self.input_csv['filepath']):
            pysam.index(filename)
            bam = bamHandler.openBam(filename)
            bam_handles.append(bam)
        
        chrom_names_and_size, non_common = getCommonChrNames(bam_handles, verbose=False)
        assert len(chrom_names_and_size) > 0, "No common chromosomes found. Check `bam` files"
        self.chrom_sizes = pd.DataFrame(chrom_names_and_size, columns=['chr_name', 'length'])
        self.chrom_sizes['start'] = 1

        self.chrom_sizes_path = join(self.out_dir, 'chrom_sizes.bed')
        self.chrom_sizes[['chr_name', 'start', 'length']].to_csv(self.chrom_sizes_path, sep='\t', header=False, index=False)

        logger.info('Completed initialisation of chrom.sizes')

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

        input_dirname = os.path.dirname(self.input_csv_filepath)
        input_csv['filepath'] = input_csv['filepath'].apply(lambda name: name if isabs(name) else join(input_dirname, name))
        for filename in list(input_csv['filepath']):
            assert os.path.exists(filename), f"File {filename} not found"
        
        self.input_csv = input_csv

    def count_reads(self, list_of_filepaths, fragment_length, is_control):
        logger.info('Starting to count reads..')
        readcount_object = crpb.CountReadsPerBin(list_of_filepaths, binLength=self.bin_size, bedFile=self.chrom_sizes_path, 
                                                 stepSize=self.bin_size, bed_and_bin=True, ignoreDuplicates=is_control, 
                                                 extendReads=fragment_length, verbose=True, center_read=is_control,
                                                 numberOfProcessors=self.num_jobs
                                                 )
        res = []
        for i in range(len(self.chrom_sizes)):
            row = self.chrom_sizes.iloc[i]
            chr_, start, end = row.chr_name, row.start, row.length
            res.append(readcount_object.count_reads_in_region(chr_, start, end))
        processed_reads = np.concatenate([r[0] for r in res])

        logger.info('Read coverage computed!')
        return processed_reads
    
    def normalise_library_size(self, list_of_filepaths, processed_reads, target_library_size=2.5e7):
        # normalise to fixed library size, choose conservative estimate to prevent
        # amplification of noise. 25 million is default
        for i, bam_file in list_of_filepaths:
            library_size = reduce(lambda x, y: x + y, [
    int(line.split('\t')[-2]) if len(line) > 0 else 0 for line in pysam.idxstats(bam_file).split('\n')])
            processed_reads[i, :] *= (target_library_size/library_size)
        return processed_reads

    def preprocess_single(self, condition, group):
        # estimate fragment length for condition
        logger.info(f'Preprocessing files for {condition}')
        list_of_filepaths = list(group['filepath'].unique())

        # if control, extend reads in both directions
        assert len(group.is_control.unique()) == 1, 'BAM files from same condition cannot have mixed `is_control` column'
        is_control = True if 1 in group.is_control.unique() else False
        
        # fragment length
        fragment_length = None
        if not is_control:
            fragment_length = get_fragment_length(list_of_filepaths, self.out_dir, self.genome_size)
            self.fragment_lengths[condition] = fragment_length
        else:
            fragment_length = np.median(
                [self.fragment_lengths[a] for a in self.fragment_lengths]
                )
            self.fragment_lengths[condition] = fragment_length

        # count reads and aggregate
        processed_reads = self.count_reads(list_of_filepaths, fragment_length, is_control)
            
        # normalise to fixed library size
        processed_reads = self.normalise_library_size(list_of_filepaths, processed_reads)

        # save results
        save_path = join(self.out_dir, 'data', f'{condition}_reads.npy')
        logger.info(f'Saving results to {save_path}')
        np.save(save_path, processed_reads)

        # update experiment_conditions
        if "sample_label" in group.columns:
            sample_names = list(group.apply(lambda row: row['sample_label'] + "_rep" + str(row['replicate']), axis=1))
        else:
            sample_names = list(group.apply(lambda row: row['exp_name'] + "_"+ row['cell_type'] + "_rep" + str(row['replicate']), axis=1))
        self.experiment_conditions[save_path] = {
            'condition': condition,
            'sample_names': sample_names,
            'filenames': list(group['filepath']),
            'bin_size': self.bin_size,
            'fragment_length': self.fragment_lengths[condition]
            }
        
    def write_experiment_conditions(self):
        # reorder experiment conditions because control condition has to be the first
        file_names = [*self.experiment_conditions]
        file_names = [file_names[-1]] + file_names[:-1]
        experiment_conditions = {k: self.experiment_conditions[k] for k in file_names}
        json.dump(experiment_conditions, open(join(self.out_dir, 'experiment_conditions.json'), 'w'))
        logger.info('Experiment conditions written to file')
 
    def run(self):
        self.init_chrom_sizes()

        Path(self.out_dir, 'data').mkdir(parents=True, exist_ok=True)

        grouped = self.input_csv.groupby('exp_name')
        control_group_name = None
        for condition, group in grouped:
            if 1 not in group['is_control'].unique():
                self.preprocess_single(condition, group)
            else:
                control_group_name = condition
        control = grouped.get_group(control_group_name)
        self.preprocess_single(control_group_name, control)

        logger.info('Successfully completed computing read coverage for all conditions')
        self.write_experiment_conditions()

    def check_preprocessed(self):
        # check if .npy files are present in save directory
        for condition in self.input_csv['exp_name'].unique():
            save_path = f'{condition}_reads.npy'
            try:
                if save_path not in os.listdir(join(self.out_dir, 'data')):
                    return False
            except FileNotFoundError:
                return False
        return True

def run_preprocessing(input_csv, bin_size, num_jobs, out_dir, genome_size):
    """Run DecoDen for all samples 

        Args:
            input_csv (string): path to CSV with details about experiment conditions and files. CSV should contain `filepath`, `exp_name` and `is_control` columns. Control should be the first condition.
            bin_size (int): width of bin for tiling. Recommended to choose a bin width from 10 - 200 bp. Smaller bin width increases run time.
            num_jobs (int): Number of parallel jobs
            out_dir (string): path to output directory 

        Returns:
            list: list of tuples (tiled_filepath, name). `tiled_filepath` is the path to the processed file.
    """
    preprocess_object = Preprocessor(input_csv, bin_size, num_jobs, out_dir, genome_size)
    if not preprocess_object.check_preprocessed():
        preprocess_object.run()
    else:
        logger.info('Exisiting preprocessed files found. Proceeding with existing files.')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', "--input_csv", required=True, help='path to CSV file with information about experimental conditions. Must contain `filepath`, `exp_name` and `is_control` columns. Control/input should be the first condition. Input files can be in BED/BAM format.')
    parser.add_argument('-bs', "--bin_size", default=200, type=int, help='size of genomic bin for tiling. Recommended value is 10-200. Smaller bin size increases space and runtime, larger binsizes may occlude small variations. Default: 200')
    parser.add_argument('-n', "--num_jobs", default=1, type=int, help='Number of parallel jobs for preprocessing. Default: 1')
    parser.add_argument('-o', "--out_dir", required=True, help='path to directory where all output files will be written')
    parser.add_argument('-a', "--assembly", required=True, help='genome assembly')
    
    logger.info('Parsing arguments...')
    args = parser.parse_args()
    
    _ = run_preprocessing(args.input_csv, args.bin_size, args.num_jobs, args.out_dir, args.assembly)
