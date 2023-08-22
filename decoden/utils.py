import pandas as pd
import os
from tqdm import tqdm
import numpy as np
from pathlib import Path
import json
import os
from os.path import join
from decoden.preprocessing.logger import logger
import random
from itertools import groupby
from collections import defaultdict
from decoden.constants import *

dtype_mapping = defaultdict(lambda: float)
dtype_mapping["seqnames"] = str
dtype_mapping["start"] = int
dtype_mapping["end"] = int

def print_message():
    """Print cool opening message when DecoDen in run :) 
    """
    dirname = Path(__file__).parent.parent
    with open(os.path.join(dirname, 'utils/message.txt')) as fp:
        print(fp.read())

def get_blacklisted_regions_mask(df, bl_regions):
    """get boolean mask for genmic bins belonging to blacklisted regions

    Args:
        df (pandas.DataFrame): Data frame with data
        bl_regions: blacklisted regions

    Returns:
        np.array: Boolean mask with 0s for bins belonging to blacklist regions
    """
    if len(bl_regions)<1:
        return np.ones(len(df)).astype(bool)
    reg = bl_regions[bl_regions["seqnames"].isin(df.index.get_level_values("seqnames").unique())]
    filter_ixs = None
    index_start = df.index.get_level_values("start")
    index_end = df.index.get_level_values("end")
    for i, r in tqdm(reg.iterrows()):
        ixs = ((index_start>=r["start"]) & (index_start<=r["end"]) | (index_end>=r["start"]) & (index_end<=r["end"]))
        if filter_ixs is None:
            filter_ixs = ixs
        else:
            filter_ixs = filter_ixs | ixs
    return ~filter_ixs


def load_files(files_ref, data_folder, sample_conditions):
    """Load all files into data matrix

    Args:
        files_ref (dict): dict containing filepath and experiment type
        data_folder (string): path to data folder
        sample_conditions (list): list of conditions

    Returns:
        tuple: (data, conditions_counts)
    """

    conditions_counts = {c: 0 for c in sample_conditions}
    
    data = None
    
    
    for fname, (c, rep, label) in tqdm(files_ref.items()):
        conditions_counts[c] += 1    
        colname = c+"_"+str(rep)
        df = pd.read_csv(os.path.join(data_folder, fname), sep="\t", 
                         names=["seqnames", "start", "end", colname], dtype=dtype_mapping)
        df.set_index(["seqnames", "start", "end"], inplace=True)
        if data is None:
            data = df
        else:
            data = pd.merge(data, df, left_index=True, right_index=True)

    return data, conditions_counts


def extract_conditions(json_file, control_label="control"):
    """Extract list of different experimental conditions, to pass to run_decoden.py.

    Args:
        json_file (string): path to `experiment_conditions.json` generated from run_preprocess.py
        control_label (string): the label of the control/input condition 
    Returns:
        list: list of different experimental conditions
    """
    json_object = json.load(open(json_file))
    
    conditions = []
    for k, v in json_object.items():
        if v[0] not in conditions:
            conditions.append(v[0])
            
    if not control_label in conditions:
        raise Exception("Invalid label for control condition")
    conditions = [control_label] + [c for c in conditions if c!=control_label]
    logger.info(conditions)
    return conditions 

def extract_control_condition(input_csv_filepath):
    input_csv = pd.read_csv(input_csv_filepath)
    control_label = input_csv[input_csv["is_control"]==1]["exp_name"].iloc[0]
    return control_label

def compress_bdg_df(df):
    if not "seqnames" in df.columns:
        df = df.reset_index()
    vals = df.iloc[:,[0,3]].values.tolist()
    count_dups = [[sum(1 for v in group), v] for v, group in groupby(vals)]
    idx = 0
    compressed_df = []
    for interval, (chrom, val) in tqdm(count_dups):
        row = {"seqnames": chrom, "start": df.iat[idx, 1], "end": df.iat[idx+interval-1, 2],
               "Value": val}
        compressed_df.append(row)
        idx += interval

    compressed_df = pd.DataFrame(compressed_df)
    return compressed_df

def save_hsr_output(hsr_df, out_dir, replicate_specific=False, files_ref=None):
    print("\nSaving HSR output")
    # label="_replicates" if replicate_specific else "_consolidated"
    if replicate_specific:
        with open(files_ref, "r") as f:
            files_mapping = json.load(f)
            
        label_mapping = {}
        for v in files_mapping.values():
            label_mapping[(v[0], v[1])] = v[2]
    
    # hsr_df.reset_index().to_feather(join(out_dir, f"HSR_results{label}.ftr"))
    
    
    bedgraph_dir = join(out_dir, BEDGRAPH_FOLDER)
    os.makedirs(bedgraph_dir, exist_ok=True)
    cols = [c for c in hsr_df.columns if c.endswith("HSR Value")]
    
    r = lambda: random.randint(0, 255) 
    for c in tqdm(cols):
        condition = c.split(" HSR")[0]
        if replicate_specific:
            condition, replicate = condition.rsplit("_", 1)
            replicate = int(replicate)
        
        random.seed(condition)
        color1 = '{},{},{}'.format(r(),r(),r())
        color2 = '{},{},{}'.format(r(),r(),r())
        
        if replicate_specific:
            label = label_mapping[(condition, replicate)]
            fname = join(bedgraph_dir, f"{label}_{condition}_DecoDen.bdg")
        else:
            fname = join(bedgraph_dir, f"{condition}_DecoDen.bdg")
        with open(fname, 'w') as f:
            f.write(f'track type=bedGraph name="{condition}" description="{condition}" visibility=full color={color1} altColor={color2} priority=20\n')
        bdg_df = hsr_df[[c]].fillna(0.0)
        bdg_df = compress_bdg_df(bdg_df)
        bdg_df.to_csv(fname, header=False, sep="\t", mode='a', index=False)
    
    
    