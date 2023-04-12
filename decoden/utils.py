import pandas as pd
import os
from tqdm import tqdm
import numpy as np
from pathlib import Path

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
    for fname, c in tqdm(files_ref.items()):
        conditions_counts[c] += 1    
        colname = c+"_"+str(conditions_counts[c])
        df = pd.read_csv(os.path.join(data_folder, fname), sep="\t", names=["seqnames", "start", "end", colname])
        df.set_index(["seqnames", "start", "end"], inplace=True)
        if data is None:
            data = df
        else:
            data = pd.merge(data, df, left_index=True, right_index=True)

    return data, conditions_counts