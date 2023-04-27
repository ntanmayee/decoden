from scipy import stats
import pandas as pd
from tqdm import tqdm
import os
import numpy as np
import json

def detect_ttest(hsr_rep_output, conditions, threshold = 0.01):
    cols = []
    for cond in conditions:
        cols.extend([f"{cond} Statistic", f"{cond} Pvalue", f"{cond} Peak"])
    out_df = pd.DataFrame(None, columns=cols, index=hsr_rep_output.index)

    condition_columns = {cond: [c for c in hsr_rep_output if cond in c and c.endswith("HSR Value")] for cond in conditions}


    for cond in conditions:
        vals = hsr_rep_output[condition_columns[cond]]
        test = stats.ttest_1samp(vals, 0.0, axis=1, alternative="greater")
        
        out_df.loc[:, [f"{cond} {colname}" for colname in ["Statistic", "Pvalue", "Peak"]]] = np.vstack([test.statistic, test.pvalue, test.pvalue<threshold]).T
    return out_df


def run_peak_calling(files_reference, out_dir):
    bedgraph_dir = os.path.join(out_dir, "output_bedgraph_files")
    peaks_output_dir = os.path.join(out_dir, "called_peaks")
    os.makedirs(peaks_output_dir, exist_ok=True)
    
    with open(files_reference, "r") as f:
        files_mapping = json.load(f)
        
    label_mapping = {}
    
    for v in files_mapping.values():
        condition, replicate, label = v
        if condition not in label_mapping:
            label_mapping[condition] = []
        label_mapping[condition].append(f"{label}_DecoDen.bdg")
        
    # Load the bdg
    
    
    # Run the hypothesis testing
        
    