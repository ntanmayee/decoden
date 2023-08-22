from scipy import stats
import pandas as pd
from tqdm import tqdm
import os
from os.path import join
import numpy as np
import json
from itertools import groupby

from statsmodels.stats.multitest import multipletests
from decoden.constants import *


def load_hsr_results(condition_files, out_dir):
    df = None
    for fname in condition_files:
        label = condition_files[0].rsplit("_", 1)[0]
        sample_df = pd.read_csv(join(out_dir, BEDGRAPH_FOLDER, fname), skiprows=1, sep="\t", names=["seqnames", "start", "end", label])
        sample_df.set_index(["seqnames", "start", "end"], inplace=True)

        if df is None:
            df = sample_df
        else:
            df = pd.merge(df, sample_df, left_index=True, right_index=True)
    return df

# def detect_condition_peaks_ttest(condition_df, condition, threshold = 0.01):
#     cols = [f"{condition} Peak", f"{condition} Statistic", f"{condition} Pvalue", ]
        
#     test = stats.ttest_1samp(condition_df.values, 0.0, axis=1, alternative="greater")
    
#     # TODO: improve selection criterion
#     test_results = np.vstack([test.pvalue<threshold, test.statistic, test.pvalue]).T
#     out_df = pd.DataFrame(test_results, columns=cols, index=condition_df.index)
        
    return out_df

# def compress_bed_peak_intervals(df, condition, min_pval = 1e-10):
#     if "seqnames" not in df.columns:
#         df = df.reset_index()
#     # ix_peak_col, peak_col = [(i, c) for i, c in enumerate(df.columns) if c.endswith(" Peak")][0]
#     ix_pval_col, pval_col = [(i, c) for i, c in enumerate(df.columns) if "Pvalue" in c][0]
#     vals = df[["seqnames", peak_col]].values.tolist()
#     count_dups = [[sum(1 for v in group), v] for v, group in groupby(vals)]

#     idx = 0
#     compressed_df = []

#     peak_counter = 1
#     for interval, (chrom, val) in tqdm(count_dups):
#         # aggregate_pval =  df.iloc[idx:idx+interval, ix_pval_col].min()
#         if val>0.0:
#             aggregate_pval =  np.min(df.iloc[idx:idx+interval, ix_pval_col].values)
#             score = -np.log10(np.max([aggregate_pval, min_pval]))
#             label = f"{condition}_n{peak_counter}_w{interval}"
#             row = {"seqnames": chrom, "start": df.iat[idx, 1], "end": df.iat[idx+interval-1, 2],
#                    "name": label, "score": score}#, peak_col: val, "width": interval, pval_col:aggregate_pval}
#             compressed_df.append(row)
#             peak_counter += 1 
#         idx += interval

#     compressed_df = pd.DataFrame(compressed_df)
#     return compressed_df


def compact_df(df):
    compacted_df = []
    
    peak_val = None
    chromosome = None
    lowest_pval = 1.0
    start = None
    end = None
    prev_row = None
    
    peak_colname = [c for c in df.columns if "Peak" in c][0]
    pval_colname = [c for c in df.columns if "Pvalue" in c][0]
    for i, row in df.iterrows():
        
        if peak_val==row[peak_colname] and i[0]==chromosome:
            end = i[2]
            if row[pval_colname]<lowest_pval:
                lowest_pval = row[pval_colname]
        else:
            if chromosome is not None:
                compacted_df.append({"seqnames": chromosome, "start": start, "end": end,
                                    peak_colname: peak_val, pval_colname: lowest_pval})
            peak_val = int(row[peak_colname])
            chromosome, start, end = i
            lowest_pval = row[pval_colname]
    
    compacted_df.append({"seqnames": chromosome, "start": start, "end": end,
                                    peak_colname: peak_val, pval_colname: lowest_pval})
    compacted_df = pd.DataFrame(compacted_df)
    # compacted_df = compacted_df.set_index(["seqnames", "start", "end"])
    return compacted_df


def detect_condition_peaks_ttest(condition_df, condition, peak_threshold = 0.01, alpha=0.05, min_width=150, eps=1e-3, min_pval=1e-8):
    cols = [f"{condition} Peak", f"{condition} Statistic", f"{condition} Pvalue"]

    noise_adj = np.tile(np.arange(condition_df.values.shape[1]), (len(condition_df), 1))
    test = stats.ttest_1samp(condition_df.values+noise_adj*eps, 0.0, axis=1, alternative="greater")

    # TODO: improve selection criterion
    test_results = np.vstack([(test.pvalue<peak_threshold).astype(int), test.statistic, np.maximum(test.pvalue, 1e-8)]).T
    out_df = pd.DataFrame(test_results, columns=cols, index=condition_df.index)
    peak_colname = [c for c in out_df.columns if "Peak" in c][0]
    out_df[peak_colname] = out_df[peak_colname].astype(int)
    pval_colname = [c for c in out_df.columns if "Pvalue" in c][0]

    out_df = compact_df(out_df)
    # TODO CHECK BEHAVIOUR IF YOU HAVE GAPS
    # Merge Peaks
    width = out_df["end"] - out_df["start"]
    out_df.loc[(out_df[peak_colname]==0)&(width<min_width), peak_colname] = 1
    out_df = compact_df(out_df.set_index(["seqnames", "start", "end"]))
    # return out_df

    # Remove short peaks
    width = out_df["end"] - out_df["start"]
    out_df.loc[(out_df[peak_colname]==1)&(width<min_width), peak_colname] = 0
    out_df = compact_df(out_df.set_index(["seqnames", "start", "end"]))

    
    corrected_pvalues = multipletests(out_df[pval_colname], alpha=alpha, method="fdr_bh")
    out_df[peak_colname+" Corrected"] = corrected_pvalues[0].astype(int)
    out_df[pval_colname+" Corrected"] = corrected_pvalues[1]
    if out_df[peak_colname+" Corrected"].sum()<1:
        print(f"\n\nNo peaks found for condition {condition} after BH correction\n\n")
    
    out_df = out_df[out_df[peak_colname+" Corrected"]==1]

    peak_counter = 1
    results = []
    for i, row in out_df.iterrows():
        interval = row.end-row.start
        score = -np.log10(np.max([row[pval_colname], min_pval]))
        label = f"{condition}_n{peak_counter}_w{interval}"
        peak_counter += 1
        results.append({"seqnames": row.seqnames, "start": row.start, "end": row.end,
                       "name": label, "score": score})
    results = pd.DataFrame(results)
    return results



def run_peak_calling(files_reference, out_dir, control_label, pval_alpha=None,
                      peak_threshold=None,
                      min_width=None):
    bedgraph_dir = os.path.join(out_dir, BEDGRAPH_FOLDER)
    
    peaks_output_dir = os.path.join(out_dir, CALLED_PEAKS_FOLDER)
    os.makedirs(peaks_output_dir, exist_ok=True)
    
    with open(files_reference, "r") as f:
        files_mapping = json.load(f)
        
    label_mapping = {}
    
    for v in files_mapping.values():
        condition, replicate, label = v
        if condition==control_label:
            continue
        if condition not in label_mapping:
            label_mapping[condition] = []
        label_mapping[condition].append(f"{label}_{condition}_DecoDen.bdg")
        
    for cond, files in label_mapping.items():
        assert len(files)>1, "Hypothesis testing requires multiple replicates"
        cond_df = load_hsr_results(files, out_dir=out_dir)
        pc_results = detect_condition_peaks_ttest(cond_df, cond, alpha=pval_alpha, peak_threshold=peak_threshold, min_width=min_width)
        pc_results.to_csv(join(peaks_output_dir, f"{cond}_peaks.bed"), sep="\t", header=False, index=False)
    
        
