from scipy import stats
import pandas as pd
from tqdm import tqdm

def detect_ttest(hsr_rep_output, conditions, threshold = 0.01):
    cols = []
    for cond in conditions:
        cols.extend([f"{cond} Statistic", f"{cond} Pvalue", f"{cond} Peak"])
    out_df = pd.DataFrame(None, columns=cols, index=hsr_rep_output.index)

    condition_columns = {cond: [c for c in hsr_rep_output if cond in c and c.endswith("HSR Value")] for cond in conditions}
    for i, row in tqdm(hsr_rep_output.iterrows()):
        
        for cond in conditions:
            vals = row[condition_columns[cond]].values
            test = stats.ttest_1samp(vals, 0.0, alternative="greater")
            out_df.loc[i, f"{cond} Statistic"] = test.statistic
            out_df.loc[i, f"{cond} Pvalue"] = test.pvalue
            out_df.loc[i, f"{cond} Peak"] = int(test.pvalue < threshold)
    return out_df