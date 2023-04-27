from scipy import stats
import pandas as pd
from tqdm import tqdm

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