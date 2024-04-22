import numpy as np
from sklearn.linear_model import LinearRegression
from tqdm import tqdm
import warnings
from decoden.constants import UNSPECIFIC_SIGNAL_LABEL


def run_HSR(wmat, bl_mask, conditions_list, eps=1e-20, signal_clip=0.1):
    """Run HSR step of DecoDen

    Args:
        wmat (np.array): signal matrix from NMF step
        bl_mask (np.array): mask of genomic bins that belong to blacklist regions
        conditions_list (list): list of experimental conditions
        eps (_type_, optional): minimum value threshold. Defaults to 1e-20.
    """
    # control_label = conditions_list[0]
    out_df = wmat.loc[:, []]

    control_transf = wmat.loc[:, UNSPECIFIC_SIGNAL_LABEL].apply(
        lambda x: np.maximum(eps, x)).apply(np.log)

#     control_transf -= np.mean(control_transf)

    for i, treatment_cond in tqdm(enumerate(conditions_list[1:])):
        # Select only values above the median for the fit, to reduce the contribution of noise
        treatment_transf = wmat.loc[:, treatment_cond].apply(
            lambda x: np.maximum(eps, x)).apply(np.log)
        mean_treatment_transf = np.mean(treatment_transf)
#         treatment_transf -= mean_treatment_transf

        # print warning message if median is too low
        if (np.median(treatment_transf[bl_mask]) < 1) or (np.median(treatment_transf[bl_mask]) < 1):
            warnings.warn("Treatment/Control coverage might be too low for half-sibling regression to be effective!")

        fit_ixs = np.where((control_transf[bl_mask] > np.median(control_transf[bl_mask])) & (
            treatment_transf[bl_mask] > np.median(treatment_transf[bl_mask])))[0]
        reg = LinearRegression(fit_intercept=False).fit(
            control_transf[bl_mask].values[fit_ixs].reshape(-1, 1), treatment_transf[bl_mask][fit_ixs])

        log_pred = np.maximum(reg.predict(
            control_transf.values.reshape(-1, 1)), np.log(signal_clip))
        pred = np.exp(log_pred)
#         track = np.exp(treatment_transf+mean_treatment_transf-log_pred)
        track = np.exp(treatment_transf-log_pred)
        out_df[treatment_cond+" HSR Value"] = track

    return out_df
