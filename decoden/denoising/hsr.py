import numpy as np
from sklearn.linear_model import LinearRegression
from tqdm import tqdm
import warnings



def run_HSR(wmat, bl_mask, conditions_list, eps=1e-20):
    """Run HSR step of DecoDen

    Args:
        wmat (np.array): signal matrix from NMF step
        bl_mask (np.array): mask of genomic bins that belong to blacklist regions
        conditions_list (list): list of experimental conditions
        eps (_type_, optional): minimum value threshold. Defaults to 1e-20.
    """
    control_condition = conditions_list[0]
    out_df = wmat.loc[:, []]

    control_transf = wmat.loc[:, control_condition].apply(
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
            control_transf.values.reshape(-1, 1)), np.log(0.5))
        pred = np.exp(log_pred)
#         track = np.exp(treatment_transf+mean_treatment_transf-log_pred)
        track = np.exp(treatment_transf-log_pred)
        out_df[treatment_cond+" HSR Value"] = track
        out_df[treatment_cond+" fit"] = pred

    return out_df


def run_HSR_replicates(replicates, wmat, mmat, bl_mask, conditions_list, conditions_counts_ref, eps=1e-20):
    """Run HSR step of DecoDen to adjust individual replicates

    Args:
        replicates (DataFrame): the original samples before the NMF step
        wmat (np.array): signal matrix from NMF step
        mmat (np.array): mixing matrix from NMF step
        bl_mask (np.array): mask of genomic bins that belong to blacklist regions
        conditions_list (list): list of experimental conditions
        conditions_counts_ref (dict): the reference of how many samples are given for each condition
        eps (_type_, optional): minimum value threshold. Defaults to 1e-20.
    """

    control_condition = conditions_list[0]
    out_df = wmat.loc[:, []]
    n_control_samples = conditions_counts_ref[control_condition]


    control_transf = wmat.loc[:, control_condition].apply(
        lambda x: np.maximum(eps, x)).apply(np.log)

    #     control_transf -= np.mean(control_transf)

    treatment_columns = [c for c in replicates.columns if not c.startswith(control_condition)]
    # samples_conditions = [c.split("_")[0] for c in treatment_columns]
    tissue_signal = wmat.loc[:, control_condition]
    for i, treatment_sample in tqdm(enumerate(treatment_columns)):
        # Select only values above the median for the fit, to reduce the contribution of noise
        sample_condition = treatment_sample.split("_")[0]
        sample_data = replicates.loc[:, treatment_sample]

        tissue_coef = mmat.iloc[0, i+n_control_samples]

        sample_specific_signal = (sample_data - tissue_coef*tissue_signal).clip(eps, None)
        treatment_transf = sample_specific_signal.apply(np.log)

        # mean_treatment_transf = np.mean(treatment_transf)
        # treatment_transf -= mean_treatment_transf
        fit_ixs = np.where((control_transf[bl_mask] > np.median(control_transf[bl_mask])) & (
            treatment_transf[bl_mask] > np.median(treatment_transf[bl_mask])))[0]
        reg = LinearRegression(fit_intercept=False).fit(
            control_transf[bl_mask].values[fit_ixs].reshape(-1, 1), treatment_transf[bl_mask][fit_ixs])

        log_pred = np.maximum(reg.predict(
            control_transf.values.reshape(-1, 1)), np.log(0.5))
        pred = np.exp(log_pred)
    #         track = np.exp(treatment_transf+mean_treatment_transf-log_pred)
        track = np.exp(treatment_transf-log_pred)
        out_df[treatment_sample+" HSR Value"] = track
        out_df[treatment_sample+" fit"] = pred

    
    return out_df
