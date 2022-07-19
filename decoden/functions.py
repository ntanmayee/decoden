import pandas as pd
import numpy as np
from sklearn.decomposition import NMF, non_negative_factorization
from sklearn.linear_model import LinearRegression
from tqdm import tqdm

def extract_mixing_matrix(data_df, conditions_list, conditions_counts_ref, alpha_W=0.01, alpha_H=0.001,
                        control_cov_threshold=1.0, n_train_bins=300000, seed=42):

    # Filter data to have sufficient control coverage
    dsel = data_df[(data_df[[c for c in data_df.columns if c.startswith(conditions_list[0])]]>control_cov_threshold).any(axis=1)]
    train_data = dsel.sample(n_train_bins, random_state=seed)
    train_data

    # Extract unspecific signal from control samples
    model = NMF(n_components=1, init='random', random_state=0, beta_loss=1, solver="mu", alpha_W=alpha_W, alpha_H=alpha_H)
    control_cols = [c for c in data_df.columns if c.startswith(conditions_list[0])]
    treatment_cols = [c for c in data_df.columns if c not in control_cols]
    W_unspec = model.fit_transform(train_data.loc[:,control_cols])
    H_unspec = model.components_

    # Calculate unspecific signal coefficients for treatment
    # I swapped and transposed the matrices to make use of the update_H parameter
    W_treat_uns, H_treat_uns, n_iter = non_negative_factorization(train_data.loc[:,treatment_cols].T, 
                                                                  n_components=1, init="custom",
                                                                  H=W_unspec.reshape(1,-1), update_H=False, alpha_W=alpha_W,
                                                                  beta_loss=1, solver="mu")
    H_unspec_coefs = np.concatenate((H_unspec.flatten(), W_treat_uns.T.flatten()))

    # Subtract the unspecific contribution from the data
    # Cap the minimum value to 0 to account for predictions higher than the signal
    specific_data_mat = np.maximum(train_data.values - W_unspec.dot(H_unspec_coefs.reshape(1,-1)), 0)
    specific_data_mat.shape

    # For each modification, extract the specific signal components
    treatment_conditions = conditions_list[1:]
    n_replicates = [conditions_counts_ref[c] for c in conditions_list]
    n_control_replicates = n_replicates[0]
    histone_n_replicates = n_replicates[1:]
    signal_matrix = [W_unspec]
    mixing_matrix = np.zeros((len(conditions_list), np.sum(n_replicates)))
    mixing_matrix[0,:] = H_unspec_coefs

    ix = n_control_replicates
    for i, modif in tqdm(enumerate(treatment_conditions)):

        c = histone_n_replicates[i]
        # print(modif, str(ix), str(ix+c))
        W_spec, H_spec, n_iter = non_negative_factorization(specific_data_mat[:, ix:ix+c], 
                                                                  n_components=1, beta_loss=1, solver="mu", alpha_W=alpha_W, alpha_H=alpha_H)

        mixing_matrix[i+1, ix:ix+c] = H_spec
        signal_matrix.append(W_spec)
        ix += c
    mm = pd.DataFrame(mixing_matrix, index=["unspecific"]+treatment_conditions, columns=train_data.columns)
    return mm


def extract_signal(data_df, mmatrix, conditions_list, chunk_size=100000, alpha_W=0.01, seed=42):
    # Use the mixing matrix to extract the signals. Can be done in chunks to fit in memory
    processed_W = []
    for ck in tqdm(np.split(data_df.values, range(chunk_size, len(data_df)+chunk_size, chunk_size))):
        if len(ck)<1:
            # print(ck)
            break
        W = np.random.uniform(0, 0.1, size=(len(ck), len(mmatrix)))
        ck_W, H, n_iter = non_negative_factorization(ck, W, mmatrix.values, 
                                                   n_components=len(mmatrix), init='custom', random_state=seed, update_H=False, 
                                                     beta_loss="kullback-leibler", solver="mu", alpha_W=alpha_W, max_iter=500
                                                         )
        processed_W.append(ck_W)
    processed_W = np.vstack(processed_W)
    processed_W = pd.DataFrame(processed_W, index=data_df.index, columns=conditions_list)
    return processed_W


def run_HSR(wmat, bl_mask, conditions_list, eps=1e-20):
    control_condition = conditions_list[0]
    out_df = wmat.loc[:,[]]
    
    control_transf = wmat.loc[:,control_condition].apply(lambda x: np.maximum(eps, x)).apply(np.log)
    
    for i, treatment_cond in tqdm(enumerate(conditions_list[1:])):
        # Select only values above the median for the fit, to reduce the contribution of noise
        treatment_transf = wmat.loc[:,treatment_cond].apply(lambda x: np.maximum(eps, x)).apply(np.log)
        
        fit_ixs = np.where((control_transf[bl_mask]>np.median(control_transf[bl_mask])) & (treatment_transf[bl_mask]>np.median(treatment_transf[bl_mask])))[0]
        reg = LinearRegression(fit_intercept=False).fit(control_transf[bl_mask].values.reshape(-1,1), treatment_transf[bl_mask])
        
        log_pred = np.maximum(reg.predict(control_transf.values.reshape(-1,1)), np.log(0.5))
        pred = np.exp(log_pred)
        track = np.exp(treatment_transf-log_pred)
        out_df[treatment_cond+" HSR Value"] = track
        out_df[treatment_cond+" fit"] = pred

    return out_df