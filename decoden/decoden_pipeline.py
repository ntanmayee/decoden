"""
This module contains a unified function to run the DecoDen pipeline, to centralize all the different CLI calls.
The main control parameter is a list of the steps to perform in the pipeline. 
The possible steps are:
- preprocess
- nmf
- hsr_consolidate (HSR), or alternatively
- hsr_replicates (Replicate-specific HSR)
- detect (only for replicate-specific hsr)
"""

from os.path import join
from pathlib import Path

from decoden.utils import extract_conditions, save_hsr_output
from decoden.preprocessing.pipeline import run_preprocessing
from decoden.denoising.nmf import run_NMF
from decoden.denoising.hsr import run_HSR, run_HSR_replicates

def _decoden_pipeline(pipeline_steps,
                      
                      # Preprocessing arguments
                      input_csv=None,
                      out_dir=None,
                      bin_size=None, 
                      num_jobs=None,
                      
                      # NMF + HSR arguments
                      files_reference=None,
                      control_label=None,
                      blacklist_file=None,
                      alpha_W=None,
                      alpha_H=None,
                      control_cov_threshold=None,
                      n_train_bins=None,
                      chunk_size=None,
                      seed=None,
                      plotting=None
                      ):
    
    if "preprocess" in pipeline_steps:
        assert input_csv is not None, "Required input csv"
        assert out_dir is not None, "Required output directory"
        assert bin_size>0, "Invalid bin size"
        assert num_jobs>0, "Invalid number of jobs"
    
        Path(out_dir).mkdir(parents=True, exist_ok=True)

        tiled_files = run_preprocessing(input_csv, bin_size, num_jobs, out_dir)
        files_reference = join(out_dir, 'experiment_conditions.json')
        
        
    if "nmf" in pipeline_steps:        
        conditions = extract_conditions(files_reference, control_label=control_label)
        wmatrix, mmatrix, data_noBL, mask, conditions_counts = run_NMF(files_reference, 
                                                                        conditions, 
                                                                        out_dir, 
                                                                        blacklist_file=blacklist_file,
                                                                        alpha_W=alpha_W, 
                                                                        alpha_H=alpha_H, 
                                                                        control_cov_threshold=control_cov_threshold, 
                                                                        n_train_bins=n_train_bins, 
                                                                        chunk_size=chunk_size, 
                                                                        seed=seed,
                                                                        plotting=plotting)
        
    if "hsr_consolidate" in pipeline_steps:
        assert "nmf" in pipeline_steps, "HSR requires NMF as a starting step"
        
        hsr_df = run_HSR(wmatrix, mask, conditions)
        save_hsr_output(hsr_df, out_dir, label="_consolidate")
    
    if "hsr_replicates" in pipeline_steps:
        assert "nmf" in pipeline_steps, "HSR requires NMF as a starting step"
        
        hsr_df = run_HSR_replicates(data_noBL, wmatrix, mmatrix, mask, conditions, conditions_counts)
        save_hsr_output(hsr_df, out_dir, label="_replicates")
        
        
    if "detect" in pipeline_steps:
        pass
    
    