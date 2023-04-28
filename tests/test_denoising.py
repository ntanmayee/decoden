import pytest
from os.path import join, exists
from decoden.main import *
from decoden.decoden_pipeline import _decoden_pipeline
from decoden.constants import *


@pytest.fixture(autouse=True)
def prepare_preprocessed_data(tmp_session_directory, correct_csv):
    bin_size = 200
    num_jobs = 1
    out_dir = join(tmp_session_directory, "denoise")
    
    preprocess(correct_csv, bin_size, num_jobs, out_dir)
    
    
def test_nmf_results_saved(tmp_session_directory, bl_file):    
    
    out_dir = join(tmp_session_directory, "denoise")
    control_label = "control"
    files_reference = join(out_dir, "experiment_conditions.json")
    alpha_W = 0.01
    alpha_H = 0.001
    control_cov_threshold = 0.5
    n_train_bins = 50000
    chunk_size = 50000
    seed = 0
    plotting = False
    
     
    _decoden_pipeline(["nmf"], 
                      files_reference=files_reference, 
                      control_label=control_label, 
                      out_dir=out_dir, 
                      blacklist_file=bl_file,
                      alpha_W=alpha_W,
                      alpha_H=alpha_H,
                      control_cov_threshold=control_cov_threshold,
                      n_train_bins=n_train_bins,
                      chunk_size=chunk_size,
                      seed=seed,
                      plotting=plotting                      
                      )
    
    
    nmf_folder = join(out_dir, "NMF")
    assert exists(nmf_folder)
    assert exists(join(nmf_folder, "mixing_matrix.csv"))
    assert exists(join(nmf_folder, "signal_matrix.ftr"))
    
    



def test_hsr_consolidated_results_saved(tmp_session_directory, bl_file):    
    
    out_dir = join(tmp_session_directory, "denoise")
    control_label = "control"
    files_reference = join(out_dir, "experiment_conditions.json")
    alpha_W = 0.01
    alpha_H = 0.001
    control_cov_threshold = 0.5
    n_train_bins = 50000
    chunk_size = 50000
    seed = 0
    plotting = False
    
    
    denoise_consolidate(
        files_reference=files_reference,
        control_label=control_label,
        out_dir=out_dir, 
        blacklist_file=bl_file, 
        alpha_W=alpha_W, 
        alpha_H=alpha_H, 
        control_cov_threshold=control_cov_threshold, 
        n_train_bins=n_train_bins, 
        chunk_size=chunk_size, 
        seed=seed, 
        plotting=plotting
    )

    bdg_folder = join(out_dir, BEDGRAPH_FOLDER)
    assert exists(bdg_folder)
    assert exists(join(out_dir, "HSR_results_consolidated.ftr"))

    
def test_nmf_hsr_replicates_results_saved(tmp_session_directory, bl_file):    
    
    out_dir = join(tmp_session_directory, "denoise")
    control_label = "control"
    files_reference = join(out_dir, "experiment_conditions.json")
    alpha_W = 0.01
    alpha_H = 0.001
    control_cov_threshold = 0.5
    n_train_bins = 50000
    chunk_size = 50000
    seed = 0
    plotting = False
    
    denoise_replicates(
        files_reference=files_reference,
        control_label=control_label,
        out_dir=out_dir, 
        blacklist_file=bl_file, 
        alpha_W=alpha_W, 
        alpha_H=alpha_H, 
        control_cov_threshold=control_cov_threshold, 
        n_train_bins=n_train_bins, 
        chunk_size=chunk_size, 
        seed=seed, 
        plotting=plotting
    )
        

    bdg_folder = join(out_dir, BEDGRAPH_FOLDER)
    assert exists(bdg_folder)
    assert exists(join(out_dir, "HSR_results_replicates.ftr"))