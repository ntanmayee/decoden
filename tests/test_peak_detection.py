import pytest
from os.path import join, exists
from decoden.main import *
from decoden.constants import *

@pytest.fixture(autouse=True)
def prepare_denoised_data(tmp_session_directory, correct_csv, bl_file):
    bin_size = 200
    num_jobs = 1
    out_dir = join(tmp_session_directory, "detect")
    
    control_label = "control"
    alpha_W = 0.01
    alpha_H = 0.001
    control_cov_threshold = 0.5
    n_train_bins = 50000
    chunk_size = 50000
    seed = 0
    plotting = False
    
    # if not exists(join(out_dir, "HSR_results_replicates.ftr")):
    run_replicates(correct_csv, 
                    bin_size, 
                    num_jobs, 
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
    


def test_peak_detection(tmp_session_directory):
    out_dir = join(tmp_session_directory, "detect")
    files_reference = join(out_dir, "experiment_conditions.json")
    control_label="control"
    conditions = extract_conditions(os.path.join(out_dir, 'experiment_conditions.json'), control_label=control_label)
    detect(files_reference=files_reference, out_dir=out_dir, control_label=control_label)

    assert exists(join(out_dir, CALLED_PEAKS_FOLDER))
    for c in conditions:
        if c==control_label:
            continue
        assert exists(join(out_dir, CALLED_PEAKS_FOLDER, f"{c}_peaks.bed"))