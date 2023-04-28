import pytest
from os.path import join, exists
from decoden.main import *
from decoden.decoden_pipeline import _decoden_pipeline
from decoden.constants import *


def test_run_replicates_results_saved(tmp_session_directory, bl_file, correct_csv):
    bin_size = 200
    num_jobs = 1
    out_dir = join(tmp_session_directory, "pipeline")
    control_label = "control"
    alpha_W = 0.01
    alpha_H = 0.001
    control_cov_threshold = 0.5
    n_train_bins = 50000
    chunk_size = 50000
    seed = 0
    plotting = False

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
    assert exists(out_dir)
    assert exists(join(out_dir, "experiment_conditions.json"))
    nmf_folder = join(out_dir, "NMF")
    assert exists(nmf_folder)
    assert exists(join(nmf_folder, "mixing_matrix.csv"))
    assert exists(join(nmf_folder, "signal_matrix.ftr"))

    bdg_folder = join(out_dir, BEDGRAPH_FOLDER)
    assert exists(bdg_folder)

    files_reference = os.path.join(out_dir, 'experiment_conditions.json')
    with open(files_reference, "r") as f:
        files_mapping = json.load(f)
    labels = []
    for v in files_mapping.values():
        if v[0]==control_label:
            continue
        labels.append(v[2])

    for label in labels:
        assert exists(join(bdg_folder, f"{label}_DecoDen.bdg"))


def test_run_consolidated_results_saved(tmp_session_directory, bl_file, correct_csv):
    bin_size = 200
    num_jobs = 1
    out_dir = join(tmp_session_directory, "pipeline")
    control_label = "control"
    alpha_W = 0.01
    alpha_H = 0.001
    control_cov_threshold = 0.5
    n_train_bins = 50000
    chunk_size = 50000
    seed = 0
    plotting = False

    run_consolidate(correct_csv,
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
    assert exists(out_dir)
    assert exists(join(out_dir, "experiment_conditions.json"))
    nmf_folder = join(out_dir, "NMF")
    assert exists(nmf_folder)
    assert exists(join(nmf_folder, "mixing_matrix.csv"))
    assert exists(join(nmf_folder, "signal_matrix.ftr"))

    bdg_folder = join(out_dir, BEDGRAPH_FOLDER)
    assert exists(bdg_folder)

    files_reference = os.path.join(out_dir, 'experiment_conditions.json')
    conditions = extract_conditions(
        files_reference, control_label=control_label)
    for cond in conditions[1:]:
        assert exists(join(bdg_folder, f"{cond}_DecoDen.bdg"))
