import pytest
import json
from os.path import join, exists

from decoden.main import *


def test_output_files_created(tmp_session_directory, correct_csv):
    bin_size = 200
    num_jobs = 1
    out_dir = tmp_session_directory
    assert exists(correct_csv)
    
    preprocess(correct_csv, bin_size, num_jobs, out_dir)
    
    assert exists(out_dir)
    assert exists(join(out_dir, "experiment_conditions.json"))
    
    with open(join(out_dir, "experiment_conditions.json"), "r") as f:
        experiment_conditions = json.load(f)
    for k in experiment_conditions.keys():
        assert exists(join(out_dir, k))
