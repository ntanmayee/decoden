import pytest

from os.path import join, exists

from decoden.main import *

correct_csv = join("tests", "sample_data", "samples.csv")
output_directory = "decoden_output"


@pytest.fixture(scope="session")
def tmp_session_directory(tmp_path_factory):
    tmp_folder = tmp_path_factory.mktemp("tmp_test_folder")
    return tmp_folder


def test_output_files_created(tmp_session_directory):
    bin_size = 200
    num_jobs = 1
    out_dir = join(tmp_session_directory, output_directory)
    assert exists(correct_csv)
    
    preprocess(correct_csv, bin_size, num_jobs, out_dir)
    
    assert exists(out_dir)
