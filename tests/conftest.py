import pytest

from os.path import join

# correct_csv = join("sample_data", "samples.csv")
# output_directory = "decoden_output"


@pytest.fixture(autouse=True)
def change_test_dir(request, monkeypatch):
    monkeypatch.chdir(request.fspath.dirname)
    
# @pytest.fixture(scope="session")
# def tmp_session_directory(tmp_path_factory):
#     tmp_folder = tmp_path_factory.mktemp("tmp_test_folder")
#     return tmp_folder

@pytest.fixture(scope="session")
def tmp_session_directory():
    return "decoden_output_devel"



@pytest.fixture(scope="session")
def correct_csv():
    return join("sample_data", "samples.csv")

@pytest.fixture(scope="session")
def bl_file():
    return join("sample_data", "annotations", "hg19-blacklist.v2.bed")

