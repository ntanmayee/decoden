import typer
import os
import json
from os.path import join
from pathlib import Path
from typing import Optional, List
from decoden.preprocessing.pipeline import run_preprocessing
from decoden.utils import print_message, extract_conditions
from decoden.denoising.nmf import run_NMF
from decoden.denoising.hsr import run_HSR, run_HSR_replicates

from decoden.apps.denoise_app import denoise_app
from decoden.apps.run_app import run_app

app = typer.Typer()

app.add_typer(denoise_app, name="denoise")
app.add_typer(run_app, name="run")

@app.callback()
def callback():
    """
    Multi-condition ChIP-Seq Analysis with DecoDen
    """
    print_message()
    

@app.command("preprocess")
def preprocess(
    input_csv: Optional[Path] = typer.Option(None, "--input_csv", "-i", help="""Path to CSV file with information about 
                                            experimental conditions. Must contain `filepath`, `exp_name` and `is_control` columns. 
                                            Control/input should be the first condition. Input files can be in BED/BAM format."""), 
    bin_size: int = typer.Option(200, "--bin_size", "-bs", help="""Aize of genomic bin for tiling. 
                                Recommended value is 10-200. Smaller bin size increases space and runtime, larger binsizes may occlude small variations. 
                                """), 
    num_jobs: int = typer.Option(1, "--num_jobs", "-n", help="Number of parallel jobs for preprocessing."), 
    out_dir: Optional[Path] = typer.Option(None, "--out_dir", "-o", help="Path to directory where all output files will be written"), 
    ):
    """
    Preprocess data to be in the correct format for DecoDen
    """
    assert input_csv is not None, "Required input csv"
    assert out_dir is not None, "Required output directory"
    assert bin_size>0, "Invalid bin size"
    assert num_jobs>0, "Invalid number of jobs"
    
    typer.echo("Preprocessing data")
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    _ = run_preprocessing(input_csv, bin_size, num_jobs, out_dir)






@app.command()
def detect():
    """
    Detect peaks
    """
    typer.echo("Detecting peaks")

