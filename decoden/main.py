import typer
import os
from pathlib import Path
from typing import Optional
from decoden.preprocessing.pipeline import run_preprocessing

app = typer.Typer()


@app.callback()
def callback():
    """
    Multi-condition ChIP-Seq Analysis with DecoDen
    """

@app.command("preprocess")
def preprocess(input_csv: Optional[Path] = typer.Option(None, "--input_csv", "-i", help="""Path to CSV file with information about 
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
def denoise():
    """
    Run decoden to denoise and pool your data
    """
    typer.echo("Running DecoDen")


@app.command()
def denoise_replicates():
    """
    Run decoden to denoise your replicates individually
    """
    typer.echo("Running DecoDen on individual replicates")


@app.command()
def detect():
    """
    Detect peaks
    """
    typer.echo("Detecting peaks")

