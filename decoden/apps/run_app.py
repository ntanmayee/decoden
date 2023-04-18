import typer
import json
from os.path import join
from pathlib import Path
from typing import Optional, List
from decoden.utils import extract_conditions
from decoden.preprocessing.pipeline import run_preprocessing
from decoden.denoising.nmf import run_NMF
from decoden.denoising.hsr import run_HSR, run_HSR_replicates



run_app = typer.Typer()



@run_app.command("consolidate")
def run_consolidate(
    input_csv: Optional[Path] = typer.Option(None, "--input_csv", "-i", help="""Path to CSV file with information about 
                                            experimental conditions. Must contain `filepath`, `exp_name` and `is_control` columns. 
                                            Control/input should be the first condition. Input files can be in BED/BAM format."""), 
    bin_size: int = typer.Option(200, "--bin_size", "-bs", help="""Aize of genomic bin for tiling. 
                                Recommended value is 10-200. Smaller bin size increases space and runtime, larger binsizes may occlude small variations. 
                                """), 
    num_jobs: int = typer.Option(1, "--num_jobs", "-n", help="Number of parallel jobs for preprocessing."), 
    out_dir: Optional[Path] = typer.Option(None, "--out_dir", "-o", help="Path to directory where all output files will be written"), 
    
    
    control_condition: str  = typer.Option("control", "--control_condition", "-con", help="The label for the control/input samples."), 
    
    blacklist_file: Optional[Path] = typer.Option(None, "--blacklist_file", "-bl", help="Path to blacklist file. Make sure to use the blacklist that is appropriate for the genome assembly/organism."), 
    alpha_W: float  = typer.Option(0.01, "--alpha_W", "-aW", help="Regularisation for the signal matrix."), 
    alpha_H: float  = typer.Option(0.001, "--alpha_H", "-aH", help="Regularisation for the mixing matrix."), 
    control_cov_threshold: float  = typer.Option(0.1, "--control_cov_threshold", "-cov", help="""Threshold for coverage in control samples. Only genomic bins above this threshold will be used. 
                                It is recommended to choose a value larger than 1/bin_size."""), 
    n_train_bins: int  = typer.Option(50000, "--n_train_bins", "-nt", help="Number of genomic bins to be used for training."), 
    chunk_size: int  = typer.Option(50000, "--chunk_size", "-ch", help="Chunk size for processing the signal matrix. Should be smaller than `n_train_bins`."), 
    seed: int  = typer.Option(0, "--seed", "-s", help="Random state."), 
    plotting: bool  = typer.Option(False, "--plotting", "-p", help="Plot sanity checks for extracted matrices."), 

):
    """
    Preprocess and denoise data with pooling
    """
    typer.echo("Running DecoDen")
    
    assert input_csv is not None, "Required input csv"
    assert out_dir is not None, "Required output directory"
    assert bin_size>0, "Invalid bin size"
    assert num_jobs>0, "Invalid number of jobs"
    
    typer.echo("Preprocessing data")
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    tiled_files = run_preprocessing(input_csv, bin_size, num_jobs, out_dir)
    files_reference = join(out_dir, 'experiment_conditions.json')
    conditions = extract_conditions(files_reference, control_condition=control_condition)
    wmatrix, mmatrix, data_noBL, mask, conditions_counts = run_NMF(out_dir, 
                                                                    files_reference, 
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
    # Perform HSR to remove multiplicative noise
    hsr_df = run_HSR(wmatrix, mask, conditions)
    hsr_df.reset_index().to_feather(join(out_dir, "HSR_results.ftr"))
    
    typer.echo("\nDecoDen complete!")
    


@run_app.command("replicates")
def run_replicates(    input_csv: Optional[Path] = typer.Option(None, "--input_csv", "-i", help="""Path to CSV file with information about 
                                            experimental conditions. Must contain `filepath`, `exp_name` and `is_control` columns. 
                                            Control/input should be the first condition. Input files can be in BED/BAM format."""), 
    bin_size: int = typer.Option(200, "--bin_size", "-bs", help="""Aize of genomic bin for tiling. 
                                Recommended value is 10-200. Smaller bin size increases space and runtime, larger binsizes may occlude small variations. 
                                """), 
    num_jobs: int = typer.Option(1, "--num_jobs", "-n", help="Number of parallel jobs for preprocessing."), 
    out_dir: Optional[Path] = typer.Option(None, "--out_dir", "-o", help="Path to directory where all output files will be written"), 
    
    
    control_condition: str  = typer.Option("control", "--control_condition", "-con", help="The label for the control/input samples."), 
    
    blacklist_file: Optional[Path] = typer.Option(None, "--blacklist_file", "-bl", help="Path to blacklist file. Make sure to use the blacklist that is appropriate for the genome assembly/organism."), 
    alpha_W: float  = typer.Option(0.01, "--alpha_W", "-aW", help="Regularisation for the signal matrix."), 
    alpha_H: float  = typer.Option(0.001, "--alpha_H", "-aH", help="Regularisation for the mixing matrix."), 
    control_cov_threshold: float  = typer.Option(0.1, "--control_cov_threshold", "-cov", help="""Threshold for coverage in control samples. Only genomic bins above this threshold will be used. 
                                It is recommended to choose a value larger than 1/bin_size."""), 
    n_train_bins: int  = typer.Option(50000, "--n_train_bins", "-nt", help="Number of genomic bins to be used for training."), 
    chunk_size: int  = typer.Option(50000, "--chunk_size", "-ch", help="Chunk size for processing the signal matrix. Should be smaller than `n_train_bins`."), 
    seed: int  = typer.Option(0, "--seed", "-s", help="Random state."), 
    plotting: bool  = typer.Option(False, "--plotting", "-p", help="Plot sanity checks for extracted matrices."), 

):
    """
    Preprocess and denoise individual replicates
    """
    typer.echo("Running DecoDen (replicates-specific)")

    
    assert input_csv is not None, "Required input csv"
    assert out_dir is not None, "Required output directory"
    assert bin_size>0, "Invalid bin size"
    assert num_jobs>0, "Invalid number of jobs"
    
    typer.echo("Preprocessing data")
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    tiled_files = run_preprocessing(input_csv, bin_size, num_jobs, out_dir)
    files_reference = join(out_dir, 'experiment_conditions.json')
    conditions = extract_conditions(files_reference, control_condition=control_condition)
    wmatrix, mmatrix, data_noBL, mask, conditions_counts = run_NMF(out_dir, 
                                                                    files_reference, 
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

    # Perform HSR to remove multiplicative noise
    hsr_df = run_HSR_replicates(data_noBL, wmatrix, mmatrix, mask, conditions, conditions_counts)
    hsr_df.reset_index().to_feather(join(out_dir, "HSR_results_replicates.ftr"))
    
    typer.echo("\nDecoDen (replicate specific) complete!")
