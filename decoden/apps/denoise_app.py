import typer
import json
from os.path import join
from pathlib import Path
from typing import Optional, List
from decoden.denoising.nmf import run_NMF
from decoden.denoising.hsr import run_HSR, run_HSR_replicates


denoise_app = typer.Typer()




@denoise_app.command("consolidate")
def denoise_consolidated(
    data_folder:  Optional[Path] = typer.Option(None, "--data_folder", "-df", help="Path to preprocessed data files in BED format"), 
    files_reference: Optional[Path] = typer.Option(None, "--files_reference", "-f", help="""Path to JSON file with experiment conditions. 
                        If you used DecoDen for pre-processing, use the `experiment_conditions.json` file"""), 
    conditions: List[str] = typer.Option(None, "--conditions", "-c", help="List of experimental conditions. First condition MUST correspond to the control/input samples."), 
    out_dir: Optional[Path] = typer.Option(None, "--out_dir", "-o", help="Path to directory where all output files will be written"), 
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
    Run decoden to denoise and pool your data
    """

    wmatrix, mmatrix, data_noBL, mask, conditions_counts = run_NMF(data_folder, 
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
    
    

@denoise_app.command("replicates")
def denoise_replicates(
    data_folder:  Optional[Path] = typer.Option(None, "--data_folder", "-df", help="Path to preprocessed data files in BED format"), 
    files_reference: Optional[Path] = typer.Option(None, "--files_reference", "-f", help="""Path to JSON file with experiment conditions. 
                        If you used DecoDen for pre-processing, use the `experiment_conditions.json` file"""), 
    conditions: List[str] = typer.Option(None, "--conditions", "-c", help="List of experimental conditions. First condition MUST correspond to the control/input samples."), 
    out_dir: Optional[Path] = typer.Option(None, "--out_dir", "-o", help="Path to directory where all output files will be written"), 
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
    Run decoden to denoise your replicates individually
    """
    typer.echo("Running DecoDen on individual replicates")

    
    wmatrix, mmatrix, data_noBL, mask, conditions_counts = run_NMF(data_folder, 
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
