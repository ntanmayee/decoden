import typer
import os
import pandas as pd
from pathlib import Path
from typing import Optional

from decoden.decoden_pipeline import _decoden_pipeline
from decoden.utils import print_message


__version__ = "0.1.0"

app = typer.Typer()
# denoise_app = typer.Typer()
# run_app = typer.Typer()


# app.add_typer(denoise_app, name="denoise", help="""Run the denoising step of DecoDen on suitably preprocessed data.
#               Use `decoden denoise --help` for more info""")
# app.add_typer(run_app, name="run", help="""Run the full decoden pipeline to preprocess end denoise BAM/BED files.
#               Use `decoden run --help` for more info""")


def version_callback(value: bool):
    if value:
        typer.echo(f"DecoDen Version: {__version__}")
        raise typer.Exit()


@app.callback()
def callback(version: bool = typer.Option(
        None, "--version", callback=version_callback, is_eager=True, help="Display DecoDen version"),):
    """
    Multi-condition ChIP-Seq Analysis with DecoDen
    """
    print_message()

@app.command("create_csv")
def create_csv(
    out_dir: Optional[Path] = typer.Option(
        None, "--out_dir", "-o", help="Path to directory where all output files will be written"),
    sample_label: bool = typer.Option(
        False, "--sample_label", "-sl", help="Flag to add the optional column `sample_label`.")
):
    """
    Create a sample .csv file with the correct fields for annotating the data for DecoDen.
    Edit the file with the correct paths before executing the pipeline.
    """
    cols = ["filepath", "exp_name", "is_control", "replicate", "cell_type"]
    if sample_label:
        cols.append("sample_label")
    default_vals = ["path/to/sample", "control", 1, 1, "myTissue", "mySampleLabel"]
    df = pd.DataFrame([default_vals[:len(cols)]]*2, columns=cols)
    
    fname = "samples.csv"
    if out_dir is not None:
        fname = os.path.join(out_dir, fname)
    df.to_csv(fname, index=False)
    


@app.command("preprocess")
def preprocess(
    input_csv: Optional[Path] = typer.Option(None, "--input_csv", "-i", help="""Path to CSV file with information about 
                                            experimental conditions. Must contain `filepath`, `exp_name` and `is_control` columns. 
                                            Control/input should be the first condition. Input files can be in BED/BAM format."""),
    bin_size: int = typer.Option(200, "--bin_size", "-bs", help="""Aize of genomic bin for tiling. 
                                Recommended value is 10-200. Smaller bin size increases space and runtime, larger binsizes may occlude small variations. 
                                """),
    num_jobs: int = typer.Option(
        1, "--num_jobs", "-n", help="Number of parallel jobs for preprocessing."),
    out_dir: Optional[Path] = typer.Option(
        None, "--out_dir", "-o", help="Path to directory where all output files will be written"),
    genome_size: str = typer.Option(
        None, "--genome_size", "-gs", help="Size of genome. Give an integer value. Alternatively, `hs` for Homo sapien and `mm` for Mus musculus are also supported."),
):
    """
    Preprocess BAM/BED data to be in the correct format for running DecoDen
    """

    typer.echo("Preprocessing data")

    _decoden_pipeline(["preprocess"],
                      input_csv=input_csv,
                      bin_size=bin_size,
                      num_jobs=num_jobs,
                      out_dir=out_dir,
                      genome_size=genome_size
                      )


@app.command("denoise")
def denoise(
    # data_folder:  Optional[Path] = typer.Option(None, "--data_folder", "-df", help="Path to preprocessed data files in BED format"),
    files_reference: Optional[Path] = typer.Option(None, "--files_reference", "-f", help="""Path to JSON file with experiment conditions. 
                        If you used DecoDen for pre-processing, use the `experiment_conditions.json` file"""),
    control_label: str = typer.Option(
        "control", "--control_label", "-con", help="The label for the control/input samples."),
    genome_size: str = typer.Option(
        None, "--genome_size", "-gs", help="Size of genome. Give an integer value. Alternatively, `hs` for Homo sapien and `mm` for Mus musculus are also supported."),
    # conditions: List[str] = typer.Option(None, "--conditions", "-c", help="List of experimental conditions. First condition MUST correspond to the control/input samples."),
    out_dir: Optional[Path] = typer.Option(
        None, "--out_dir", "-o", help="Path to directory where all output files will be written"),
    blacklist_file: Optional[Path] = typer.Option(
        None, "--blacklist_file", "-bl", help="Path to blacklist file. Make sure to use the blacklist that is appropriate for the genome assembly/organism."),
    alpha_W: float = typer.Option(
        0.01, "--alpha_W", "-aW", help="Regularisation for the signal matrix."),
    alpha_H: float = typer.Option(
        0.001, "--alpha_H", "-aH", help="Regularisation for the mixing matrix."),
    control_cov_threshold: float = typer.Option(0.1, "--control_cov_threshold", "-cov", help="""Threshold for coverage in control samples. Only genomic bins above this threshold will be used. 
                                It is recommended to choose a value larger than 1/bin_size."""),
    n_train_bins: int = typer.Option(
        50000, "--n_train_bins", "-nt", help="Number of genomic bins to be used for training."),
    chunk_size: int = typer.Option(
        50000, "--chunk_size", "-ch", help="Chunk size for processing the signal matrix. Should be smaller than `n_train_bins`."),
    seed: int = typer.Option(0, "--seed", "-s", help="Random state."),
    plotting: bool = typer.Option(
        False, "--plotting", "-p", help="Plot sanity checks for extracted matrices."),
):
    """
    Run the denoising step of DecoDen in pooling mode.
    """

    _decoden_pipeline(["nmf", "hsr"],
                      files_reference=files_reference,
                      control_label=control_label,
                      out_dir=out_dir,
                      blacklist_file=blacklist_file,
                      alpha_W=alpha_W,
                      alpha_H=alpha_H,
                      control_cov_threshold=control_cov_threshold,
                      n_train_bins=n_train_bins,
                      chunk_size=chunk_size,
                      seed=seed,
                      plotting=plotting,
                      genome_size=genome_size
                      )

    typer.echo("\nDecoDen complete!")




@app.command("run")
def run(
    input_csv: Optional[Path] = typer.Option(None, "--input_csv", "-i", help="""Path to CSV file with information about 
                                            experimental conditions. Must contain `filepath`, `exp_name` and `is_control` columns. 
                                            Control/input should be the first condition. Input files can be in BED/BAM format."""),
    bin_size: int = typer.Option(200, "--bin_size", "-bs", help="""Aize of genomic bin for tiling. 
                                Recommended value is 10-200. Smaller bin size increases space and runtime, larger binsizes may occlude small variations. 
                                """),
    num_jobs: int = typer.Option(
        1, "--num_jobs", "-n", help="Number of parallel jobs for preprocessing."),
    out_dir: Optional[Path] = typer.Option(
        None, "--out_dir", "-o", help="Path to directory where all output files will be written"),
    genome_size: str = typer.Option(
        None, "--genome_size", "-gs", help="Size of genome. Give an integer value. Alternatively, `hs` for Homo sapien and `mm` for Mus musculus are also supported."),

    # control_label: str = typer.Option(
    #     "control", "--control_label", "-con", help="The label for the control/input samples."),

    blacklist_file: Optional[Path] = typer.Option(
        None, "--blacklist_file", "-bl", help="Path to blacklist file. Make sure to use the blacklist that is appropriate for the genome assembly/organism."),
    alpha_W: float = typer.Option(
        0.01, "--alpha_W", "-aW", help="Regularisation for the signal matrix."),
    alpha_H: float = typer.Option(
        0.001, "--alpha_H", "-aH", help="Regularisation for the mixing matrix."),
    control_cov_threshold: float = typer.Option(0.1, "--control_cov_threshold", "-cov", help="""Threshold for coverage in control samples. Only genomic bins above this threshold will be used. 
                                It is recommended to choose a value larger than 1/bin_size."""),
    n_train_bins: int = typer.Option(
        50000, "--n_train_bins", "-nt", help="Number of genomic bins to be used for training."),
    chunk_size: int = typer.Option(
        50000, "--chunk_size", "-ch", help="Chunk size for processing the signal matrix. Should be smaller than `n_train_bins`."),
    seed: int = typer.Option(0, "--seed", "-s", help="Random state."),
    plotting: bool = typer.Option(
        False, "--plotting", "-p", help="Plot sanity checks for extracted matrices."),

):
    """
    Preprocess and denoise data in pooling mode.
    """
    typer.echo("Running DecoDen")
    _decoden_pipeline(["preprocess", "nmf", "hsr"],
                      input_csv=input_csv,
                      bin_size=bin_size,
                      num_jobs=num_jobs,
                      out_dir=out_dir,
                      genome_size=genome_size,
                      blacklist_file=blacklist_file,
                      alpha_W=alpha_W,
                      alpha_H=alpha_H,
                      control_cov_threshold=control_cov_threshold,
                      n_train_bins=n_train_bins,
                      chunk_size=chunk_size,
                      seed=seed,
                      plotting=plotting
                      )

    typer.echo("\nDecoDen complete!")



@app.command()
def detect(files_reference: Optional[Path] = typer.Option(None, "--files_reference", "-f", help="""Path to JSON file with experiment conditions. 
                        If you used DecoDen for pre-processing, use the `experiment_conditions.json` file"""),
           out_dir: Optional[Path] = typer.Option(
                       None, "--out_dir", "-o", help="Path to directory where all output files will be written"),
           control_label: str = typer.Option(
                       "control", "--control_label", "-con", help="The label for the control/input samples."),
            pval_alpha: float = typer.Option(
                       0.05, "--pval_alpha", help="Value for multiple hypthesis test correction (BH)."),
            peak_threshold: float = typer.Option(
                0.01, "--peak_threshold", help="Pvalue threshold for the initial identification of peak regions."),
            min_width: int = typer.Option(
                150, "--min_width", help="Minimum width in bp for merging and discarding peaks."),    
    ):
    """
    Detect peaks in the processed DecoDen signals. Must be used with the `replicates` mode.
    """
    typer.echo("Detecting peaks")
    
    _decoden_pipeline(["detect"],
                      files_reference=files_reference,
                      out_dir=out_dir,
                      control_label=control_label,
                      pval_alpha=pval_alpha,
                      peak_threshold=peak_threshold,
                      min_width=min_width
                      )
    
    typer.echo("\nDecoDen Peak Detection complete!")
