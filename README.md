<p align="center">
    <img src="utils/logo.png" alt="decoden logo">
</p>

# Multi-assay ChIP-Seq Analysis with DecoDen
[![DOI:10.1101/2022.10.18.512665](https://img.shields.io/badge/DOI-10.1101/2022.10.18.512665-B31B1B.svg)](https://doi.org/10.1101/2022.10.18.512665)
![GPLv3 license](https://img.shields.io/github/license/ntanmayee/DecoDen)

DecoDen uses replicates and multi-histone ChIP-Seq experiments for a target cell type to learn and remove shared biases from fragmentation, PCR amplification and sequence mappability.

## Installation
1. Install [Poetry](https://python-poetry.org/)
2. Clone the repository and install with poetry
```sh
# Clone the repository
git clone git@github.com:ntanmayee/DecoDen.git
cd decoden

# Install the external dependencies and DecoDen
conda install pyarrow poetry
poetry install
```

## Quick Start

### Input data

Running decoden requires two inputs:
1. Aligned reads in `.bam` format from ChIP-Seq experiments
2. Sample annotation file in `.csv` format

#### Auto-generate a sample annotation file
To generate a skeleton sample annotation file, run -
```sh
decoden create_csv 
```
This will create `samples.csv` in your current directory. Edit this file and fill in the columns with appropriate information. There are more details [here](https://github.com/ntanmayee/decoden/wiki/Preparing-a-sample-annotation-file).

### Run DecoDen
Run the DecoDen pipeline with default parameters -
```sh
decoden run replicates -i samples.csv -o output_directory
```

## Detailed Usage Guidelines
The following commands are available in DecoDen. Please click on the links to know more about them.
| Command | Description |
|---|---|
| [`create_csv`](https://github.com/ntanmayee/decoden/wiki/Preparing-a-sample-annotation-file) | Create a skeleton sample annotation file |
| [`run`](https://github.com/ntanmayee/decoden/wiki/Run-the-DecoDen-pipeline) | Run the full DecoDen pipeline to preprocess end denoise BAM/BED files |
| [`preprocess`](https://github.com/ntanmayee/decoden/wiki/Preprocess-alignment-files) | Pre-process BAM/BED data to be in the correct format for running DecoDen |
| [`denoise`](https://github.com/ntanmayee/decoden/wiki/Run-denoising) | Run the denoising step of DecoDen on suitably preprocessed data |
| [`detect`](https://github.com/ntanmayee/decoden/wiki/Detect-peaks) | Detect peaks in the processed DecoDen signals |

There is more helpful information in the [wiki](https://github.com/ntanmayee/DecoDen/wiki).


## Bug Reports and Suggestions for Improvement
Please [raise an issue](https://github.com/ntanmayee/DecoDen/issues/new) if you find bugs or if you have any suggestions for improvement.
