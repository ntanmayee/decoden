# Multi-condition ChIP-Seq Analysis with DecoDen
[![DOI:10.1101/2022.10.18.512665](https://img.shields.io/badge/DOI-10.1101/2022.10.18.512665-B31B1B.svg)](https://doi.org/10.1101/2022.10.18.512665)
![GPLv3 license](https://img.shields.io/github/license/ntanmayee/DecoDen)

![DecoDen Schematic](utils/decoden_schematic.png "DecoDen")

DecoDen uses replicates and multi-histone ChIP-Seq experiments for a target cell type to learn and remove shared biases from fragmentation, PCR amplification and seqeunce mappability. 

Details about DecoDen are in the paper [**Multi-histone ChIP-Seq Analysis with DecoDen**](https://www.biorxiv.org/content/10.1101/2022.10.18.512665v1).


## Installation

DecoDen is available as a python package on PyPi and Bioconda (**SOON**). To ensure the dependencies are satisfied with the correct package version, we recommend the use of [Anaconda](https://www.anaconda.com/) to create suitable environment. Alternative solutions like [https://virtualenv.pypa.io/en/latest/](Virtualenv) and [Poetry](https://python-poetry.org/) work as well, however you will need to ensure the installation of the external dependencies [BEDOPS](https://github.com/bedops/bedops) and [bedtools](https://github.com/arq5x/bedtools2).

### Temporary development instructions

Section to remove once the package is published. 
Clone the github repository and run the following commands

```sh
# Create and activate the environment
conda create -n dden python=3.9
conda activate dden

# Install the external dependencies and DecoDen
conda install -c bioconda bedops bedtools

# From the cloned folder
poetry install
```


### Setup with Conda

To setup a suitable environment with Conda, use the following commands:
```sh
# Create and activate the environment
conda create -n dden python=3.9
conda activate dden

# Install the external dependencies and DecoDen
conda install -c bioconda bedops bedtools decoden

# Verify the correct installation of DecoDen
decoden --version
decoden --help
```


## Decoden Pipeline

DecoDen is employed to denoise ChIP-seq experiments, starting from *.bam* files.
The pipeline involves three main steps, which can be run separately or all at once: preprocessing, denoising, and peak calling.

In the *preprocessing* setp, the *.bam* files are processed to filted duplicate reads, extend reads and finally tile the measured counts to obtain a binned signal. DecoDen makes use of the external tools BEDOPS, bedtools, and some preprocessing functionalities from MACS2.

This functionality is covered by the `decoden preprocess` command:

```sh
# Print the possible arguments for the preprocess command
decoden preprocess --help

# Example of a call to run the preprocessing step
decoden preprocess -i "samples.csv" -o "output_directory" -bs 200 -n 1
```


The *denoising* procedure is the core of the DecoDen algorithm. It involves a sequential use of Non-negative Matrix Factorization (NMF) and Half-Sibling Regression (HSR).
The denoising step is available with two distinct functionalities: **denoise consolidate** performs the two aforementioned steps to generate a consolidated signal for each experimental condition (i.e. histone modification), while **denoise replicates** results in an adjusted signal for each input replicate. 

This process is available through the commands `decoden denoise consolidate` and `decoden denoise replicates`:

```sh
# Print the possible arguments for the denoising commands
decoden denoise conoslidate --help
decoden denoise replicates --help


# Example of calls to denoise the preprocessed data
decoden denoise consolidate --files_reference "output_directory/experiment_conditions.json" \
        --control_label "control" \
        --out_dir "output_directory" \
        --blacklist_file "data/annotations/hg19-blacklist.v2.bed" \
        --plotting

decoden denoise replicates --files_reference "output_directory/experiment_conditions.json" \
        --control_label "control" \
        --out_dir "output_directory" \
        --blacklist_file "data/annotations/hg19-blacklist.v2.bed" \
        --plotting
```

Finally, the peak calling step is available only for signals processed using the **denoise replicates** command, as it requires multiple replicates for statistical hypothesis testing. 

To run the detection process use the command `decoden detect`:

```sh
# Print the possible arguments for the preprocess command
decoden detect --help

#TODO example
```


The full DecoDen pipeline can be run using the `decoden run` command, which is available in the two modalities for consolidation or adjusting individual replicates:


```sh
# Print the possible arguments for the complete pipeline commands
decoden run conoslidate --help
decoden run replicates --help


# Example of calls to run the full DecoDen pipeline
decoden run consolidate -i "samples.csv" -o "output_directory" -bs 200 -n 2 \
    --control_label "control" \
    --blacklist_file "hg19-blacklist.v2.bed" \
    --plotting

decoden run replicates -i "samples.csv" -o "output_directory" -bs 200 -n 2 \
    --control_label "control" \
    --blacklist_file "hg19-blacklist.v2.bed" \
    --plotting
```



## Quick Start

### Input data


Running decoden requires three inputs:
- The *.bam* or *.bed* files for your ChIP-seq experiments
- A file with the annotations for the ChIP-seq experiments files in *.csv* format
- (Optional) a *.bed* file with the blacklisted regions to exclude for the target genome alignment

The annotations .csv must contain the following columns:

| Column | Description |
|---|---|
| `filepath` | the path pointing to the .bam file for the sample. It can be either an absolute path, or a relative path from the directory where the .csv file is saved |
| `exp_name` | the label for the histone modification corresponding to the sample (e.g. "control", "H3K4me3") |
| `is_control` | a binary indicator to mark the control/input samples. It should be 1 for the control samples and 0 for the treatment samples |
| `replicate` | the numbering of the replicate corresponding to the samples. For N replicates corresponding to the same condition (`exp_name`) it should not present repeated values.
| `cell_type` | the cell type of the sample assayed. Used to generate the output file names if `sample_label` is not provided.
| `sample_label` (OPTIONAL) | if provided, it is used to name the output files corresponding to the sample.


An example of a suitable annotation .csv is as follows:

|filepath                      |exp_name|is_control|replicate|cell_type       |sample_label      |
|------------------------------|--------|----------|---------|----------------|------------------|
|wce/ENCFF234NNJ_chr21.bam     |control |1         |1        |transverse colon|MySample1_WCE     |
|wce/ENCFF346JZT_chr21.bam     |control |1         |2        |transverse colon|MySample1_WCE     |
|h3k4me3/ENCFF404DOT_chr21.bam |h3k4me3 |0         |1        |transverse colon|MySample2_H3K4me3 |
|h3k4me3/ENCFF779XRN_chr21.bam |h3k4me3 |0         |2        |transverse colon|MySample2_H3K4me3 |
|h3k27me3/ENCFF623DRR_chr21.bam|h3k27me3|0         |1        |transverse colon|MySample3_H3K27me3|
|h3k27me3/ENCFF228ABC_chr21.bam|h3k27me3|0         |2        |transverse colon|MySample3_H3K27me3|


A skeleton file for the sample annotation, which can then be edited, can be generated using the command `decoden create_csv`:

```sh
# Create an annotation csv with the correct structure
decoden create_csv --sample_label
```


### DecoDen outputs

All the results produced by DecoDen are saved in the output directory specified through the `-o` option.

The *NMF* folder contains the intermediate results produced by the factorization step. Optionally, with the `--plotting` argument several plots for the extracted matrices are produced. This step is a valuable sanity check to verify the choice of regularization parameters.

The *output_bedgraph_files* directory contains the outputs of DecoDen for each histone modification (if run in `consolidate` mode) or for each replicate (if run in `replicates` mode) in .bdg files. These files can directly be uploaded to the [UCSC Genome Browser](https://genome.ucsc.edu/) for visualization.

If the DecoDen pipeline was run in `replicates` mode, the output folder will also contain a directory `called_peaks` that contains the peaks detected from the adjusted replicates in *.bed* format, where the score represent -log10(pval).

## Example: how to run DecoDen on your data

Here is a brief example of the workflow to run DecoDen on your data:

```sh
# Create the annotation csv
decoden create_csv --sample_label

# Edit the generated file `samples.csv` with the information for your ChIP-seq samples

# Run the DecoDen pipeline
decoden run replicates -i "samples.csv" -o "output_directory" -bs 200
```





## Bug Reports and Suggestions for Improvement
Please [raise an issue](https://github.com/ntanmayee/DecoDen/issues/new) if you find bugs or if you have any suggestions for improvement.
