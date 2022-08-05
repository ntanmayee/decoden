# DecoDen:  Using replicates to remove cell-type specific bias in multi-histone ChIP-Seq

This is the accompanying code for the paper 
``` 
Using replicates to remove cell-type specific bias in multi-histone ChIP-Seq
``` 

## Dependencies
DecoDen depends on MACS2, BEDOPS and BEDTools.

```sh
conda install -c bioconda macs2 bedops bedtools
```

## Usage

### Preprocessing


### Running DecoDen

To run DecoDen, the data must be preprocessed into bedGraph format and binned correctly. The algorithm requires a small .json file that indicates the correspondance of each file to the experimental condition.

```javascript
{
    "control_1.bdg": "control",
    "control_2.bdg": "control",
    "h3k27me3_1.bdg": "H3K27me3",
    "h3k27me3_2.bdg": "H3K27me3",
    "h3k27me3_3.bdg": "H3K27me3",
    "h3k4me3_1.bdg": "H3K4me3",
    "h3k4me3_2.bdg": "H3K4me3",
    "h3k4me3_3.bdg": "H3K4me3"
}
```

To run DecoDen, an example command would be:
```
python run_decoden --data_folder "data/my_experiment" \ # where the .bdg files are saved
                   --output_folder "outputs/my_experiment_results" \ # where to save the results
                   --files_reference "data/my_experiment_files.json" \ # the aforementioned mapping
                   --blacklist_file "data/annotations/hg19-blacklist.v2.bed" \
                   --conditions "control" "H3K27me3" "H3K4me3" \ # Ordering of the experimental conditions. The first one must be the control.
                   --control_cov_threshold 1.0 \ # Minimum coverage for the training data for the NMF
                   --n_train_bins 300000 \ # Number of training bins for the extraction of the mixing matrix
                   --chunk_size 100000 \ # For the processing of the signal matrix
                   --seed 42 \ # Random state for reproductibility
                   --alpha_W 0.01 \ # Regularisation for the signal matrix
                   --alpha_H 0.001 \ # Regularisation for the mixing matrix

```