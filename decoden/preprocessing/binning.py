import numpy as np
import pandas as pd
from tqdm import tqdm
from collections import namedtuple
import os
import csv                                

def preprocess_bedgraph(fname, bin_size=25, output_fname=None, output_dir=None):
    """Preprocess tiled bedgraph files

    Args:
        fname (string): path to input file
        bin_size (int, optional): width of genomic bin. Defaults to 25.
        output_fname (string, optional): name of the output file. Defaults to None.
        output_dir (string, optional): path to output directory. Defaults to None.
    """
    with open(fname, 'r') as fp:
        for count, line in enumerate(fp):
            pass
    nlines = count+1
    
    Row = namedtuple("Row", ["seqname", "start", "end", "value"])
    if output_fname is None:
        output_fname = os.path.split(fname)[-1].rsplit(".", 1)[0]
    st_bp = 1
    bin_data = []
    
    with open(fname, 'r') as f_in:          # Read lines separately
        reader = csv.reader(f_in, delimiter='\t')
        # Select first row
        line = next(reader)
        row = Row(line[0], int(line[1]), int(line[2]), float(line[3]))
        current_chr = row.seqname
        print(f"Processing chromosome {current_chr}")
        while st_bp<row.start:
            st_bp += bin_size

        # Start looping through rows
        outpath = output_fname + "_binned_py.tsv"
        if output_dir is not None:
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            outpath = os.path.join(output_dir, outpath)

        pbar = tqdm(total=nlines)
        with open(outpath, "w") as outfile:
            while line is not None:
                end_bp = st_bp + bin_size - 1
                for ix, bp in enumerate(range(st_bp, end_bp+1)):
                    while bp>row.end:
                        try:
                            line = next(reader)
                            row = Row(line[0], int(line[1]), int(line[2]), float(line[3]))
                            pbar.update(1)
                        except StopIteration:
                            line = None
                            break
                        if row.seqname != current_chr:
                            break

                    if row.seqname == current_chr:
                        bin_data.append(row.value)
                    else:
                        break
                    
                outfile.write("\t".join([current_chr, str(st_bp), str(end_bp), str(np.mean(bin_data))])+"\n")
                bin_data = []
                if row.seqname != current_chr:
                    st_bp = 1
                    current_chr = row.seqname
                    print(f"Processing chromosome {current_chr}")
                else:
                    st_bp += bin_size

                while st_bp>row.end:
                    try:
                        line = next(reader)
                        row = Row(line[0], int(line[1]), int(line[2]), float(line[3]))
                        pbar.update(1)
                    except StopIteration:
                        line = None
                        break
         