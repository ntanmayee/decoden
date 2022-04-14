import pandas as pd
from joblib import Parallel, delayed
import subprocess
import os


def compute_coverage(filename):
    subprocess.run(
        f'bedtools genomecov -i {base_save_path}/{filename} -g ../sh/hg19.chrom.sizes -bga > {base_save_path}/{filename}_coverage.bdg', 
        shell=True)

def tile(filename):
    # tile only chr1
    # filename = f'{base_save_path}/{filename}_coverage.bdg'

    subprocess.run("""sort-bed %s | awk -vOFS="\t" '{ print $1, $2, $3, ".", $4 }' - > %s_signal.bed""" % (filename, filename), shell=True)
    subprocess.run(f"""bedops --chrom {chrom} --merge {filename}_signal.bed | bedops --chop {bin_size} - > {filename}_bins.bed""", shell=True)
    subprocess.run(f'bedmap --echo --mean --prec 3 {filename}_bins.bed {filename}_signal.bed > {filename}_tiled.bed', shell=True)
    os.remove(f'{filename}_signal.bed')
    os.remove(f'{filename}_bins.bed')