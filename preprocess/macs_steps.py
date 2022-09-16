import argparse
import subprocess
import os
import uuid
from preprocess.logger import logger 
from pathlib import Path

def filter_duplicate_reads(filename, out_dir, name):
    logger.info(f'Filtering duplicate reads for {filename}')

    out_filepath = os.path.join(out_dir, 'data', os.path.splitext(os.path.basename(filename))[0] + '_filterdup.bed')
    logger.info(f'Filtered filepath will be {out_filepath}')

    subprocess.run(f"macs2 filterdup -i {filename} --keep-dup=1 -o {out_filepath}", shell=True)
    return out_filepath

def extend_reads(filepath, extend_length, is_control=0):
    logger.info(f'Extending reads for {filepath}')
    out_filepath = os.path.join(os.path.splitext(filepath)[0] + '_pileup.bdg')
    logger.info(f'Extended reads filepath will be {out_filepath}')

    if is_control == 0:
        subprocess.run(f'macs2 pileup -i {filepath} -o {out_filepath} --extsize {extend_length} -f BED', shell=True)
    else:
        logger.info(f'{filepath} is a control file. Extending reads in both directions')
        extend_length = extend_length // 2
        subprocess.run(f'macs2 pileup -i {filepath} -o {out_filepath} --extsize {extend_length} -B  -f BED', shell=True)
    os.remove(filepath)
    return out_filepath

def get_fragment_length(filtered_filepath):
    logger.info(f'Getting fragment length for {filtered_filepath}')

    result = subprocess.run(f'macs2 predictd -i {filtered_filepath} -g hs -m 5 50', capture_output=True, text=True, shell=True)
    try:
        fragment_length = int([s for s in result.stderr.split('\n') if 'tag size is' in s][0].split()[-2])
    except:
        logger.error(f'Unable to compute fragment length for {filtered_filepath}')
        raise Exception('Unable to compute fragment length')
    return fragment_length

def tile(extended_filepath, bin_size):
    logger.info(f'Tiling {extended_filepath} into bins of size {bin_size}')
    out_filepath = os.path.splitext(extended_filepath)[0]+ '_tiled.bed'
    name = str(uuid.uuid4())

    command = f'''sort-bed {extended_filepath} | awk -vOFS="\t" '{{ print $1, $2, $3, ".", $4 }}' - > signal_{name}.bed;'''

    subprocess.run(command, shell=True)
    subprocess.run(f'bedops --merge signal_{name}.bed | bedops --chop {bin_size} - > bins_{name}.bed;', shell=True)
    subprocess.run(f"bedmap --echo --delim '\t' --mean --prec 3 bins_{name}.bed signal_{name}.bed > {out_filepath};", shell=True)

    os.remove(f'signal_{name}.bed')
    os.remove(f'bins_{name}.bed')
    os.remove(extended_filepath)
    logger.info(f'Tiled files in {out_filepath}')
    
    return out_filepath

def run_pipeline(input_filepath, name, out_dir, is_control, bin_size):
    logger.info(f'Running preprocessing pipeline for {input_filepath}')

    filterdup_filepath = filter_duplicate_reads(input_filepath, out_dir, name)
    fragment_length = get_fragment_length(filterdup_filepath)
    extended_filepath = extend_reads(filterdup_filepath, fragment_length, is_control)
    tiled_filepath = tile(extended_filepath, bin_size)

    return tiled_filepath

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', "--input_filepath", required=True)
    parser.add_argument('-n', "--name", required=True)
    parser.add_argument('-o', "--out_dir", required=True)
    parser.add_argument('-c', "--is_control", default=0)
    parser.add_argument('-bs', "--bin_size", default=200, type=int)

    args = parser.parse_args()

    tiled_filepath = run_pipeline(args.input_filepath, args.name, args.out_dir, args.is_control, args.bin_size)
#     tile('temp/e1_filterdup_pileup.bdg', 200)