import argparse
import subprocess
import os

def filter_duplicate_reads(filename):
    out_filepath = os.path.join(args.out_dir, args.name + '_filterdup.bed')
    print(f'Filtered filepath will be {out_filepath}')

    subprocess.run(f"macs2 filterdup -i {filename} --keep-dup=1 -o {out_filepath}", shell=True)
    return out_filepath

def extend_reads(filepath, extend_length, is_control=0):
    out_filepath = os.path.join(os.path.splitext(filepath)[0] + '_pileup.bdg')
    print(f'Extended reads filepath will be {out_filepath}')

    if is_control == 0:
        subprocess.run(f'macs2 pileup -i {filepath} -o {out_filepath} --extsize {extend_length} -f BED', shell=True)
    else:
        extend_length = extend_length // 2
        subprocess.run(f'macs2 pileup -i {filepath} -o {out_filepath} --extsize {extend_length} -B  -f BED', shell=True)
    os.remove(filepath)
    return out_filepath

def get_fragment_length(filtered_filepath):
    print(f'Getting fragment length for {filtered_filepath}')

    result = subprocess.run(f'macs2 predictd -i {filtered_filepath} -g hs -m 5 50', capture_output=True, text=True, shell=True)
    try:
        fragment_length = int([s for s in result.stderr.split('\n') if 'tag size is' in s][0].split()[-2])
    except:
        raise Exception('Unable to compute fragment length')
    return fragment_length


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', "--input_filepath", required=True)
    parser.add_argument('-n', "--name", required=True)
    parser.add_argument('-o', "--out_dir", required=True)
    parser.add_argument('-c', "--is_control", default=0)

    args = parser.parse_args()
    print(args)

    filterdup_filepath = filter_duplicate_reads(args.input_filepath)
    fragment_length = get_fragment_length(filterdup_filepath)
    extended_filepath = extend_reads(filterdup_filepath, fragment_length, args.is_control)
