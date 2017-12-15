import csv
import os
import re
from collections import defaultdict

ref_dir="/c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy"

def get_sequences(filename):
    files = defaultdict(list)
    with open(filename, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            sample = row[1].split('-')[0] # remove rep number from sample name
            file = os.path.join(row[4], row[0])
            files[sample].append(file)
    for sample in files.keys():
        files[sample].sort() # need to make sure that R reads follow F reads
    return files

if __name__ == '__main__':
    files = get_sequences("Data/sequences.txt")
    for sample in files.keys():
        print sample, files[sample]
