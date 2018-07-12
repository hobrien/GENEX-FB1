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
            if 
        samples = defaultdict(list)
        filenames = defaultdict(list)
        folders = defaultdict(list)
        readfiles = {}
    
        brain_bank_id = sample.split('-')[0]
        samples[brain_bank_id].append(sample)
        folders[sample].append(row[4])
        filenames[sample].append(row[0])
    for brain_bank_id in samples.keys():
        sequences = []
        for sample in set(samples[brain_bank_id]):
            for read in range(2):
                folder = folders[sample][read]
                file = filenames[sample][read]
                sequences = sequences + [f for f in os.listdir(folder) if re.match(re.escape(file) + r'.*f(ast)?q(.gz)?$', f)]
        readfiles[brain_bank_id] = " ".join(sorted(sequences))

if __name__ == '__main__':
    print(get_sequences("Data/sequences.txt"))