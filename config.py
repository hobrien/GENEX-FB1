import csv
import os
import re
from collections import defaultdict


with open("Data/sequences.txt", "r") as f:
    reader = csv.reader(f, delimiter="\t")
    samples = defaultdict(list)
    filenames = defaultdict(list)
    folders = defaultdict(list)
    for row in reader:
        sample = row[1]
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
        print(": ".join([brain_bank_id, ",".join(sorted(sequences))]))
