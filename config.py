import csv
import os
from collections import defaultdict


with open("Data/sequences.txt", "r") as f:
    reader = csv.reader(f, delimiter="\t")
    sequence_folders = defaultdict(list)
    for row in reader:
        sequence_folders[row[1].split('-')[0]].append(row[4])
    for sample in sequence_folders.keys():
        for folder in set(sequence_folders[seq]):
            print(os.listdir(folder))
        