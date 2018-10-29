import csv
import numpy
from skbio import DistanceMatrix as dm
from skbio.stats.distance import anosim
from skbio.stats.distance import permanova
from skbio.stats.distance import permdisp
from scipy.spatial.distance import braycurtis

fn = '18S_controls.csv'
with open(fn, 'rU') as f:
    rows = []
    for row in csv.reader(f):
        rows.append(row)

sample_names = rows[0][1:]
samples = []
for a in range(len(rows[0])):
    if a > 0:
        this_sample = []
        for b in range(len(rows)):
            if b > 0:
                this_sample.append(float(rows[b][a]))
        samples.append(this_sample)
"""
only_samples = ['LR', 'SR']
new_samples, new_names = [], []
for a in range(len(sample_names)):
    for b in range(len(only_samples)):
        if sample_names[a] == only_samples[b]:
            new_samples.append(samples[a])
            new_names.append(sample_names[a])
samples = new_samples
sample_names = new_names
print(len(samples), len(sample_names))
"""
        
sam_dm = dm.from_iterable(samples, metric=braycurtis)
pdisp = permdisp(sam_dm, sample_names, column=None, test='median', permutations=999)
print(pdisp)
asim = anosim(sam_dm, sample_names, column=None, permutations=999)
print(asim)
perm = permanova(sam_dm, sample_names, column=None, permutations=999)
print(perm)