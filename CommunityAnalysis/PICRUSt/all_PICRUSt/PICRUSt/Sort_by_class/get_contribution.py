import csv
import os

basic_path = '/Users/u1560915/Documents/OneDrive/PhD_Plastic_Oceans/Experiments/MiSeq_Dada/PICRUSt/Sort_by_class/'
folders = sorted(os.listdir(basic_path))
not_this = ['get_contribution.py', '.DS_Store', 'sort_by_class.py', 'taxa_and_seqs.csv', 'Woesearchaeia_None', 'Other_None', 'Babeliae_None', 'Chitinivibrionia_none', 'new_plot']

new_folders = []
for a in range(len(folders)):
    adding = True
    for b in range(len(not_this)):
        if folders[a] == not_this[b]:
            adding = False
    if adding:
        new_folders.append(folders[a])
print new_folders

all_otus, old_otus, prop = [], [], []
for a in range(len(new_folders)):
    os.chdir(basic_path+new_folders[a])
    this_new, this_old = [], []
    with open('otu_table.csv', 'rU') as f:
        rows = []
        for row in csv.reader(f):
            rows.append(row)
    for b in range(len(rows[1])):
        if b > 0:
            all_otus.append(rows[1][b])
            this_new.append(rows[1][b])
    with open('16S_shared.csv', 'rU') as f:
        rows = []
        for row in csv.reader(f):
            rows.append(row)
    for c in range(len(rows)):
        if c > 0:
            this_old.append(rows[c][0])
            old_otus.append(rows[c][0])
    p = str(len(this_new))+'/'+str(len(this_old))
    prop.append(p)

os.chdir(basic_path)
with open('taxa_and_seqs.csv', 'rU') as f:
    rows = []
    for row in csv.reader(f):
        rows.append(row)
new_file = [rows[0]]
for a in range(len(rows)):
    for b in range(len(all_otus)):
        if rows[a][0] == all_otus[b]:
            new_file.append(rows[a])
print len(new_file)

with open('new_taxa_and_seqs.csv', 'w') as f:
    writer= csv.writer(f)
    for b in range(len(new_file)):
        writer.writerow(new_file[b])
    
