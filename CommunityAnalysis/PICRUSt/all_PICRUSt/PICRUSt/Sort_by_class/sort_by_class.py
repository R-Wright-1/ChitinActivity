import csv
import os

fn = 'taxa_and_seqs.csv'
pc = 4 #if 4, this will sort by class, if 3, by phyla
path = '/Users/u1560915/Documents/OneDrive/PhD_Plastic_Oceans/Experiments/MiSeq_Dada/PICRUSt/Sort_by_class/'
os.chdir(path)

with open(fn, 'rU') as f:
    rows = []
    for row in csv.reader(f):
        rows.append(row)

classes = []
numbers = []
for a in range(len(rows)):
    if a > 0:
        c = rows[a][pc]
        adding = True
        for b in range(len(classes)):
            if classes[b] == c:
                adding = False
                numbers[b] += 1
        if adding:
            classes.append(c)
            numbers.append(1)

other = []
other_seqs = []
for a in range(len(classes)):
    c = classes[a]
    new_rows = []
    sequences = []
    names = []
    for b in range(len(rows)):
        if b > 0:
            if rows[b][pc] == c:
                new_rows.append(rows[b])
                sequences.append(rows[b][1])
                names.append(rows[b][0])
    if len(new_rows) > 11:
        os.mkdir(path+c)
        os.chdir(path+c)
        with open('16S_shared.csv', 'w') as f:
            writer = csv.writer(f)
            writer.writerow([rows[0][0]]+rows[0][9:])
            for row in new_rows:
                writer.writerow([row[0]]+row[9:])
        with open('seqs.fna', 'w') as f:
            for s in range(len(sequences)):
                f.write('>'+names[s]+'\n')
                f.write(sequences[s]+'\n')
    else:
        for d in range(len(new_rows)):
            other.append(new_rows[d])
        for e in range(len(sequences)):
            name = '>'+names[e]+'\n'
            name += sequences[e]+'\n'
if len(other) > 0:
    os.mkdir(path+'Other')
    os.chdir(path+'Other')
    with open('16S_shared.csv', 'w') as f:
        writer = csv.writer(f)
        writer.writerow([rows[0][0]]+rows[0][9:])
        for row in other:
            writer.writerow([row[0]]+row[9:])
    with open('seqs.fna', 'w') as f:
        for f in range(len(other_seqs)):
            f.write(other_seqs[f])
        