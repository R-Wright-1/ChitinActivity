######
#Get sequences for tree
######
import csv
import os

file_path = '/Users/u1560915/Documents/GitHub/CommunityAnalysis/Fig4_FigS3/16S/'
os.chdir(file_path)
is16S = True #set this to true so that the new file will also contain the sequences of the isolates
limit = 0.5

with open('Grouped.csv', 'rU') as f:
    rows = []
    for row in csv.reader(f):
        rows.append(row)

adding, otus = [], []
new_grouped = False

for a in range(len(rows)):
    if a == 0:
        adding.append(rows[a])
    else:
        if rows[a][0][:3] == 'Otu':
            rows[a][0] = 'ASV'+rows[a][0][5:]
            new_grouped = True
        for b in range(len(rows[a])):
            if b > 0:
                rows[a][b] = float(rows[a][b])
        if max(rows[a][1:]) > limit:
            adding.append(rows[a])
            otus.append(rows[a][0])

with open('Sequences.fasta', 'rU') as f:
    sequences = []
    for row in f:
        sequences.append(row)

new_sequences = []
count = 0
for a in range(len(sequences)):
    if sequences[a][0] == '>':
        count += 1
        new_sequences.append('>'+str(count)+'\n')
    else:
        new_sequences.append(sequences[a])

count = 0
new_sequences = []
for a in range(len(sequences)):
    if sequences[a][0] == '>':
        ASV = int(sequences[a][4:])
        for b in range(len(otus)):
            otu = int(otus[b][3:])
            if ASV == otu:
                count += 1
                #new_sequences.append(sequences[a])
                new_sequences.append('>'+sequences[a][1:8]+'\n')
                if sequences[a+1][0] != '>':
                    new_sequences.append(sequences[a+1])

                
with open('Sequences'+str(limit)+'%.fasta', 'w') as f:
    for row in new_sequences:
        f.write(str(row))

with open('Taxonomy.csv', 'rU') as f:
    taxo = []
    for row in csv.reader(f):
        taxo.append(row)

new_tax = [taxo[0]]
for a in range(len(taxo)):
    if a > 0:
        taxo[a][0] = 'ASV'+taxo[a][0][4:]
        for b in range(len(otus)):
            if int(taxo[a][0][3:]) == int(otus[b][3:]):
                new_tax.append(taxo[a])

with open('New_grouped.csv', 'w') as f:
    writer = csv.writer(f)
    for a in range(len(new_tax)):
        writer.writerow(new_tax[a]+adding[a][1:])

####FOR 16S ONLY###
if is16S:
    with open('Isolates.fasta', 'rU') as f:
        rows = []
        for row in f:
            rows.append(row)
        
    with open('Sequences'+str(limit)+'% and isolates.fasta', 'w') as f:
        for row in new_sequences:
            f.write(str(row))
        for row in rows:
            f.write(str(row))

