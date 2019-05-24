import csv

fn = 'gamma_metagenome_contributions.csv'
otu = 'otu_table.csv'

with open(fn, 'rU') as f:
    rows = []
    for row in csv.reader(f):
        rows.append(row)

totu = []
for a in range(len(rows)):
    totu.append(rows[a][2])

with open(otu, 'rU') as f:
    otu_rows = []
    for row in csv.reader(f):
        otu_rows.append(row)

new_otus = totu

for a in range(len(otu_rows)):
    for b in range(len(totu)):
        if otu_rows[a][0] == totu[b]:
            this_one = [otu_rows[a][0]]
            for c in range(len(otu_rows[a])):
                if otu_rows[a][c] == '1':
                    this_one.append(otu_rows[1][c])
            if len(this_one) == 0:
                this_one = ['']
            new_otus[b] = this_one

with open('chitin_conts_gamma.csv', 'w') as f:
    writer = csv.writer(f)
    for a in range(len(rows)):
        if a > 0:
            writer.writerow(rows[a]+new_otus[a])
        else:
            writer.writerow(rows[a])