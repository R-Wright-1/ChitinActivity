import csv

fn = 'chitin_conts_overall.csv'
tax_fn = 'Taxonomy.csv'

with open(fn, 'rU') as f:
    rows = []
    for row in csv.reader(f):
        rows.append(row)

with open(tax_fn, 'rU') as f:
    f_tax = []
    for row in csv.reader(f):
        f_tax.append(row)

asv, kegg, otu, count, abun, cont, tax = [], [], [], [], [], [], []
rest = []
col_names = ['GG OTU', 'KEGG', 'ASV', 'Gene copies', 'Abundance (%)', 'Contribution (%)', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']

for a in range(len(rows)):
    for b in range(len(rows[a])):
        if rows[a][b][:3] == 'ASV':
            asv.append(rows[a][b])
            kegg.append(rows[a][0])
            otu.append(rows[a][2])
            count.append(rows[a][3])
            abun.append(rows[a][4])
            cont.append(rows[a][6])
            rest.append([rows[a][2], rows[a][0], rows[a][b], rows[a][3], rows[a][4], rows[a][6]])
            for c in range(len(f_tax)):
                if rows[a][b] == f_tax[c][0]:
                    tax.append(f_tax[c][1:])

with open('Chitin_conts_overall_ASV.csv', 'w') as f:
    writer = csv.writer(f)
    writer.writerow(col_names)
    for a in range(len(rest)):
        writer.writerow(rest[a]+tax[a])

new = sorted(zip(asv, rest))
new2 = sorted(zip(asv, tax))

new_rest, new_tax = [], []
for a in range(len(new)):
    new_rest.append(new[a][1])
    new_tax.append(new2[a][1])

with open('Chitin_conts_overall_ASV_sorted.csv', 'w') as f:
    writer = csv.writer(f)
    writer.writerow(col_names)
    for a in range(len(new_rest)):
        writer.writerow(new_rest[a]+new_tax[a])

            
