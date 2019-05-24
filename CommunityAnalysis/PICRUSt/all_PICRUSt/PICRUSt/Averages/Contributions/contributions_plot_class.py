import csv
import numpy
import matplotlib.pyplot as plt
from colorsys import hls_to_rgb
import random

fn = 'chitin_conts_four.csv'
tax_fn = 'Taxonomy.csv'

with open(fn, 'rU') as f:
    rows = []
    for row in csv.reader(f):
        rows.append(row)

with open(tax_fn, 'rU') as f:
    f_tax = []
    for row in csv.reader(f):
        f_tax.append(row)

ko, otu, count_cont, perc, ASV = [], [], [], [], []
for a in range(len(rows)):
    ko.append(rows[a][0])
    otu.append(rows[a][2])
    count_cont.append(rows[a][5])
    perc.append(rows[a][6])
    ASV.append(rows[a][16])

tax = ASV

for a in range(len(f_tax)):
    for b in range(len(tax)):
        if f_tax[a][0] == tax[b]:
            tax[b] = f_tax[a][3]

unique_ko = []
for a in range(len(ko)):
    adding = True
    for b in range(len(unique_ko)):
        if ko[a] == unique_ko[b]:
            adding = False
    if adding:
        unique_ko.append(ko[a])
del unique_ko[0]

unique_tax = []
for a in range(len(tax)):
    adding = True
    for b in range(len(unique_tax)):
        if tax[a] == unique_tax[b]:
            adding = False
    if adding:
        unique_tax.append(tax[a])
print unique_tax

new_otu, new_count_cont, new_perc, new_tax = [], [], [], []
for c in range(len(unique_ko)):
    this_otu, this_count_cont, this_perc, this_tax = [], [], [], []
    for d in range(len(ko)):
        if ko[d] == unique_ko[c]:
            this_otu.append(otu[d])
            this_count_cont.append(float(count_cont[d])/100)
            this_perc.append(perc[d])
            this_tax.append(tax[d])
    new_otu.append(this_otu)
    new_count_cont.append(this_count_cont)
    new_perc.append(this_perc)
    new_tax.append(this_tax)

def get_distinct_colors(n):
    colors = []
    for i in numpy.arange(0., 360., 360. / n):
        h = i / 360.
        l = (50 + numpy.random.rand() * 10) / 100.
        s = (90 + numpy.random.rand() * 10) / 100.
        colors.append(hls_to_rgb(h, l, s))
    random.shuffle(colors)
    return colors

"""
ax1 = plt.subplot(111)
plot = new_count_cont
for a in range(len(plot)):
    colors = get_distinct_colors(len(plot[a]))
    x = a+1.1
    ax1.bar([x], [plot[a][0]], color=[colors[0]], label=new_tax[a][0])
    for j in xrange(1, len(plot[a])):
        ax1.bar([x], [plot[a][j]], color=[colors[j]], bottom=[plot[a][j-1]], label=new_tax[a][j])
ax1.set_yscale('log')
ax1.set_ylabel('Gene copies')
plt.xticks([1.5, 2.5, 3.5, 4.5], unique_ko)
ax1.legend(bbox_to_anchor=(3, 1.03), ncol=3)
plt.savefig(fn[:-4]+'.png', bbox_inches='tight', dpi=600)
"""    
    




"""
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
col_names = ['GG OTU', 'KEGG', 'Sample', 'ASV', 'Gene copies', 'Abundance (%)', 'Contribution (%)', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']

for a in range(len(rows)):
    for b in range(len(rows[a])):
        if rows[a][b][:3] == 'ASV':
            asv.append(rows[a][b])
            kegg.append(rows[a][0])
            otu.append(rows[a][2])
            count.append(rows[a][3])
            abun.append(rows[a][4])
            cont.append(rows[a][6])
            rest.append([rows[a][2], rows[a][0], rows[a][1], rows[a][b], rows[a][3], rows[a][4], rows[a][6]])
            for c in range(len(f_tax)):
                if rows[a][b] == f_tax[c][0]:
                    tax.append(f_tax[c][1:])

with open('Chitin_conts_gamma_ASV.csv', 'w') as f:
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

with open('Chitin_conts_gamma_ASV_sorted.csv', 'w') as f:
    writer = csv.writer(f)
    writer.writerow(col_names)
    for a in range(len(new_rest)):
        writer.writerow(new_rest[a]+new_tax[a])
"""
            
