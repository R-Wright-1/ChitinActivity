import csv
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy
from colorsys import hls_to_rgb
import random
import matplotlib.patches as patches

otu = [540, 43, 132, 767, 218, 351, 274, 140, 56, 15, 325, 124, 177, 515, 52, 10, 55, 208, 592, 19, 226, 617, 78, 40, 54, 328, 105, 163, 1700, 255, 470, 188, 9, 114, 45, 72, 62, 447, 99, 91, 338, 22, 202, 46, 316, 279, 311, 13, 178, 173, 721, 160, 93, 148, 190, 119, 262, 133, 11, 7, 94, 654, 268, 254, 69, 6, 329, 64, 95, 92, 246, 276, 294, 621, 366, 3, 304, 154, 141, 67, 18, 2, 37, 1, 199, 81, 134, 138, 8, 252, 5]
"""
for a in range(len(otu)):
    otu[a] = int(otu[a])
otu = sorted(otu)
print otu

"""
for a in range(len(otu)):
    otu[a] = str(otu[a])
    while len(otu[a]) < 5:
        otu[a] = '0'+otu[a]
    otu[a] = 'Otu'+otu[a]

with open ('3_daily_percent_grouped.csv', 'rU') as f:
    rows = []
    for row in csv.reader(f):
        rows.append(row)

new_rows = []

for o in otu:
    for row in rows:
        if o == row[0]:
            new_rows.append(row[1:])

taxonomy = []
rest_of_taxonomy = []
with open('18S_taxonomy.csv', 'rU') as f:
    rows = []
    for row in csv.reader(f):
        rows.append(row)
for o in otu:
    adding = False
    for row in rows:
        if o == row[0]:
            taxonomy.append('OTU'+o[-3:]+' '+row[6])
            rest_of_taxonomy.append(row[3])
            #'OTU'+o[-3:]+' '+row[6]
            adding = True
    if adding == False:
        print o

list.reverse(taxonomy)
list.reverse(new_rows)
list.reverse(otu)
list.reverse(rest_of_taxonomy)

all_class = []

for a in rest_of_taxonomy:
    adding = True
    for b in all_class:
        if a == b:
            adding = False
    if adding == True:
        all_class.append(a)

def get_distinct_colors(n):
    colors = []
    for i in numpy.arange(0., 360., 360. / n):
        h = i / 360.
        l = (50 + numpy.random.rand() * 10) / 100.
        s = (90 + numpy.random.rand() * 10) / 100.
        colors.append(hls_to_rgb(h, l, s))
    random.shuffle(colors)
    return colors

cols = get_distinct_colors(len(all_class))  
fontcols = []
for a in rest_of_taxonomy:
    for b in range(len(all_class)):
        if a == all_class[b]:
            col = cols[b]
            if all_class[b] == 'Eukaryota_unclassified':
                col = 'k'
            fontcols.append(col)
print len(fontcols)


for a in range(len(new_rows)):
    for b in range(len(new_rows[a])):
        if new_rows[a][b] == '':
            new_rows[a][b] = 0.0
        else:
            new_rows[a][b] = float(new_rows[a][b])

cmap = 'YlOrRd_r'
norm = mpl.colors.Normalize(vmin=0, vmax=1)
colormap = mpl.cm.get_cmap(cmap, 256)
m = mpl.cm.ScalarMappable(norm=norm, cmap=colormap)
colors = []
x = []
y = []
for a in range(len(new_rows)):
    y.append(a+0.5)
    maximum = max(new_rows[a])
    this_cols = []
    for b in range(len(new_rows[a])):
        if a == 0:
            x.append(b)
        if new_rows[a][b] != 0:
            new_rows[a][b] = new_rows[a][b]/maximum
            if maximum == 0:
                new_rows[a][b] = maximum
        this_cols.append(m.to_rgba(new_rows[a][b]))
        new_rows[a][b] = 1
    colors.append(this_cols)

    
fig = plt.figure(figsize=(1.5, 30)) 
data = numpy.array(new_rows)
bottom = numpy.cumsum(data, axis=0)
ax = plt.subplot(111)
ax.bar(x, data[0], color=colors[0], width=1.0)
for j in xrange(1, data.shape[0]):
    ax.bar(x, data[j], color=colors[j], bottom=bottom[j-1], width=1.0)
plt.xlim([0,len(x)])
plt.ylim([0,len(new_rows)])
plt.yticks(y, taxonomy)
plt.tick_params(axis='y', which='both', labelleft='off', labelright='on')
"""
for l,i in zip(ax.yaxis.get_ticklabels(),fontcols):
    l.set_color(i)
    l.set_fontsize(11)
"""
plt.xticks([])
plt.savefig('18S heatmap OTU daily.pdf', bbox_inches='tight')  
plt.close()

all_classes = ["Tremellomycetes", "Agaricomycetes", "Exobasidiomycetes", "Dothideomycetes", "Saccharomycetes", "Incertae Sedis", "Mammalia", "Intramacronucleata", "Ulvophyceae", "Embryophyta", "Diatomea", "Chrysophyceae"]
r_cols = ['#%02x%02x%02x' % (255, 127, 0), '#%02x%02x%02x' % (152, 245, 255), '#%02x%02x%02x' % (138, 43, 226), '#%02x%02x%02x' % (255, 215, 0), '#%02x%02x%02x' % (0, 0, 255), '#%02x%02x%02x' % (0, 205, 0), '#%02x%02x%02x' % (255, 0, 0), '#%02x%02x%02x' % (139, 136, 120), '#%02x%02x%02x' % (202, 255, 112), '#%02x%02x%02x' % (139, 0, 139), '#%02x%02x%02x' % (139, 0, 0)]

new_list = sorted(zip(all_class, cols))
#new_list = zip(all_class, cols)
new_list = sorted(zip(all_classes, r_cols))

pos = [0, 20]
ax1 = plt.subplot(111)
for i in new_list:
    ax1.add_patch(patches.Rectangle((pos[0], pos[1]), 0.05, 0.5, facecolor=i[1], alpha=0.6))
    ax1.text(pos[0]+0.1, pos[1], i[0])
    pos[1] -= 1
plt.ylim([0, 21])
plt.axis('off')
plt.savefig('List of colors.pdf')