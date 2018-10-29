import csv
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy
from colorsys import hls_to_rgb
import random
import matplotlib.patches as patches
import os
import math

order_18S = ["ASV000003", "ASV000047", "ASV000004", "ASV000001", "ASV000002",
             "ASV000013","ASV000011","ASV000039", "ASV000058", "ASV000675", "ASV000026",
             "ASV000040","ASV000050","ASV000177",
             "ASV000016", "ASV000006", "ASV000092","ASV000024","ASV000078","ASV000014","ASV000017",
             "ASV000049","ASV000045","ASV000038","ASV000282","ASV000028",
             "ASV000030","ASV000005","ASV000530","ASV000010","ASV000025","ASV000048","ASV000084","ASV000021","ASV000007",
             "ASV000015","ASV000191","ASV000033","ASV000221","ASV000098","ASV000070","ASV000133","ASV000184"]

file_path = '/Users/u1560915/Documents/GitHub/CommunityAnalysis/Fig4_FigS3/18S/'
os.chdir(file_path)
order = order_18S

with open('New_grouped.csv', 'rU') as f:
    rows = []
    for row in csv.reader(f):
        rows.append(row)
ordered_rows = [rows[0]]

count = 0
for a in range(len(order)):
    for b in range(len(rows)):
        if b == 0:
            continue
        if order[a][0:3] != 'ASV':
            if order[a] == rows[b][0]:
                ordered_rows.append(rows[b])
        elif rows[b][0][0:3] != 'ASV':
            if order[a] == rows[b][0]:
                ordered_rows.append(rows[b])
        elif int(order[a][3:]) == int(rows[b][0][3:]):
            ordered_rows.append(rows[b])
rows = ordered_rows

ASV, cls, gen_sp = [], [], []
colors = []
"""
for a in range(len(rows)):
    color = 'k'
    if a > 0:
        ASV.append(rows[a][0])
        if rows[a][7] != 'NA':
            if rows[a][0][0:3] != 'ASV':
                new_otu = 'Isolate: '+r'$'+rows[a][6]+'$'
                new_otu2 = r'$'+rows[a][7]+'$ '+rows[a][3]
                gen_sp.append(new_otu+' '+new_otu2)
            else:
                new_otu = 'ASV'+str(int(rows[a][0][3:]))+': '+r'$'+rows[a][6]+'$'
                new_otu2 = r'$'+rows[a][7]+'$ '+rows[a][3]
                gen_sp.append(new_otu+' '+new_otu2)
        elif rows[a][6] != 'NA':
            gen_sp.append('ASV'+str(int(rows[a][0][3:]))+': '+r'$'+rows[a][6]+'$ '+rows[a][3])
        elif rows[a][5] != 'NA':
            gen_sp.append('ASV'+str(int(rows[a][0][3:]))+': '+rows[a][5]+' '+rows[a][3])
        elif rows[a][4] != 'NA':
            gen_sp.append('ASV'+str(int(rows[a][0][3:]))+': '+rows[a][4]+' '+rows[a][3])
        elif rows[a][3] != 'NA':
            gen_sp.append('ASV'+str(int(rows[a][0][3:]))+': '+rows[a][3]+' '+rows[a][3])
        elif rows[a][2] != 'NA':
            gen_sp.append('ASV'+str(int(rows[a][0][3:]))+': '+rows[a][2]+' '+rows[a][3])
        elif rows[a][1] != 'NA':
            gen_sp.append('ASV'+str(int(rows[a][0][3:]))+': '+rows[a][1]+' '+rows[a][3])
    if rows[a][3] != 'NA':
        cls.append(rows[a][3])
    elif rows[a][2] != 'NA':
        cls.append(rows[a][2])
    elif rows[a][1] != 'NA':
        cls.append(rows[a][1])
    colors.append(color)
"""    
for a in range(len(rows)):
    color = 'k'
    if a > 0:
        ASV.append(rows[a][0])
        if rows[a][7] != 'NA':
            if rows[a][0][0:3] != 'ASV':
                new_otu = 'Isolate: '+r'$'+rows[a][6]+'$'
                new_otu2 = r'$'+rows[a][7]+'$ '
                gen_sp.append(new_otu+' '+new_otu2)
            else:
                new_otu = 'ASV'+str(int(rows[a][0][3:]))+': '+r'$'+rows[a][6]+'$'
                new_otu2 = r'$'+rows[a][7]+'$ '
                gen_sp.append(new_otu+' '+new_otu2)
        elif rows[a][6] != 'NA':
            gen_sp.append('ASV'+str(int(rows[a][0][3:]))+': '+r'$'+rows[a][6]+'$ ')
        elif rows[a][5] != 'NA':
            gen_sp.append('ASV'+str(int(rows[a][0][3:]))+': '+rows[a][5])
        elif rows[a][4] != 'NA':
            gen_sp.append('ASV'+str(int(rows[a][0][3:]))+': '+rows[a][4])
        elif rows[a][3] != 'NA':
            gen_sp.append('ASV'+str(int(rows[a][0][3:]))+': '+rows[a][3])
        elif rows[a][2] != 'NA':
            gen_sp.append('ASV'+str(int(rows[a][0][3:]))+': '+rows[a][2])
        elif rows[a][1] != 'NA':
            gen_sp.append('ASV'+str(int(rows[a][0][3:]))+': '+rows[a][1])
    if rows[a][3] != 'NA':
        cls.append(rows[a][3])
    elif rows[a][2] != 'NA':
        cls.append(rows[a][2])
    elif rows[a][1] != 'NA':
        cls.append(rows[a][1])
    colors.append(color)


def get_distinct_colors(n):
    colors = []
    for i in numpy.arange(0., 360., 360. / n):
        h = i / 360.
        l = (50 + numpy.random.rand() * 10) / 100.
        s = (90 + numpy.random.rand() * 10) / 100.
        colors.append(hls_to_rgb(h, l, s))
    random.shuffle(colors)
    return colors
   
unique = []
for d in range(len(cls)):
    adding = True
    for e in range(len(unique)):
        if cls[d] == unique[e]:
            adding = False
    if adding:
        unique.append(cls[d])
cols = get_distinct_colors(len(unique))
text_cols = []
for f in range(len(cls)):
    for g in range(len(unique)):
        if cls[f] == unique[g]:
            text_cols.append(cols[g])

new_rows, new_text_cols = [], []
for a in range(len(cls)):
    for b in range(len(unique)):
        if unique[b] == cls[a]:
            new_rows.append(rows[a])
            new_text_cols.append(text_cols[b])
del new_text_cols[0]

rows = new_rows
text_cols = new_text_cols

max_abun = []
new_rows, plotting = [], []
del rows[0]
for a in range(len(rows)):
    for b in range(len(rows[a])):
        if b == 0:
            continue
        if b > 7:
            rows[a][b] = float(rows[a][b])
for a in range(len(rows)):
    new_row = []
    this_plot = []
    mi, ma = min(rows[a][8:]), max(rows[a][8:])
    max_abun.append(ma)
    cmap = 'YlOrRd'
    norm = mpl.colors.Normalize(vmin=0, vmax=1)
    colormap = mpl.cm.get_cmap(cmap, 256)
    m = mpl.cm.ScalarMappable(norm=norm, cmap=colormap)
    for c in range(len(rows[a])):
        if c > 7:
            if rows[a][c] == 0:
                num = 0
            else:
                num = rows[a][c]/ma
            rows[a][c] = m.to_rgba(num)
            new_row.append(m.to_rgba(num))
            this_plot.append(1)
    new_rows.append(new_row)
    plotting.append(this_plot)

fig = plt.figure(figsize=(10, 25))
data = numpy.array(plotting)
bottom = numpy.cumsum(data, axis=0)
ax = plt.subplot2grid((43, 20), (0, 0), rowspan=43, colspan=5)
x = [1, 2, 3, 4]
y = [0.5]
x2 = [1.5, 2.5, 3.5, 4.5]
ax.bar(x, data[0], color=new_rows[0], width=1.0)
for j in xrange(1, data.shape[0]):
    y.append(j+0.5)
    ax.bar(x, data[j], color=new_rows[j], bottom=bottom[j-1], width=1.0)
del cls[0]
plt.sca(ax)

gen_sp[0] = r'ASV3: ($Developayella$ $elegans$)'
gen_sp[3] = r'ASV1: ($Cafeteria$ sp.)'
gen_sp[4] = r'ASV2: ($Cafeteria$ sp.)'
gen_sp[21] = r'ASV49: (Eustigmatophyceae)'
gen_sp[22] = r'ASV45: (Cercomonas)'
gen_sp[24] = r'ASV282: (Thaumatomonadida)'


plt.yticks(y, gen_sp)
plt.tick_params(axis='y', which='both', labelleft='off', labelright='on')
plt.tick_params(axis='x', which='both', bottom='on', top='off')
plt.xticks(x2, x, fontsize=16)
if order == order_16S:
    plt.ylim([0, 56])
else:
    plt.ylim([0, 43])
ax.set_xlabel('Days        ', fontsize=16)
for l,i in zip(ax.yaxis.get_ticklabels(),text_cols):
    #l.set_color(i)
    l.set_fontsize(16)
ax.set_xlim([1, 6])
for a in range(len(y)):
    ax.scatter(5.5, y[a], marker='o', s=int(15*max_abun[a]), color='k')
    ax.plot([5, 6], [a+1, a+1], 'k-')

ax2 = plt.subplot2grid((43, 20), (0, 5), rowspan=43, colspan=15, frameon=False)
plt.sca(ax2)
plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax2.get_yticklabels(), visible=False)
plt.tick_params(axis='y', which='both', left='off', right='off')
plt.tick_params(axis='x', which='both', bottom='off', top='off')
plt.ylim([0, 43])
plt.xlim([1, 6])

xc = 4.65
xt = 4.75
plt.plot([xc, xc], [42.9, 41.1], 'k-')
plt.text(xt, 41.65+0.7, 'Tremellomycetes', va='center', fontsize=20, weight='bold')
plt.text(xt, 41.65, '  1        2        3        4\n1.2%  0.99%  0%    0%', va='center', fontsize=14)
plt.plot([1.25, xc], [42.9, 42.9], 'k-')
plt.plot([1.25, xc], [41.1, 41.1], 'k-')

plt.plot([xc, xc], [40.9, 38.1], 'k-')
plt.text(xt, 39.35+0.7, 'Agaricomycetes', va='center', fontsize=20, weight='bold')
plt.text(xt, 39.35, '  1        2        3        4\n2.1% 4.25%   0%    0%', va='center', fontsize=14)
plt.plot([1.25, xc], [40.9, 40.9], 'k-')
plt.plot([1.25, xc], [38.1, 38.1], 'k-')

plt.plot([xc, xc], [37.9, 37.1], 'k-')
plt.text(xt, 37.5, 'Microbotryomycetes', va='center', fontsize=20, weight='bold')
plt.plot([1.25, xc], [37.9, 37.9], 'k-')
plt.plot([1.25, xc], [37.1, 37.1], 'k-')

plt.plot([xc, xc], [36.9, 33.1], 'k-')
plt.text(xt, 35+0.7, 'Malasseziomycetes', va='center', fontsize=20, weight='bold')
plt.text(xt, 35, '  1        2        3        4\n4.35%  9%    0%   0.15%', va='center', fontsize=14)
plt.plot([1.25, xc], [36.9, 36.9], 'k-')
plt.plot([1.25, xc], [33.1, 33.1], 'k-')

plt.plot([xc, xc], [32.9, 32.1], 'k-')
plt.text(xt, 32.5, 'Aphelidea', va='center', fontsize=20, weight='bold')
plt.plot([1.25, xc], [32.9, 32.9], 'k-')
plt.plot([1.25, xc], [32.1, 32.1], 'k-')

plt.plot([xc, xc], [31.9, 30.1], 'k-')
plt.text(xt, 31+0.7, 'Saccharomycetes', va='center', fontsize=20, weight='bold')
plt.text(xt, 31, '  1        2        3        4\n3.65% 0.05% 0%   0.05%', va='center', fontsize=14)
plt.plot([1.25, xc], [31.9, 31.9], 'k-')
plt.plot([1.25, xc], [30.1, 30.1], 'k-')

plt.plot([xc, xc], [29.9, 26.1], 'k-')
plt.text(xt, 28+0.7, 'Dothideomycetes', va='center', fontsize=20, weight='bold')
plt.text(xt, 28, '  1        2        3        4\n18%  3.35% 0.05%  0.1%', va='center', fontsize=14)
plt.plot([1.25, xc], [29.9, 29.9], 'k-')
plt.plot([1.25, xc], [26.1, 26.1], 'k-')

plt.plot([xc, xc], [25.9, 25.1], 'k-')
plt.text(xt, 25.5, 'Insecta', va='center', fontsize=20, weight='bold')
plt.plot([1.25, xc], [25.9, 25.9], 'k-')
plt.plot([1.25, xc], [25.1, 25.1], 'k-')

plt.plot([xc, xc], [24.9, 24.1], 'k-')
plt.text(xt, 24.5, '(Imbricatea)', va='center', fontsize=20, weight='bold')
plt.plot([1.25, xc], [24.9, 24.9], 'k-')
plt.plot([1.25, xc], [24.1, 24.1], 'k-')


plt.plot([xc, xc], [23.9, 23.1], 'k-')
plt.text(xt, 23.5, 'Glissomonadida', va='center', fontsize=20, weight='bold')
plt.plot([1.25, xc], [23.9, 23.9], 'k-')
plt.plot([1.25, xc], [23.1, 23.1], 'k-')

plt.plot([xc, xc], [22.9, 22.1], 'k-')
plt.text(xt, 22.5, '(Sarcomonadea)', va='center', fontsize=20, weight='bold')
plt.plot([1.25, xc], [22.9, 22.9], 'k-')
plt.plot([1.25, xc], [22.1, 22.1], 'k-')

plt.plot([xc, xc], [21.9, 21.1], 'k-')
plt.text(xt, 21.5, '(Eustigmatophyceae)', va='center', fontsize=20, weight='bold')
plt.plot([1.25, xc], [21.9, 21.9], 'k-')
plt.plot([1.25, xc], [21.1, 21.1], 'k-')

plt.plot([xc, xc], [20.9, 13.1], 'k-')
plt.text(xt, 17+0.7, 'Intramacaonucleata', va='center', fontsize=20, weight='bold')
plt.text(xt, 17, '  1        2        3        4\n35%    9%    0.55% 0.8%', va='center', fontsize=14)
plt.plot([1.25, xc], [20.9, 20.9], 'k-')
plt.plot([1.25, xc], [13.1, 13.1], 'k-')

plt.plot([xc, xc], [12.9, 12.1], 'k-')
plt.text(xt, 12.5, 'Perkinsidae', va='center', fontsize=20, weight='bold')
plt.plot([1.25, xc], [12.9, 12.9], 'k-')
plt.plot([1.25, xc], [12.1, 12.1], 'k-')

plt.plot([xc, xc], [11.9, 11.1], 'k-')
plt.text(xt, 11.5, 'Ulvophyceae', va='center', fontsize=20, weight='bold')
plt.plot([1.25, xc], [11.9, 11.9], 'k-')
plt.plot([1.25, xc], [11.1, 11.1], 'k-')

plt.plot([xc, xc], [10.9, 10.1], 'k-')
plt.text(xt, 10.5, 'Chlorophyceae', va='center', fontsize=20, weight='bold')
plt.plot([1.25, xc], [10.9, 10.9], 'k-')
plt.plot([1.25, xc], [10.1, 10.1], 'k-')

plt.plot([xc, xc], [9.9, 8.1], 'k-')
plt.text(xt, 9+0.7, 'Embryophyta', va='center', fontsize=20, weight='bold')
plt.text(xt, 9, '  1        2        3        4\n0.85% 0.7%   0%    0%', va='center', fontsize=14)
plt.plot([1.25, xc], [9.9, 9.9], 'k-')
plt.plot([1.25, xc], [8.1, 8.1], 'k-')

plt.plot([xc, xc], [7.9, 5.1], 'k-')
plt.text(xt, 6.5+0.7, 'Diatomea', va='center', fontsize=20, weight='bold')
plt.text(xt, 6.5, '  1        2        3        4\n6%      9%    0.05%  0.1%', va='center', fontsize=14)
plt.plot([1.25, xc], [7.9, 7.9], 'k-')
plt.plot([1.25, xc], [5.1, 5.1], 'k-')

plt.plot([xc, xc], [4.9, 3.1], 'k-')
plt.text(xt, 4+0.7, '(Bicosoecophyceae)', va='center', fontsize=20, weight='bold')
plt.text(xt, 4, '  1        2        3        4\n2.6%   33%   76%   89%', va='center', fontsize=14)
plt.plot([1.25, xc], [4.9, 4.9], 'k-')
plt.plot([1.25, xc], [3.1, 3.1], 'k-')

plt.plot([xc, xc], [2.9, 1.1], 'k-')
plt.text(xt, 2+0.7, 'Chrysophyceae', va='center', fontsize=20, weight='bold')
plt.text(xt, 2, '  1        2        3        4\n0.45% 11%   14%   3.4%', va='center', fontsize=14)
plt.plot([1.25, xc], [2.9, 2.9], 'k-')
plt.plot([1.25, xc], [1.1, 1.1], 'k-')

plt.plot([xc, xc], [0.9, 0.1], 'k-')
plt.text(xt, 0.5, '(Bigyromonadea)', va='center', fontsize=20, weight='bold')
plt.plot([1.25, xc], [0.9, 0.9], 'k-')
plt.plot([1.25, xc], [0.1, 0.1], 'k-')

plt.text(xt+1.25, -1, 'Totals', va='center', fontsize=20, ha='center', weight='bold')
plt.text(xt+1.25, -2, '1        2        3        4\n87%   93%   99%   99%\nDays', fontsize=14, va='center', ha='center')

#plt.show()
    
plt.savefig('Heatmap.png', bbox_inches='tight', dpi=300)
plt.close()

a = numpy.array([[0,100]])
plt.figure(figsize=(9, 1.5))
img = plt.imshow(a, cmap="YlOrRd")
plt.gca().set_visible(False)
cax = plt.axes([0.1, 0.2, 0.8, 0.6])
plt.colorbar(orientation='horizontal',cax=cax)
plt.xticks([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100], [0, '', '', '', '', '', '', '', '', '', 1], fontsize=22)
cax.set_xlabel('Proportion of maximum \nrelative abundance within OTU', fontsize=22)
plt.savefig("colorbar.png", bbox_inches='tight')
plt.close()

