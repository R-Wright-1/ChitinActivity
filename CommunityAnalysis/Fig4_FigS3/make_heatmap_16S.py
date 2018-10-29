import csv
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy
from colorsys import hls_to_rgb
import random
import matplotlib.patches as patches
import os
import math

order_16S = ['11_A_macleodii', 'ASV0037', '17_A_australica', '19_P_shioyasakaiensis', '18_V_tubiashii', 'ASV0024', 'ASV0021', 'ASV0003', 
               'ASV0025', '12_H_saccharevitans', '13_H_aestuarii', '09_H_campaniensis', 'ASV0042', 'ASV0040', 'ASV0004', 'ASV0035',
               'ASV0019', 'ASV0036', 'ASV0023', 'ASV0014', '03_P_gallaeciensis', '01_R_mobilis', '06_D_eberneus', 'ASV0001', 
               '04_R_calidilacus', '16_R_aestuari', '14_C_halophilus', '15_T_dalianensis', '10_N_aquibiodomus', 'ASV0076',
               'ASV0026', 'ASV0073', 'ASV0007', 'ASV0029', 'ASV0022', 'ASV0031', 'ASV0047', 'ASV0080', 'ASV0012', 'ASV0034', 'ASV0103',
               'ASV0125', '20_G_terrae', '08_L_aquatica', 'ASV0030', 'ASV0011', '07_J_marina', 'ASV0005', '02_M_ruestringensis', '05_M_antarctica', 'ASV0018', 'ASV0057', 
               'ASV0039', 'ASV0131', 'ASV0146', 'ASV0002']

file_path = '/Users/u1560915/Documents/GitHub/CommunityAnalysis/Fig4_FigS3/16S/'
os.chdir(file_path)
order = order_16S
order.reverse()

with open('New_grouped_isolates.csv', 'rU') as f:
    rows = []
    for row in csv.reader(f):
        rows.append(row)
ordered_rows = [rows[0]]

re, ora, gr = 'r', '#FFA500', '#32CD32'

chitin = [ora, ora, ora, gr, re, re, ora, ora, gr, gr, gr, re, re, gr, gr, gr, gr, gr, gr, gr]
chitin = [re, re, gr, gr, gr, ora, ora, ora, gr, re, ora, ora, re, gr, gr, gr, gr, gr, gr, gr]

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

with open('New_grouped_order.csv', 'w') as f:
    writer = csv.writer(f)
    for a in range(len(rows)):
        writer.writerow(rows[a])

ASV, cls, gen_sp = [], [], []
colors = []
count = 0

for a in range(len(rows)):
    color = 'k'
    if a > 0:
        sim = ''
        ASV.append(rows[a][0])
        add_cls = ''#' '+rows[a][3]
        if rows[a][7] != 'NA':
            if rows[a][0][0:3] != 'ASV':
                new_otu = 'Isolate: '+r'$'+rows[a][6]+'$'
                new_otu2 = r'$'+rows[a][7]+'$'+add_cls
                gen_sp.append(new_otu+' '+new_otu2)
                color = chitin[count]
                count += 1
            else:
                new_otu = 'ASV'+str(int(rows[a][0][3:]))+': '+r'$'+rows[a][6]+'$'
                new_otu2 = r'$'+rows[a][7]+'$'+add_cls
                gen_sp.append(new_otu+' '+new_otu2+sim)
        elif rows[a][6] != 'NA':
            gen_sp.append('ASV'+str(int(rows[a][0][3:]))+': '+r'$'+rows[a][6]+'$'+sim+add_cls)
        elif rows[a][5] != 'NA':
            gen_sp.append('ASV'+str(int(rows[a][0][3:]))+': '+rows[a][5]+sim+add_cls)
        elif rows[a][4] != 'NA':
            gen_sp.append('ASV'+str(int(rows[a][0][3:]))+': '+rows[a][4]+sim+add_cls)
        elif rows[a][3] != 'NA':
            gen_sp.append('ASV'+str(int(rows[a][0][3:]))+': '+rows[a][3]+sim+add_cls)
        elif rows[a][2] != 'NA':
            gen_sp.append('ASV'+str(int(rows[a][0][3:]))+': '+rows[a][2]+sim+add_cls)
        elif rows[a][1] != 'NA':
            gen_sp.append('ASV'+str(int(rows[a][0][3:]))+': '+rows[a][1]+sim+add_cls)
    if rows[a][3] != 'NA':
        cls.append(rows[a][3])
    elif rows[a][2] != 'NA':
        cls.append(rows[a][2])
    elif rows[a][1] != 'NA':
        cls.append(rows[a][1])
    colors.append(color)
is_colors = colors
cls_cols = colors

simper = ["ASV0002", "ASV0003", "ASV0004", "ASV0005", "ASV0007"]

for a in range(len(gen_sp)):
    if gen_sp[a][:15] == 'Isolate: $Phaeo':
        gen_sp[a] = r'Isolate: $Phaeobacter$ $gallaeciensis$'
    elif gen_sp[a][:6] == 'ASV125':
        gen_sp[a] = r'ASV125: ($Thermobrachium$)' #85%, Clostridia
    elif gen_sp[a][:5] == 'ASV2:':
        gen_sp[a] = r'ASV2: ($Spirochaeta$)*' #82%, Spirochaeta
    elif gen_sp[a][:5] == 'ASV3:':
        gen_sp[a] = gen_sp[a]+'*'
    elif gen_sp[a][:5] == 'ASV4:':
        gen_sp[a] = gen_sp[a]+'*'
    elif gen_sp[a][:5] == 'ASV5:':
        gen_sp[a] = gen_sp[a]+'*'
    elif gen_sp[a][:5] == 'ASV7:':
        gen_sp[a] = gen_sp[a]+'*'


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

fs = 22
x = [1, 2, 3, 4]
fig = plt.figure(figsize=(20, 30))
plt.subplots_adjust(wspace=0, hspace=0)
for a in range(len(plotting)):
    if gen_sp[a][:3] == 'ASV':
        ax = plt.subplot2grid((56, 20), (56-(a+1), 0), colspan=3)
        ax.bar(x, plotting[a], color=new_rows[a], width=1.0)
        ax.scatter(5.5, 0.5, marker='o', s=int(15*max_abun[a]), color='k')
        ax.set_xlim([1, 6])
        ax.set_ylim([0, 1])
        plt.setp(ax.get_yticklabels(), visible=False)
        ax.tick_params(axis='y',which='both',right='off', left='off')
        if a == 0:
            x1 = [1.5, 2.5, 3.5, 4.5]
            plt.xticks(x1, x, fontsize=fs)
            ax.set_xlabel('Days      ', fontsize=fs)
            ax.tick_params(axis='x',which='both',bottom='on', top='off')
        else:
            plt.setp(ax.get_xticklabels(), visible=False)
            ax.tick_params(axis='x',which='both',bottom='off', top='off')
        

xc = 3.5
xt = 3.6
xc1 = 0.5
ax = plt.subplot2grid((56, 20), (0, 3), colspan=15, rowspan=56, frameon=False)
ax.tick_params(axis='x',which='both',bottom='off', top='off')
ax.tick_params(axis='y',which='both',right='off', left='off')
plt.setp(ax.get_xticklabels(), visible=False)
plt.setp(ax.get_yticklabels(), visible=False)
ax.set_ylim([0, 56])
ax.set_xlim([0, 6])

ax.plot([xc1, xc], [55.9, 55.9], 'k-')
ax.plot([xc, xc], [55.9, 37.1], 'k-')
ax.text(xt, 46.5+0.5, 'Gammaproteobacteria', fontsize=fs, weight='bold')
ax.text(xt, 46.5, '\n  1        2        3        4\n73%   73%   46%   21%', va='center', fontsize=fs)
ax.plot([xc1, xc], [37.1, 37.1], 'k-')

ax.plot([xc1, xc], [36.9, 36.9], 'k-')
ax.plot([xc, xc], [36.9, 22.1], 'k-')
ax.text(xt, 29.5+0.5, 'Alphaproteobacteria', fontsize=fs, weight='bold')
ax.text(xt, 29.5, '\n  1        2        3        4\n12%   13%    13%   21%', va='center', fontsize=fs)
ax.plot([xc1, xc], [22.1, 22.1], 'k-')

ax.plot([xc1, xc], [21.9, 21.9], 'k-')
ax.plot([xc, xc], [21.9, 19.1], 'k-')
ax.text(xt, 20.5+0.5, 'Deltaproteobacteria', fontsize=fs, weight='bold')
ax.text(xt, 20.5, '\n  1        2        3        4\n0.3%  1.5%  4.5%   3.5%', va='center', fontsize=fs)
ax.plot([xc1, xc], [19.1, 19.1], 'k-')

ax.plot([xc1, xc], [18.9, 18.9], 'k-')
ax.plot([xc, xc], [18.9, 14.1], 'k-')
ax.text(xt, 16.5+0.5, 'Clostridia', fontsize=fs, weight='bold')
ax.text(xt, 16.5, '\n  1        2        3        4\n0.1%  0.3%   2%    13.5%', va='center', fontsize=fs)
ax.plot([xc1, xc], [14.1, 14.1], 'k-')

ax.plot([xc1, xc], [13.9, 13.9], 'k-')
ax.plot([xc, xc], [13.9, 12.1], 'k-')
ax.text(xt, 13, 'Actinobacteria', va='center', fontsize=fs, weight='bold')
ax.plot([xc1, xc], [12.1, 12.1], 'k-')

ax.plot([xc1, xc], [11.9, 11.9], 'k-')
ax.plot([xc, xc], [11.9, 4.1], 'k-')
ax.text(xt, 8+0.5, 'Bacteroidia', fontsize=fs, weight='bold')
ax.text(xt, 8, '\n  1        2        3        4\n3.8%   6%     19%    22%', va='center', fontsize=fs)
ax.plot([xc1, xc], [4.1, 4.1], 'k-')

ax.plot([xc1, xc], [3.9, 3.9], 'k-')
ax.plot([xc, xc], [3.9, 2.1], 'k-')
ax.text(xt, 3+0.5, 'Chitinivibrionia', fontsize=fs, weight='bold')
ax.text(xt, 3, '\n  1        2        3        4\n0.1%  0.6%   2.9%  0.4%', va='center', fontsize=fs)
ax.plot([xc1, xc], [2.1, 2.1], 'k-')

ax.plot([xc1, xc], [1.9, 1.9], 'k-')
ax.plot([xc, xc], [1.9, 1.1], 'k-')
ax.text(xt, 1.5, 'Cloacimonadia', va='center', fontsize=fs, weight='bold')
ax.plot([xc1, xc], [1.1, 1.1], 'k-')

ax.plot([xc1, xc], [0.9, 0.9], 'k-')
ax.plot([xc, xc], [0.9, 0.1], 'k-')
ax.text(xt, 0.5, 'Spirochaetes', va='center', fontsize=fs, weight='bold')
ax.plot([xc1, xc], [0.1, 0.1], 'k-')

ax.text(xt+xc1*2, -1, 'Totals', va='center', ha='center', fontsize=fs, weight='bold')
ax.text(xt+xc1*2, -3, '1        2        3        4\n89%   95%    93%   91%', fontsize=fs, ha='center')
ax.text(xc, -3, 'Days\nRelative abundance', ha='right', fontsize=fs+4, weight='bold')
#ax.text(xt+xc1*2, -3, '\n\n\nDays', ha='center', fontsize=fs+4, weight='bold')

n = 0.01384
n1 = 0.1295

for a in range(len(gen_sp)):
    fig.text(0.25, n1+n*a, gen_sp[a], color=is_colors[a+1], fontsize=fs)
    #fig.text(0.95, n1+n*a, gen_sp[a], color=text_cols[a], fontsize=16)

#plt.show()    
plt.savefig('Heatmap 16S.png', bbox_inches='tight', dpi=600)
plt.close()
"""
fs=22

a = numpy.array([[0,100]])
plt.figure(figsize=(5, 1.5))
img = plt.imshow(a, cmap="YlOrRd")
plt.gca().set_visible(False)
cax = plt.axes([0.1, 0.2, 0.8, 0.6])
plt.colorbar(orientation='horizontal',cax=cax)
plt.xticks([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100], [0, '', '', '', '', '', '', '', '', '', 1], fontsize=22)
cax.set_xlabel('Proportion of maximum \nrelative abundance within ASV', fontsize=fs)
plt.savefig("colorbar.png", bbox_inches='tight', dpi=600)
plt.close()
"""