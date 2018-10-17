import csv
import numpy
from scipy import stats
import pandas
import matplotlib.pyplot as plt
import classify_file as cf
from colorsys import hls_to_rgb
import random
import pandas as pd
from sklearn.metrics import jaccard_similarity_score
import matplotlib
import matplotlib.patches as mpatches
import os

#run this after analysis scripts 1 and 2!
#names, indices, all_samples, all_treatments = analysis_script_1.get_meta_treatments('16S_metadata.csv') 
files = ['0_percent_grouped.csv', '1_percent_grouped.csv', '2_daily_r_percent_grouped.csv', '3_daily_percent_grouped.csv', '4_LG_percent_grouped.csv', '5_RG_percent_grouped.csv', '6_beg_end_percent_grouped.csv']
simper_means_files = ['', '', '2_daily_r_percent_simper_mean.csv', '3_daily_percent_simper_mean.csv', '4_LG_percent_simper_mean.csv', '5_RG_percent_simper_mean.csv', '']
taxonomy = 'Taxonomy.csv'

def get_OTU_taxonomy(level, OTU):
    #Genus = 6, Family = 5, Order = 4, Class = 3, Phylum = 2
    with open(taxonomy, 'rU') as f:
        reader = csv.reader(f)
        rows = []
        for row in reader:
            rows.append(row)
    for a in range(len(rows)):
        if rows[a][0] == OTU:
            returning = rows[a]
    this_level = returning[level]
    return this_level
    
def get_distinct_colors(n):
    colors = []
    for i in numpy.arange(0., 360., 360. / n):
        h = i / 360.
        l = (50 + numpy.random.rand() * 10) / 100.
        s = (90 + numpy.random.rand() * 10) / 100.
        colors.append(hls_to_rgb(h, l, s))
    random.shuffle(colors)
    return colors
    
plt1otu, plt1tax = ['OTU1', 'OTU2', 'OTU5', 'OTU8', 'OTU13'], ['Eukaryota', 'Eukaryota', r'$Developayella$', r'$Paraphysomonas$', 'Intramacronucleata']
plt2otu, plt2tax = ['OTU1', 'OTU5', 'OTU2', 'OTU8', 'OTU13'], ['Eukaryota', r'$Developayella$', 'Eukaryota', r'$Paraphysomonas$', 'Intramacronucleata']

def plot_simper(sf, axis, colors, labels, plc, text, ylim, ot, ta):
    with open(sf, 'rU') as f:
        reader = csv.reader(f)
        rows = []
        for row in reader:
            rows.append(row)
    a, b = 0, 1
    lines = []
    while a < 5:
        lines.append(rows[b])
        a += 1
        b += 1
    OTUs = []
    for a in range(len(lines)):
        for b in range(len(lines[a])):
            if b == 0:
                OTUs.append(lines[a][b])
                lines[a][b] = get_OTU_taxonomy(6, lines[a][b])
            else:
                lines[a][b] = float(lines[a][b])
    #nums = ['', 1, '', 2, '', 3, '', 4]
    means, stds = [[], [], [], [], [], [], [], [], [], []], [[], [], [], [], [], [], [], [], [], []]
    p_vals, names = [], []
    for c in range(2):
        c = (c*2)+1
        e = c+1
        for d in range(5):
            means[d].append(lines[d][c]*100)
            stds[d].append(lines[d][e]*100)
            p_vals.append(lines[d][-1])
            names.append(lines[d][0])
    for i in range(5):
        number = 0
        new_OTU = 'OTU'
        for a in range(len(OTUs[i])):
            a = OTUs[i][a]
            if a != 'O' and a!= 't' and a!= 'u':
                number += float(a)
            if number != 0:
                new_OTU += a
        OTUs[i] = new_OTU
        name = names[i]
        new_name = ''
        count = 0
        new_line = False
        for j in name:
            count += 1
            if j == '_':
                new_line = True
                new_name = '${'+new_name+'}$'
                new_name += '\n'
            else:
                new_name += j
            if count == len(name) and not new_line:
                new_name = r'${'+new_name+'}$'
        new_name = OTUs[i]+'\n'+new_name        
        names[i] = new_name
    nums2 = [0.5, 1.5]
    #ax3.set_xlabel('Day')
    for e in range(5):
        axis[e].bar(nums2, means[e], yerr=stds[e], color=colors, error_kw=dict(ecolor='gray', lw=1, capsize=3, capthick=1))
        axis[e].text(1.4, text, '${'+r'p = '+'}$'+'%.3f'%float(p_vals[e]), va='bottom', ha='left', color='#bd0303', fontsize=8)
        #axis[e].set_xticklabels(nums)
        axis[e].set_title(ot[e]+'\n'+ta[e], fontsize=10)
        axis[e].set_ylim([0, ylim])
        if e == 9:
            axis[e].bar(nums2, means[e], yerr=stds[e], color=colors, error_kw=dict(ecolor='gray', lw=1, capsize=3, capthick=1))
        axis[e].set_xlim([0.25, 2.5])
    patch0 = mpatches.Patch(color=colors[0], label=labels[0])
    patch1 = mpatches.Patch(color=colors[1], label=labels[1])
    axis[4].legend(handles=[patch0, patch1], bbox_to_anchor=(plc,1.05), fontsize=10)
    #ax1.set_yticks([0, 10, 20, 30, 40])
    axis[0].tick_params(axis='y', which='both', left='on', right='off', labelbottom='off')
    axis[0].tick_params(axis='x', which='both', bottom='off', top='off')
    for a in axis:
        plt.setp(a.get_xticklabels(), visible=False)
        a.tick_params(axis='x', which='both', bottom='off', top='off')
    axis2 = [axis[1], axis[2], axis[3], axis[4]]
    for ax in axis2:
        ax.tick_params(axis='y', which='both', left='off', right='off')
        plt.setp(ax.get_yticklabels(), visible=False)
    #axes = [ax6, ax7, ax8, ax9, ax10]
    #for a in axes:
    #    plt.setp(a.get_xticklabels(), visible=False)
    #    a.tick_params(axis='x', which='both', bottom='off', top='off')
    plt.tight_layout()
    #ax6.set_ylabel('Relative abundance (%)')
    
    return
fig = plt.figure(figsize=(8.27, 5))    
colors = get_distinct_colors(15)
labels1, labels2 = labels = ['Beginning', 'End'], ['Selection \n (9 day)', 'Selection \n (4 day)']
colors1, colors2 = [colors[0], colors[1]], [colors[2], colors[3]]
ax1, ax6 = plt.subplot2grid((2,5), (0,0)), plt.subplot2grid((2,5), (1,0))
ax2, ax3, ax4, ax5, ax7, ax8, ax9, ax10 = plt.subplot2grid((2,5), (0,1), sharey=ax1), plt.subplot2grid((2,5), (0,2), sharey=ax1), plt.subplot2grid((2,5), (0,3), sharey=ax1), plt.subplot2grid((2,5), (0,4), sharey=ax1), plt.subplot2grid((2,5), (1,1), sharey=ax6), plt.subplot2grid((2,5), (1,2), sharey=ax6), plt.subplot2grid((2,5), (1,3), sharey=ax6), plt.subplot2grid((2,5), (1,4), sharey=ax6)
axis1 = [ax1, ax2, ax3, ax4, ax5]
axis2 = [ax6, ax7, ax8, ax9, ax10]
ax1.set_ylabel('Simper (a) \n Relative abundance (%)')
ax6.set_ylabel('Simper (b) \n Relative abundance (%)')
plot_simper('11_beg_end_S_percent_simper_mean.csv', axis1, colors1, labels1, plc=2.15, text=90, ylim=100, ot=plt1otu, ta=plt1tax)
plot_simper('12_LS_percent_simper_mean.csv', axis2, colors2, labels2, plc=2.15, text=90, ylim=100, ot=plt2otu, ta=plt2tax)
os.chdir('/Users/u1560915/Documents/GitHub/ChitinActivity/MiSeq/18S/Figures/')
plt.savefig('Fig supp 14.pdf', bbox_inches='tight')

