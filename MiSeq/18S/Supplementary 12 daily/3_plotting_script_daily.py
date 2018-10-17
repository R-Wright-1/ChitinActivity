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
import os

#run this after analysis scripts 1 and 2!
#names, indices, all_samples, all_treatments = analysis_script_1.get_meta_treatments('16S_metadata.csv') 
files = ['0_percent_grouped.csv', '1_percent_grouped.csv', '2_daily_r_percent_grouped.csv', '3_daily_percent_grouped.csv', '4_LG_percent_grouped.csv', '5_RG_percent_grouped.csv', '6_beg_end_percent_grouped.csv']
simper_means_files = ['', '', '2_daily_r_percent_simper_mean.csv', '3_daily_percent_simper_mean.csv', '4_LG_percent_simper_mean.csv', '5_RG_percent_simper_mean.csv', '']
taxonomy = 'Taxonomy.csv'

def get_OTU_taxonomy(level, OTU):
    os.chdir('/Users/u1560915/Documents/GitHub/ChitinActivity/MiSeq/18S/Supplementary 12 daily/')
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
   
def stacked_barchart(fn, level=6):
    os.chdir('/Users/u1560915/Documents/GitHub/ChitinActivity/MiSeq/18S/Supplementary 12 daily/')
    new_fn = cf.group_file(fn, True)
    new_fn = cf.classify_file(new_fn, level)
    count = 0
    fig = plt.figure(figsize=(8.27, 8.27))
    ax7 = plt.subplot2grid((3,10), (0,0), colspan=5, rowspan=2)
    OTUs = ['OTU1', 'OTU2', 'OTU3', 'OTU5', 'OTU8', 'OTU9', 'OTU10', 'OTU13', 'OTU15', 'OTU18', 'OTU22', 'OTU45', 'OTU46', 'OTU55', '']
    with open(new_fn, 'rU') as f:
        reader = csv.reader(f)
        rows = []
        for row in reader:
            rows.append(row)
    OTU, numbers = [], []
    for a in range(len(rows)):
        if a > 0:
            count += 1
            this_row = []
            for b in range(len(rows[a])):
                if b == 0:
                    OTU.append(rows[a][b])
                else:
                    this_row.append(float(rows[a][b]))
            numbers.append(this_row)
    print len(OTUs)
    print len(OTU)
    new_OTU = []
    for o in range(len(OTU)):
        string = ''
        new_line = False
        count1 = 0
        for p in range(len(OTU[o])):
            count1 += 1
            if OTU[o][p] == '_':
                new_line = True
                string = OTUs[o]+'\n'+'${'+string+'}$'
                string+=' \n '
            else:
                string+=OTU[o][p]
            if count1 == len(OTU[o]) and not new_line:
                string = OTUs[o]+'\n'+r'${'+string+'}$'
        new_OTU.append(string)
    OTU = new_OTU
    OTU[-1] = 'Other'
    del numbers[-1][-1]
    data = numpy.array(numbers)
    bottom = numpy.cumsum(data, axis=0)
    nums = [1, 2, 3, 4]
    colors = get_distinct_colors(count)
    OTUs = ['OTU1', 'OTU2', 'OTU3', 'OTU5', 'OTU8', 'OTU9', 'OTU10', 'OTU13', 'OTU15', 'OTU18', 'OTU22', 'OTU45', 'OTU46', 'OTU55', '']
    nline = '\n'
    tax = ['Eukaryota', 'Eukaryota', r'$Amphora$', r'$Developayella$', r'$Paraphysomonas$', 'Mammalia', 'Capnodiales', 'Intramacronucleata', r'$Malassezia$', r'$Navicula$', 'Choreotrichia', 'Thysanoptera', 'Choreotrichia', 'Dothideales', 'Other']
    ax7.bar(nums, data[0], color=colors[0], label=OTUs[0]+nline+tax[0])
    for j in xrange(1, data.shape[0]):
        ax7.bar(nums, data[j], color=colors[j], bottom=bottom[j-1], label=OTUs[j]+nline+tax[j])
    ax7.tick_params(axis='x', which='both', bottom='off', top='off')
    ax7.tick_params(axis='y', which='both', left='on', right='off')
    ax7.legend(loc='upper left', bbox_to_anchor=(1.085, 1.015), ncol=2, fontsize=10)
    ax7.set_ylabel('Relative abundance (%)', fontsize=12)
    labels = ['Day 1', 'Day 2', 'Day 3', 'Day 4']
    nums = [1.4, 2.4, 3.4, 4.4]
    plt.setp(ax7, xticks=nums, xticklabels=labels)
    ax7.set_ylim([0, 100])
    ax7.set_xlim([0.8, 5])
    #ax7.set_title('A', loc='left')
    #ax1.savefig('Daily_stacked_bar_'+str(level)+'.pdf', bbox_inches='tight')
    return
stacked_barchart('3_daily_percent_grouped.csv', 6)

def plot_simper(sf):
    os.chdir('/Users/u1560915/Documents/GitHub/ChitinActivity/MiSeq/18S/Supplementary 12 daily/')
    with open(sf, 'rU') as f:
        reader = csv.reader(f)
        rows = []
        for row in reader:
            rows.append(row)
    a, b = 0, 1
    lines = []
    while a < 5:
        """
        p_val = float(rows[b][-1])
        if p_val <= 0.05:
            lines.append(rows[b])
            a += 1
        """
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
    ax1 = plt.subplot2grid((3,10), (2,0), colspan=2)
    ax2, ax3, ax4, ax5 = plt.subplot2grid((3,10), (2,2), sharey=ax1, colspan=2), plt.subplot2grid((3,10), (2,4), sharey=ax1, colspan=2), plt.subplot2grid((3,10), (2,6), sharey=ax1, colspan=2), plt.subplot2grid((3,10), (2,8), sharey=ax1, colspan=2)
    axis = [ax1, ax2, ax3, ax4, ax5]
    colors = ['#33FFFF', '#33CCFF', '#3399FF', '#3366FF']
    nums = ['', 1, '', 2, '', 3, '', 4]
    means, stds = [[], [], [], [], []], [[], [], [], [], []]
    p_vals, names = [], []
    for c in range(4):
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
        print name
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
    nums2 = [0.5, 1.5, 2.5, 3.5]
    nums3 = [1, 2, 3, 4]
    nums = ['1', '2', '3', '4']
    ax3.set_xlabel('Day')
    OTUs = ['OTU2', 'OTU55', 'OTU13', 'OTU8', 'OTU3']
    nline = '\n'
    tax = ['Eukaryota', 'Dothideales', 'Intramacronucleata', r'$Paraphysomonas$', r'$Amphora$']
    for e in range(5):
        axis[e].bar(nums2, means[e], yerr=stds[e], color=colors, error_kw=dict(ecolor='gray', lw=1, capsize=3, capthick=1))
        plt.setp(axis[e], xticks=nums3, xticklabels=nums)
        axis[e].set_title(OTUs[e]+nline+tax[e], fontsize=10)
        axis[e].text(2.75, 85, '${'+r'p = '+'}$'+'%.3f'%float(p_vals[e]), va='bottom', ha='left', color='#bd0303', fontsize=8)
        axis[e].set_xlim([0.25, 4.5])
    ax1.set_yticks([0, 20, 40, 60, 80, 100])
    ax1.set_ylim([0, 100])
    ax1.tick_params(axis='y', which='both', left='on', right='off', labelbottom='off')
    ax1.tick_params(axis='x', which='both', bottom='off', top='off')
    axes = [ax2, ax3, ax4, ax5]
    for a in axes:
        plt.setp(a.get_yticklabels(), visible=False)
        a.tick_params(axis='y', which='both', left='off', right='off', labelbottom='off')
        a.tick_params(axis='x', which='both', bottom='off', top='off')
    plt.tight_layout()
    #ax1.set_title('B', loc='left'), ax2.set_title('C', loc='left'), ax3.set_title('D', loc='left'), ax4.set_title('E', loc='left'), ax5.set_title('F', loc='left')
    ax1.set_ylabel('Relative abundance (%)')
    os.chdir('/Users/u1560915/Documents/GitHub/ChitinActivity/MiSeq/18S/Figures/')
    plt.savefig('Fig supp 12.pdf', bbox_inches='tight')
    return 
plot_simper('3_daily_percent_simper_mean.csv')

