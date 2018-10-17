import csv
import numpy
import matplotlib.pyplot as plt
import classify_file as cf
from colorsys import hls_to_rgb
import random
import matplotlib
import numpy as np
import matplotlib as mpl
import os
#import skbio.math.diversity.alpha

#run this after analysis scripts 1 and 2!
#names, indices, all_samples, all_treatments = analysis_script_1.get_meta_treatments('16S_metadata.csv') 
files = ['0_percent_grouped.csv', '1_percent_grouped.csv', '2_daily_r_percent_grouped.csv', '3_daily_percent_grouped.csv', '4_LG_percent_grouped.csv', '5_RG_percent_grouped.csv', '6_beg_end_percent_grouped.csv']
simper_means_files = ['', '', '2_daily_r_percent_simper_mean.csv', '3_daily_percent_simper_mean.csv', '4_LG_percent_simper_mean.csv', '5_RG_percent_simper_mean.csv', '']
taxonomy = 'Taxonomy.csv'

def get_OTU_taxonomy(level, OTU):
    #levels: Genus = 6, Family = 5, Order = 4, Class = 3, Phylum = 2
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
    new_fn = cf.classify_file(fn, level)
    new_fn = cf.group_file(new_fn, True)
    #new_fn = '0_percent_grouped_Genus_grouped.csv'
    new_fn = fn
    count = 0
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
                elif b < 71:
                    this_row.append(float(rows[a][b]))
            numbers.append(this_row)
    new_rows = rows[0]
    del new_rows[0]
    gens = []
    for a in range(len(new_rows)):
        if new_rows[a][2] == '_':
            gens.append(int(new_rows[a][1]))
        else:
            gens.append(int(new_rows[a][1]+new_rows[a][2]))
    count = 0
    for a in range(len(numbers)):
        count += numbers[a][5]
    G, R, L, SG, SR = [], [], [], [], []
    for c in range(len(numbers)):
        G.append([])
        R.append([])
        L.append([])
        SG.append([])
        SR.append([])
    begs, ends, mids, daily = [], [], [], []
    shorts, new = [], []
    gens_G, gens_R, gens_L, gens_SG, gens_SR = [], [], [], [], []
    colors = get_distinct_colors(96)
    for z in range(len(new_rows)):
        begs.append(new_rows[z][0])
        ends.append(new_rows[z][-1])
        if len(new_rows[z]) > 5:
            shorts.append(new_rows[z])
            mids.append(int(new_rows[z][4]))
            daily.append(1)
        else:
            mids.append(0)
            daily.append(0)
    for d in range(len(numbers)):
        for e in range(len(new_rows)):
            this_number = numbers[d][e]
            beg, end, gen, mid = begs[e], ends[e], int(gens[e]), mids[e]
            if beg == 'L':
                if end == 'G':
                    G[d].append(this_number)
                    if d == 0:
                        gens_G.append(gen)
                if end == 'R':
                    R[d].append(this_number)
                    if d == 0:
                        gens_R.append(gen)
                if end == 'L':
                    L[d].append(this_number)
                    if d == 0:
                        gens_L.append(gen)
            if beg == 'S':
                if mid == 0 or mid == 2:
                    if end == 'G':
                        SG[d].append(this_number)
                        if d == 0:
                            gens_SG.append(gen)
                    if end == 'R':
                        SR[d].append(this_number)
                        if d == 0:
                            gens_SR.append(gen)
    new_OTU = []
    for o in range(len(OTU)):
        string = ''
        count1 = 0
        new_line = False
        for p in range(len(OTU[o])):
            count1 += 1
            if OTU[o][p] == '_':
                if not new_line:
                    string = r'${'+string+'}$'
                new_line = True
                string+=' \n '
            else:
                string+=OTU[o][p]
            if p == len(OTU[o])-1 and not new_line:
                string = r'${'+string+'}$'
        new_OTU.append(string)
    
    OTU = new_OTU
    ax1 = plt.subplot2grid((3, 4), (0, 0), colspan=3)    
    ax2 = plt.subplot2grid((3, 4), (1, 0), colspan=3, sharex=ax1)
    ax3 = plt.subplot2grid((3, 4), (2, 0), colspan=3, sharex=ax1)
    ax4 = plt.subplot2grid((3, 4), (0, 3), colspan=1)
    ax5 = plt.subplot2grid((3, 4), (1, 3), colspan=1, sharex=ax4)
    ax1.set_title('Long')
    ax4.set_title('Short')
    colors = get_distinct_colors(len(numbers))
    
    data = numpy.array(G)
    bottom = numpy.cumsum(data, axis=0)
    ax1.bar(gens_G, data[0], color=colors[0], label=OTU[0])
    for j in xrange(1, data.shape[0]):
        ax1.bar(gens_G, data[j], color=colors[j], bottom=bottom[j-1], label=OTU[j])
        
    data = numpy.array(R)
    bottom = numpy.cumsum(data, axis=0)
    ax2.bar(gens_R, data[0], color=colors[0], label=OTU[0])
    for j in xrange(1, data.shape[0]):
        ax2.bar(gens_R, data[j], color=colors[j], bottom=bottom[j-1], label=OTU[j])
        
    data = numpy.array(L)
    bottom = numpy.cumsum(data, axis=0)
    ax3.bar(gens_L, data[0], color=colors[0], label=OTU[0])
    for j in xrange(1, data.shape[0]):
        ax3.bar(gens_L, data[j], color=colors[j], bottom=bottom[j-1], label=OTU[j])
    
    data = numpy.array(SG)
    bottom = numpy.cumsum(data, axis=0)
    ax4.bar(gens_SG, data[0], color=colors[0], label=OTU[0])
    for j in xrange(1, data.shape[0]):
        ax4.bar(gens_SG, data[j], color=colors[j], bottom=bottom[j-1], label=OTU[j])
        
    data = numpy.array(SR)
    bottom = numpy.cumsum(data, axis=0)
    ax5.bar(gens_SR, data[0], color=colors[0], label=OTU[0])
    for j in xrange(1, data.shape[0]):
        ax5.bar(gens_SR, data[j], color=colors[j], bottom=bottom[j-1], label=OTU[j])
    
    ax4.set_xlim([15.8, 21])
    ax5.set_xticks([16, 17, 18, 19, 20])
    ax1.set_ylim([0, 100])
    ax2.set_ylim([0, 100])
    ax3.set_ylim([0, 100])
    ax4.set_ylim([0, 100])
    ax5.set_ylim([0, 100])
    ax1.set_xlim([-0.2, 20])
    plt.tight_layout()
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax4.get_xticklabels(), visible=False)
    plt.setp(ax4.get_yticklabels(), visible=False)
    plt.setp(ax5.get_yticklabels(), visible=False)
    #ax1.set_ylabel('Relative abundance (%)')
    ax2.set_ylabel('Relative abundance (%)')
    #ax3.set_ylabel('Relavive abundance (%)')
    plt.subplots_adjust(left=0.09,right=0.97,top=0.96,bottom=0.08,wspace=0.1, hspace=0.1)
    ax4.legend(loc='upper left', bbox_to_anchor=(1, 1.1), ncol=3, fontsize=11)
    os.chdir('/Users/u1560915/Documents/GitHub/ChitinActivity/MiSeq/16S/Figures/')
    plt.savefig('Supplementary fig 10.pdf', bbox_inches='tight')
    plt.close()
    return G, R, L, SG, SR, gens_G, gens_R, gens_L, gens_SG, gens_SR
os.chdir('/Users/u1560915/Documents/GitHub/ChitinActivity/MiSeq/16S/Supplementary 6 diversity 10 abundance/')
stacked_barchart('0_percent_grouped.csv', 6)
#levels: Genus = 6, Family = 5, Order = 4, Class = 3, Phylum = 2

def plot_heatmap():
    os.chdir('/Users/u1560915/Documents/GitHub/ChitinActivity/MiSeq/16S/Supplementary 6 diversity 10 abundance/')
    #plot heatmap with diversity and richness represented as a color for each generation
    #You can run them with Bray-Curtis, Jaccard, weighted or unweighted UniFrac to answer different questions. For example, if your variable is significant for Bray-Curtis/weighted UniFrac but not Jaccard/unweighted UniFrac, this means your groups tend to have the same OTUs (richness) but different abundances of those OTUs (diversity). When variables are signficant for Bray-Curtis/Jaccard but not UniFrac, this indicates that your samples have different specific OTUs but similar taxa.
    #High Bray-Curtis = diversity and abundance
    with open('0_diversity.csv', 'rU') as f:
        reader = csv.reader(f)
        rows = []
        for row in reader:
            rows.append(row)            
    coverage = [[], [], [], [], []]
    bergerparker = [[], [], [], [], []]
    chao = [[], [], [], [], []]
    shannon = [[], [], [], [], []]
    simpsons = [[], [], [], [], []]
    for a in range(len(rows[0])):
        if a > 0:
            name = rows[0][a]
            beg, end = name[0], name[-1]
            if beg == 'L' and end == 'G':
                coverage[0].append(rows[1][a])
                bergerparker[0].append(rows[2][a])
                chao[0].append(rows[3][a])
                shannon[0].append(rows[4][a])
                simpsons[0].append(rows[5][a])
            elif beg == 'L' and end == 'R':
                coverage[1].append(rows[1][a])
                bergerparker[1].append(rows[2][a])
                chao[1].append(rows[3][a])
                shannon[1].append(rows[4][a])
                simpsons[1].append(rows[5][a])
            elif beg == 'L' and end == 'L':
                coverage[2].append(rows[1][a])
                bergerparker[2].append(rows[2][a])
                chao[2].append(rows[3][a])
                shannon[2].append(rows[4][a])
                simpsons[2].append(rows[5][a])
            elif beg == 'S' and end == 'G':
                coverage[3].append(rows[1][a])
                bergerparker[3].append(rows[2][a])
                chao[3].append(rows[3][a])
                shannon[3].append(rows[4][a])
                simpsons[3].append(rows[5][a])
            elif beg == 'S' and end == 'R':
                coverage[4].append(rows[1][a])
                bergerparker[4].append(rows[2][a])
                chao[4].append(rows[3][a])
                shannon[4].append(rows[4][a])
                simpsons[4].append(rows[5][a])
    all_diversity = [coverage, bergerparker, chao, shannon, simpsons]
    all_max, all_min = [], []
    for b in range(len(all_diversity)):
        m, s = 0, 1
        maxs, mins = [], []
        for c in range(len(all_diversity[b])):
            maxs.append(max(all_diversity[b][c]))
            mins.append(min(all_diversity[b][c]))
            for d in range(len(all_diversity[b][c])):
                all_diversity[b][c][d] = float(all_diversity[b][c][d])
                if all_diversity[b][c][d] > m:
                    m = all_diversity[b][c][d]
                if all_diversity[b][c][d] < s:
                    s = all_diversity[b][c][d]
        ma, mi = max(maxs), min(mins)
        all_max.append(float(ma))
        all_min.append(float(mi))
    for e in range(len(all_diversity)):
        m = all_max[e]
        for f in range(len(all_diversity[e])):
            for g in range(len(all_diversity[e][f])):
                num = all_diversity[e][f][g]/all_max[e]
                all_diversity[e][f][g] = num
    all_max, all_min = [], []
    for z in range(len(all_diversity)):
        maxs, mins = [], []
        for y in range(len(all_diversity[z])):
            maxs.append(max(all_diversity[z][y]))
            mins.append(min(all_diversity[z][y]))
        all_max.append(max(maxs))
        all_min.append(min(mins))
    ax1 = plt.subplot2grid((3, 4), (0, 0), colspan=3)    
    ax2 = plt.subplot2grid((3, 4), (1, 0), colspan=3, sharex=ax1, sharey=ax1)
    ax3 = plt.subplot2grid((3, 4), (2, 0), colspan=3, sharex=ax1, sharey=ax1)
    ax4 = plt.subplot2grid((3, 4), (0, 3), colspan=1, sharey=ax1)
    ax5 = plt.subplot2grid((3, 4), (1, 3), colspan=1, sharex=ax4, sharey=ax1)
    x1 = [0, 1, 3, 4, 6, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20]
    x2 = [0, 1, 2, 3, 5, 6, 8, 9, 10, 11, 12, 13, 14, 16, 18, 19, 20]
    x3 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20]
    x4 = [16, 17, 18, 19, 20]
    x5 = [16, 17, 18, 19, 20]
    for a in range(21):
        if a != x1[a]:
            x1.insert(a, a)
            all_diversity[0][0].insert(a, 0)
            all_diversity[1][0].insert(a, 0)
            all_diversity[2][0].insert(a, 0)
            all_diversity[3][0].insert(a, 0)
            all_diversity[4][0].insert(a, 0)
    for b in range(21):
        if b != x2[b]:
            x2.insert(b, b)
            all_diversity[0][1].insert(b, 0)
            all_diversity[1][1].insert(b, 0)
            all_diversity[2][1].insert(b, 0)
            all_diversity[3][1].insert(b, 0)
            all_diversity[4][1].insert(b, 0)
    for c in range(21):
        if c != x3[c]:
            x3.insert(c, c)
            all_diversity[0][2].insert(c, 0)
            all_diversity[1][2].insert(c, 0)
            all_diversity[2][2].insert(c, 0)
            all_diversity[3][2].insert(c, 0)
            all_diversity[4][2].insert(c, 0)
    G, R, L, GS, RS = [all_diversity[0][0], all_diversity[1][0], all_diversity[2][0], all_diversity[3][0], all_diversity[4][0]], [all_diversity[0][1], all_diversity[1][1], all_diversity[2][1], all_diversity[3][1], all_diversity[4][1]], [all_diversity[0][2], all_diversity[1][2], all_diversity[2][2], all_diversity[3][2], all_diversity[4][2]], [all_diversity[0][3], all_diversity[1][3], all_diversity[2][3], all_diversity[3][3], all_diversity[4][3]], [all_diversity[0][4], all_diversity[1][4], all_diversity[2][4], all_diversity[3][4], all_diversity[4][4]]
    cmap = 'Blues'
    norm = matplotlib.colors.Normalize(vmin=0, vmax=1)
    colormap = matplotlib.cm.get_cmap(cmap, 256)
    m = matplotlib.cm.ScalarMappable(norm=norm, cmap=colormap)
    colors1 = []
    for a in range(len(G)):
        this_row = []
        for b in range(len(G[a])):
            if G[a][b] == 0:
                color = (0, 0, 0, 0)
            else:
                color = m.to_rgba(G[a][b])
            this_row.append(color)
        colors1.append(this_row)
    colors = colors1
    l = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    plot = [l, l, l, l]
    data = numpy.array(plot)
    bottom = numpy.cumsum(data, axis=0)
    ax1.bar(x1, data[0], color=colors[0], width=1.0)
    for j in xrange(1, data.shape[0]):
        ax1.bar(x1, data[j], color=colors[j], bottom=bottom[j-1], width=1.0)
    colors2 = []
    for a in range(len(R)):
        this_row = []
        for b in range(len(R[a])):
            if R[a][b] == 0:
                color = (0, 0, 0, 0)
            else:
                color = m.to_rgba(R[a][b])
            this_row.append(color)
        colors2.append(this_row)
    colors = colors2
    ax2.bar(x1, data[0], color=colors[0], width=1.0)
    for j in xrange(1, data.shape[0]):
        ax2.bar(x1, data[j], color=colors[j], bottom=bottom[j-1], width=1.0)
    colors3 = []
    for a in range(len(L)):
        this_row = []
        for b in range(len(L[a])):
            if L[a][b] == 0:
                color = (0, 0, 0, 0)
            else:
                color = m.to_rgba(L[a][b])
            this_row.append(color)
        colors3.append(this_row)
    colors = colors3
    ax3.bar(x1, data[0], color=colors[0], width=1.0)
    for j in xrange(1, data.shape[0]):
        ax3.bar(x1, data[j], color=colors[j], bottom=bottom[j-1], width=1.0)
        
    s = [1, 1, 1, 1, 1]
    plot = [s, s, s, s]
    data = numpy.array(plot)
    bottom = numpy.cumsum(data, axis=0)
    colors4 = []
    for a in range(len(GS)):
        this_row = []
        for b in range(len(GS[a])):
            color = m.to_rgba(GS[a][b])
            this_row.append(color)
        colors4.append(this_row)
    colors = colors4
    ax4.bar(x4, data[0], color=colors[0], width=1.0)
    for j in xrange(1, data.shape[0]):
        ax4.bar(x4, data[j], color=colors[j], bottom=bottom[j-1], width=1.0)
        
    colors5 = []
    for a in range(len(RS)):
        this_row = []
        for b in range(len(RS[a])):
            color = m.to_rgba(RS[a][b])
            this_row.append(color)
        colors5.append(this_row)
    colors = colors5
    ax5.bar(x4, data[0], color=colors[0], width=1.0)
    for j in xrange(1, data.shape[0]):
        ax5.bar(x4, data[j], color=colors[j], bottom=bottom[j-1], width=1.0)
    ax1.set_xlim([0, 21])
    plt.setp(ax4.get_yticklabels(), visible=False)
    plt.setp(ax5.get_yticklabels(), visible=False)
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax4.get_xticklabels(), visible=False)
    ax3.set_xticks(x1)
    os.chdir('/Users/u1560915/Documents/GitHub/ChitinActivity/MiSeq/16S/Figures/')
    #plt.savefig('All_diversity.pdf', bbox_inches='tight')
    plt.close()    
    
    names = ['Coverage', 'Berger-Parker', 'Chao', 'Shannon', 'Simpsons']
    x = [x1, x2, x3, x4, x5]
    all_colors = [colors1, colors2, colors3, colors4, colors5]
    y = [l, l, l, s, s]
    for a in range(5):
        diversity = all_diversity[a]
        ax1 = plt.subplot2grid((10, 4), (0, 0), colspan=3)    
        ax2 = plt.subplot2grid((10, 4), (1, 0), colspan=3, sharex=ax1, sharey=ax1)
        ax3 = plt.subplot2grid((10, 4), (2, 0), colspan=3, sharex=ax1, sharey=ax1)
        ax4 = plt.subplot2grid((10, 4), (0, 3), colspan=1, sharey=ax1)
        ax5 = plt.subplot2grid((10, 4), (1, 3), colspan=1, sharex=ax4, sharey=ax1)
        ax6 = plt.subplot2grid((10, 4), (2, 3), colspan=1)
        ax6.set_ylim([0.1,1])
        maxs, mins = [], []
        ma, mi = 0, 1
        for c in range(len(diversity)):
            maxs.append(max(diversity[c]))
            mins.append(min(diversity[c]))
            for d in range(len(diversity[c])):
                if diversity[c][d] > ma:
                    ma = diversity[c][d]
                if diversity[c][d] < mi:
                    mi = diversity[c][d]
        axes = [ax1, ax2, ax3, ax4, ax5]
        for b in range(5):
            ma, mi = all_max[b], all_min[b]
            cmap = mpl.cm.Blues
            norm = matplotlib.colors.Normalize(vmin=mi, vmax=ma)
            colormap = matplotlib.cm.get_cmap(cmap, 256)
            m = matplotlib.cm.ScalarMappable(norm=norm, cmap=colormap)
            colors = []
            for c in range(len(diversity[b])):
                if diversity[b][c] == 0:
                    color = (0, 0, 0, 0)
                else:
                    color = m.to_rgba(diversity[b][c])
                colors.append(color)
            axes[b].bar(x[b], y[b], color=colors, width=1.0)
        plt.setp(ax1.get_yticklabels(), visible=False)
        plt.setp(ax1.get_xticklabels(), visible=False)
        plt.setp(ax2.get_yticklabels(), visible=False)
        plt.setp(ax3.get_yticklabels(), visible=False)
        plt.setp(ax3.get_xticklabels(), visible=False)
        plt.setp(ax4.get_yticklabels(), visible=False)
        plt.setp(ax5.get_yticklabels(), visible=False)
        plt.setp(ax1.get_xticklabels(), visible=False)
        plt.setp(ax2.get_xticklabels(), visible=False)
        plt.setp(ax4.get_xticklabels(), visible=False)
        plt.setp(ax5.get_xticklabels(), visible=False)
        ax1.tick_params(axis='y',which='both',left='off',right='off')
        ax2.tick_params(axis='y',which='both',left='off',right='off')
        ax3.tick_params(axis='y',which='both',left='off',right='off')
        ax4.tick_params(axis='y',which='both',left='off',right='off')
        ax5.tick_params(axis='y',which='both',left='off',right='off')
        ax1.set_title('Long')
        ax4.set_title('Short')
        ax1.text(-2.5, 0.25, 'Good')
        ax2.text(-3.5, 0.25, 'Random')
        ax3.text(-2.5, 0.25, 'Light')
        ax1.plot([2,3], [0,1], 'k')
        ax1.plot([5,6], [0,1], 'k')
        ax1.plot([7,8], [0,1], 'k')
        ax1.plot([15,16], [0,1], 'k')
        ax2.plot([4,5], [0,1], 'k')
        ax2.plot([7,8], [0,1], 'k')
        ax2.plot([15,16], [0,1], 'k')
        ax2.plot([17,18], [0,1], 'k')
        ax3.plot([15,16], [0,1], 'k')
        cmap = mpl.cm.Blues
        norm = mpl.colors.Normalize(vmin=mi, vmax=ma)
        cb1 = mpl.colorbar.ColorbarBase(ax6, cmap=cmap, norm=norm, orientation='horizontal')
        cb1.set_ticks([])
        #cb1.set_ticklabels(['Low', 'High'])
        #cb1.ax.tick_params(labelsize=8) 
        cb1.ax.text(-0.1, -1, 'Low')
        cb1.ax.text(0.9, -1, 'High')
        ax1.set_xlim([0, 21])
        plt.subplots_adjust(left=0.09,right=0.97,top=0.96,bottom=0.08,wspace=0.15, hspace=0.5)
        os.chdir('/Users/u1560915/Documents/GitHub/ChitinActivity/MiSeq/16S/Figures/')
        #plt.savefig('Diversity'+names[a]+'.pdf', bbox_inches='tight')
        if names[a] == 'Simpsons':
            plt.savefig('Supplementary fig 6.pdf', bbox_inches='tight')
        plt.close()      
    return
plot_heatmap()