#Daily 16S and 18S
import matplotlib.pyplot as plt
import csv
from colorsys import hls_to_rgb
import numpy
import random
import statsmodels.stats.multitest as smm
import group_file as gf
import matplotlib
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import matplotlib.patches as mpatches

#these take files that are given from the initial analysis (for each of 16S and 18S):
#all = treatment1_daily_percent_grouped.csv
#simper = treatment1_daily_simper_means.csv
#no_percent = treatment1_daily_samples_grouped.csv
fn_18S_DADA = '18S_all_DADA_bar.csv' 
sim_18S_DADA = '18S_simper_DADA.csv'
tax_18S_DADA = '18S_taxonomy_DADA.csv'

def get_tax(otus, tax):
    with open(tax, 'rU') as f:
        rows = []
        for r in csv.reader(f):
            rows.append(r)
    new_tax = []
    for a in range(len(otus)):
        for b in range(len(rows)):
            if otus[a] == rows[b][0]:
                new_tax.append(rows[b])
    new_new_tax = []
    for c in range(len(new_tax)):
        if new_tax[c][7] != 'NA':
            new_new_tax.append(r'$'+new_tax[c][6]+' '+new_tax[c][7]+'$')
        elif new_tax[c][6] != 'NA':
            new_new_tax.append(r'$'+new_tax[c][6]+'$')
        elif new_tax[c][5] != 'NA':
            new_new_tax.append(new_tax[c][5])
        elif new_tax[c][4] != 'NA':
            new_new_tax.append(new_tax[c][4])
        elif new_tax[c][3] != 'NA':
            new_new_tax.append(new_tax[c][3])
        elif new_tax[c][2] != 'NA':
            new_new_tax.append(new_tax[c][2])
        elif new_tax[c][1] != 'NA':
            new_new_tax.append(new_tax[c][1])
        else:
            new_new_tax.append('NA')
    return new_new_tax
    
def get_old_otus(otus):
    old_otus = []
    prev_otus = otus
    for d in range(len(prev_otus)):
        if prev_otus[d][:3] != 'ASV':
            old_otus.append(prev_otus[d])
            continue
        otu = int(prev_otus[d][3:])
        otu = 'ASV'+str(otu)
        old_otus.append(otu)
    return old_otus

def simper(fn, tax, rng, axis, x, xlim, ylim, colors, xtxt, ytxt, ylab):
    with open(fn, 'rU') as f:
        rows = []
        for r in csv.reader(f):
            rows.append(r)
    otus, rest_of_row = [], []
    cont, treat_mean, treat_sd, krusk, krusk_p = [], [], [], [], []
    for a in range(len(rows)):
        for b in range(len(rows[a])):
            if a > 0 and b > 0:
                rows[a][b] = float(rows[a][b])
    for a in range(rng):
        a += 1
        otus.append(rows[a][0])
        rest_of_row.append(rows[a][1:])
        cont.append(rows[a][1]*100)
        this_mean, this_sd = [rows[a][2], rows[a][4], rows[a][6], rows[a][8]], [rows[a][3], rows[a][5], rows[a][7], rows[a][9]]
        treat_mean.append(this_mean)
        treat_sd.append(this_sd)
        krusk.append(rows[a][10])
        krusk_p.append(rows[a][11])
    old = get_old_otus(otus)
    otus = get_tax(otus, tax)
    for a in range(len(axis)):
        ax = axis[a]
        if a > 0:
            plt.setp(ax.get_yticklabels(), visible=False)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.bar(x, treat_mean[a], yerr=treat_sd[a], color=colors, error_kw=dict(ecolor='gray', lw=1, capsize=3, capthick=1, alpha=0.5), edgecolor='k')
        title = old[a]+'\n'
        if otus[a] == 'Bacteria':
             otus[a] = r'($Spirochaeta$)'
        if otus[a] == 'Eukaryota':
            otus[a] = r'($Cafeteria$)'
        for b in otus[a]:
            if b != '_':
                title += b
            else:
                title += ' '
        title += '\nSIMPER: %.0f'%cont[a]+'%'
        ax.set_title(title, fontsize=8)
        h, p = 'H=%.2f'%krusk[a], '${'+r'p = '+'}$'+'%.3f'%float(krusk_p[a])
        ax.text(xtxt, ytxt, h+', '+p, va='bottom', ha='left', color='#bd0303', fontsize=8)
        ax.tick_params(axis='y',which='both',left='on',right='off')
        ax.tick_params(axis='x',which='both',top='off',bottom='on')
        plt.sca(ax)
        plt.xticks(x)
        #plt.setp(ax, xticks=[1.4, 2.4, 3.4, 4.4], xticklabels=[1, 2, 3, 4])
        ax.set_xlim(xlim)
    axis[0].set_ylabel(ylab+'\n', fontsize=12)
    ylab = '\n Relative abundance (%)'
    axis[0].set_ylabel('\n'+ylab, fontsize=8)
    return
    
def get_distinct_colors(n):
    colors = []
    for i in numpy.arange(0., 360., 360. / n):
        h = i / 360.
        l = (50 + numpy.random.rand() * 10) / 100.
        s = (90 + numpy.random.rand() * 10) / 100.
        colors.append(hls_to_rgb(h, l, s))
    random.shuffle(colors)
    return colors

def barplot(fn, tax, lim, ax, alpha):
    new_rows, old_otus = gf.gf(fn, lim, tax, return_old_otus=True)
    taxa = []
    for a in range(len(old_otus)):
        if a > 0:
            if new_rows[a][0] == 'Clostridiales_vadinBB60_group':
                new_rows[a][0] = 'Clostridiales'
            if new_rows[a][0] == 'Incertae_Sedis':
                new_rows[a][0] = 'Incertae Sedis'
            if new_rows[a][0] == '$Kurtzmaniella-Candida_clade$':
                new_rows[a][0] = '$Kurtzmaniella$'
            if new_rows[a][0] == '$Thalassolituus$ $oleivorans$':
                new_rows[a][0] = '$Thalassolituus$\n          $oleivorans$'
            name = old_otus[a]+' '+new_rows[a][0]
            if name == 'ASV2 Bacteria':
                name = 'ASV2 $(Spirochaeta)$'
            if name == 'ASV125 Bacteria':
                name = 'ASV125 $(Thermobrachium)$'
            if name == 'ASV1 Eukaryota':
                name = 'ASV1 $(Cafeteria)$'
            if name == 'ASV2 Eukaryota':
                name = 'ASV2 $(Cafeteria)$'
            if name == 'ASV35 Eukaryota':
                name = 'ASV35 $(Oxyrrhis$ $marina)$'
            taxa.append(name)
    taxa.append(new_rows[-1][0])
    samples = []
    for b in range(len(new_rows)):
        if b > 0:
            samples.append(new_rows[b][1:])
    data = numpy.array(samples)
    bottom = numpy.cumsum(data, axis=0)
    x = [1, 2, 3, 4]
    colors = get_distinct_colors(len(samples))
    ax.bar(x, data[0], color=colors[0], label=taxa[0], alpha=alpha, edgecolor='k')
    for j in range(1, data.shape[0]):
        ax.bar(x, data[j], color=colors[j], bottom=bottom[j-1], label=taxa[j], alpha=alpha, edgecolor='k')
    ax.set_ylim([0, 100])
    if ax == l16S:
        ax.legend(bbox_to_anchor=(2.23, 1.01), fontsize=7)
    elif ax == l18S:
        ax.legend(bbox_to_anchor=(2.12, 1.01), fontsize=7)
    ax.set_xlim([0.5, 4.5])
    plt.sca(ax)
    plt.xticks(x)
    #plt.setp(ax, xticks=[1.4, 2.4, 3.4, 4.4], xticklabels=[1, 2, 3, 4])
    return
    
def get_diversity(diversity, sample):
    for a in range(len(sample)):
        sample[a] = float(sample[a])
    total = sum(sample)
    if diversity == 'Simpsons':
        for b in range(len(sample)):
            sample[b] = (sample[b]/total)**2
        simpsons = 1-(sum(sample))
        return simpsons
    elif diversity == 'Shannon':
        for b in range(len(sample)):
            sample[b] = (sample[b]/total)
            if sample[b] != 0:
                sample[b] = -(sample[b] * (numpy.log(sample[b])))
        shannon = sum(sample)
        return shannon
    elif diversity == 'Richness':
        rich = 0
        for b in range(len(sample)):
            if sample[b] != 0:
                rich += 1
        return rich
    elif diversity == 'Evenness':
        for b in range(len(sample)):
            sample[b] = (sample[b]/total)
            if sample[b] != 0:
                sample[b] = -(sample[b] * (numpy.log(sample[b])))
        shannon = sum(sample)
        rich = 0
        for b in range(len(sample)):
            if sample[b] != 0:
                rich += 1
        even = shannon/(numpy.log(rich))
        return even
    return

def get_diversity_plot(fn, ax,fn2, ax2, diversity='Simpsons'):
    with open(fn, 'rU') as f:
        rows = []
        for row in csv.reader(f):
            rows.append(row)
    transpose = []
    for a in range(len(rows[0])):
        if a == 0:
            continue
        new_col = []
        for b in range(len(rows)):
            if a > 0 and b > 0:
                new_col.append(rows[b][a])
        transpose.append(new_col)
    d = [[], []]
    for c in range(len(transpose)):
        simpsons = get_diversity(diversity, transpose[c])
        d[0].append(simpsons)
    with open(fn2, 'rU') as f:
        rows = []
        for row in csv.reader(f):
            rows.append(row)
    transpose = []
    for a in range(len(rows[0])):
        if a == 0:
            continue
        new_col = []
        for b in range(len(rows)):
            if a > 0 and b > 0:
                new_col.append(rows[b][a])
        transpose.append(new_col)
    for c in range(len(transpose)):
        simpsons = get_diversity(diversity, transpose[c])
        d[1].append(simpsons)
    print(d[0])
    print(d[1])
    cmap = 'Blues'
    colors1 = []
    mi = 2
    ma = 0
    for x in range(len(d)):
        if min(d[x]) < mi:
            mi = min(d[x])
        if max(d[x]) > ma:
            ma = max(d[x])
    for a in range(len(d)):
        norm = matplotlib.colors.Normalize(vmin=mi, vmax=ma)
        colormap = matplotlib.cm.get_cmap(cmap, 256)
        m = matplotlib.cm.ScalarMappable(norm=norm, cmap=colormap)
        this_row = []
        for b in range(len(d[a])):
            if d[a][b] == 0:
                color = (0, 0, 0, 0)
            else:
                color = m.to_rgba(d[a][b])
            this_row.append(color)
        colors1.append(this_row)
    x = [1, 2, 3, 4]
    y = [1, 1, 1, 1]
    ax.bar(x, y, color=colors1[0], width=1.0, edgecolor='k')
    #plt.setp(ax, xticks=[1.4, 2.4, 3.4, 4.4], xticklabels=[1, 2, 3, 4])
    plt.setp(ax.get_yticklabels(), visible=False)
    plt.setp(ax.get_xticklabels(), visible=False)
    ax.set_xlim([0.5, 4.5])
    ax2.tick_params(axis='x',which='both',top='off', bottom='off')
    ax2.bar(x, y, color=colors1[1], width=1.0, edgecolor='k')
    ax.tick_params(axis='x',which='both',top='off', bottom='off')
    #plt.setp(ax2, xticks=[1.4, 2.4, 3.4, 4.4], xticklabels=[1, 2, 3, 4])
    plt.setp(ax2.get_yticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)
    ax.set_ylim([0,1])
    ax2.set_ylim([0,1])
    ax2.set_xlim([0.5, 4.5])
    print(mi, ma)
    return mi, ma, diversity

fig = plt.figure(figsize=(8.27, 12))
h1, w, ss1, ss2 = 17, 10, 23, 30
h = 35
rsb, rss, rssimp = 10, 5, 1
l16S = plt.subplot2grid((h1,12), (0,0), colspan=3, rowspan=rsb)
l18S = plt.subplot2grid((h1,12), (0,6), colspan=3, rowspan=rsb, sharey=l16S, sharex=l16S)
s1_16S = plt.subplot2grid((h,w), (ss1,0), colspan=2, rowspan=rss)
s2_16S, s3_16S, s4_16S, s5_16S = plt.subplot2grid((h,w), (ss1,2), sharey=s1_16S, colspan=2, rowspan=rss), plt.subplot2grid((h,w), (ss1,4), rowspan=rss, sharey=s1_16S, colspan=2), plt.subplot2grid((h,w), (ss1,6), rowspan=rss, sharey=s1_16S, colspan=2), plt.subplot2grid((h,w), (ss1,8), rowspan=rss, sharey=s1_16S, colspan=2)
s1_18S = plt.subplot2grid((h,w), (ss2,0), colspan=2, rowspan=rss)
s2_18S, s3_18S, s4_18S, s5_18S = plt.subplot2grid((h,w), (ss2,2), rowspan=rss, sharey=s1_18S, colspan=2), plt.subplot2grid((h,w), (ss2,4), rowspan=rss, sharey=s1_18S, colspan=2), plt.subplot2grid((h,w), (ss2,6), rowspan=rss, sharey=s1_18S, colspan=2), plt.subplot2grid((h,w), (ss2,8), rowspan=rss, sharey=s1_18S, colspan=2)

lblank = plt.subplot2grid((h1,12), (0,9), colspan=2, rowspan=rsb, frameon=False)

s16S, s18S = [s1_16S, s2_16S, s3_16S, s4_16S, s5_16S], [s1_18S, s2_18S, s3_18S, s4_18S, s5_18S]
removey = [s2_16S, s3_16S, s4_16S, s5_16S, s2_18S, s3_18S, s4_18S, s5_18S, lblank]
removex = [s1_16S, s2_16S, s3_16S, s4_16S, s5_16S]
fsl, fst, fsst, fspv = 6, 9, 8, 8
cols_16S, cols_18S = ['#33FFFF', '#33CCFF', '#3399FF', '#3366FF'], ['#CCFF99', '#66FF66', '#009900', '#006400']

for a in removey:
    plt.setp(a.get_yticklabels(), visible=False)
plt.setp(lblank.get_xticklabels(), visible=False)
lblank.tick_params(axis='x',which='both',top='off', bottom='off')
lblank.tick_params(axis='y',which='both',left='off',right='off')

simper(sim_18S_DADA, tax_18S_DADA, 5, s16S, [1, 2, 3, 4], [0.5, 4.5], [0, 100], cols_18S, 0.6, 85, '18S')

barplot(fn_18S_DADA, tax_18S_DADA, 1, l16S, 0.8)

l16S.text(-0.98, 99, 'A', fontsize=16, fontweight='bold', color='k')#, bbox=dict(facecolor='white', edgecolor='white'))
l16S.text(-0.98, -16, 'B', fontsize=16, fontweight='bold', color='k')
l16S.set_title('DADA2 analysis', fontsize=12)
l18S.set_title('Mothur analysis', fontsize=12)

l16S.set_ylabel('Relative abundance (%)')
l16S.set_xlabel('Days')
l18S.set_xlabel('Days')
s3_18S.set_xlabel('Days')
s1_16S.text(-1.6, 50, 'DADA2 analysis', fontsize=10, ha='center', va='center', rotation=90)

fn_18S = '18S_daily_percent_grouped.csv'
sim_18S = '18S_daily_simper.csv'
ord_18S = '18S_3_order.csv'
tax_18S = '18S_taxonomy_mothur.csv'
meta_18S = '18S_daily_meta.csv'

def simper_mothur(fn, order, meta, tax, rng):
    with open(order, 'rU') as f:
        rows = []
        for r in csv.reader(f):
            rows.append(r)
    simp_order, cont = [], []
    for row in range(len(rows)):
        if row > 0 and row < 6:
            simp_order.append(rows[row][0])
            cont.append(float(rows[row][1])*100)
    with open(tax, 'rU') as f:
        rows = []
        for row in csv.reader(f):
            rows.append(row)
    tax = []
    for a in range(len(simp_order)):
        for b in range(len(rows)):
            if simp_order[a] == rows[b][0]:
                phylo = [rows[b][2], rows[b][3], rows[b][4], rows[b][5], rows[b][6]]
                if phylo[4][-12:] != 'unclassified':
                    this_tax = r'$'+str(phylo[4])+'$'
                else:
                    this_tax = phylo[4][:-13]
                tax.append(this_tax)          
    print_otus = []
    for c in simp_order:
        totu = 'OTU'
        d = 0
        while d < len(c):
            if d > 2 and c[d] != '0':
                totu += c[d:]
                d = len(c)
            d += 1
        totu += '\n'
        print_otus.append(totu)
    with open(fn, 'rU') as f:
        rows = []
        for row in csv.reader(f):
            rows.append(row)
    simp_rows = []
    for e in range(len(simp_order)):
        for f in range(len(rows)):
            if rows[f][0] == simp_order[e]:
                simp_rows.append(rows[f])
    krusk, krusk_p, treat_mean, treat_sd = [], [], [], []
    for g in simp_rows:
        krusk.append(float(g[-2]))
        krusk_p.append(float(g[-1]))
        this_mean, this_sd = [], []
        for h in range(rng):
            h += 1
            if h % 2 != 0:
                this_mean.append(float(g[h])*100)
            else:
                this_sd.append(float(g[h])*100)
        treat_mean.append(this_mean)
        treat_sd.append(this_sd)
    krusk_p = smm.fdrcorrection(krusk_p)[1]
    return krusk, krusk_p, treat_mean, treat_sd, cont, print_otus, tax

def simper_plot_mothur(krusk, krusk_p, treat_mean, treat_sd, cont, print_otus, tax, axis, x, xlim, ylim, colors, xtxt, ytxt, ylab):
    for a in range(len(axis)):
        ax = axis[a]
        if a > 0:
            plt.setp(ax.get_yticklabels(), visible=False)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.bar(x, treat_mean[a], yerr=treat_sd[a], color=colors, error_kw=dict(ecolor='gray', lw=1, capsize=3, capthick=1, alpha=0.5), edgecolor='k')
        if tax[a] == 'Eukaryota':
            tax[a] = r'$(Cafeteria)$'
        title = print_otus[a]+tax[a]+'\n'
        title += 'SIMPER: %.0f'%cont[a]+'%'
        ax.set_title(title, fontsize=8)
        h, p = 'H=%.2f'%krusk[a], '${'+r'p = '+'}$'+'%.3f'%float(krusk_p[a])
        ax.text(xtxt, ytxt, h+', '+p, va='bottom', ha='left', color='#bd0303', fontsize=8)
        ax.tick_params(axis='y',which='both',left='on',right='off')
        ax.tick_params(axis='x',which='both',top='off',bottom='on')
        plt.setp(ax, xticks=[1, 2, 3, 4], xticklabels=[1, 2, 3, 4])
        ax.set_xlim([0.5, 4.5])
    ylab += '\n Relative abundance (%)'
    axis[0].set_ylabel(ylab, fontsize=8)
    return
    
def get_tax_mothur(otus, tax):
    with open(tax, 'rU') as f:
        rows = []
        for row in csv.reader(f):
            rows.append(row)
    tax = []
    for a in range(len(otus)):
        for b in range(len(rows)):
            if otus[a] == rows[b][0]:
                phylo = [rows[b][2], rows[b][3], rows[b][4], rows[b][5], rows[b][6]]
                if phylo[4][-12:] != 'unclassified' or phylo[4][-12:] == '':
                    this_tax = r'$'+str(phylo[4])+'$'
                else:
                    this_tax = phylo[4][:-13]
                for c in range(len(this_tax)):
                    if this_tax[c] == '_':
                        this_tax = this_tax[:c]+' '+this_tax[c+1:]
                tax.append(this_tax)
    for c in range(len(otus)):
        totu = 'OTU'
        d = 3
        while otus[c][d] == '0':
            d += 1
        totu += otus[c][d:]+' '
        otus[c] = totu
    for e in range(len(otus)):
        otus[e] += tax[e]
    return otus
    
def barplot_mothur(fn, tax, lim, ax, alpha):
    with open(fn, 'rU') as f:
        rows = []
        for row in csv.reader(f):
            rows.append(row)
    otus, samples = [], []
    for r in range(len(rows)):
        if r > 0:
            this_row = []
            for c in range(len(rows[r])):
                if c == 0:
                    otus.append(rows[r][c])
                else:
                    this_row.append(float(rows[r][c])*100)
            samples.append(this_row)
    other, new_samples, new_otus = [], [], []
    for a in range(len(samples[0])):
        other.append(0)
    for b in range(len(samples)):
        if max(samples[b]) > lim:
            new_samples.append(samples[b])
            new_otus.append(otus[b])
        else:
            for c in range(len(samples[b])):
                other[c] += samples[b][c]
    otus = get_tax_mothur(new_otus, tax)
    otus.append('Other')
    for c in range(len(otus)):
        if otus[c] == 'OTU2 Bacteria':
            otus[c] = 'OTU2 '+r'$(Spirochaeta)$'
        if otus[c] == 'OTU190 Bacteria':
            otus[c] = 'OTU190 '+r'$(Thermobrachium)$'
        if otus[c] == 'OTU1 Eukaryota':
            otus[c] = 'OTU1 '+r'$(Cafeteria)$'
        if otus[c] == 'OTU2 Eukaryota':
            otus[c] = 'OTU2 '+r'$(Cafeteria)$'
    new_samples.append(other)
    data = numpy.array(new_samples)
    bottom = numpy.cumsum(data, axis=0)
    x = [1, 2, 3, 4]
    colors = get_distinct_colors(len(new_samples))
    ax.bar(x, data[0], color=colors[0], label=otus[0], alpha=alpha, edgecolor='k')
    for j in range(1, data.shape[0]):
        ax.bar(x, data[j], color=colors[j], bottom=bottom[j-1], label=otus[j], alpha=alpha, edgecolor='k')
    ax.set_ylim([0, 100])
    ax.legend(bbox_to_anchor=(1.0, 1.01), fontsize=7)
    plt.setp(ax, xticks=[1, 2, 3, 4], xticklabels=[1, 2, 3, 4])
    ax.set_xlim([0.5, 4.5])
    return

krusk, krusk_p, treat_mean, treat_sd, cont, print_otus, tax = simper_mothur(sim_18S, ord_18S, meta_18S, tax_18S, 8)
simper_plot_mothur(krusk, krusk_p, treat_mean, treat_sd, cont, print_otus, tax, s18S, [1, 2, 3, 4], [0.75, 5], [0, 100], cols_18S, 0.6, 85, '')

barplot_mothur(fn_18S, tax_18S, 1, l18S, 0.8)



s1_18S.text(-1.6, 50, 'Mothur analysis', fontsize=10, ha='center', va='center', rotation=90)
plt.setp(l18S.get_yticklabels(), visible=False)
fig.subplots_adjust(hspace=200, wspace=0.4)
plt.savefig('18S daily comparison.png', bbox_inches='tight', dpi=600)

