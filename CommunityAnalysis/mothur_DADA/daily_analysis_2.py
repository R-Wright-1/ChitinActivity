#Daily 16S and 18S
import matplotlib.pyplot as plt
import csv
from colorsys import hls_to_rgb
import numpy
import random
import statsmodels.stats.multitest as smm

fn_16S, fn_18S = '16S_daily_percent_grouped.csv', '18S_daily_percent_grouped.csv'
sim_16S, sim_18S = '16S_daily_simper.csv', '18S_daily_simper.csv'
ord_16S, ord_18S = '16S_3_order.csv', '18S_3_order.csv'
tax_16S, tax_18S = '16S_taxonomy.csv', '18S_taxonomy.csv'
meta_16S, meta_18S = '16S_daily_meta.csv', '18S_daily_meta.csv'

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
    
def get_distinct_colors(n):
    colors = []
    for i in numpy.arange(0., 360., 360. / n):
        h = i / 360.
        l = (50 + numpy.random.rand() * 10) / 100.
        s = (90 + numpy.random.rand() * 10) / 100.
        colors.append(hls_to_rgb(h, l, s))
    random.shuffle(colors)
    return colors
    
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
    otus = get_tax(new_otus, tax)
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


fig = plt.figure(figsize=(8.27, 12))
h, w, ss1, ss2 = 5, 10, 3, 4
l16S = plt.subplot2grid((h,12), (0,0), colspan=3, rowspan=3)
l18S = plt.subplot2grid((h,12), (0,6), colspan=3, rowspan=3, sharey=l16S, sharex=l16S)
s1_16S = plt.subplot2grid((h,w), (ss1,0), colspan=2)
s2_16S, s3_16S, s4_16S, s5_16S = plt.subplot2grid((h,w), (ss1,2), sharey=s1_16S, colspan=2), plt.subplot2grid((h,w), (ss1,4), sharey=s1_16S, colspan=2), plt.subplot2grid((h,w), (ss1,6), sharey=s1_16S, colspan=2), plt.subplot2grid((h,w), (ss1,8), sharey=s1_16S, colspan=2)
s1_18S = plt.subplot2grid((h,w), (ss2,0), colspan=2)
s2_18S, s3_18S, s4_18S, s5_18S = plt.subplot2grid((h,w), (ss2,2), sharey=s1_18S, colspan=2), plt.subplot2grid((h,w), (ss2,4), sharey=s1_18S, colspan=2), plt.subplot2grid((h,w), (ss2,6), sharey=s1_18S, colspan=2), plt.subplot2grid((h,w), (ss2,8), sharey=s1_18S, colspan=2)
s16S, s18S = [s1_16S, s2_16S, s3_16S, s4_16S, s5_16S], [s1_18S, s2_18S, s3_18S, s4_18S, s5_18S]
removey = [s2_16S, s3_16S, s4_16S, s5_16S, s2_18S, s3_18S, s4_18S, s5_18S]
removex = [s1_16S, s2_16S, s3_16S, s4_16S, s5_16S]
fsl, fst, fsst, fspv = 7, 10, 8, 8
cols_16S, cols_18S = ['#33FFFF', '#33CCFF', '#3399FF', '#3366FF'], ['#CCFF99', '#66FF66', '#009900', '#006400']

for a in removey:
    plt.setp(a.get_yticklabels(), visible=False)

krusk, krusk_p, treat_mean, treat_sd, cont, print_otus, tax = simper(sim_16S, ord_16S, meta_16S, tax_16S, 8)
simper_plot(krusk, krusk_p, treat_mean, treat_sd, cont, print_otus, tax, s16S, [1, 2, 3, 4], [0.75, 5], [0, 40], cols_16S, 0.6, 35, '')
krusk, krusk_p, treat_mean, treat_sd, cont, print_otus, tax = simper(sim_18S, ord_18S, meta_18S, tax_18S, 8)
simper_plot(krusk, krusk_p, treat_mean, treat_sd, cont, print_otus, tax, s18S, [1, 2, 3, 4], [0.75, 5], [0, 90], cols_18S, 0.6, 79, '')

barplot(fn_16S, tax_16S, 1, l16S, 0.8)
barplot(fn_18S, tax_18S, 1, l18S, 0.5)


l16S.set_title('16S rRNA gene', fontsize=14)
l18S.set_title('18S rRNA gene', fontsize=14)
l16S.text(-0.98, 99, 'A', fontsize=16, weight='bold')
l16S.text(-0.98, -16, 'B', fontsize=16, weight='bold')


l16S.set_ylabel('Relative abundance (%)')
l16S.set_xlabel('Days')
l18S.set_xlabel('Days')
s3_18S.set_xlabel('Days')
s1_16S.text(-1.4, 20, '16S rRNA gene', fontsize=10, ha='center', va='center', rotation=90)
s1_18S.text(-1.6, 50, '18S rRNA gene', fontsize=10, ha='center', va='center', rotation=90)
plt.setp(l18S.get_yticklabels(), visible=False)
fig.subplots_adjust(hspace=1, wspace=0.4)
plt.savefig('16S and 18S daily.png', bbox_inches='tight', dpi=600)




