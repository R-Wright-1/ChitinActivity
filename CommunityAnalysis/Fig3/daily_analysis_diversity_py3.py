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
fn_16S, fn_18S = '16S_all.csv', '18S_all.csv'
sim_16S, sim_18S = '16S_simper.csv', '18S_simper.csv'
tax_16S, tax_18S = '16S_taxonomy.csv', '18S_taxonomy.csv'
simp_16S, simp_18S = '16S_no_percent.csv', '18S_no_percent.csv'

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

def barplot(fn, tax, lim, ax, alpha, colors):
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
    #colors = get_distinct_colors(len(samples))
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
rsb, rss, rssimp = 9, 5, 1
l16S = plt.subplot2grid((h1,12), (1,0), colspan=3, rowspan=rsb)
l18S = plt.subplot2grid((h1,12), (1,6), colspan=3, rowspan=rsb, sharey=l16S, sharex=l16S)
s1_16S = plt.subplot2grid((h,w), (ss1,0), colspan=2, rowspan=rss)
s2_16S, s3_16S, s4_16S, s5_16S = plt.subplot2grid((h,w), (ss1,2), sharey=s1_16S, colspan=2, rowspan=rss), plt.subplot2grid((h,w), (ss1,4), rowspan=rss, sharey=s1_16S, colspan=2), plt.subplot2grid((h,w), (ss1,6), rowspan=rss, sharey=s1_16S, colspan=2), plt.subplot2grid((h,w), (ss1,8), rowspan=rss, sharey=s1_16S, colspan=2)
s1_18S = plt.subplot2grid((h,w), (ss2,0), colspan=2, rowspan=rss)
s2_18S, s3_18S, s4_18S, s5_18S = plt.subplot2grid((h,w), (ss2,2), rowspan=rss, sharey=s1_18S, colspan=2), plt.subplot2grid((h,w), (ss2,4), rowspan=rss, sharey=s1_18S, colspan=2), plt.subplot2grid((h,w), (ss2,6), rowspan=rss, sharey=s1_18S, colspan=2), plt.subplot2grid((h,w), (ss2,8), rowspan=rss, sharey=s1_18S, colspan=2)

colors_16S = [(0.9641834781270211, 0.15251184998598855, 0.3398206872493026), (0.07936810482295042, 0.9547788788709628, 0.5507431370026498), (0.6695520108781028, 0.9671073697556801, 0.19346343667397947), (0.6650056719776007, 0.17693929562032362, 0.9700471572009001), (0.11962506085563329, 0.9632832564974338, 0.3792121979761872), (0.8396842995562442, 0.9832236409716707, 0.05021792177139728), (0.9583864993175466, 0.077270282764492, 0.077270282764492), (0.174811339798026, 0.41859784199869243, 0.9671174719501932), (0.9736738684794177, 0.046334735073729094, 0.6883387505084362), (0.19506354058937758, 0.958039822707026, 0.9580398227070263), (0.9700812808792767, 0.0952440986914207, 0.9027861130186713), (0.9800216148402335, 0.7426103532548091, 0.20843501468760361), (0.20125061364286978, 0.058901855425394056, 0.9841687838389895), (0.9757051028270705, 0.9054562590683243, 0.06247013396336909), (0.040180968086908075, 0.7567656720614202, 0.9717410832537747), (0.9699665571401991, 0.18057241108798727, 0.5449081708043919), (0.1806756403420965, 0.9946432508641482, 0.24328853345917745), (0.48051420626970676, 0.17512682594624385, 0.9691340147872475), (0.04092623822488106, 0.9751092949163293, 0.7595285895259953), (0.12715303798476396, 0.19175980882978919, 0.967041058970095), (0.9682015785830554, 0.37629040270014685, 0.19871704993527428), (0.9547810109210887, 0.47695963594877117, 0.06739845740107031), (0.39223849092991037, 0.979596740817666, 0.02513958475006417), (0.057475826073571934, 0.559921195280077, 0.9905886545999394), (0.8379266654690272, 0.11793289054653644, 0.9688346245458439), (0.31914190281009075, 0.9873667064494908, 0.19764648396656337)]
colors_18S = [(0.11485896106865778, 0.993404093226795, 0.8957879674314467), (0.5513770486254422, 0.21907372141532422, 0.9667562076380887), (0.9875951951552691, 0.07092754627059095, 0.6820393121937093), (0.05173475942148553, 0.9937995437969563, 0.6797779490051327), (0.248073886500512, 0.95632713442833, 0.045715815663992965), (0.966922745003283, 0.1608769715353252, 0.1608769715353252), (0.8769101958295187, 0.9827835454318841, 0.029923399010596707), (0.7216943514799523, 0.2293452267778503, 0.9678689138310029), (0.9991486334463409, 0.19301415443741687, 0.5512961451080487), (0.9792110661501915, 0.18338565312009247, 0.8907860202579577), (0.5576206886062759, 0.9609278631687087, 0.23497494895632964), (0.11077851453765442, 0.11077851453765442, 0.9583711458955805), (0.19380255199668805, 0.5423930710824111, 0.9781312199395661), (0.03184113462580318, 0.8891938596057141, 0.9963629502282036), (0.8856167255204997, 0.1646679459730228, 0.9757353229639344), (0.9671273979159809, 0.14323233315331474, 0.3263201253227953), (0.04472396738130713, 0.661380354802624, 0.969708548513283), (0.05884502108045253, 0.2640836067636504, 0.9824186566548445), (0.243765447736264, 0.035002670847713246, 0.9744351668461931), (0.9713732045481511, 0.5549394809048407, 0.2217925019901923), (0.17913783460386778, 0.9805912013166912, 0.17913783460386778), (0.20156954874985233, 0.9662510394347255, 0.5414279890542406), (0.96829352082177, 0.6823994101842377, 0.11061118890917299), (0.14291686560043826, 0.9969595475082869, 0.3327041282466271), (0.9808262603802382, 0.3297928068929199, 0.1437832487536861), (0.9657878715691559, 0.8758723798464674, 0.15654844606495777), (0.6793775810358893, 0.9914421799709516, 0.055248383165765746)]
barplot(fn_16S, tax_16S, 1, l16S, 0.8, colors_16S)
barplot(fn_18S, tax_18S, 1, l18S, 0.5, colors_18S)

ax_16S = plt.subplot2grid((h1, 12), (0,0), colspan=3, rowspan=1)
ax_18S = plt.subplot2grid((h1, 12), (0,6), colspan=3, rowspan=1)
axcol = plt.subplot2grid((h1*2, 12), (1,9), colspan=2, rowspan=1)

#ax_16S_divider = make_axes_locatable(l16S)
#ax_16S = ax_16S_divider.append_axes("top", size="7%", pad="4%")
#ax_18S_divider = make_axes_locatable(l18S)
#ax_18S = ax_18S_divider.append_axes("top", size="7%", pad="4%")
#lblank = plt.subplot2grid((h1,12), (0,9), colspan=2, rowspan=rsb, frameon=False)
#axcol_divider = make_axes_locatable(lblank)
#axcol = axcol_divider.append_axes("top", size="7%", pad="60%")

s16S, s18S = [s1_16S, s2_16S, s3_16S, s4_16S, s5_16S], [s1_18S, s2_18S, s3_18S, s4_18S, s5_18S]
removey = [s2_16S, s3_16S, s4_16S, s5_16S, s2_18S, s3_18S, s4_18S, s5_18S, ax_16S, ax_18S]
removex = [s1_16S, s2_16S, s3_16S, s4_16S, s5_16S, ax_16S, ax_18S]
fsl, fst, fsst, fspv = 6, 9, 8, 8
cols_16S, cols_18S = ['#33FFFF', '#33CCFF', '#3399FF', '#3366FF'], ['#CCFF99', '#66FF66', '#009900', '#006400']
                      

for a in removey:
    plt.setp(a.get_yticklabels(), visible=False)
#plt.setp(lblank.get_xticklabels(), visible=False)
#lblank.tick_params(axis='x',which='both',top='off', bottom='off')
#lblank.tick_params(axis='y',which='both',left='off',right='off')

mi, ma, diversity = get_diversity_plot(simp_16S, ax_16S, simp_18S, ax_18S, 'Simpsons')

simper(sim_16S, tax_16S, 5, s16S, [1, 2, 3, 4], [0.5, 4.5], [0, 40], cols_16S, 0.6, 35, '16S')
simper(sim_18S, tax_18S, 5, s18S, [1, 2, 3, 4], [0.5, 4.5], [0, 100], cols_18S, 0.6, 86, '18S')

ax_16S.set_title('16S rRNA gene')
ax_18S.set_title('18S rRNA gene')
lA = l16S.text(-0.98, 109, 'A', fontsize=16, fontweight='bold', color='k')#, bbox=dict(facecolor='white', edgecolor='white'))
lB = l16S.text(-0.98, 99, 'B', fontsize=16, fontweight='bold', color='k')#, bbox=dict(facecolor='white', edgecolor='white'))
lC = l16S.text(-0.98, -18, 'C', fontsize=16, fontweight='bold', color='k')#, bbox=dict(facecolor='white', edgecolor='white'))


l16S.set_ylabel('Relative abundance (%)')
l16S.set_xlabel('Days')
l18S.set_xlabel('Days')
s3_18S.set_xlabel('Days')
s1_16S.text(-1.4, 20, '16S rRNA gene', fontsize=10, ha='center', va='center', rotation=90)
s1_18S.text(-1.6, 50, '18S rRNA gene', fontsize=10, ha='center', va='center', rotation=90)
plt.setp(l18S.get_yticklabels(), visible=False)

axcol.tick_params(axis='x',which='both',top='off', bottom='off')
axcol.tick_params(axis='y',which='both',left='off',right='off')
ax_16S.tick_params(axis='y',which='both',left='off',right='off')
ax_18S.tick_params(axis='y',which='both',left='off',right='off')
plt.setp(axcol.get_xticklabels(), visible=False)
plt.setp(axcol.get_yticklabels(), visible=False)
cmap = matplotlib.cm.Blues
norm = matplotlib.colors.Normalize(vmin=mi, vmax=ma)
cb1 = matplotlib.colorbar.ColorbarBase(axcol, cmap=cmap, norm=norm, orientation='horizontal')
cb1.set_ticks([])
cb1.ax.text(0.05, 0.35, 'Low', fontsize=8, color='gray')
cb1.ax.text(0.7, 0.35, 'High', fontsize=8, color='w')

fig.subplots_adjust(hspace=200, wspace=0.4)
plt.savefig('16S and 18S daily '+diversity+'p3.png', bbox_inches='tight', dpi=600)

#plt.close()

