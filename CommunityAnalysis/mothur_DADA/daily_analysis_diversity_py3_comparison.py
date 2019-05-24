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

fn_16S, fn_18S = '16S_all_pres.csv', '18S_all_pres.csv'
tax_16S, tax_18S = '16S_taxonomy_DADA.csv', '18S_taxonomy_DADA.csv'

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
    x = [1, 2, 3, 4, 5]
    if ax == l16S:
        colors = [(0.9784561615118298, 0.8544526316405848, 0.11043145241311547), (0.20015147663099214, 0.9752845415843673, 0.1405258562499635), (0.2355334788858654, 0.9602365332487572, 0.9602365332487575), (0.5692453105690995, 0.9795050492777015, 0.1589855718604971), (0.5824628057375845, 0.18474425574558218, 0.9801813557295874), (0.1533182103913031, 0.6371021121284499, 0.99994003843131), (0.1512473190618882, 0.08807742030502586, 0.9724560029011022), (0.9669797640117426, 0.5223094007170633, 0.18880662824605388), (0.1633168997160368, 0.820010933188215, 0.9991093059533547), (0.14608558739885424, 0.4403676619617287, 0.9700753961749042), (0.21710576449418473, 0.9626100823361772, 0.4833573065806109), (0.11414594323754845, 0.9990088230852705, 0.8093953488321874), (0.7489378889588101, 0.12358724124585263, 0.9990781480439929), (0.9703040330592666, 0.17675506882904968, 0.8569398953120929), (0.9209638926431851, 0.12839534652503937, 0.9819307038830429), (0.9704780729631548, 0.6960071852119865, 0.2019595872598834), (0.9968577762986314, 0.06298771658430391, 0.26310272938023116), (0.2918368133300881, 0.9857160759840947, 0.014285108268485947), (0.99807399645539, 0.02623618072940448, 0.02623618072940448), (0.8978297835866152, 0.9620647138802005, 0.06277568977000492), (0.9579720543054757, 0.3027723777165785, 0.12408155682869748), (0.05840245005662492, 0.9786746791606609, 0.18986991135720163), (0.40706716820400585, 0.1861479381611485, 0.9593652433111485), (0.08413585895436615, 0.21115811584195723, 0.973291657167507), (0.9941953083479688, 0.18928435172215308, 0.7067271095530344), (0.9632776829039245, 0.05428066332342807, 0.4438508145722114), (0.0644938369930973, 0.9829129205955643, 0.5893047419087928), (0.7469534608075091, 0.9633040903295119, 0.20607688700250215)]
        alpha=0.1
        ax.bar(x, data[0], color=colors[0], label=taxa[0], alpha=alpha, edgecolor='k')
        for j in range(1, data.shape[0]):
            ax.bar(x, data[j], color=colors[j], bottom=bottom[j-1], label=taxa[j], alpha=alpha, edgecolor='k')
            """
            if j == 2 or j == 3:
                ax.bar(x, data[j], color=colors[j], bottom=bottom[j-1], label=taxa[j], alpha=1, edgecolor='k', linewidth=2)
                #ax.bar(x[:3], data[j][:3], color=colors[j][:3], bottom=bottom[j-1][:3], label=taxa[j], alpha=1, edgecolor='k', linewidth=2)
                #ax.bar(x[3:], data[j][3:], color=colors[j][3:], bottom=bottom[j-1][3:], alpha=alpha, edgecolor='k')
            else:
                ax.bar(x, data[j], color=colors[j], bottom=bottom[j-1], label=taxa[j], alpha=alpha, edgecolor='k')
            """
    elif ax == l18S:
        colors = [(0.8046628587291401, 0.9847340964016936, 0.08437790803892675), (0.16068100826509213, 0.812725598821964, 0.9757367464611822), (0.6436875010120449, 0.14402242251126474, 0.9767975533458975), (0.1711921518161551, 0.9803335451227438, 0.3330204304774728), (0.15210023365541248, 0.9654221473432107, 0.6400933818680916), (0.4831155384359089, 0.16674931260219694, 0.957664877186477), (0.9807780937856913, 0.07931416136055469, 0.25960694784558197), (0.9837397533803838, 0.09670477839123093, 0.45151876838689137), (0.13182563895300092, 0.4646501037103113, 0.963886800846277), (0.9890448336366396, 0.3479571614890893, 0.18768524345220172), (0.9833610439316122, 0.8165161217478611, 0.14913643301285695), (0.15788061367972028, 0.15788061367972028, 0.9684884132083601), (0.22189305336408016, 0.9721974637542361, 0.5220148175201427), (0.05962894416054154, 0.9749015092853106, 0.05962894416054154), (0.9607155787073374, 0.9607155787073376, 0.17402712858047287), (0.20871927804490387, 0.9655241879001555, 0.9655241879001557), (0.9969160025600691, 0.16784646939971049, 0.996916002560069), (0.9783492220943032, 0.6016522496924305, 0.03660679108962117), (0.16896336774581666, 0.6425806369800509, 0.9583254831362071), (0.983223988424602, 0.07534102345091298, 0.07534102345091298), (0.634346909674582, 0.9700991345369391, 0.13071857238104656), (0.0925429104376172, 0.2697573840925528, 0.9786152787122955), (0.9750248939667724, 0.17051115981659548, 0.8141221471367364), (0.9713266768301306, 0.520553637887629, 0.22003827859262792), (0.07733748452818678, 0.9838852945756328, 0.8025757325661438), (0.8034228927487032, 0.11352326761401266, 0.975897799032376), (0.9880639417025762, 0.09028477322054818, 0.6289522743097649), (0.23962110378416382, 0.9564484109039381, 0.06041427700422031), (0.45009313460979883, 0.9676734612144107, 0.10503958354005771), (0.3258338839894387, 0.16119823022031277, 0.9843764990659452)]
        alpha=0.1
        ax.bar(x, data[0], color=colors[0], label=taxa[0], alpha=0.5, edgecolor='k', linewidth=2)
        for j in range(1, data.shape[0]):
            if j == 0 or j == 1:
                ax.bar(x, data[j], color=colors[j], bottom=bottom[j-1], label=taxa[j], alpha=0.5, edgecolor='k', linewidth=2)
            else:
                ax.bar(x, data[j], color=colors[j], bottom=bottom[j-1], label=taxa[j], alpha=alpha, edgecolor='k')
            #ax.bar(x, data[j], color=colors[j], bottom=bottom[j-1], label=taxa[j], alpha=alpha, edgecolor='k')
        ax.set_ylim([0, 100])
    if ax == l16S:
        ax.legend(bbox_to_anchor=(2.13, 1.01), fontsize=7)
    elif ax == l18S:
        ax.legend(bbox_to_anchor=(1.1, 1.01), fontsize=7)
    ax.set_xlim([0.5, 5.5])
    plt.sca(ax)
    #plt.xticks(x)
    plt.setp(ax, xticks=[1, 2, 3, 4, 5], xticklabels=[0, 1, 2, 3, 4])
    return


fig = plt.figure(figsize=(10, 7))
h1, w, ss1, ss2 = 10, 10, 23, 30
h = 35
rsb, rss, rssimp = 10, 5, 1
l16S = plt.subplot2grid((h1,12), (0,0), colspan=3, rowspan=rsb)
l18S = plt.subplot2grid((h1,12), (0,6), colspan=3, rowspan=rsb, sharey=l16S, sharex=l16S)

fsl, fst, fsst, fspv = 6, 9, 8, 8
cols_16S, cols_18S = ['#33FFFF', '#33CCFF', '#3399FF', '#3366FF'], ['#CCFF99', '#66FF66', '#009900', '#006400']

barplot(fn_16S, tax_16S, 1, l16S, 0.8)
barplot(fn_18S, tax_18S, 1, l18S, 0.5)

l16S.set_title('16S rRNA gene', fontsize=16)
l18S.set_title('18S rRNA gene', fontsize=16)

l16S.set_ylabel('Relative abundance (%)', fontsize=12)
l16S.set_xlabel('Days', fontsize=12)
l18S.set_xlabel('Days', fontsize=12)
plt.setp(l18S.get_yticklabels(), visible=False)

fig.subplots_adjust(hspace=200, wspace=0.4)
plt.savefig('16S and 18S daily presentation 4.png', bbox_inches='tight', dpi=600)