import csv
import numpy
import matplotlib.pyplot as plt
import scipy.stats as stats
from colorsys import hls_to_rgb
import random

fn = 'New_file_chitobiose_GlcNAc.csv'

with open(fn, 'rU') as f:
    rows = []
    for row in csv.reader(f):
        rows.append(row)

labels = []
for a in range(len(rows)):
    if a > 0:
        labels.append(rows[a][0])

for a in range(len(rows[0])):
    if a == 0:
        continue
    if rows[0][a][-1] == '7' or rows[0][a][-1] == '8' or rows[0][a][-1] == '9' or rows[0][a][-1] == '4' or rows[0][a][-1] == '5' or rows[0][a][-1] == '6':
        rows[0][a] = rows[0][a][:-1]+'G'
    elif rows[0][a][-1] == 'B' or rows[0][a][-1] == 'A':
        rows[0][a] = rows[0][a][:-1]+'R'
        
ASVs = ['11/14', '120/130', '361/1387', '15/28', '94/100', '19/440', '504/994', '15/23', '2/12', '35/170', '232/517', '163/607', '10/13', '478/1280', '3/37', '3/41', '10/71', '5/124', '21/23', '1/119', '10/23', '12/22', '13/30', '26/68', '34/51']

tota, totb = 0, 0
for a in range(len(ASVs)):
    n1, n2 = '', ''
    before = True
    for b in ASVs[a]:
        if b == '/':
            before = False
        elif before:
            n1 += b
        else:
            n2 += b
    tota += int(n1)
    totb += int(n2)
#print('Proportion of these ASVs = ', (float(tota)/float(totb))*100)


gens = [[], [], [], []]
gens_nums = [[], [], [], []]
for a in range(len(rows[0])):
    if a == 0 or a == len(rows[0])-1:
        continue
    if rows[0][a][2] != '_':
        gen = int(rows[0][a][1:3])
    else:
        gen = int(rows[0][a][1])
    ls = rows[0][a][0]
    rg = rows[0][a][-1]
    if ls == 'L' and rg == 'G':
        gens[0].append(rows[0][a])
        gens_nums[0].append(gen)
    elif ls == 'S' and rg == 'G':
        gens[1].append(rows[0][a])
        gens_nums[1].append(gen)
    elif ls == 'L' and rg == 'R':
        gens[2].append(rows[0][a])
        gens_nums[2].append(gen)
    elif ls == 'S' and rg == 'R':
        gens[3].append(rows[0][a])
        gens_nums[3].append(gen)

sorted_gens = [[], [], [], []]
for a in range(len(gens)):
    this_sort = sorted(zip(gens_nums[a], gens[a]))
    sort = []
    new_nums = []
    for b in range(len(this_sort)):
        sort.append(this_sort[b][1])
        new_nums.append(this_sort[b][0])
    gens_nums[a] = []
    for c in range(len(sort)):
        adding = True
        for d in range(len(sorted_gens[a])):
            if sort[c] == sorted_gens[a][d]:
                adding = False
        if adding:
            sorted_gens[a].append(sort[c])
            gens_nums[a].append(new_nums[c])

new_rows = []
for e in range(len(rows)):
    del rows[e][0]
    if e > 0:
        new_rows.append(rows[e])
row0 = rows[0] 

means, stds = [[], [], [], [], []], [[], [], [], [], []]
for a in range(len(new_rows)):
    for g in range(len(sorted_gens)):
        this_mean, this_std = [], []
        for h in range(len(sorted_gens[g])):
            this_gen = []
            for i in range(len(row0)):
                if sorted_gens[g][h] == row0[i]:
                    this_gen.append(float(new_rows[a][i]))
            this_mean.append(numpy.mean(this_gen))
            this_std.append(numpy.std(this_gen))
        if g == 1:
            means[g].append(this_mean[0:5]+this_mean[-1])
            means[4].append(this_mean[4:])
            stds[g].append(this_std[0:5]+this_std[-1])
            stds[4].append(this_std[4:])
        else:
            means[g].append(this_mean)
            stds[g].append(this_std)


fig = plt.figure(figsize=(8.27, 6))
ax1 = plt.subplot2grid((3, 10), (0, 0), colspan=10)
ax2 = plt.subplot2grid((3, 10), (1, 7), colspan=3)
ax3 = plt.subplot2grid((3, 10), (2, 7), colspan=3)
ax4, ax5, ax6 = ax1.twinx(), ax2.twinx(), ax3.twinx()
plt.setp(ax4.get_xticklabels(), visible=False)
plt.setp(ax4.get_yticklabels(), visible=False)
plt.setp(ax5.get_xticklabels(), visible=False)
plt.setp(ax5.get_yticklabels(), visible=False)
plt.setp(ax6.get_xticklabels(), visible=False)

def get_distinct_colors(n):
    colors = []
    for i in numpy.arange(0., 360., 360. / n):
        h = i / 360.
        l = (50 + numpy.random.rand() * 10) / 100.
        s = (90 + numpy.random.rand() * 10) / 100.
        colors.append(hls_to_rgb(h, l, s))
    random.shuffle(colors)
    return colors


colors = get_distinct_colors(len(means[0]))

sums = [0, 0, 0, 0, 0]
for a in range(len(means[0])):
    for b in range(len(means[0][a])):
        means[0][a][b] = means[0][a][b]/100
        if b == 17:
            sums[4] += means[0][a][b]
        elif b == 16:
            sums[3] += means[0][a][b]
        elif b == 15:
            sums[2] += means[0][a][b]
        elif b == 14:
            sums[1] += means[0][a][b]
        elif b == 13:
            sums[0] += means[0][a][b]
fdm = numpy.mean(sums)

alpha = [] #2 = alpha, 6 = bacteroidia
col_tot = []
for a in range(len(means[0][0])):
    col = []
    for b in range(len(means[0])):
        col.append(means[0][b][a])
    col_tot.append(col)
for c in range(len(means[0][2])):
    alpha.append(means[0][2][c])
for d in range(len(col_tot)):
    col_tot[d] = sum(col_tot[d])
#print('Alpha % = ', (sum(alpha)/sum(col_tot)) *100)

colors = [(0.9593533900362557, 0.2217496362084278, 0.2217496362084278), (0.3677314152473781, 0.07257738045673634, 0.9949337391774936), (0.2910712046189042, 0.9663161833145264, 0.23235424994971976), (0.970173021651961, 0.9326790193567375, 0.03282296427137166), (0.14345130483543167, 0.9647338305650909, 0.866179927477532), (0.15613635978414375, 0.867808363750623, 0.9648545461096886), (0.07310998002153213, 0.43611722096414357, 0.980628082378061), (0.6329096696683079, 0.17603663887121235, 0.9918813367231692), (0.8184033310250391, 0.9789477059179617, 0.1762258314533498), (0.9941115348280554, 0.7588102045408464, 0.15374964094516574), (0.06626494648907855, 0.9861582191429079, 0.4342222555506105), (0.9756225564280319, 0.07357440741814925, 0.2900659631805207), (0.976242562847172, 0.5167252919629874, 0.09255550345450936), (0.9683787217219271, 0.4053171009589734, 0.22750816808646168), (0.6172688207901094, 0.9664472586192219, 0.17285989991669415), (0.02430300223087567, 0.9789663219941769, 0.6352875268793885), (0.16339823415468147, 0.674111567308188, 0.9613878172070353), (0.34161143266306737, 0.9813837230751649, 0.040542119527962606), (0.2211817907838447, 0.3401855482541831, 0.9649552749734613), (0.26825638944403984, 0.20717148422463894, 0.9707327994671454), (0.9586218644695517, 0.05116984768354205, 0.9223237837981103), (0.9721567816162855, 0.061654180431966354, 0.4986954290004393), (0.9849353619149365, 0.21317537631841443, 0.7688425659479102), (0.17931912357067814, 0.9781093080767069, 0.3071255530916428), (0.8084475222438796, 0.04573310422980392, 0.9991261267473988)]

for a in range(len(labels)):
    labels[a] = labels[a]+' ('+ASVs[a]+')'
data = numpy.array(means[0])
bottom = numpy.cumsum(data, axis=0)
ax1.bar(gens_nums[0], data[1], color=colors[1], width=1.0, label=labels[1], edgecolor='k')
#print(data[0], labels[0])
for j in range(1, data.shape[0]):
    #print(data[j], labels[j])
    if j == 2 or j == 5 or j == 6 or j == 10 or j == 11 or j == 12 or j == 13 or j == 22 or j == 23 or j == 24:
        ax1.bar(gens_nums[0], data[j], color=colors[j], bottom=bottom[j-1], width=1.0, label=labels[j], edgecolor='k')
ax1.plot([15.5, 22], [fdm, fdm], 'k--')
print('Mean gens 16-20 nine-day = ', fdm)
new_plc = []
for a in gens_nums[0]:
    new_plc.append(a)
plt.setp(ax1, xticks=new_plc+[2, 5, 7, 15], xticklabels=gens_nums[0]+['*', '*', '*', '*'])
#plt.setp(ax1, yticks=[0, 0.05, 0.1])
ax1.set_xlim([-0.75, 20.75])
ax1.set_ylim([0, 0.5])
ax1.legend(bbox_to_anchor=(1.05, 1.03), fontsize=8)
"""
data = numpy.array(means[0])
bottom = numpy.cumsum(data, axis=0)
ax1.bar(gens_nums[0], data[2], color=colors[2], width=1.0, label=labels[2])
for j in xrange(3, data.shape[0]):
    if j == 5 or j == 6 or j == 13:
        ax1.bar(gens_nums[0], data[j], color=colors[j], bottom=bottom[j-1], width=1.0, label=labels[j])
ax1.plot([16, 22], [fdm, fdm], 'k--')
print('Mean gens 16-20 nine-day = ', fdm)
new_plc = []
for a in gens_nums[0]:
    new_plc.append(a+0.5)
plt.setp(ax1, xticks=new_plc+[2.5, 5.5, 7.5, 15.5], xticklabels=gens_nums[0]+['*', '*', '*', '*'])
good = [0.032303538824836796, 0.058417920889193974, 0.49210619391892269, 0.2772596471402996, 0.35195691524440936, 0.40619853651538845, 1.0727540795422981, 1.0842478080218705, 0.83568821770169899, 0.75373907371974136, 0.83913402368103318, 0.11779409114726323, 0.10533969173919007, 0.22796557505651693, 0.46866529829456766, 3.0545284350774673, 0.20010810714532598, 0.56149340784676616, 0.37405815803627901, 0.16648247933766447, 0.90354763218748213]
new_good = []
for a in range(len(good)):
    if a != 2 and a != 5 and a!= 7 and a != 15:
        new_good.append(good[a])
good = new_good
gradient, intercept, r_value, p_value, std_err = stats.linregress(good, means[0][13])
print('Long gamma r2 = ', r_value**2)
gradient, intercept, r_value, p_value, std_err = stats.linregress(good, means[0][2])
print('Long apha r2 = ', r_value**2)
ax1.set_xlim([-0.2, 21])
ax1.set_ylim([0, 0.4])
ax1.legend(bbox_to_anchor=(1.45, 1.05), fontsize=8)
"""

gens_n = [16, 17, 18, 19, 20]

sums = [0, 0, 0, 0, 0]
gamma = [] #13 == gamma, 6 = bacteroidia
col_tot = []
for a in range(len(means[1][0])):
    col = []
    for b in range(len(means[1])):
        col.append(means[1][b][a])
    col_tot.append(col)
for c in range(len(means[1][13])):
    gamma.append(means[1][13][c])
for d in range(len(col_tot)):
    col_tot[d] = sum(col_tot[d])
#print('Gamma % = ', (sum(gamma)/sum(col_tot)) *100)

for a in range(len(means[1])):
    means[1][a][-1] = means[4][a][-1]
    for b in range(len(means[1][a])):
        means[1][a][b] = means[1][a][b]/100
        sums[b] += means[1][a][b]
fdm = numpy.mean(sums)

data = numpy.array(means[1])
bottom = numpy.cumsum(data, axis=0)
ax2.bar(gens_n, data[6], color=colors[6], width=1.0, edgecolor='k', label=labels[6])
#print(data[0], labels[0])
for j in range(1, data.shape[0]):
    #print(data[j], labels[j])
    ax2.bar(gens_n, data[j], color=colors[j], bottom=bottom[j-1], width=1.0, edgecolor='k', label=labels[j])
short = [0.45579855376400991, 0.90921753435498598, 0.3131633845834822, 0.19110119009477697, 0.78238427403644228]
gradient, intercept, r_value, p_value, std_err = stats.linregress(short, means[1][13])
#print('Short gamma r2 = ', r_value**2)
gradient, intercept, r_value, p_value, std_err = stats.linregress(short, means[1][2])
#print('Short apha r2 = ', r_value**2)
ax2.plot([14.5, 22], [fdm, fdm], 'k--')
print('Mean gens 16-20 four-day = ', fdm)
new_plc = []
for a in gens_n:
    new_plc.append(a)
plt.setp(ax2, xticks=new_plc, xticklabels=gens_n)
#ax2.legend(bbox_to_anchor=(2.5, 3.1), fontsize=8)
plt.setp(ax2, yticks=[0, 0.5, 1, 1.5])
ax2.set_ylim([0, 1.5])


daily = [1, 2, 3, 4]

sums = [0, 0, 0, 0]
for a in range(len(means[4])):
    for b in range(len(means[4][a])):
        means[4][a][b] = means[4][a][b]/100
        sums[b] += means[4][a][b]

data = numpy.array(means[4])
bottom = numpy.cumsum(data, axis=0)
ax3.bar(daily, data[0], color=colors[0], width=1.0, edgecolor='k')
#print(data[0], labels[0])
for j in range(1, data.shape[0]):
    #print(data[j], labels[j])
    ax3.bar(daily, data[j], color=colors[j], bottom=bottom[j-1], width=1.0, edgecolor='k')
plt.setp(ax3, xticks=[1, 2, 3, 4], xticklabels=[1, 2, 3, 4])
#plt.setp(ax3, yticks=[0, 0.01, 0.02, 0.03, 0.04])
chi_daily = [0.237026614506512, 0.8669166986037341, 0.80789838336188657, 0.10367476695375623]
gradient, intercept, r_value, p_value, std_err = stats.linregress(sums, chi_daily)
#print('Daily r2 = ', r_value**2)
ax6.plot([1, 2, 3, 4], chi_daily, 'o-', color='blue')
#ax6.set_ylim([0, 1])
#ax6.set_yticks([])
ax3.set_ylim([0, 1.5])
ax6.set_ylim([0, 1.2])


ax2.set_xlim([15.25, 20.75])
ax3.set_xlim([0.25, 4.75])


ax1.set_ylabel('Gene copies\n(per bacterium)')
ax2.set_ylabel('Gene copies\n(per bacterium)')
ax3.set_ylabel('Gene copies\n(per bacterium)')
ax4.set_ylabel('Nine-day')
ax5.set_ylabel('Four-day')
ax1.set_xlabel('Generation')
ax2.set_xlabel('Generation')
ax3.set_xlabel('Day')
ax6.set_ylabel(r'Chitinase $\mu$M day$^{-1}$'+'\nDaily')
ax1.set_title('A', loc='left', weight='bold')
ax2.set_title('B', loc='left', weight='bold')
ax3.set_title('C', loc='left', weight='bold')
ax1.set_title(r'$\beta$-N-acetyl-hexosaminidase'+'\nK01207 and K12373')
#ax1.set_title('K01452')

#plt.tight_layout()
plt.subplots_adjust(hspace=1)
plt.savefig('Good '+fn[9:-4]+'.png', bbox_inches='tight', dpi=600)
plt.close()
