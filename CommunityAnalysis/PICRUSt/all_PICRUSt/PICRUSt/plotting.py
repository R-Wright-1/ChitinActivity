import csv
import numpy
import matplotlib.pyplot as plt
import scipy.stats as stats

fn = 'chitinase_abs.csv'
fn = 'chitinase_percent.csv'

with open(fn, 'rU') as f:
    rows = []
    for row in csv.reader(f):
        rows.append(row)

for a in range(len(rows[0])):
    if a == 0:
        continue
    if rows[0][a][-1] == '7' or rows[0][a][-1] == '8' or rows[0][a][-1] == '9' or rows[0][a][-1] == '4' or rows[0][a][-1] == '5' or rows[0][a][-1] == '6':
        rows[0][a] = rows[0][a][:-1]+'G'
    elif rows[0][a][-1] == 'B' or rows[0][a][-1] == 'A':
        rows[0][a] = rows[0][a][:-1]+'R'

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

added_row2, row0 = [], []
for e in range(len(rows[0])):
    if e == 0 or e == len(rows[0])-1:
        continue
    this_col = []
    total = 0
    """
    for f in range(len(rows)):
        if f > 0 and f < len(rows)-1:
            this_col.append(float(rows[f][e]))
        elif f > 0:
            total = float(rows[f][e])
    perc = (sum(this_col)/total)*100
    added_row2.append(perc)
    """
    for f in range(len(rows)):
        if f > 0:
            this_col.append(float(rows[f][e]))
    added_row2.append(sum(this_col))
    row0.append(rows[0][e])

means, stds = [[], [], [], []], [[], [], [], []]
for g in range(len(sorted_gens)):
    for h in range(len(sorted_gens[g])):
        this_gen = []
        for i in range(len(row0)):
            if sorted_gens[g][h] == row0[i]:
                this_gen.append(added_row2[i])
        means[g].append(numpy.mean(this_gen))
        stds[g].append(numpy.std(this_gen))

"""""""""""""""""""""""
fig = plt.figure(figsize=(10, 6))
ax1 = plt.subplot2grid((2,5), (0,0), colspan=3)
ax2 = plt.subplot2grid((2,5), (0,3))
ax3 = plt.subplot2grid((2,5), (1,0), colspan=3, sharex=ax1)
ax4 = plt.subplot2grid((2,5), (1,3), sharex=ax2)
ax5 = plt.subplot2grid((2,5), (0,4))
ax6 = plt.subplot2grid((2,5), (1,4), sharex=ax5)

ax1.bar(gens_nums[0], means[0], yerr=stds[0], error_kw=dict(ecolor='gray', lw=1, capsize=3, capthick=1, alpha=0.5))
new_plc = []
for a in gens_nums[0]:
    new_plc.append(a+0.5)
plt.setp(ax1, xticks=new_plc+[2.5, 5.5, 7.5, 15.5], xticklabels=gens_nums[0]+['*', '*', '*', '*'])
ax1.set_xlim([-0.25, 21])
ax3.bar(gens_nums[2], means[2], yerr=stds[2], error_kw=dict(ecolor='gray', lw=1, capsize=3, capthick=1, alpha=0.5))
new_plc = []
for a in gens_nums[2]:
    new_plc.append(a+0.5)
plt.setp(ax3, xticks=new_plc+[4.5, 7.5, 15.5, 17.5], xticklabels=gens_nums[2]+['*', '*', '*', '*'])
short_gens, daily = [[], []], [[], []]
short_stds, daily_stds = [[], []], [[], []]
short_names, daily_names = [[], []], [[], []]
for j in range(len(sorted_gens[1])):
    if len(sorted_gens[1][j]) > 5:
        daily[0].append(means[1][j])
        daily_stds[0].append(means[1][j])
        daily_names[0].append(gens_nums[1][j])
    else:
        short_gens[0].append(means[1][j])
        short_stds[0].append(means[1][j])
        short_names[0].append(gens_nums[1][j])
    if sorted_gens[1][j][:5] == 'S20_4':
        short_gens[0].append(means[1][j])
        short_stds[0].append(means[1][j])
        short_names[0].append(gens_nums[1][j])

for j in range(len(sorted_gens[3])):
    if len(sorted_gens[3][j]) > 5:
        daily[1].append(means[3][j])
        daily_stds[1].append(means[3][j])
        daily_names[1].append(gens_nums[3][j])
    else:
        short_gens[1].append(means[3][j])
        short_stds[1].append(means[3][j])
        short_names[1].append(gens_nums[3][j])
    if sorted_gens[3][j][:5] == 'S20_4':
        short_gens[1].append(means[3][j])
        short_stds[1].append(means[3][j])
        short_names[1].append(gens_nums[3][j])

ax2.bar(short_names[0], short_gens[0], yerr=short_stds[0], color='g', error_kw=dict(ecolor='gray', lw=1, capsize=3, capthick=1, alpha=0.5))
new_plc = []
for a in short_names[0]:
    new_plc.append(a+0.5)
plt.setp(ax2, xticks=new_plc, xticklabels=short_names[0])
ax4.bar(short_names[1], short_gens[1], yerr=short_stds[1], color='g', error_kw=dict(ecolor='gray', lw=1, capsize=3, capthick=1, alpha=0.5))
plt.setp(ax4, xticks=new_plc, xticklabels=short_names[0])
ax2.set_xlim([15.75, 21])
ax5.bar([1, 2, 3, 4], daily[0], yerr=daily_stds[0], color='r', error_kw=dict(ecolor='gray', lw=1, capsize=3, capthick=1, alpha=0.5))
plt.setp(ax5, xticks=[1.5, 2.5, 3.5, 4.5], xticklabels=[1, 2, 3, 4])
ax6.bar([1, 2, 3, 4], daily[1], yerr=daily_stds[1], color='r', error_kw=dict(ecolor='gray', lw=1, capsize=3, capthick=1, alpha=0.5))
plt.setp(ax5, xticks=[1.5, 2.5, 3.5, 4.5], xticklabels=[1, 2, 3, 4])
ax5.set_xlim([0.75, 5])

ax1.set_ylabel('Chitinase gene copies')
ax3.set_ylabel('Chitinase gene copies')
ax3.set_xlabel('Generation')
ax4.set_xlabel('Generation')
ax6.set_xlabel('Day')
ax1.set_title('Nine-day')
ax2.set_title('Four-day')
ax5.set_title('Generation 20 \ndaily')

random = [0.091551722965788063, 0.05578774456113262, 0.19773343859446771, 0.4028983968586955, 0.11320979376189184, 0.0834252413280238, 0.82712022135273211, 0.35369461961865095, 0.13607574184682639, 0.58081616109583611, 0.67705765994884948, 0.15636653214928983, 0.14229152585045449, 0.12877847544165252, 0.13825991526193623, 0.9362907669500059, 0.18202097226894032, 0.30241137752641739, 0.15505295713128878, 0.13059046071159483, 0.92241478375389363]
good = [0.032303538824836796, 0.058417920889193974, 0.49210619391892269, 0.2772596471402996, 0.35195691524440936, 0.40619853651538845, 1.0727540795422981, 1.0842478080218705, 0.83568821770169899, 0.75373907371974136, 0.83913402368103318, 0.11779409114726323, 0.10533969173919007, 0.22796557505651693, 0.46866529829456766, 3.0545284350774673, 0.20010810714532598, 0.56149340784676616, 0.37405815803627901, 0.16648247933766447, 0.90354763218748213]
short_random = [0.9362907669500059, 0.12497142569080828, 0.076422349721072078, 0.047963580738016376, 0.050683583615973103, 0.22039840838767263]
short = [3.0545284350774673, 0.45579855376400991, 0.90921753435498598, 0.3131633845834822, 0.19110119009477697, 0.78238427403644228]
del short[0]
del short_random[0]
chi_daily = [0.237026614506512, 0.8669166986037341, 0.80789838336188657, 0.10367476695375623]
#chi_daily_error = [0.059736179903026806, 0.22382194238678188, 0.1198038210570813, 0.042979078330859642]
chi_daily_r = [0.060815516547624314, 0.20574092459971904, 0.28381263064226919, 0.10472441680031566]
#chi_daily_r_error = [0.023864665241773816, 0.0903933660773827, 0.075946470829352422, 0.042281604809644227]

gens, gens_s, gens_daily = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20], [16, 17, 18, 19, 20], [1, 2, 3, 4]
a1, a2, a3, a4, a5, a6 = ax1.twinx(), ax2.twinx(), ax3.twinx(), ax4.twinx(), ax5.twinx(), ax6.twinx()
a1.plot(gens, good, 'o-', color='gray')
a2.plot(gens_s, short, 'o-', color='gray')
a3.plot(gens, random, 'o-', color='gray')
a4.plot(gens_s, short_random, 'o-', color='gray')
a5.plot(gens_daily, chi_daily, 'o-', color='gray')
a6.plot(gens_daily, chi_daily_r, 'o-', color='gray')

plt.tight_layout()

plt.savefig(fn[:-4]+'.png', bbox_inches='tight', dpi=600)
plt.close()
"""""""""""""""""""""""""""""
fig = plt.figure(figsize=(8.27, 6))
ax1 = plt.subplot2grid((3, 10), (0, 0), colspan=10)
ax2 = plt.subplot2grid((3, 10), (1, 7), colspan=3)
ax3 = plt.subplot2grid((3, 10), (2, 7), colspan=3)
ax4, ax5, ax6 = ax1.twinx(), ax2.twinx(), ax3.twinx()
plt.setp(ax4.get_xticklabels(), visible=False)
plt.setp(ax4.get_yticklabels(), visible=False)
plt.setp(ax5.get_xticklabels(), visible=False)
plt.setp(ax5.get_yticklabels(), visible=False)
#plt.setp(ax6.get_yticklabels(), visible=False)
plt.setp(ax6.get_xticklabels(), visible=False)
ax1.bar(gens_nums[0], means[0], yerr=stds[0], error_kw=dict(ecolor='gray', lw=1, capsize=3, capthick=1, alpha=0.5))
new_plc = []
for a in gens_nums[0]:
    new_plc.append(a+0.4)
plt.setp(ax1, xticks=new_plc+[2.5, 5.5, 7.5, 15.5], xticklabels=gens_nums[0]+['*', '*', '*', '*'])
ax1.set_xlim([-0.2, 21])
ax1.plot([16, 21], [numpy.mean(means[0][-5:]), numpy.mean(means[0][-5:])], '--', color='gray')
print numpy.mean(means[0][-5:])

short_gens, daily = [[], []], [[], []]
short_stds, daily_stds = [[], []], [[], []]
short_names, daily_names = [[], []], [[], []]
for j in range(len(sorted_gens[1])):
    if len(sorted_gens[1][j]) > 5:
        daily[0].append(means[1][j])
        daily_stds[0].append(means[1][j])
        daily_names[0].append(gens_nums[1][j])
    else:
        short_gens[0].append(means[1][j])
        short_stds[0].append(means[1][j])
        short_names[0].append(gens_nums[1][j])
    if sorted_gens[1][j][:5] == 'S20_4':
        short_gens[0].append(means[1][j])
        short_stds[0].append(means[1][j])
        short_names[0].append(gens_nums[1][j])

for j in range(len(sorted_gens[3])):
    if len(sorted_gens[3][j]) > 5:
        daily[1].append(means[3][j])
        daily_stds[1].append(means[3][j])
        daily_names[1].append(gens_nums[3][j])
    else:
        short_gens[1].append(means[3][j])
        short_stds[1].append(means[3][j])
        short_names[1].append(gens_nums[3][j])
    if sorted_gens[3][j][:5] == 'S20_4':
        short_gens[1].append(means[3][j])
        short_stds[1].append(means[3][j])
        short_names[1].append(gens_nums[3][j])


ax2.bar(short_names[0], short_gens[0], yerr=short_stds[0], color='g', error_kw=dict(ecolor='gray', lw=1, capsize=3, capthick=1, alpha=0.5))
new_plc = []
for a in short_names[0]:
    new_plc.append(a+0.4)
plt.setp(ax2, xticks=new_plc, xticklabels=short_names[0])
ax2.plot([16, 21], [numpy.mean(short_gens[0][-5:]), numpy.mean(short_gens[0][-5:])], '--', color='gray')
print numpy.mean(short_gens[0][-5:])
ax3.bar([1, 2, 3, 4], daily[0], yerr=daily_stds[0], color='r', error_kw=dict(ecolor='gray', lw=1, capsize=3, capthick=1, alpha=0.5))
plt.setp(ax3, xticks=[1.4, 2.4, 3.4, 4.4], xticklabels=[1, 2, 3, 4])
chi_daily = [0.237026614506512, 0.8669166986037341, 0.80789838336188657, 0.10367476695375623]
ax6.plot([1.4, 2.4, 3.4, 4.4], chi_daily, 'o-', color='blue')
print numpy.mean(daily[0])
gradient, intercept, r_value, p_value, std_err = stats.linregress(daily[0], chi_daily)
print r_value**2

ax2.set_xlim([15.8, 21])
ax3.set_xlim([0.8, 5])


ax1.set_ylabel('Chitinase gene copies')
ax2.set_ylabel('Chitinase gene copies')
ax3.set_ylabel('Chitinase gene copies')
ax4.set_ylabel('Nine-day')
ax5.set_ylabel('Four-day')
ax1.set_xlabel('Generation')
ax2.set_xlabel('Generation')
ax3.set_xlabel('Day')
ax6.set_ylabel(r'Chitinase $\mu$M day$^{-1}$'+'\nDaily')
ax1.set_title('A', loc='left')
ax2.set_title('B', loc='left')
ax3.set_title('C', loc='left')

plt.tight_layout()

plt.savefig('Good only.png', bbox_inches='tight', dpi=600)
plt.close()

plt.plot([1, 2, 3, 4], chi_daily, 'o-', color='blue')
plt.ylabel(r'Chitinase $\mu$M day$^{-1}$')
plt.xlabel('Day')
plt.xlim([0, 5])
plt.ylim([0, 1.2])
plt.savefig('Daily.png', dpi=600)