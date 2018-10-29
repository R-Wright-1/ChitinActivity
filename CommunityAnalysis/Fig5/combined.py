import csv
import numpy
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from matplotlib.patches import Ellipse
from sklearn import manifold
import matplotlib as mpl
from scipy.spatial import distance
import matplotlib.patches as mpatches
import statsmodels.stats.multitest as smm
from pylab import *

def get_data(fn, tax, rng):
    krusk, krusk_p, treat_mean, treat_sd, cont, print_otus, tax_name = [], [], [], [], [], [], []
    with open(fn, 'rU') as f:
        rows = []
        for row in csv.reader(f):
            rows.append(row)
    others = []
    for a in range(rng+1):
        a += 1
        krusk.append(float(rows[a][-2]))
        krusk_p.append(float(rows[a][-1]))
        print_otus.append(rows[a][0])
        cont.append(float(rows[a][1]))
        others.append(rows[a][2:-2])
    for a in range(len(others)):
        this_mean, this_sd = [], []
        for b in range(len(others[a])):
            if (b+1) % 2 == 0:
                this_sd.append(float(others[a][b]))
            else:
                this_mean.append(float(others[a][b]))
        treat_mean.append(this_mean)
        treat_sd.append(this_sd)
    with open(tax, 'rU') as f:
        rows = []
        for row in csv.reader(f):
            rows.append(row)
    otus = []
    for a in range(len(print_otus)):
        for b in range(len(rows)):
            if print_otus[a] == rows[b][0]:
                otus.append(rows[b][1:])
    for c in range(len(otus)):
        if otus[c][6] != 'NA':
            new_otu = r'$'+otus[c][5]+'$'
            new_otu2 = r'$'+otus[c][6]+'$'
            otus[c] = new_otu+' '+new_otu2
        elif otus[c][5] != 'NA':
            otus[c] = r'$'+otus[c][5]+'$'
        elif otus[c][4] != 'NA':
            otus[c] = otus[c][4]
        elif otus[c][3] != 'NA':
            otus[c] = otus[c][3]
        elif otus[c][2] != 'NA':
            otus[c] = otus[c][2]
        elif otus[c][1] != 'NA':
            otus[c] = otus[c][1]
        elif otus[c][0] != 'NA':
            otus[c] = otus[c][0]
    for a in range(len(otus)):
        changing = False
        for b in range(len(otus[a])):
            if otus[a][b] == '_':
                changing = True
        new_str, count = '', 0
        if changing:
            for c in range(len(otus[a])):
                if otus[a][c] == '_':
                    count += 1
                    if count > 1:
                        otus[a] = new_str
                        break
                    else:
                        new_str += ' '
                else:
                    new_str += otus[a][c]
            otus[a] = new_str
    for d in range(len(print_otus)):
        print_otus[d] = print_otus[d][3:]
        otu = int(print_otus[d])
        otu = 'ASV'+str(otu)+'\n'+otus[d]+'\n'
        tax_name.append(otu)
    return krusk, krusk_p, treat_mean, treat_sd, cont, tax_name
    
def plot_simper(fn, tax, rng, axis, x, lxim, ylim, colors, treats, xtxt, ytxt, legend, leg_coords, ylab):
    krusk, krusk_p, treat_mean, treat_sd, cont, otu = get_data(fn, tax, rng)
    for a in range(len(axis)):
        ax = axis[a]
        if a > 0:
            plt.setp(ax.get_yticklabels(), visible=False)
        #plt.setp(ax.get_xticklabels(), visible=False)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.bar(x, treat_mean[a], yerr=treat_sd[a], color=colors, error_kw=dict(ecolor='gray', lw=1, capsize=3, capthick=1, alpha=0.5),  label=treats)
        x1 = [1.4, 2.4, 3.4, 4.4]
        plt.sca(ax)
        plt.xticks(x1, ['+\n        9-day', 'R', '+\n        4-day', 'R'], fontsize=8)
        if ax == ax4:
            otu[a] = 'ASV2\n'+r'$(Spirochaeta)$'+'\n'
        if ax == ax8:
            otu[a] = 'ASV1\n'+r'$(Cafeteria)$'+'\n'
        if ax == ax9:
            otu[a] = 'ASV2\n'+r'$(Cafeteria)$'+'\n'
        title = otu[a]
        title += 'SIMPER: %.0f'%(cont[a]*100)+'%'
        ax.set_title(title, fontsize=8)
        h, p = 'H=%.2f'%krusk[a], '${'+r'p = '+'}$'+'%.3f'%float(krusk_p[a])
        ax.text(xtxt, ytxt, h+', '+p, va='bottom', ha='left', color='#bd0303', fontsize=8)
        ax.tick_params(axis='y',which='both',left='on',right='off')
        ax.tick_params(axis='x',which='both',top='off',bottom='off')
    ylab += '\n Relative abundance (%)'
    axis[0].set_ylabel(ylab, fontsize=8)
    return

def get_gens_and_samples(fn):
    with open(fn, 'rU') as f:
        rows = []
        for row in csv.reader(f):
            rows.append(row)
    gens, lg, ind, samples = [[], [], [], []], 0, 0, []
    for a in range(len(rows[0])):
        if a > 0:
            gen = float(rows[0][a][1:-2])
            if gen < lg:
                ind += 1
            lg = gen
            gens[ind].append(gen)
    for b in range(len(rows)):
        if b > 0:
            this_row = []
            for c in range(len(rows[b])):
                if c > 0:
                    this_row.append(float(rows[b][c]))
            samples.append(this_row)
    return gens, samples
    
def transform_for_NMDS(df):
    X = df.iloc[0:].values
    y = df.iloc[:,0].values
    n_samples = 14
    seed = np.random.RandomState(seed=3)
    X_true = seed.randint(0, 20, 2 * n_samples).astype(np.float)
    X_true = X_true.reshape((n_samples, 2))
    X_true = X
    similarities = distance.cdist(X_true, X_true, 'braycurtis')
    mds = manifold.MDS(n_components=3, max_iter=3000, eps=1e-9, random_state=seed,
                   dissimilarity="precomputed", n_jobs=1)
    pos = mds.fit(similarities).embedding_
    
    nmds = manifold.MDS(n_components=3, metric=False, max_iter=3000, eps=1e-12,
                        dissimilarity="precomputed", random_state=seed, n_jobs=1,
                        n_init=1)
    npos = nmds.fit_transform(similarities, init=pos)
    #stress =  nmds.stress_
    #print stress
    # Rescale the data
    pos *= np.sqrt((X_true ** 2).sum()) / np.sqrt((pos ** 2).sum())
    npos *= np.sqrt((X_true ** 2).sum()) / np.sqrt((npos ** 2).sum())
    # Rotate the data
    clf = PCA(n_components=2)
    X_true = clf.fit_transform(X_true)
    pos = clf.fit_transform(pos)
    npos = clf.fit_transform(npos)
    return pos, npos
    
def eigsorted(cov):
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    return vals[order], vecs[:,order]

def plot_nmds(fn, ax):
    gens, samples = get_gens_and_samples(fn)
    s = pd.DataFrame(samples)
    s = s.transpose()
    pos, npos = transform_for_NMDS(s)
    colors = ['r', 'b', 'g', 'm']
    shapes = ['o', 's', '*', '^']
    ind, ind1 = 0, 0
    alpha = [1./(len(gens[0])*1.25), 1./(len(gens[1])*1.25), 1./(len(gens[2])*1.25), 1./(len(gens[3])*1.25)]
    all_len = len(gens[0])+len(gens[1])+len(gens[2])+len(gens[3])
    x, y = [[], [], [], []], [[], [], [], []]
    al = 0.15
    labels = ['Positive selection \n nine-day', 'Positive selection \n four-day', 'Random selection \n nine-day', 'Random selection \n four-day']
    for a in range(len(samples[0])):
        num = int(gens[ind][ind1])
        color = colors[ind]
        shape = shapes[ind]
        x[ind].append(npos[a,0])
        y[ind].append(npos[a,1])
        al += alpha[ind]
        ax.scatter(npos[a,0], npos[a,1], marker=shape, color=color, alpha=al)
        if gens[ind][ind1] == gens[ind][-1]:
            ax.scatter(npos[a,0], npos[a,1], marker=shape, color=color, alpha=al, label=labels[ind])
            ind += 1
            ind1 = -1
            al = 0.15
        ind1 += 1
    for z in range(len(gens)):
        nstd = 2
        x1, y1 = x[z], y[z]
        cov = np.cov(x1, y1)
        vals, vecs = eigsorted(cov)
        theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))
        w, h = 2 * nstd * np.sqrt(vals)
        ell = Ellipse(xy=(np.mean(x1), np.mean(y1)),width=w, height=h,angle=theta, color=colors[z])     
        ell.set_facecolor('none')
        ax.add_artist(ell)
    ax.set_xlabel('nMDS 1')
    ax.set_ylabel('nMDS 2')
    if ax == ax2:
        ax.legend(bbox_to_anchor=(0, 1.03), fontsize=10, scatterpoints=1)
        """
        patch0 = mpatches.Patch(color=colors[0], label=labels[0])
        patch1 = mpatches.Patch(color=colors[1], label=labels[1])
        patch3 = mpatches.Patch(color=colors[2], label=labels[2])
        patch4 = mpatches.Patch(color=colors[3], label=labels[3])
        ax2.legend(handles=[patch0, patch1, patch3, patch4], bbox_to_anchor=(1.7, 1.03), fontsize=10)
        """
    return colors, labels

##############
#Long vs short vs random vs good last 5
##############
fig = plt.figure(figsize=(12, 8.27))   
ax1, ax2 = plt.subplot2grid((4, 5), (0, 0), rowspan=2, colspan=2), plt.subplot2grid((4, 5), (0, 2), rowspan=2, colspan=2)
ax3, ax8 = plt.subplot2grid((4,7), (2,0)), plt.subplot2grid((4,7), (3,0))
ax4, ax5, ax6, ax7 = plt.subplot2grid((4,7), (2,1), sharey=ax3), plt.subplot2grid((4,7), (2,2), sharey=ax3), plt.subplot2grid((4,7), (2,3), sharey=ax3), plt.subplot2grid((4,7), (2,4), sharey=ax3)
ax9, ax10, ax11, ax12 = plt.subplot2grid((4,7), (3,1), sharey=ax8), plt.subplot2grid((4,7), (3,2), sharey=ax8), plt.subplot2grid((4,7), (3,3), sharey=ax8), plt.subplot2grid((4,7), (3,4), sharey=ax8)

ax13, ax14 = plt.subplot2grid((4,7), (2,5), sharey=ax3), plt.subplot2grid((4,7), (2,6), sharey=ax3)
ax15, ax16 = plt.subplot2grid((4,7), (3,5), sharey=ax8), plt.subplot2grid((4,7), (3,6), sharey=ax8)

axis_16S = [ax3, ax4, ax5, ax6, ax7, ax13, ax14]
axis_18S = [ax8, ax9, ax10, ax11, ax12, ax15, ax16]

#for ax in axis_16S:
#    plt.setp(ax.get_xticklabels(), visible=False)

fn, order, meta, tax, x, xlim, ylim, colors, treats, xtxt, ytxt, legend, leg_coords, ylab = '16S_all_rg_ls.csv', '16S_13_order.csv', '16S_13_LS_last5_meta.csv', '16S_taxonomy.csv', [1, 2, 3, 4], [0.75, 5.25], [0, 40], ['r', 'g', 'b', 'm'], ['Positive \n nine-day', 'Random \n nine-day', 'Positive \n four-day', 'Random \n four-day'], 0.8, 33, True, [1.05, 1.05], '16S'
plot_simper(fn, tax, 6, axis_16S, x, xlim, ylim, colors, treats, xtxt, ytxt, legend, leg_coords, ylab)
fn, order, meta, tax, x, xlim, ylim, colors, treats, xtxt, ytxt, legend, leg_coords, ylab = '18S_all_rg_ls.csv', '18S_13_order.csv', '18S_13_LS_last5_meta.csv', '18S_taxonomy.csv', [1, 2, 3, 4], [0.75, 5.25], [0, 100], ['r', 'g', 'b', 'm'], ['Positive \n nine-day', 'Random \n nine-day', 'Positive \n four-day', 'Random \n four-day'], 0.8, 82, True, [1.05, 1.05], '18S'
plot_simper(fn, tax, 6, axis_18S, x, xlim, ylim, colors, treats, xtxt, ytxt, legend, leg_coords, ylab)

ax1.text(-90, 90, 'A', fontsize=14)
ax1.text(-90, -100, 'B', fontsize=14)

colors, labels = plot_nmds('16S_all.csv', ax1)
colors, labels = plot_nmds('18S_all.csv', ax2)
ax2.legend(bbox_to_anchor=(1.65, 1.03), scatterpoints=1)

ax2.set_ylabel('')
ax1.set_title('16S')
ax2.set_title('18S')

fig.subplots_adjust(hspace=0.9, wspace=0.22)
#plt.savefig('LS_RG_all.pdf', bbox_inches='tight')
plt.savefig('LS_RG_all.png', bbox_inches='tight', dpi=600)
#plt.close()

"""
##############
#Long vs short last 5
##############
fig = plt.figure(figsize=(8.27, 8.27))   
ax1, ax2 = plt.subplot2grid((4, 2), (0, 0), rowspan=2), plt.subplot2grid((4, 2), (0, 1), rowspan=2)
ax3, ax8 = plt.subplot2grid((4,5), (2,0)), plt.subplot2grid((4,5), (3,0))
ax4, ax5, ax6, ax7 = plt.subplot2grid((4,5), (2,1), sharey=ax3), plt.subplot2grid((4,5), (2,2), sharey=ax3), plt.subplot2grid((4,5), (2,3), sharey=ax3), plt.subplot2grid((4,5), (2,4), sharey=ax3)
ax9, ax10, ax11, ax12 = plt.subplot2grid((4,5), (3,1), sharey=ax8), plt.subplot2grid((4,5), (3,2), sharey=ax8), plt.subplot2grid((4,5), (3,3), sharey=ax8), plt.subplot2grid((4,5), (3,4), sharey=ax8)

axis_16S = [ax3, ax4, ax5, ax6, ax7]
axis_18S = [ax8, ax9, ax10, ax11, ax12]

fn, order, meta, tax, x, xlim, ylim, colors, treats, xtxt, ytxt, legend, leg_coords, ylab = '16S_last5_g_ls.csv', '16S_13_order.csv', '16S_13_LS_last5_meta.csv', '16S_taxonomy.csv', [1, 2], [0.75, 3.25], [0, 45], ['r', 'b'], ['Positive \n nine-day', 'Positive \n four-day'], 0.8, 37, True, [1.05, 1.05], '16S'
plot_simper(fn, tax, 4, axis_16S, x, xlim, ylim, colors, treats, xtxt, ytxt, legend, leg_coords, ylab)
fn, order, meta, tax, x, xlim, ylim, colors, treats, xtxt, ytxt, legend, leg_coords, ylab = '18S_last5_g_ls.csv', '18S_13_order.csv', '18S_13_LS_last5_meta.csv', '18S_taxonomy.csv', [1, 2], [0.75, 3.25], [0, 55], ['r', 'b'], ['Positive \n nine-day', 'Positive \n four-day'], 0.8, 45, True, [1.05, 1.05], '18S'
plot_simper(fn, tax, 4, axis_18S, x, xlim, ylim, colors, treats, xtxt, ytxt, legend, leg_coords, ylab)

ax1.text(-90, 90, 'A', fontsize=14)
ax1.text(-90, -100, 'B', fontsize=14)

colors, labels = plot_nmds('16S_all.csv', ax1)
colors, labels = plot_nmds('18S_all.csv', ax2)
ax2.legend(bbox_to_anchor=(1.8, 1.03), scatterpoints=1)

ax2.set_ylabel('')
ax1.set_title('16S')
ax2.set_title('18S')

fig.subplots_adjust(hspace=0.9, wspace=0.22)
plt.savefig('LS_last5.png', bbox_inches='tight', dpi=600)
plt.close()


##############
#Long vs short vs random vs good last 5
##############
fig = plt.figure(figsize=(8.27, 8.27))   
ax1, ax2 = plt.subplot2grid((4, 2), (0, 0), rowspan=2), plt.subplot2grid((4, 2), (0, 1), rowspan=2)
ax3, ax8 = plt.subplot2grid((4,5), (2,0)), plt.subplot2grid((4,5), (3,0))
ax4, ax5, ax6, ax7 = plt.subplot2grid((4,5), (2,1), sharey=ax3), plt.subplot2grid((4,5), (2,2), sharey=ax3), plt.subplot2grid((4,5), (2,3), sharey=ax3), plt.subplot2grid((4,5), (2,4), sharey=ax3)
ax9, ax10, ax11, ax12 = plt.subplot2grid((4,5), (3,1), sharey=ax8), plt.subplot2grid((4,5), (3,2), sharey=ax8), plt.subplot2grid((4,5), (3,3), sharey=ax8), plt.subplot2grid((4,5), (3,4), sharey=ax8)

axis_16S = [ax3, ax4, ax5, ax6, ax7]
axis_18S = [ax8, ax9, ax10, ax11, ax12]

fn, order, meta, tax, x, xlim, ylim, colors, treats, xtxt, ytxt, legend, leg_coords, ylab = '16S_last5_rg_ls.csv', '16S_13_order.csv', '16S_13_LS_last5_meta.csv', '16S_taxonomy.csv', [1, 2, 3, 4], [0.75, 5.25], [0, 45], ['r', 'g', 'b', 'm'], ['Positive \n nine-day', 'Random \n nine-day', 'Positive \n four-day', 'Random \n four-day'], 0.8, 37, True, [1.05, 1.05], '16S'
plot_simper(fn, tax, 4, axis_16S, x, xlim, ylim, colors, treats, xtxt, ytxt, legend, leg_coords, ylab)
fn, order, meta, tax, x, xlim, ylim, colors, treats, xtxt, ytxt, legend, leg_coords, ylab = '18S_last5_rg_ls.csv', '18S_13_order.csv', '18S_13_LS_last5_meta.csv', '18S_taxonomy.csv', [1, 2, 3, 4], [0.75, 5.25], [0, 75], ['r', 'g', 'b', 'm'], ['Positive \n nine-day', 'Random \n nine-day', 'Positive \n four-day', 'Random \n four-day'], 0.8, 60, True, [1.05, 1.05], '18S'
plot_simper(fn, tax, 4, axis_18S, x, xlim, ylim, colors, treats, xtxt, ytxt, legend, leg_coords, ylab)

ax1.text(-90, 90, 'A', fontsize=14)
ax1.text(-90, -100, 'B', fontsize=14)

colors, labels = plot_nmds('16S_all.csv', ax1)
colors, labels = plot_nmds('18S_all.csv', ax2)
ax2.legend(bbox_to_anchor=(1.8, 1.03), scatterpoints=1)

ax2.set_ylabel('')
ax1.set_title('16S')
ax2.set_title('18S')

fig.subplots_adjust(hspace=0.9, wspace=0.22)
#plt.savefig('LS_RG_last5.pdf', bbox_inches='tight')
plt.savefig('LS_RG_last5.png', bbox_inches='tight', dpi=600)
plt.close()

##############
#Long vs short
##############
fig = plt.figure(figsize=(8.27, 8.27))   
ax1, ax2 = plt.subplot2grid((4, 3), (0, 0), rowspan=2), plt.subplot2grid((4, 2), (0, 1), rowspan=2)
ax3, ax8 = plt.subplot2grid((4,7), (2,0)), plt.subplot2grid((4,7), (3,0))
ax4, ax5, ax6, ax7, ax13, ax14 = plt.subplot2grid((4,7), (2,1), sharey=ax3), plt.subplot2grid((4,7), (2,2), sharey=ax3), plt.subplot2grid((4,7), (2,3), sharey=ax3), plt.subplot2grid((4,7), (2,4), sharey=ax3), plt.subplot2grid((4,7), (2,5), sharey=ax3), plt.subplot2grid((4,7), (2,6), sharey=ax3)
ax9, ax10, ax11, ax12, ax15, ax16 = plt.subplot2grid((4,7), (3,1), sharey=ax8), plt.subplot2grid((4,7), (3,2), sharey=ax8), plt.subplot2grid((7,5), (3,3), sharey=ax8), plt.subplot2grid((4,7), (3,4), sharey=ax8), plt.subplot2grid((4,7), (3,5), sharey=ax8), plt.subplot2grid((4,7), (3,6), sharey=ax8)

axis_16S = [ax3, ax4, ax5, ax6, ax7, ax13, ax14]
axis_18S = [ax8, ax9, ax10, ax11, ax12, ax15, ax16]

fn, order, meta, tax, x, xlim, ylim, colors, treats, xtxt, ytxt, legend, leg_coords, ylab = '16S_all_g_ls.csv', '16S_13_order.csv', '16S_13_LS_last5_meta.csv', '16S_taxonomy.csv', [1, 2], [0.75, 3.25], [0, 40], ['r', 'b'], ['Positive \n nine-day', 'Positive \n four-day'], 0.8, 32, True, [1.05, 1.05], '16S'
plot_simper(fn, tax, 4, axis_16S, x, xlim, ylim, colors, treats, xtxt, ytxt, legend, leg_coords, ylab)
fn, order, meta, tax, x, xlim, ylim, colors, treats, xtxt, ytxt, legend, leg_coords, ylab = '18S_all_g_ls.csv', '18S_13_order.csv', '18S_13_LS_last5_meta.csv', '18S_taxonomy.csv', [1, 2], [0.75, 3.25], [0, 85], ['r', 'b'], ['Positive \n nine-day', 'Positive \n four-day'], 0.8, 65, True, [1.05, 1.05], '18S'
plot_simper(fn, tax, 4, axis_18S, x, xlim, ylim, colors, treats, xtxt, ytxt, legend, leg_coords, ylab)

ax1.text(-90, 90, 'A', fontsize=14)
ax1.text(-90, -100, 'B', fontsize=14)

colors, labels = plot_nmds('16S_all.csv', ax1)
colors, labels = plot_nmds('18S_all.csv', ax2)
ax2.legend(bbox_to_anchor=(1.8, 1.03), scatterpoints=1)

ax2.set_ylabel('')
ax1.set_title('16S')
ax2.set_title('18S')

fig.subplots_adjust(hspace=0.9, wspace=0.22)
#plt.savefig('LS_all.pdf', bbox_inches='tight')
plt.savefig('LS_all.png', bbox_inches='tight', dpi=600)
plt.close()
"""


