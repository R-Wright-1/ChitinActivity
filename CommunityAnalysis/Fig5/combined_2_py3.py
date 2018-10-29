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
        ax.bar(x, treat_mean[a], yerr=treat_sd[a], color=colors, error_kw=dict(ecolor='gray', lw=1, capsize=3, capthick=1, alpha=0.5),  label=treats, edgecolor='k')
        #x1 = [1.4, 2.4, 3.4, 4.4]
        plt.sca(ax)
        plt.xticks(x, ['+\n        9-day', 'R', '+\n        4-day', 'R'], fontsize=12)
        if ax == ax4:
            otu[a] = 'ASV2\n'+r'$(Spirochaeta)$'+'\n'
        if ax == ax8:
            otu[a] = 'ASV1\n'+r'$(Cafeteria)$'+'\n'
        if ax == ax9:
            otu[a] = 'ASV2\n'+r'$(Cafeteria)$'+'\n'
        title = otu[a]
        title += 'SIMPER: %.0f'%(cont[a]*100)+'%'
        ax.set_title(title, fontsize=14)
        h, p = 'H=%.2f'%krusk[a], '${'+r'p = '+'}$'+'%.3f'%float(krusk_p[a])
        ax.text(xtxt, ytxt, h+', '+p, va='bottom', ha='left', color='#bd0303', fontsize=12)
        ax.tick_params(axis='y',which='both',left='on',right='off')
        ax.tick_params(axis='x',which='both',top='off',bottom='off')
    ylab += '\n Relative abundance (%)'
    axis[0].set_ylabel(ylab, fontsize=14)
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
    
def get_control_samples(fn):
    with open(fn, 'rU') as f:
        rows = []
        for row in csv.reader(f):
            rows.append(row)
    samples = []
    for a in range(len(rows)):
        if a > 0:
            this_row = []
            for b in range(len(rows[a])):
                if b > 0:
                    this_row.append(float(rows[a][b]))
            samples.append(this_row)
    return samples
    
def transform_for_NMDS(df):
    X = df.iloc[0:].values
    y = df.iloc[:,0].values
    n_samples = 14
    seed = np.random.RandomState(seed=3)
    X_true = seed.randint(0, 20, 2 * n_samples).astype(np.float)
    X_true = X_true.reshape((n_samples, 2))
    X_true = X
    similarities = distance.cdist(X_true, X_true, 'braycurtis')
    mds = manifold.MDS(n_components=2, max_iter=3000, eps=1e-9, random_state=seed,
                   dissimilarity="precomputed", n_jobs=1)
    pos = mds.fit(similarities).embedding_
    
    nmds = manifold.MDS(n_components=2, metric=False, max_iter=3000, eps=1e-12,
                        dissimilarity="precomputed", random_state=seed, n_jobs=1,
                        n_init=1)
    npos = nmds.fit_transform(similarities, init=pos)
    print(nmds.stress_)
    # Rescale the data
    pos *= np.sqrt((X_true ** 2).sum()) / np.sqrt((pos ** 2).sum())
    npos *= np.sqrt((X_true ** 2).sum()) / np.sqrt((npos ** 2).sum())
    # Rotate the data
    clf = PCA()
    X_true = clf.fit_transform(X_true)
    pos = clf.fit_transform(pos)
    npos = clf.fit_transform(npos)
    return pos, npos
    
def eigsorted(cov):
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    return vals[order], vecs[:,order]

def plot_nmds(fn, fnc, ax):
    gens, samples = get_gens_and_samples(fn)
    controls = get_control_samples(fnc)
    s = pd.DataFrame(samples)
    #c = pd.DataFrame(controls)
    s = s.transpose()
    #c = c.transpose()
    #s = s.append(c)
    pos, npos = transform_for_NMDS(s)
    colors = ['r', 'b', 'g', 'm']
    shapes = ['o', 's', '*', '^']
    ind, ind1 = 0, 0
    alpha = [1./(len(gens[0])*1.25), 1./(len(gens[1])*1.25), 1./(len(gens[2])*1.25), 1./(len(gens[3])*1.25)]
    all_len = len(gens[0])+len(gens[1])+len(gens[2])+len(gens[3])
    x, y = [[], [], [], []], [[], [], [], []]
    al = 0.15
    labels = ['Positive selection nine-day', 'Positive selection four-day', 'Random selection nine-day', 'Random selection four-day']
    count = 0
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
        count += 1
    """for b in range(len(controls[0])):
        ax.scatter(npos[b+count,0], npos[b+count,1], marker='1', color='k')
        if b == 0:
            ax.scatter(npos[b+count,0], npos[b+count,1], marker='1', color='k', label='Controls')"""
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
    ax.set_xlabel('nMDS 1', fontsize=14)
    ax.set_ylabel('nMDS 2', fontsize=14)
    return colors, labels

##############
#Long vs short vs random vs good last 5
##############
fig = plt.figure(figsize=(14, 12))   
ax1, ax2 = plt.subplot2grid((5, 5), (0, 0), rowspan=2, colspan=2), plt.subplot2grid((5, 5), (0, 2), rowspan=2, colspan=2)
ax3, ax8 = plt.subplot2grid((10,5), (5,0), rowspan=2), plt.subplot2grid((10,5), (8,0), rowspan=2)
ax4, ax5, ax6, ax7 = plt.subplot2grid((10,5), (5,1), sharey=ax3, rowspan=2), plt.subplot2grid((10,5), (5,2), sharey=ax3, rowspan=2), plt.subplot2grid((10,5), (5,3), sharey=ax3, rowspan=2), plt.subplot2grid((10,5), (5,4), sharey=ax3, rowspan=2)
ax9, ax10, ax11, ax12 = plt.subplot2grid((10,5), (8,1), sharey=ax8, rowspan=2), plt.subplot2grid((10,5), (8,2), sharey=ax8, rowspan=2), plt.subplot2grid((10,5), (8,3), sharey=ax8, rowspan=2), plt.subplot2grid((10,5), (8,4), sharey=ax8, rowspan=2)

axis_16S = [ax3, ax4, ax5, ax6, ax7]
axis_18S = [ax8, ax9, ax10, ax11, ax12]

#for ax in axis_16S:
#    plt.setp(ax.get_xticklabels(), visible=False)

fn, order, meta, tax, x, xlim, ylim, colors, treats, xtxt, ytxt, legend, leg_coords, ylab = '16S_all_rg_ls.csv', '16S_13_order.csv', '16S_13_LS_last5_meta.csv', '16S_taxonomy.csv', [1, 2, 3, 4], [0.5, 4.5], [0, 40], ['r', 'g', 'b', 'm'], ['Positive \n nine-day', 'Random \n nine-day', 'Positive \n four-day', 'Random \n four-day'], 0.5, 33, True, [1.05, 1.05], '16S rRNA gene'
plot_simper(fn, tax, 4, axis_16S, x, xlim, ylim, colors, treats, xtxt, ytxt, legend, leg_coords, ylab)
fn, order, meta, tax, x, xlim, ylim, colors, treats, xtxt, ytxt, legend, leg_coords, ylab = '18S_all_rg_ls.csv', '18S_13_order.csv', '18S_13_LS_last5_meta.csv', '18S_taxonomy.csv', [1, 2, 3, 4], [0.5, 4.5], [0, 100], ['r', 'g', 'b', 'm'], ['Positive \n nine-day', 'Random \n nine-day', 'Positive \n four-day', 'Random \n four-day'], 0.5, 82, True, [1.05, 1.05], '18S rRNA gene'
plot_simper(fn, tax, 4, axis_18S, x, xlim, ylim, colors, treats, xtxt, ytxt, legend, leg_coords, ylab)

ax1.text(-80,70, 'A', fontsize=18, weight='bold')
ax1.text(-80, -95, 'B', fontsize=18, weight='bold')

colors, labels = plot_nmds('16S_all.csv', '16S_controls.csv', ax1)
colors, labels = plot_nmds('18S_all.csv', '18S_controls.csv', ax2)
ax2.legend(bbox_to_anchor=(1.9, 1.03), scatterpoints=1, fontsize=14)

ax2.set_ylabel('')
ax1.set_title('16S rRNA gene', fontsize=16)
ax2.set_title('18S rRNA gene', fontsize=16)

plt.subplots_adjust(hspace=0.5, wspace=0.4)
plt.savefig('LS_RG_all_py3.png', bbox_inches='tight', dpi=600)
#plt.close()


