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
from matplotlib import pyplot

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
    return pos, npos, nmds.stress_
    
def eigsorted(cov):
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    return vals[order], vecs[:,order]

def plot_nmds(fn, ax):
    gens, samples = get_gens_and_samples(fn)
    s = pd.DataFrame(samples)
    s = s.transpose()
    pos, npos, stress = transform_for_NMDS(s)
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
    linestyles = ['-', '--', '-.', ':', '']
    for z in range(len(gens)):
        nstd = 2
        x1, y1 = x[z], y[z]
        cov = np.cov(x1, y1)
        vals, vecs = eigsorted(cov)
        theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))
        w, h = 2 * nstd * np.sqrt(vals)
        ell = Ellipse(xy=(np.mean(x1), np.mean(y1)),width=w, height=h,angle=theta, color=colors[z], linestyle=linestyles[z], label=labels[z])     
        ell.set_facecolor('none')
        ax.add_artist(ell)
    ax.set_xlabel('nMDS 1')
    ax.set_ylabel('nMDS 2')
    plt.sca(ax)
    plt.annotate('Stress = %.3f'%stress, xy=(0.05, 0.95), xycoords='axes fraction', fontsize=10)
    return colors, labels

##############
#Long vs short vs random vs good last 5
##############
colors = ['r', 'b', 'g', 'm']
linestyles = ['-', '--', '-.', ':', '']
labels = ['Positive selection nine-day', 'Positive selection four-day', 'Random selection nine-day', 'Random selection four-day']

fig = plt.figure(figsize=(14, 10))   
ax1 = plt.subplot(221)
ax2 = plt.subplot(222)
ax3 = plt.subplot(223)
ax4 = plt.subplot(224)

colors, labels = plot_nmds('16S_all_DADA.csv', ax1)
colors, labels = plot_nmds('18S_all_DADA.csv', ax2)
colors, labels = plot_nmds('16S_all_mothur.csv', ax3)
colors, labels = plot_nmds('18S_all_mothur.csv', ax4)
l1 = ax2.plot([0,1], [0,1], color=colors[0], linestyle=linestyles[0], label=labels[0])
l2 = ax2.plot([0,1], [0,1], color=colors[1], linestyle=linestyles[1], label=labels[1])
l3 = ax2.plot([0,1], [0,1], color=colors[2], linestyle=linestyles[2], label=labels[2])
l4 = ax2.plot([0,1], [0,1], color=colors[3], linestyle=linestyles[3], label=labels[3])
l5 = ax2.plot([0,1], [0,1], color='white')

handles, labels = ax2.get_legend_handles_labels()
ax2.legend([handles[4], handles[0], handles[5], handles[1], handles[6], handles[2], handles[7], handles[3]], ['Positive selection', 'nine-day', 'Positive selection', 'four-day', 'Random selection', 'nine-day', 'Random selection', 'four-day'], bbox_to_anchor=(1.05, 1.03), scatterpoints=1, fontsize=14)

ax1.text(-70, 0, 'DADA2 analysis', fontsize=16, rotation=90, ha='center', va='center')
ax3.text(-0.65, 0, 'Mothur analysis', fontsize=16, rotation=90, ha='center', va='center')
ax1.set_title('A', loc='left', fontsize=14, weight='bold')
ax2.set_title('B', loc='left', fontsize=14, weight='bold')
ax3.set_title('C', loc='left', fontsize=14, weight='bold')
ax4.set_title('D', loc='left', fontsize=14, weight='bold')

#ax2.set_ylabel('')
ax1.set_title('16S rRNA gene', fontsize=16)
ax2.set_title('18S rRNA gene', fontsize=16)

#plt.subplots_adjust(hspace=1, wspace=0.4)
plt.savefig('Mothur DADA nMDS.png', bbox_inches='tight', dpi=600)
#plt.close()


