import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.metrics import euclidean_distances
from sklearn.decomposition import PCA
from matplotlib.patches import Ellipse
import csv
from sklearn import manifold
import matplotlib as mpl
from scipy.spatial import distance
import numpy
import os


fn = '16S_percent_grouped_new.csv'
fn2 = '18S_percent_grouped_new.csv'

def get_groups(fn):
    count = 0
    with open(fn, 'rU') as f:
        rows = []
        for row in csv.reader(f):
            rows.append(row)
    OTU, numbers = [], []
    for a in range(len(rows)):
        if a > 0:
            count += 1
            this_row = []
            for b in range(len(rows[a])):
                if b == 0:
                    OTU.append(rows[a][b])
                elif b < 23:
                    this_row.append(float(rows[a][b]))
            numbers.append(this_row)
    new_rows = rows[0]
    del new_rows[0]
    G, SG = [], []
    if fn == '16S_percent_grouped_new.csv':
        end = 17
    elif fn == '18S_percent_grouped_new.csv':
        end = 11
    for a in range(len(numbers)):
        G.append(numbers[a][0:end])
        SG.append(numbers[a][end:])
    return G, SG

    
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
    stress =  nmds.stress_
    print stress
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
    
def NMDS_plot(fn, ax, legend, title, lp):
    os.chdir('/Users/u1560915/Documents/GitHub/ChitinActivity/MiSeq/16S/Fig 4/')
    s, maps, el_cols, labels, shape, lw = 25, ['Reds', 'Purples', 'Blues'], ['#B31204', '#7D3C98', '#036197'], ['Selection \n (9 day)', 'Shortened \n selection \n (4 day)'], ['o', '^', 's'], 1
    #16S
    G, SG = get_groups(fn)
    G, SG = pd.DataFrame(G), pd.DataFrame(SG)
    G, SG = G.transpose(), SG.transpose()
    lengths = [len(G), len(SG)]
    result = pd.concat([G, SG])
    pos, npos = transform_for_NMDS(result)
    x, y, this_x, this_y = [], [], [], []
    S16 = False
    if fn == '16S_percent_grouped_new.csv':
        S16 = True
    if S16:
        gens = [0, 1, 3, 4, 6, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20]
        gens_s = [16, 17, 18, 19, 20]
    else:
        gens = [0, 1, 8, 9, 10, 12, 13, 14, 17, 19, 20]
        gens_s = [16, 17, 18, 19, 20]
    for a in range(len(G)):
        cmap = maps[0]
        norm = mpl.colors.Normalize(vmin=-3, vmax=lengths[0])
        colormap = mpl.cm.get_cmap(cmap, 256)
        m = mpl.cm.ScalarMappable(norm=norm, cmap=colormap)
        color = m.to_rgba(a)
        #16S 15, 18S 10
        if S16:
            tx = 15
        else:
            tx = 10
        if a == tx:
            ax.scatter(npos[a,0], npos[a,1], color=color, s=s, lw=lw, marker=shape[0], label=labels[0])
            ax.text((npos[a,0]+0.02), (npos[a,1]-0.01), gens[a], color=color)
        else:
            ax.text((npos[a,0]+0.02), (npos[a,1]-0.01), gens[a], color=color)
            ax.scatter(npos[a,0], npos[a,1], color=color, s=s, lw=lw, marker=shape[0])
        this_x.append(npos[a,0]), this_y.append(npos[a,1])
    x.append(this_x), y.append(this_y)
    this_x, this_y = [], []
    for b in range(len(SG)):
        c = b
        b += a+1
        cmap = maps[1]
        if S16:
            mn, mx = 14, 21
        else:
            mn, mx = 9, 15
        norm = mpl.colors.Normalize(vmin=mn, vmax=mx)
        colormap = mpl.cm.get_cmap(cmap, 256)
        m = mpl.cm.ScalarMappable(norm=norm, cmap=colormap)
        color = m.to_rgba(b)
        ax.scatter(npos[b,0], npos[b,1], color=color, s=s, lw=lw, marker=shape[1])
        ax.text(npos[b,0]+0.02, npos[b,1]-0.01, gens_s[c], color=color)
        if S16:
            tx = 20
        else:
            tx = 14
        if b == tx:
            ax.scatter(npos[b,0], npos[b,1], color=color, s=s, lw=lw, marker=shape[1], label=labels[1])
        this_x.append(npos[b,0]), this_y.append(npos[b,1])
    first_group, second_group = npos[:len(G)], npos[len(G):]
    #anova = permanova(first_group, second_group)
    x.append(this_x), y.append(this_y)
    this_x, this_y = [], []
    this_x, this_y = [], []
    for z in range(2):
        nstd = 2
        x1, y1 = x[z], y[z]
        cov = np.cov(x1, y1)
        vals, vecs = eigsorted(cov)
        theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))
        w, h = 2 * nstd * np.sqrt(vals)
        ell = Ellipse(xy=(np.mean(x1), np.mean(y1)),width=w, height=h,angle=theta, color=el_cols[z])     
        ell.set_facecolor('none')
        ax.add_artist(ell)
    if legend:
        ax.legend(scatterpoints=1, bbox_to_anchor=(lp[0],lp[1]))
    ax.set_title(title, loc='left')
    #plt.xlim([-0.3, 0.6])
    #plt.ylim([-0.3, 0.6])
    #plt.savefig('New_NMDS_OTU_18S_new.pdf', bbox_inches='tight')
    #plt.close()
    return
fig = plt.figure(figsize=(15, 5))
NMDS_plot(fn, plt.subplot(121), legend=False, title='16S', lp=[1.26, 1.025])
NMDS_plot(fn2, plt.subplot(122), legend=True, title='18S', lp=[1.26, 1.025])
plt.tight_layout()
os.chdir('/Users/u1560915/Documents/GitHub/ChitinActivity/MiSeq/16S/Figures/')
plt.savefig('Fig 4.pdf', bbox_inches='tight')
plt.close()

NMDS_plot(fn, plt.subplot(111), legend=True, title='', lp=[1.37, 1.025])
os.chdir('/Users/u1560915/Documents/GitHub/ChitinActivity/MiSeq/16S/Figures/')
plt.savefig('Fig 4 16S.pdf', bbox_inches='tight')
plt.close()

NMDS_plot(fn2, plt.subplot(111), legend=True, title='', lp=[1.37, 1.025])
os.chdir('/Users/u1560915/Documents/GitHub/ChitinActivity/MiSeq/16S/Figures/')
plt.savefig('Fig 4 18S.pdf', bbox_inches='tight')
plt.close()