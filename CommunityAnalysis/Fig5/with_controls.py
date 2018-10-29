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
    clf = PCA(n_components=2)
    X_true = clf.fit_transform(X_true)
    pos = clf.fit_transform(pos)
    npos = clf.fit_transform(npos)
    return pos, npos
    
def eigsorted(cov):
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    return vals[order], vecs[:,order]
    
def plot(fn, ax, shapes, colors):
    controls = get_control_samples(fn)
    c = pd.DataFrame(controls)
    c = c.transpose()
    pos, npos = transform_for_NMDS(c)
    for b in range(len(controls[0])):
        ax.scatter(npos[b,0], npos[b,1], marker=shapes[b], color=colors[b])
    return
    
ax = plt.subplot(211)
fn = '18S_controls.csv'

shapes_16S = ['o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 's', 's', 's', 's', 's', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '^', '^', '^', '^', '^', 'o', 'o', 'o']
colors_16S = ['r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'b', 'b', 'b', 'b', 'b', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'm', 'm', 'm', 'm', 'm', 'k', 'k', 'k']
shapes_18S = ['o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 'o', 's', 's', 's', 's', 's', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '*', '^', '^', '^', 'o', 'o', 'o']
colors_18S = ['r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 'b', 'b', 'b', 'b', 'b', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'g', 'm', 'm', 'm', 'k', 'k', 'k']

plot(fn, ax, shapes_16S, colors_16S)

