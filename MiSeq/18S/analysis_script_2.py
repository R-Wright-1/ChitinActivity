import csv
import numpy
import analysis_script_1
from scipy import stats
import pandas

#run this after the r script!

def get_most_important_simpers(fn):
    with open(fn, 'rU') as f:
        reader = csv.reader(f)
        rows = []
        for row in reader:
            rows.append(row)
    all_OTUs = []
    for a in rows[1:]:
        adding = True
        for b in all_OTUs:
            if b == a[3]:
                adding = False
        if adding:
            all_OTUs.append(a[3])
    all_means = []
    for c in range(len(all_OTUs)):
        this_mean = []
        for d in rows:
            if d[3] == all_OTUs[c]:
                this_mean.append(float(d[2]))
        all_means.append(numpy.mean(this_mean))
    all_means, all_OTUs = (list(t) for t in zip(*sorted(zip(all_means, all_OTUs))))
    all_OTUs, all_means = list(reversed(all_OTUs)), list(reversed(all_means))
    return all_OTUs, all_means

def get_values_shared(all_OTUs, fn, treatments, samples):
    with open(fn, 'rU') as f:
        reader = csv.reader(f)
        rows = []
        for row in reader:
            rows.append(row)
    new_lines = [rows[0]]
    for a in all_OTUs:
        for b in rows:
            if a == b[0]:
                new_lines.append(b)
    for c in range(len(rows[0])):
        for d in range(len(samples)):
            if rows[0][c] == samples[d]:
                rows[0][c] = treatments[d]
    with open(fn[0:-4]+'_simper.csv', 'w') as f:
        writer = csv.writer(f)
        for line in new_lines:
            writer.writerow(line)
    return fn[0:-4]+'_simper.csv'
    
def get_means(percent):
    with open(percent, 'rU') as f:
        reader = csv.reader(f)
        rows = []
        for row in reader:
            rows.append(row)
    treatments, mean_values = [], []
    for a in range(len(rows[0])):
        if a > 0:
            adding = True
            for b in range(len(treatments)):
                if rows[0][a] == treatments[b]:
                    adding = False
            if adding:
                treatments.append(rows[0][a])
    print treatments
    means, stds = [], []
    for c in range(len(rows)):
        vals = []
        for d in range(len(treatments)):
            vals.append([])
        mean_values.append(vals)
        means.append([])
        stds.append([])
    for d in range(len(rows[0])):
        if d > 0:
            name = rows[0][d]
            for e in range(len(treatments)):
                if name == treatments[e]:
                    ind = e
                    rows[0][d] = treatments[e]
            for f in range(len(rows)):
                if f > 0:
                    mean_values[f][ind].append(float(rows[f][d]))
    with open(percent[:-4]+'_treat.csv', 'w') as f:
        writer = csv.writer(f)
        for row in rows:
            writer.writerow(row)
    del mean_values[0]
    del means[0]
    del stds[0]
    krusk_stats, krusk_p = [], []
    for g in range(len(mean_values)):
        for h in range(len(mean_values[g])):
            means[g].append(numpy.mean(mean_values[g][h]))
            stds[g].append(numpy.std(mean_values[g][h]))
        print len(mean_values[g])
        kstats = stats.kruskal(numpy.array(mean_values[g][0]), numpy.array(mean_values[g][1]))
        krusk_stats.append(kstats[0])
        krusk_p.append(kstats[1])
    with open(percent[:-4]+'_mean.csv', 'w') as f:
        writer = csv.writer(f)
        writing = [rows[0][0]]
        for z in range(len(treatments)):
            writing.append(treatments[z]+'_mean')
            writing.append(treatments[z]+'_std')
        writing.append('krusk_stats')
        writing.append('krusk_p')
        writer.writerow(writing)
        for a in range(len(rows[1:])):
            writing = [rows[a+1][0]]
            for b in range(len(treatments)):
                writing.append(means[a][b])
                writing.append(stds[a][b])
            writing.append(krusk_stats[a])
            writing.append(krusk_p[a])
            writer.writerow(writing)
    return
"""
def group_percent_treat(treatments, files):
    for a in range(len(treatments)):
        treatment = treatments[a]
        new_name = files[a]
        with open(files[a], 'rU') as f:
            reader = csv.reader(f)
            rows = []
            for row in reader:
                rows.append(row)
        for a in range(len(rows[0])):
            if a > 0:
                rows[0][a] = treatment[a-1]
        all_names = []
        for b in range(len(rows[0])):
            if b > 0:
                adding = True
                for c in range(len(all_names)):
                    if rows[0][b] == all_names[c]:
                        adding = False
                if adding:
                    all_names.append(rows[0][b])
        means = []
        for d in range(len(rows)):
            this_row = []
            for e in range(len(all_names)):
                this_row.append([])
            means.append(this_row)
        for f in range(len(rows[0])):
            if f > 0:
                name = rows[0][f]
                for g in range(len(all_names)):
                    if name == all_names[g]:
                        for h in range(len(rows)):
                            if h > 0:
                                means[h][g].append(float(rows[h][f]))
        for i in range(len(means)):
            for j in range(len(means[i])):
                means[i][j] = numpy.mean(means[i][j])
        del means[0]
        with open(new_name[0:-4]+'_grouped.csv', 'w') as f:
            writer = csv.writer(f)
            writing = [rows[0][0]]
            for k in range(len(all_names)):
                writing.append(all_names[k])
            writer.writerow(writing)
            for row in range(len(means)):
                writing = []
                writing.append(rows[row+1][0])
                for l in range(len(means[row])):
                    writing.append(means[row][l])
                writer.writerow(writing)
    return    
"""    
names, indices, all_samples, all_treatments = analysis_script_1.get_meta_treatments('18S_metadata.csv') 

simper = ['12_LS_clean_simper.csv']
percent = ['12_LS_percent.csv']
all_files = ['', '', '', '', '', '', '9_L_GR_percent.csv', '', '11_beg_end_S_percent.csv', '12_LS_percent.csv']

#need to run these separately for each dataset to add the right number of groups in kruskal test
all_OTUs, all_means = get_most_important_simpers(simper[0])
name = get_values_shared(all_OTUs, percent[0], all_treatments[10], all_samples[10])
get_means(name)
#group_percent_treat(all_treatments, all_files)

#now go to plotting file!
