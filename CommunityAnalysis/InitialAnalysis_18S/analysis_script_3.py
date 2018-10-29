import csv
import numpy
from scipy import stats
import pandas
import statsmodels.stats.multitest as smm

#run this after the simper has been run in R, and after analysis script 1!
#if you have kept with my file naming, then you should only need to change name here to reflect the
#treatment name that you have used

name = 'treatment4_controls'

simper = name+'_simper.csv'
percent = name+'_percent.csv'
all_files = name+'_percent.csv'
meta = name+'_meta.csv'
no_perc = name+'_samples.csv'

def get_most_important_simpers(fn):
    #This part runs through all of the output of the simper in r and groups the contributions for each OTU
    #(from each comparison), takes a mean, and then ranks them in order of importance

    #read in the file
    with open(fn, 'rU') as f:
        reader = csv.reader(f)
        rows = []
        for row in reader:
            rows.append(row)
            
    #check how many groups (and how many comparisons) there are
    all_OTUs = []
    comps = []
    for a in range(len(rows)):
        if a > 0:
            adding = True
            for b in range(len(comps)):
                if rows[a][1] == comps[b]:
                    adding = False
            if adding:
                comps.append(rows[a][1])
    # comps should now have a list of all competitions
    
    #now add all of the unique OTU's that came from this analysis            
    for a in rows[1:]:
        adding = True
        for b in all_OTUs:
            if b == a[3]:
                adding = False
        if adding:
            all_OTUs.append(a[3])
    
    #get a mean for each OTU that was important in one comparison
    all_means = []
    for c in range(len(all_OTUs)):
        this_mean = []
        for d in rows:
            if d[3] == all_OTUs[c]:
                this_mean.append(float(d[2]))
        all_means.append(sum(this_mean)/len(comps))
    #sort them so that they're in order of importance
    all_means, all_OTUs = (list(t) for t in zip(*sorted(zip(all_means, all_OTUs))))
    #this just changes the names from OTU (which R needed for the SIMPER script) back to say ASV
    for a in range(len(all_OTUs)):
        all_OTUs[a] = 'ASV'+all_OTUs[a][3:]
    #now reverse the order, so that we have most important first rather than last
    all_OTUs, all_means = list(reversed(all_OTUs)), list(reversed(all_means))
    return all_OTUs, all_means
    
def get_kruskal(all_abun):
    zeros = True
    #this part just checks that we're not running this with something that is never there...
    #i.e. if the abundance for all replicates of all treatments is 0, then it won't continue with 
    #the stats test, causing an error
    for a in range(len(all_abun)):
        if sum(all_abun[a]) != 0:
            zeros = False
    if zeros:
        return 1, 1
    #convert numbers to an array rather than a list (needed as input to the stats test)
    for a in range(len(all_abun)):
        all_abun[a] = numpy.array(all_abun[a])
    #if you want to run this with more than 4 treatments then you will need to add more in here,
    #as currently this is only for a maximum of 4 treatments
    if len(all_abun) == 4:
        kstats = stats.kruskal(all_abun[0], all_abun[1], all_abun[2], all_abun[3])
    elif len(all_abun) == 3:
        kstats = stats.kruskal(all_abun[0], all_abun[1], all_abun[2])
    elif len(all_abun) == 2:
        kstats = stats.kruskal(all_abun[0], all_abun[1])
    krusk_p = smm.fdrcorrection(kstats[1])[1]
    return kstats[0], krusk_p[0]
    
def get_means(percent, OTUs, meta, means):
    #open file with percent abundances
    with open(percent, 'rU') as f:
        rows = []
        for row in csv.reader(f):
            rows.append(row)
    new_rows = []
    #get abundances for each ASV identified by simper
    for a in range(len(OTUs)):
        for b in range(len(rows)):
            if b > 0:
                if rows[b][0] == OTUs[a]:
                    new_rows.append(rows[b])
            else:
                if a == 0:
                    new_rows.append(rows[b])
    #rename rows according to meta file, if necessary
    with open(meta, 'rU') as f:
        rows = []
        for row in csv.reader(f):
            rows.append(row)
    for c in range(len(new_rows[0])):
        for d in range(len(rows)):
            if new_rows[0][c] == rows[d][0]:
                new_rows[0][c] = rows[d][1]
    for c in range(len(new_rows)):
        for d in range(len(new_rows[c])):
            if c > 0 and d > 0:
                new_rows[c][d] = float(new_rows[c][d])
    new_otu, new_name, new_means, new_stds, krusk, krusk_p, means = [], [], [], [], [], [], means
    #get only unique treatments/column headers
    names = []
    for a in range(len(new_rows[0])):
        adding = True
        for b in range(len(names)):
            if names[b] == new_rows[0][a]:
                adding = False
        if adding:
            names.append(new_rows[0][a])
    numbers = []
    #group into unique treatments
    for c in range(len(new_rows[0])):
        for d in range(len(names)):
            if new_rows[0][c] == names[d]:
                numbers.append(d)
    for a in range(len(new_rows)):
        row0, new_rows[a] = (list(t) for t in zip(*sorted(zip(numbers, new_rows[a]))))
    #for each ASV, get the simper contribution, abundances for each treatment, and kruskal-wallis
    #stats values
    for f in range(len(new_rows)-1):
        f += 1
        row0 = new_rows[0]
        prev = ''
        this_abun, this_row_mean, this_row_stds = [], [], []
        this_all_abun = []
        col_names = []
        new_otu.append(new_rows[f][0])
        #get abundance data, and get a mean for each treatment
        for g in range(len(row0)):
            if g == 0:
                continue
            if g == 1:
                this_abun.append(new_rows[f][g])
                prev = row0[g]
            elif prev == row0[g]:
                this_abun.append(new_rows[f][g])
                prev = row0[g]
            else:
                col_names.append(prev)
                prev = row0[g]
                this_row_mean.append(numpy.mean(this_abun))
                this_row_stds.append(numpy.std(this_abun))
                this_all_abun.append(this_abun)
                this_abun = [new_rows[f][g]]
        col_names.append(prev)
        this_row_mean.append(numpy.mean(this_abun))
        this_row_stds.append(numpy.std(this_abun))
        this_all_abun.append(this_abun)
        #get kruskal values
        k, kp = get_kruskal(this_all_abun)
        new_means.append(this_row_mean)
        new_stds.append(this_row_stds)
        krusk.append(k)
        krusk_p.append(kp)
    headers = ['OTU', 'Simper contribution']
    for a in range(len(col_names)):
        headers.append(col_names[a]+' mean')
        headers.append(col_names[a]+' std')
    headers.append('kruskal stats')
    headers.append('kruskal p')
    #save this file
    with open(percent[:-12]+'_simper_means.csv', 'w') as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        for a in range(len(means)):
            this_row = [new_otu[a], means[a]]
            for b in range(len(new_means[a])):
                this_row.append(new_means[a][b])
                this_row.append(new_stds[a][b])
            this_row.append(krusk[a])
            this_row.append(krusk_p[a])
            writer.writerow(this_row)
    return
    
def transpose(rows): #this just flips rows and columns with each other
    cols = []
    for a in range(len(rows[0])):
        col = []
        for b in range(len(rows)):
            col.append(rows[b][a])
        cols.append(col)
    return cols
    
def group_file(fn, meta):
    #open meta file
    with open(meta, 'rU') as f:
        meta = []
        for m in csv.reader(f):
            meta.append(m)
    #get old and new treatment names
    old, new = [], []
    for a in range(len(meta)):
        if a > 0:
            old.append(meta[a][0])
            new.append(meta[a][1])
    #open up file to be grouped
    with open(fn, 'rU') as f:
        rows = []
        for row in csv.reader(f):
            rows.append(row)
    rows = transpose(rows)
    #rename treatments with new name
    for a in range(len(rows)):
        for b in range(len(old)):
            if rows[a][0] == old[b]:
                rows[a][0] = new[b]
    #get unique treatment names
    unique = []
    for c in range(len(new)):
        adding = True
        for d in range(len(unique)):
            if new[c] == unique[d]:
                adding = False
        if adding:
            unique.append(new[c])
    #go through the rows and group those together that are the same
    new_rows = [rows[0]]
    for e in range(len(unique)):
        these_rows = []
        for f in range(len(rows)):
            if unique[e] == rows[f][0]:
                these_rows.append(rows[f][1:])
        this_row = [unique[e]]
        for g in range(len(these_rows[0])):
            this_val = []
            for h in range(len(these_rows)):
                this_val.append(float(these_rows[h][g]))
            this_row.append(numpy.mean(this_val))
        new_rows.append(this_row)
    new_rows = transpose(new_rows)
    #now save the new file
    with open(fn[:-4]+'_grouped.csv', 'w') as f:
        writer = csv.writer(f)
        for a in range(len(new_rows)):
            writer.writerow(new_rows[a])
    return

#all_OTUs, all_means = get_most_important_simpers(simper) #simper = name of file that r generated with simper contributions
#get_means(percent, all_OTUs, meta, all_means) #percent file = file with relative abundances, meta = file with old and new treatment names, all_OTUs and all_means are output from simper function
group_file(percent, meta) #percent = file to be grouped, meta = meta file with old and new treatment names
group_file(no_perc, meta)
