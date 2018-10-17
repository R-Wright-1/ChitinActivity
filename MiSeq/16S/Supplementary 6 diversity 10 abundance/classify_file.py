import csv
taxonomy = 'Taxonomy.csv'
import numpy
#Genus = 6, Family = 5, Order = 4, Class = 3, Phylum = 2

def classify_file(fn, level):
    with open(taxonomy, 'rU') as f:
        reader = csv.reader(f)
        taxo = []
        for row in reader:
            taxo.append(row)
    with open(fn, 'rU') as f:
        reader = csv.reader(f)
        rows = []
        for row in reader:
            rows.append(row)
    for row in rows:
        OTU = row[0]
        for a in range(len(taxo)):
            if taxo[a][0] == OTU:
                row[0] = taxo[a][level]
    with open(fn[:-4]+'_'+taxo[0][level]+'.csv', 'w') as f:
        writer = csv.writer(f)
        for row in rows:
            writer.writerow(row) 
    return fn[:-4]+'_'+taxo[0][level]+'.csv'

def group_file(fn, other=True, limit=1):
    with open(fn, 'rU') as f:
        reader = csv.reader(f)
        rows = []
        for row in reader:
            rows.append(row)
    row0 = rows[0]
    del rows[0]
    all_OTUs = []
    for a in rows:
        adding = True
        for b in all_OTUs:
            if a[0] == b:
                adding = False
        if adding:
            all_OTUs.append(a[0])
    means = []
    for c in all_OTUs:
        this_row = []
        for d in rows[0]:
            this_row.append([])
        for f in range(len(rows)):
            if rows[f][0] == c:
                for g in range(len(rows[f])):
                    if g > 0:
                        this_row[g].append(float(rows[f][g]))
        means.append(this_row)
    new_means = []
    for h in range(len(means)):
        new_mean = []
        for g in range(len(means[h])):
            if g > 0:
                new_mean.append(sum(means[h][g])*100)
        new_means.append(new_mean)
    if other:
        other = []
        for z in range(len(rows[0])):
            other.append(0)
        deleting = []
        for a in range(len(new_means)):
            if sum(new_means[a]) < limit*4:
                for b in range(len(new_means[a])):
                    other[b] += new_means[a][b]
                deleting.append(a)
        new_new_means, new_new_OTUs = [], []
        for c in range(len(new_means)):
            adding = True
            for d in deleting:
                if c == d:
                    adding = False
            if adding:
                new_new_means.append(new_means[c])
                new_new_OTUs.append(all_OTUs[c])
        new_new_means.append(other)
        new_new_OTUs.append('Other')
        new_means = new_new_means
        all_OTUs = new_new_OTUs
    new_new = []
    for i in range(len(new_means)):
        new = [all_OTUs[i]]
        for j in range(len(new_means[i])):
            new.append(new_means[i][j])
        new_new.append(new)
    new_means = new_new
    with open(fn[:-4]+'_grouped.csv', 'w') as f:
        writer = csv.writer(f)
        writer.writerow(row0)
        for row in new_means:
            writer.writerow(row)
    return fn[:-4]+'_grouped.csv'
    
#new_fn = classify_file('3_daily_percent_grouped.csv', 6)
#new_fn = group_file(new_fn)