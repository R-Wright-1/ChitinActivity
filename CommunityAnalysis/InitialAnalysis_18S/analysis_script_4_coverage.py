import csv
import numpy
import os
import random

#This script is to calculate the total number of reads from each sample, and also to calculate the
#resulting coverage of the sample
#Just change fn to be of the script that you want to calculate coverage (this should have the samples
#starting at column 2, i.e. without full taxonomy data)

fn = 'treatment3_all_samples.csv'
meta = '16S_metadata.csv'

def transpose(rows):
    transpose = []
    for a in range(len(rows[0])):
        new_col = []
        for b in range(len(rows)):
            if a > 0 and b > 0:
                new_col.append(float(rows[b][a]))
            else:
                new_col.append(rows[b][a])
        transpose.append(new_col)
    return transpose
    
def transpose_no_float(rows):
    transpose = []
    for a in range(len(rows[0])):
        new_col = []
        for b in range(len(rows)):
            new_col.append(rows[b][a])
        transpose.append(new_col)
    return transpose

def remove_singles(rows):
    new_rows = []
    for a in range(len(rows)):
        if a == 0:
            new_rows.append(rows[a])
            continue
        for b in range(len(rows[a])):
            if b == 0:
                continue
            rows[a][b] =float(rows[a][b])
        if sum(rows[a][1:]) > 1:
            new_rows.append(rows[a])
    return new_rows
    
def coverage_total(rows):
    totals = []
    all_coverage = []
    for a in range(len(rows)):
        if a == 0:
            rows[a].append('Total')
            rows[a].append('Coverage')
            continue
        singles = 0
        for b in range(len(rows[a])):
            if rows[a][b] == 1:
                singles += 1
        total = sum(rows[a][1:])
        totals.append(total)
        rows[a].append(total)
        coverage = 1-(singles/total)
        all_coverage.append(coverage)
        rows[a].append(coverage)
    minimum = min(totals)
    return rows, all_coverage, totals, minimum
    
def normalise(rows, minimum):
    for a in range(len(rows)):
        if a == 0:
            continue
        total = sum(rows[a][1:])
        divide = total/minimum
        for b in range(len(rows[a])):
            if b == 0:
                continue
            rows[a][b] = rows[a][b]/divide
            rows[a][b] = round(rows[a][b], 0)
    return rows

def rarefy_remove(rows, minimum):
    for a in range(len(rows)):
        if a == 0:
            continue
        total = sum(rows[a][1:])
        remove = total-minimum
        count = 0
        rownot0 = []
        place = []
        for b in range(len(rows[a])):
            if b == 0:
                continue
            if rows[a][b] != 0:
                rownot0.append(rows[a][b])
                place.append(b)
        l = len(rownot0)-1
        while count < remove:
            rnum = random.randint(0,l)
            if (rownot0[rnum] - 1) > -1:
                rownot0[rnum] -= 1
                count += 1
        for b in range(len(rownot0)):
            rows[a][place[b]] = rownot0[b]
    return rows
    
def rarefy_keep(rows, minimum):
    for a in range(len(rows)):
        if a == 0:
            continue
        count = 0
        rownot0 = []
        place = []
        for b in range(len(rows[a])):
            if b == 0:
                continue
            if rows[a][b] != 0:
                rownot0.append(rows[a][b])
                place.append(b)
        l = len(rownot0)-1
        new_row = []
        for c in range(len(rownot0)):
            new_row.append(0)
        while count < minimum:
            rnum = random.randint(0,l)
            if (rownot0[rnum] - 1) > -1:
                rownot0[rnum] -= 1
                new_row[rnum] += 1
                count += 1
        for d in range(len(new_row)):
            rows[a][place[d]] = new_row[d]
    return rows
    
def group_file(names, new_names, name):
    with open('Total and coverage '+fn, 'rU') as f:
        rows = []
        for row in csv.reader(f):
            rows.append(row)
    for a in range(len(rows)):
        for b in range(len(rows[a])):
            if a > 0 and b > 0:
                rows[a][b] = float(rows[a][b])
    del rows[-1]
    del rows[-1]
    rows = transpose(rows)
    new_rows = []
    for a in range(len(names)):
        for b in range(len(rows)):
            if rows[b][0] == names[a]:
                rows[b][0] = new_names[a]
                new_rows.append(rows[b])
    new_rows[0][0] = 'Group'
    unique = []
    for c in range(len(new_rows)):
        adding = True
        for d in range(len(unique)):
            if new_rows[c][0] == unique[d]:
                adding = False
        if adding:
            unique.append(new_rows[c][0])
    new_new_rows = []
    for e in range(len(unique)):
        if e == 0:
            new_new_rows.append(new_rows[e])
        else:
            this_row = []
            for f in range(len(new_rows)):
                if new_rows[f][0] == unique[e]:
                    this_row.append(new_rows[f])
            if len(this_row) < 2:
                new_new_rows.append(this_row[0])
            else:
                trow = []
                for g in range(len(this_row[0])):
                    if g == 0:
                        trow.append(this_row[0][g])
                    else:
                        this_num = []
                        for h in range(len(this_row)):
                            this_num.append(this_row[h][g])
                        trow.append(round(numpy.mean(this_num), 0))
                new_new_rows.append(trow)
    new_rows = transpose_no_float(new_new_rows)
    """
    with open(name+'_grouped_'+fn, 'w') as f:
        writer = csv.writer(f)
        for row in new_rows:
            writer.writerow(row)
    """
    rows = transpose_no_float(new_rows)
    for a in range(len(rows)):
        if a == 0:
            continue
        total = sum(rows[a][1:])
        for b in range(len(rows[a])):
            if b == 0:
                continue
            rows[a][b] = (rows[a][b]/total)*100
    rows = transpose_no_float(rows)
    """
    with open(name+'_grouped_percent_'+fn, 'w') as f:
        writer = csv.writer(f)
        for row in rows:
            writer.writerow(row)
    """
    return
    
def get_meta_treats(rows, meta):
    with open(meta, 'rU') as f:
        meta = []
        for row in csv.reader(f):
            meta.append(row)
    meta = transpose_no_float(meta)
    for a in range(len(meta)-1):
        a += 1
        this_treat, this_treat_names = [], []
        for b in range(len(meta[a])):
            if meta[a][b] != '':
                this_treat.append(meta[a][b])
                this_treat_names.append(meta[0][b])
        group_file(this_treat_names, this_treat, meta[a][0])
    return

with open(fn, 'rU') as f:
    rows = []
    for row in csv.reader(f):
        rows.append(row)
rows = remove_singles(rows)
rows = transpose(rows)
new_rows, coverage, total, minimum = coverage_total(rows)
print numpy.mean(coverage)
new_rows = transpose(new_rows)
with open(fn[:-4]+'_total_and_coverage.csv', 'w') as f:
    writer = csv.writer(f)
    for row in new_rows:
        writer.writerow(row) 
get_meta_treats(rows, meta)   
with open(fn[:-4]+'_total_and_coverage.csv', 'rU') as f:
    rows = []
    for row in csv.reader(f):
        rows.append(row)

del rows[-1]
del rows[-1]
rows = transpose(rows)
norm_rows = normalise(rows, minimum)
new_norm_rows, coverage, total, minimum = coverage_total(norm_rows)
new_norm_rows = transpose(new_norm_rows)
"""
with open('Normalised '+fn, 'w') as f:
    writer = csv.writer(f)
    for row in new_norm_rows:
        writer.writerow(row)
"""