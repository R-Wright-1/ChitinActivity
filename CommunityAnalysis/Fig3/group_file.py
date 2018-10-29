import csv

def get_tax(fn, otus):
    with open(fn, 'rU') as f:
        rows = []
        for row in csv.reader(f):
            rows.append(row)
    prev_otus = otus
    for a in range(len(otus)):
        for b in range(len(rows)):
            if otus[a] == rows[b][0]:
                otus[a] = rows[b][1:]
    for c in range(len(otus)):
        if c == 0 or otus[c] == 'Group':
            continue
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
        else:
            otus[c] = 'NA'
    return otus
    
def get_old_otus(otus):
    old_otus = []
    prev_otus = otus
    for d in range(len(prev_otus)):
        if prev_otus[d][:3] != 'ASV':
            old_otus.append(prev_otus[d])
            continue
        otu = int(prev_otus[d][3:])
        otu = 'ASV'+str(otu)
        old_otus.append(otu)
    return old_otus

def gf(fn, limit, tax, return_old_otus=False):
    with open(fn, 'rU') as f:
        rows = []
        for row in csv.reader(f):
            rows.append(row)
    unique, unique_num = [], []
    for z in range(len(rows)):
        for y in range(len(rows[z])):
            if z > 0 and y > 0:
                rows[z][y] = float(rows[z][y])
    otus = []
    for a in range(len(rows)):
        otus.append(rows[a][0])
    old = get_old_otus(otus)
    new_otus = get_tax(tax, otus)
    for a in range(len(new_otus)):
        adding = True
        for b in range(len(unique)):
            if new_otus[a] == unique[b]:
                adding = False
                unique_num[b] += 1
        if adding:
            unique.append(new_otus[a])
            unique_num.append(1)
    new_rows = []
    for c in range(len(unique)):
        for d in range(len(rows)):
            rows[d][0] = new_otus[d]
    """
            if rows[d][0] == unique[c]:
                if this_row == []:
                    new_old.append(old[d])
                    this_row = rows[d]
                    count = 1
                else:
                    count += 1
                    for e in range(len(rows[d])):
                        if e > 0:
                            this_row[e] += rows[d][e]
        this_row[0] = this_row[0]+' ('+str(count)+')'
        new_rows.append(this_row)
    """
    new_old = old
    old = []
    transpose = rows
    new_rows, other = [], []
    for h in range(len(transpose)):
        if h == 0:
            old.append(new_old[h])
            new_rows.append(transpose[h])
        else:
            total = max(transpose[h][1:])
            if total > limit:
                old.append(new_old[h])
                new_rows.append(transpose[h])
            else:
                if other == []:
                    transpose[h][0] = 'Other'
                    other = transpose[h]
                    count = 1
                else:
                    count += 1
                    for i in range(len(transpose[h])):
                        if i > 0:
                            other[i] += transpose[h][i]
    other[0] = 'Other ('+str(count)+')'
    new_rows.append(other)
    new_name = fn[:-7]+'grouped.csv'
    with open(new_name, 'w') as f:
        writer = csv.writer(f)
        for a in new_rows:
            writer.writerow(a)
    if return_old_otus:
        return new_rows, old
    return new_rows