import csv
fn = 'otu_table.csv'
shared = '16S_shared_overall_mean.csv'

with open(fn, 'rU') as f:
    rows = []
    for row in csv.reader(f):
        rows.append(row)

otus, ids = [], []
for a in range(len(rows[1])):
    if a > 0:
        otu = rows[1][a]
        for b in range(len(rows)):
            if rows[b][a] == '1':
                this_id = rows[b][0]
        otus.append(otu)
        ids.append(this_id)

with open('IDs.csv', 'w') as f:
    writer = csv.writer(f)
    writer.writerow(['OTU', 'ID'])
    for a in range(len(otus)):
        writing = [otus[a], ids[a]]
        writer.writerow(writing)

with open(shared, 'rU') as f:
    rows = []
    for row in csv.reader(f):
        rows.append(row)
row0 = rows[0]
del rows[0]

new_id, rest_of_row = [], []
for a in range(len(rows)):
    otu = rows[a][0]
    otui = int(rows[a][0][3:])
    for b in range(len(otus)):
        if otui == int(otus[b][3:]):
            new_id.append(ids[b])
            rest_of_row.append(rows[a][1:])
row0[0] = '#OTU ID'

unique_ids = []
for i in new_id:
    adding = True
    for u in unique_ids:
        if i == u:
            adding = False
    if adding == True:
        unique_ids.append(i)

all_rows = []
for u in unique_ids:
    all_rows.append([])

for u in range(len(unique_ids)):
    for a in range(len(new_id)):
        if unique_ids[u] == new_id[a]:
            all_rows[u].append(rest_of_row[a])

for b in range(len(all_rows)):
    if len(all_rows[b]) > 1:
        new_row = []
        for c in range(len(all_rows[b])):
            for d in range(len(all_rows[b][c])):
                all_rows[b][c][d] = float(all_rows[b][c][d])
        for e in range(len(all_rows[b][0])):
            for f in range(len(all_rows[b])):
                all_rows[b][0][e] += all_rows[b][f][e]
        all_rows[b] = [all_rows[b][0]]

new_id, rest_of_row = unique_ids, all_rows

new_r0 = ''
for c in range(len(row0)):
    new_r0 += row0[c]+'\t'
new_r0 += '\n'
    
for_text = [new_r0]
for a in range(len(new_id)):
    this_row = new_id[a]+'\t'
    for b in range(len(rest_of_row[a])):
        for c in range(len(rest_of_row[a][b])):
            this_row += str(rest_of_row[a][b][c])+'\t'
    this_row += '\n'
    for_text.append(this_row)

f = open('new_otu_file_overall_mean.txt', 'w')
for t in for_text:
    t = str(t)
    f.write(t)
f.close()


