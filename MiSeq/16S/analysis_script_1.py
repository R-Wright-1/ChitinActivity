import csv

#run this script first!

def convert_shared(fn, new_fn):
    with open(fn, 'rU') as f:
        rows = []
        for row in f:
            rows.append(row)
    writing, this_word = [], ''
    all_writing = []
    for a in range(len(rows)):
        for b in range(len(rows[a])):
            letter = rows[a][b]
            if letter == '\t':
                writing.append(this_word)
                this_word = ''
            else:
                this_word += letter
        all_writing.append(writing)
        writing = []
    with open(new_fn, 'w') as f:
        writer = csv.writer(f)
        for c in range(len(all_writing[0])):
            to_write = []
            if c == 0 or c == 2:
                continue
            for d in range(len(rows)):
                to_write.append(all_writing[d][c])
            writer.writerow(to_write)
    return

def get_meta_treatments(fn):
    with open(fn, 'rU') as f:
        reader = csv.reader(f)
        rows = []
        for row in reader:
            rows.append(row)
    names, indices = [], []
    for a in range(len(rows[0])):
        tr = rows[0][a]
        if tr[0:9] == 'treatment':
            name = tr[9:]
            index = a
            if len(tr) == 10:
                name = tr[-1]
            names.append(name), indices.append(index)
    all_samples, all_treatments = [], []
    for b in range(len(names)):
        ind = indices[b]
        sample, treatment = [], []
        for c in range(len(rows)):
            if c > 0 and rows[c][ind] != '':
                sample.append(rows[c][0])
                treatment.append(rows[c][ind])
        all_samples.append(sample)
        all_treatments.append(treatment)           
    return names, indices, all_samples, all_treatments
    
def treatment_shared(names, all_samples):
    with open('Shared.csv', 'rU') as f:
        reader = csv.reader(f)
        rows = []
        for row in reader:
            rows.append(row)
    files = []
    for a in range(len(names)):
        indices, new_rows = [0], []
        for b in range(len(rows[0])):
            for c in range(len(all_samples[a])):
                if rows[0][b] == all_samples[a][c]:
                    indices.append(b)
        for d in range(len(rows)):
            this_row = []
            for e in range(len(indices)):
                ind = indices[e]
                this_row.append(rows[d][ind])
            new_rows.append(this_row)
        with open(names[a]+'_shared.csv', 'w') as f:
            writer = csv.writer(f)
            for g in new_rows:
                writer.writerow(g)
        files.append(names[a]+'_shared.csv')
    return files
    
def clean_shared(files):
    new_files = []
    for fn in files:
        with open(fn, 'rU') as f:
            reader = csv.reader(f)
            rows = []
            for row in reader:
                rows.append(row)
        new_rows = [rows[0]]
        for a in rows[1:]:
            not_all_zeros = False
            for b in range(len(a)):
                if b > 0:
                    if int(a[b]) > 0:
                        not_all_zeros = True
            if not_all_zeros:
                new_rows.append(a)
        with open(fn, 'w') as f:
            writer = csv.writer(f)
            for row in new_rows:
                writer.writerow(row)
    return
    
def clean_meta(names, indices, all_samples, all_treatments, meta):
    cols = ['sample', 'treatment', 'coverage', 'nseqs', 'bergerparker', 'chao', 'chao_lci', 'chao_hci', 'shannon', 'shannon_lci', 'shannon_hci']
    all_meta = []
    with open(meta, 'rU') as f:
        reader = csv.reader(f)
        for row in reader:
            all_meta.append(row)
    #row0 = all_meta[0]
    #cov, nseqs, bp, chao, chao_l, chao_h, shannon, shannon_l, shannon_h = all_meta[0][-9], all_meta[0][-8], all_meta[0][-7], all_meta[0][-6], all_meta[0][-5], all_meta[0][-4], all_meta[0][-3], all_meta[0][-2], all_meta[0][-1]
    for i in range(len(names)):
        fn, samples, treatments = names[i], all_samples[i], all_treatments[i]
        new = [cols]
        for a in all_meta:
            this_row = []
            for b in range(len(samples)):
                if a[0] == samples[b]:
                    this_row = [samples[b], treatments[b], a[-9], a[-8], a[-7], a[-6], a[-5], a[-4], a[-3], a[-2], a[-1]]
            if len(this_row) != 0:
                new.append(this_row)
        with open(fn+'_meta.csv', 'w') as f:
            writer = csv.writer(f)
            for n in new:
                writer.writerow(n)
    return
    
def convert_for_r_file(files):
    for a in files:
        new_fn = a[0:-4]+'_r.csv'
        with open(a, 'rU') as f:
            reader = csv.reader(f)
            rows = []
            for row in reader:
                rows.append(row)
        new_rows = []
        for a in range(len(rows[0])):
            this_row = []
            for b in range(len(rows)):
                this_row.append(rows[b][a])
            new_rows.append(this_row)
        with open(new_fn, 'w') as f:
            writer = csv.writer(f)
            for row in new_rows:
                writer.writerow(row)
    return
    
def convert_percentages(files):
    new_files = []
    for a in files:
        with open(a, 'rU') as f:
            rows = []
            reader = csv.reader(f)
            for row in reader:
                rows.append(row)
        last_row = ['Totals']
        for b in range(len(rows[0])):
            total = 0
            if b > 0:
                for c in range(len(rows)):
                    if c > 0:
                        total += int(rows[c][b])
                last_row.append(total)
        for d in range(len(rows)):
            for e in range(len(rows[d])):
                if d > 0 and e > 0:
                    rows[d][e] = float(rows[d][e])
                    rows[d][e] = rows[d][e]/last_row[e]
        with open(a[:-11]+'_percent.csv', 'w') as f:
            writer = csv.writer(f)
            for row in rows:
                writer.writerow(row)
        new_files.append(a[:-11]+'_percent.csv')
    return new_files
    
convert_shared('stability.final.pick.shared', 'Shared.csv')
names, indices, all_samples, all_treatments = get_meta_treatments('16S_metadata.csv')
files = treatment_shared(names, all_samples)
#clean_shared(files)
clean_meta(names, indices, all_samples, all_treatments, '16S_metadata.csv')
new_files = convert_percentages(files)
convert_for_r_file(new_files)

#now run the r script!