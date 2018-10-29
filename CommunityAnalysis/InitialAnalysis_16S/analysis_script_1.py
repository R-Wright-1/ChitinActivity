import csv

#run this script first!
meta = '16S_metadata.csv'
tax_seqs = 'taxa_and_seqs.csv'

def get_meta_treatments(fn, table):
    #get all of the treatments, the new sample names, and the old sample names
    with open(fn, 'rU') as f:
        rows = []
        for row in csv.reader(f):
            rows.append(row)
    meta_treatments = []
    for a in range(len(rows[0])):
        if a > 0:
            treatment = rows[0][a]
            these_samples = []
            sample_names = []
            for b in range(len(rows)):
                if b > 0:
                    if rows[b][a] != '':
                        these_samples.append(rows[b][a])
                        sample_names.append(rows[b][0])
            meta_treatments.append([treatment, these_samples, sample_names])
    #open the taxa_and_seqs file
    with open(table, 'rU') as f:
        rows = []
        for row in csv.reader(f):
            rows.append(row)
    #first get all of the sequences, number these as ASV's, and save the taxonomy associated
    #with each one
    sequences, names, taxa = [], [], []
    for z in range(len(rows)):
        sequences.append(rows[z][0])
        otu_num = str(z)
        if len(otu_num) < 6:
            zeros = 6-len(otu_num)
            new_str = ''
            for a in range(zeros):
                new_str += '0'
            otu_num = new_str+otu_num
        name = 'ASV'+otu_num
        names.append(name)
        taxa.append(rows[z][1:8])
        rows[z][0] = name
        del rows[z][1:8]
    rows[0][0] = 'Group'
    #write the taxonomy file, with details of the DADA classification of each ASV
    with open('Taxonomy.csv', 'w') as f:
        writer = csv.writer(f)
        for a in range(len(taxa)):
            if a > 0:
                writer.writerow([names[a]]+taxa[a])
            else:
                writer.writerow(['Name']+taxa[a])
    #write a .fasta file with a list of all of the sequences
    with open('Sequences.fasta', 'w') as f:
        for a in range(len(sequences)):
            if a > 0:
                asv = 'ASV'+names[a][5:]
                name = '>'+asv
                name += '\n'
                f.write(name)
                f.write(sequences[a])
                f.write('\n')
    transpose = []
    for a in range(len(rows[0])):
        this_col = []
        for b in range(len(rows)):
            this_col.append(rows[b][a])
        transpose.append(this_col)
    names = []
    #now go through and save separate files for each treatment, with only those samples
    for c in range(len(meta_treatments)):
        mt = meta_treatments[c]
        #first save a new meta file, with the old and new treatment names
        with open(mt[0]+'_meta.csv', 'w') as f:
            writer = csv.writer(f)
            writer.writerow(['sample', 'treatment'])
            for d in range(len(mt[1])):
                writer.writerow([mt[2][d], mt[1][d]])
        #now get all of the samples that we need for this file
        this_meta = []
        for e in range(len(transpose)):
            if e > 0:
                for f in range(len(mt[2])):
                    if transpose[e][0] == mt[2][f]:
                        this_meta.append(transpose[e])
            else:
                this_meta.append(transpose[e])
        trans_back = []
        old_meta = this_meta
        for g in range(len(old_meta[0])):
            this_row = []
            for h in range(len(old_meta)):
                this_row.append(old_meta[h][g])
            trans_back.append(this_row)
        #and save the file (still with old names)
        with open(mt[0]+'_samples.csv', 'w') as f:
            writer = csv.writer(f)
            for row in trans_back:
                writer.writerow(row)
        names.append(mt[0]+'_samples.csv')
    return  names
    
def transpose(rows): #this just flips rows and columns with each other
    cols = []
    for a in range(len(rows[0])):
        col = []
        for b in range(len(rows)):
            col.append(rows[b][a])
        cols.append(col)
    return cols
    
def get_simper_files(names):
    #this goes through and changes the names of every 'ASV' to be 'Otu', as is needed for the R
    #simper script. This also needs the file to be transposed, and abundance to be relative
    for z in range(len(names)):
        #open previous file
        with open(names[z], 'rU') as f:
            rows = []
            for row in csv.reader(f):
                rows.append(row)
        #change the name to 'Otu'
        for a in range(len(rows)):
            if a > 0:
                rows[a][0] = 'Otu'+rows[a][0][3:]
        rows = transpose(rows)
        #convert each abundance value to a percentage
        for a in range(len(rows)):
            if a > 0:
                for b in range(len(rows[a])):
                    if b > 0:
                        rows[a][b] = float(rows[a][b])
                total = sum(rows[a][1:])
                for c in range(len(rows[a])):
                    if c > 0:
                        rows[a][c] = (rows[a][c]/total)*100
        #save the new file
        with open(names[z][:-4]+'_for_simper_r.csv', 'w') as f:
            writer = csv.writer(f)
            for a in range(len(rows)):
                writer.writerow(rows[a])
    return
    
def get_percent(names):
    #this goes through and calculates relative abundance for each ASV
    for z in range(len(names)):
        #open previous file
        with open(names[z], 'rU') as f:
            rows = []
            for row in csv.reader(f):
                rows.append(row)
        rows = transpose(rows)
        #calculate relative abundace
        for a in range(len(rows)):
            if a > 0:
                for b in range(len(rows[a])):
                    if b > 0:
                        rows[a][b] = float(rows[a][b])
                total = sum(rows[a][1:])
                for c in range(len(rows[a])):
                    if c > 0:
                        rows[a][c] = (rows[a][c]/total)*100
        rows = transpose(rows)
        #save the new file
        with open(names[z][:-11]+'percent.csv', 'w') as f:
            writer = csv.writer(f)
            for a in range(len(rows)):
                writer.writerow(rows[a])
    return
names = get_meta_treatments(meta, tax_seqs)
get_simper_files(names)
get_percent(names)