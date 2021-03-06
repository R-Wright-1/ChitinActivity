import csv
import numpy

fn_16S = '16S_shared.csv'
fasta_16S = '16S_clean_repFasta.fasta'



def get_OTU(fn):
    with open(fn, 'rU') as f:
        rows = []
        for row in csv.reader(f):
            rows.append(row)
    OTU, all_max = [], []
    for a in range(len(rows)):
        if a > 0:
            maximum = 0
            for b in range(len(rows[a])):
                if b == 0:
                    OTU.append(rows[a][b])
                if b > 0:
                    num = int(rows[a][b])
                    if num > maximum:
                        maximum = num
            all_max.append(maximum)
    new_OTU = []
    for c in range(len(all_max)):
        if all_max[c] >= 10:
            new_OTU.append(OTU[c])
    return new_OTU
    
def get_new_fasta(fn, fasta):
    OTU = get_OTU(fn)
    S = fn[0:3]
    f_asta, this_row, grouped_fasta = [], [], []
    with open(fasta, 'rU') as f:
        for row in f:
            f_asta.append(row)
    for r in range(len(f_asta)):
        row = f_asta[r]
        
        if row[0:4] != '>Otu':
            new_sequence = ''
            for s in row:
                if s != '-':
                    new_sequence += s
            this_row.append(new_sequence)
            grouped_fasta.append(this_row)
            this_row = []
        else:
            this_row.append(row)
    new_fasta, fasta_otus = [], []
    new_otus = []
    for a in grouped_fasta:
        fasta_otu = a[0][1:-1]
        fasta_otus.append(fasta_otu)
        for b in OTU:
            if fasta_otu == b:
                new_fasta.append(a)
                new_otus.append(fasta_otu)
    print new_otus

    with open(S+'_new_fasta_1.fasta', 'w') as f:
        for row in new_fasta:
            f.write(row[0]+'\n'+row[1]+'\n')
    return S+'_new_fasta_1.fasta'
#get_new_fasta(fn_18S, fasta_18S)
    
def get_both_fasta(fn_16S, fasta_16S):
    new_16S = get_new_fasta(fn_16S, fasta_16S)
    rows_16S, rows_18S = [], []
    with open(new_16S, 'rU') as f:
        for row in csv.reader(f):
            rows_16S.append(row)
    return
    
get_both_fasta(fn_16S, fasta_16S)
