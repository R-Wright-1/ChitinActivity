import csv
import numpy
import os

keggs_chitosan = ['K01233', 'K14659']
keggs_chitin = ['K01452', 'K01183', 'K03791']
keggs_GlcNAc = ['K13020', 'K13007', 'K13018', 'K13019', 'K13016', 'K13017', 'K13015', 'K02473', 'K02474', 'K02852']
keggs_cellulose = ['K01225', 'K00702', 'K01180', 'K01181', 'K01179', 'K00700', 'K01811', 'K01199', 'K05349', 'K01182', 'K01187', 'K01188', 'K01210', 'K01215', 'K01232']

keggs_chitin_chitobiose_GlcNAc = ['K01183', 'K20547', 'K13381', 'K01207', 'K12373']
keggs_chitobiose_GlcNAc = ['K01207', 'K12373']
keggs_chitin_chitosan_glucosaminide_GlcN = ['K01452', 'K01233 ', 'K15855']
keggs_GlcNAc_GlcNAc6P_GlcN6P_Fru6P = ['K00884', 'K01443', 'K00844', 'K18676', 'K02564', ]
files = ['Acidimicrobiia.csv', 'Actinobacteria.csv', 'Alphaproteobacteria.csv', 'Anaerolinea.csv', 'Bacilli.csv', 'Bacteria.csv', 'Bacteroidia.csv', 'Campylobacteria.csv', 'Chlamydiae.csv', 'Cloacimonadia.csv', 'Clostridia.csv', 'Deltaproteobacteria.csv', 'Erysipelotrichia.csv', 'Gammaproteobacteria.csv', 'Gracilibacteria.csv', 'Ignavibacteria.csv', 'Microgenomatia.csv', 'Mollicutes.csv', 'Negativicutes.csv', 'Parcubacteria.csv', 'Phycisphaerae.csv', 'Planctomycetacia.csv', 'Saccharimonadia.csv', 'Spirochaetia.csv', 'Verrucomicrobiae.csv']
keggs = ['K01207', 'K12373']

new_file = []
for fi in files:
    with open(fi, 'rU') as f:
        rows = []
        for row in csv.reader(f):
            rows.append(row)
    this_file = []
    for a in range(len(rows)):
        del rows[a][-1]
        if a > 0:
            r0 = rows[a][0]
            for b in range(len(keggs)):
                if r0 == keggs[b]:
                    this_file.append(rows[a])
    new_row = [fi[:-4]]
    for c in range(len(this_file[0])):
        if c > 0:
            count = 0
            for d in range(len(this_file)):
                count += int(this_file[d][c])
            new_row.append(count)
    new_file.append(new_row)

with open('New_file_chitobiose_GlcNAc.csv', 'w') as f:
    writer = csv.writer(f)
    writer.writerow(rows[1])
    for a in range(len(new_file)):
        writer.writerow(new_file[a])