import csv
import numpy
import os

this_dir = '/Users/u1560915/Documents/GitHub/ChitinActivity/CommunityAnalysis/PICRUSt/all_PICRUSt/PICRUSt/Sort_by_class/'
classes = ['Mollicutes', 'Chlamydiae', 'Bacteroidia', 'Actinobacteria', 'Anaerolineae', 'Parcubacteria', 'NA', 'Erysipelotrichia', 'Spirochaetia', 'Woesearchaeia_None', 'Negativicutes', 'Gammaproteobacteria', 'Other_None', 'Microgenomatia', 'Bacilli', 'Clostridia', 'Gracilibacteria', 'Campylobacteria', 'Saccharimonadia', 'Phycisphaerae', 'Babeliae_None', 'Ignavibacteria', 'Chitinivibrionia_none', 'Deltaproteobacteria', 'Verrucomicrobiae', 'Acidimicrobiia', 'Cloacimonadia', 'Planctomycetacia', 'Alphaproteobacteria']
classes = sorted(classes)
print(classes)

for a in range(len(classes)):
    os.chdir(this_dir+classes[a]+'/')
    if classes[a][-4:] != 'none' and classes[a][-4:] != 'None':
        with open('otu_table.csv', 'rU') as f:
            rows = []
            for row in csv.reader(f):
                rows.append(row)
        ASVs_predicted = rows[1][1:]
        classes[a] = [classes[a], ASVs_predicted]
    else:
        classes[a] = classes[a][:-5]
        classes[a] = [classes[a], 'None']

os.chdir(this_dir)
with open('taxa_and_seqs_LSD.csv', 'rU') as f:
    rows = []
    for row in csv.reader(f):
        rows.append(row)

classes_long = []
for a in range(len(rows)):
    if a > 0:
        adding = True
        for b in range(len(classes_long)):
            if rows[a][1] == classes_long[b]:
                adding = False
        if adding:
            classes_long.append(rows[a][1])

new_classes = []
for a in range(len(classes_long)):
    adding = False
    for b in range(len(classes)):
        if classes_long[a] == classes[b][0]:
            adding = True
    if adding:
        new_classes.append(classes_long[a])
new_classes.append('Other')

new_rows = []
count = 0
for a in range(len(new_classes)):
    new_rows.append([])
for b in range(len(rows)):
    if b > 0:
        ind = -1
        adding = True
        for c in range(len(new_classes)):
            if rows[b][1] == new_classes[c]:
                ind = c
                adding = True
        if adding != True:
            rows[b][1] = 'Other'
        new_rows[ind].append(rows[b])
       
row_info = [['Class', 'ASVs predicted', 'Total ASVs', '', '', 'Original', '', '', '', '', '', 'PICRUSt', '', '', '']]
row_info.append(['', '', '', 'Long', 'Short', 'Daily', 'Long_R', 'Short_R', 'Daily_R', 'Long', 'Short', 'Daily', 'Long_R', 'Short_R', 'Daily_R'])

print(len(new_rows), len(new_classes))

for a in range(len(classes)):
    for b in range(len(new_classes)):
        if classes[a][0] == new_classes[b]:
            this_predic = classes[a][1]
            this_abun = new_rows[b]
            old_abun = [0, 0, 0, 0, 0, 0]
            for c in range(len(this_abun)):
                for d in range(len(old_abun)):
                    this_abun[c][d+2] = float(this_abun[c][d+2])
                    old_abun[d] += this_abun[c][d+2]
            print(classes[a][0], new_classes[b])
            if this_predic == 'None':
                new_abun = [0, 0, 0, 0, 0, 0]
                row_info.append([classes[a][0], str(0), str(len(this_abun))]+old_abun+new_abun)
            else:
                predicted_abun = []
                new_abun = [0, 0, 0, 0, 0, 0]
                for e in range(len(this_abun)):
                    for f in range(len(this_predic)):
                        if this_abun[e][0] == this_predic[f]:
                            predicted_abun.append(this_abun[e])
                for g in range(len(predicted_abun)):
                    for h in range(len(new_abun)):
                        new_abun[h] += (predicted_abun[g][h+2])
                row_info.append([classes[a][0], str(len(predicted_abun)), str(len(this_abun))]+old_abun+new_abun)
                
with open('Predicted_ASV_info.csv', 'w') as f:
    writer = csv.writer(f)
    for a in range(len(row_info)):
        writer.writerow(row_info[a])
