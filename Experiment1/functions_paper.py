import csv
import os
import numpy
import matplotlib.pyplot as plt
import scipy.stats as stats
import pandas as pd
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
import statsmodels.api as sm
from sklearn.pipeline import make_pipeline

bdir = '/Users/u1560915/Documents/GitHub/ChitinActivity/Experiment1/'

def read_next_generation():
    os.chdir(bdir)
    rows, new_rows = [], []
    l, dr, dg, s6, drs, dgs = [], [], [], [], [], []
    for a in range(21):
        l.append([]), dr.append([]), dg.append([]), s6.append([]), drs.append([]), dgs.append([])
    lng, short = [l, dr, dg, s6], [drs, dgs]
    with open('Next generation.csv', 'rU') as f:
        reader = csv.reader(f)
        for row in reader:
            rows.append(row)
    for a in range(len(rows)):
        if rows[a] != [] and rows[a] != '':
            new_rows.append(rows[a])
    rows = new_rows
    del rows[0]
    for b in range(len(rows)):
        row = rows[b]
        if row[0] != '':
            gen, ls, comm, i1, i2, i3 = int(row[0]), row[1], int(row[2])-1, int(row[3])-1, int(row[4])-1, int(row[5])-1
            if ls == 'L':
                lng[comm][gen] = [i1, i2, i3]
            elif ls == 'S':
                short[comm-1][gen] = [i1, i2, i3]
    l, dr, dg, s6, drs, dgs = lng[0], lng[1], lng[2], lng[3], short[0], short[1]
    for a in range(15):
        del drs[0]
        del dgs[0]
    return l, dr, dg, s6, drs, dgs
    
def read_plate(fn):
    r, c = 27, 2
    rows, plate = [], []
    os.chdir(bdir+'data/')
    with open(fn, 'rU') as f:
        reader = csv.reader(f)
        for row in reader:
            rows.append(row)
        for i in range(8):
            for j in range(12):
                plate.append(int(rows[r+i][c+j]))
    standards = plate[0:6]
    samples = plate[6:96]
    new_s = [standards[5], standards[4], standards[3], standards[2], standards[1], standards[0]]
    return new_s, samples
    
def read_plate_txt(fn):
    r, c = 35, 1
    rows, plate = [], []
    os.chdir(bdir+'daily/')
    with open(fn, 'rU') as f:
        reader = csv.reader(f)
        for row in reader:
            rows.append(row)
        for i in range(8):
            for j in range(12):
                plate.append(int(rows[r+i][c+j]))
    standards = plate[0:6]
    samples = plate[6:96]
    new_s = [standards[5], standards[4], standards[3], standards[2], standards[1], standards[0]]
    return new_s, samples

def get_triplicate(samples):
    triplicate, plate_trip = [], []
    for a in range(90):
        triplicate.append(samples[a])
        if (a+1) % 3 == 0:
            diff_1_2, diff_2_3, diff_1_3 = abs(triplicate[0]-triplicate[1]), abs(triplicate[1]-triplicate[2]), abs(triplicate[0]-triplicate[2])
            trip_mean = (numpy.mean(triplicate))*0.05
            trip_mean_2 = (numpy.mean(triplicate))*0.10
            if diff_1_2 < trip_mean and diff_2_3 < trip_mean and diff_1_3 < trip_mean:
                vals = triplicate
            elif diff_1_2 > trip_mean_2 and diff_2_3 > trip_mean_2 and diff_1_3 > trip_mean_2:
                vals = triplicate
            elif diff_1_2 < diff_2_3 and diff_1_2 < diff_1_3:
                vals = [triplicate[0], triplicate[1]]
            elif diff_2_3 < diff_1_2 and diff_2_3 < diff_1_3:
                vals = [triplicate[1], triplicate[2]]
            elif diff_1_3 < diff_1_2 and diff_1_3 < diff_2_3:
                vals = [triplicate[0], triplicate[2]]
            else:
                vals = triplicate
            plate_trip.append(numpy.mean(vals))
            triplicate = []
    return plate_trip

def get_this_file(fn):
    standards, samples = read_plate(fn)
    community = get_triplicate(samples)
    return standards, community
    
def get_this_file_txt(fn):
    standards, samples = read_plate_txt(fn)
    community = get_triplicate(samples)
    return standards, community

count = 0 
def get_standards(standards):
    this_standard = [[], [], [], [], []]
    for a in range(len(standards)):
        for b in range(5):
            this_standard[b].append(standards[a][b])
    standards, errors = [], []
    for c in range(len(this_standard)):
        standards.append(numpy.mean(this_standard[c]))
        errors.append(numpy.std(this_standard[c]))
    return standards, errors
    
def get_standards_daily(standards):
    this_standard = [[], [], [], [], []]
    for a in range(len(standards)):
        for b in range(5):
            this_standard[b].append(standards[a][b])
    standards, errors = [], []
    for c in range(len(this_standard)):
        standards.append(numpy.mean(this_standard[c]))
        errors.append(numpy.std(this_standard[c]))
    return standards, errors
    
def normalise_community(community, standards):
    concs = [0, 0.00016, 0.0008, 0.004, 0.02]
    z = numpy.polyfit(standards, concs, 4)
    f = numpy.poly1d(z)
    gradient, intercept, r_value, p_value, std_err = stats.linregress(concs, standards)
    community = f(community)
    for a in range(len(community)):
        community[a] = (community[a]*1000*1.44)
    return community
    
def transform_data(l, dr, dg, s6, drs, dgs, standards_l, standards_s):
    lng, short = [l, dr, dg, s6], [drs, dgs]
    ls, le, ss, se = [], [], [], []
    for a in range(len(standards_l)):
        standards, errors = get_standards(standards_l[a])
        ls.append(standards)
        le.append(errors)
    for b in range(len(standards_s)):
        standards, errors = get_standards(standards_s[b])
        ss.append(standards)
        se.append(errors)
    for c in range(len(lng)):
        for d in range(len(lng[c])):
            standards, community = ls[c], lng[c][d]
            community = normalise_community(community, standards)
            lng[c][d] = community
    for e in range(len(short)):
        for f in range(len(short[e])):
            standards, community = ss[c], short[e][f]
            community = normalise_community(community, standards)
            short[e][f] = community
    l, dr, dg, s6, drs, dgs = lng[0], lng[1], lng[2], lng[3], short[0], short[1]
    return l, dr, dg, s6, drs, dgs, ls, le, ss, se

def get_all_data():
    folder = bdir+'data/'
    all_files = []
    for item in os.listdir(folder):
        all_files.append(item)
    all_files = sorted(all_files)
    l, dr, dg, s6, drs, dgs, standards_l, standards_s = [], [], [], [], [], [], [], []
    for a in range(21):
        l.append([]), dr.append([]), dg.append([]), s6.append([]), drs.append([]), dgs.append([]), standards_l.append([]), standards_s.append([])
    lng, short = [l, dr, dg, s6], [drs, dgs]
    for b in range(len(all_files)):
        if all_files[b].endswith('.csv'):
            gen = ''
            gen += all_files[b][1]
            if all_files[b][2] != '_': 
                gen += all_files[b][2]
            if all_files[b][3] != '_': 
                comm = all_files[b][3]
            else: 
                comm = all_files[b][4]
            gen, comm = int(gen), int(comm)
            if all_files[b][6] == 'S':
                standards, community = get_this_file(all_files[b])
                comm -= 2
                standards_s[gen].append(standards)
                short[comm][gen] = community
            else:
                standards, community = get_this_file(all_files[b])
                comm -= 1
                standards_l[gen].append(standards)
                lng[comm][gen] = community
    for c in range(15):
        del drs[0]
        del dgs[0]
        del standards_s[0]
    l, dr, dg, s6, drs, dgs = lng[0], lng[1], lng[2], lng[3], short[0], short[1]
    l, dr, dg, s6, drs, dgs, ls, le, ss, se = transform_data(l, dr, dg, s6, drs, dgs, standards_l, standards_s)
    drs[0], dgs[0], ss[0], se[0] = dr[15], dg[15], ls[15], le[15]
    return l, dr, dg, s6, drs, dgs, ls, le, ss, se
    
def get_means(l, dr, dg, s6, drs, dgs):
    os.chdir(bdir+'figures/')
    communities = [l, dr, dg, s6, drs, dgs]
    new_means = [[], [], [], [], [], []]
    new_errors = [[], [], [], [], [], []]
    for a in range(len(communities)):
        for b in range(len(communities[a])):
            new_means[a].append(numpy.mean(communities[a][b]))
            new_errors[a].append(numpy.std(communities[a][b]))
        for c in range(len(new_means[a])):
            if a == 4 or a == 5:
                new_means[a][c] -= new_means[a-3][0]
            else:
                new_means[a][c] -= new_means[a][0]
    return new_means, new_errors
    
def plot_standards(ls, le, ss, se):
    os.chdir(bdir+'figures/')
    colors = ['#330000', '#990000', '#FF0000', '#FF0066', '#990066', '#330066', '#003399', '#003366', '#0099CC', '#00CCCC', '#009999',
          '#006666', '#006633', '#006600', '#339900', '#66CC00', '#CCCC00', '#FFFF00', '#FFCC00', '#FF9900', '#FF6600']
    concs = [0, 0.00016, 0.0008, 0.004, 0.02]
    ax1 = plt.subplot(111)
    for a in range(len(ls)):
        this_standard, this_error = ls[a], le[a]
        label = 'CA'+str(a)
        gradient, intercept, r_value, p_value, std_err = stats.linregress(concs, ls[a])
        r2 = r_value**2
        r2 = str("%.3f" % r2)
        label += r', $r^2$='
        label += r2
        ax1.errorbar(concs, this_standard, yerr=this_error, marker='o', color=colors[a], label=label)
    ax1.set_ylabel('Fluorescence')
    ax1.set_xlabel(r'Chitinase U $ml^{-1}$')
    plt.legend(bbox_to_anchor=(1.4, 1.05), fontsize=8)
    plt.tight_layout()
    ax1.set_xlim([0, 0.0045])
    plt.savefig('Standards long.pdf', bbox_inches='tight')
    plt.close()
    ax1 = plt.subplot(111)
    for a in range(len(ss)):
        this_standard, this_error = ss[a], se[a]
        label = 'CA'+str(a+15)
        gradient, intercept, r_value, p_value, std_err = stats.linregress(concs, ss[a])
        r2 = r_value**2
        r2 = str("%.3f" % r2)
        label += r', $r^2$='
        label += r2
        ax1.errorbar(concs, this_standard, yerr=this_error, marker='o', color=colors[a+15], label=label)
    ax1.set_ylabel('Fluorescence')
    ax1.set_xlabel(r'Chitinase U $ml^{-1}$')
    plt.legend(bbox_to_anchor=(1.4, 1.05), fontsize=8)
    plt.tight_layout()
    ax1.set_xlim([0, 0.0045])
    plt.savefig('Standards short.pdf', bbox_inches='tight')
    plt.close()
    for a in range(len(ls)):
        ax1 = plt.subplot(111)
        this_standard, this_error = ls[a], le[a]
        label = 'CA'+str(a)
        z = numpy.polyfit(this_standard, concs, 4)
        f = numpy.poly1d(z)
        x_new = f(this_standard)
        ax1.plot(concs, this_standard, 'o-')
        ax1.plot(x_new, this_standard, 'r')
        ax1.set_ylabel('Fluorescence')
        ax1.set_xlabel(r'Chitinase U $ml^{-1}$')
        plt.savefig(label+' long.pdf', bbox_inches='tight')
        plt.close()  
    for a in range(len(ss)):
        ax1 = plt.subplot(111)
        this_standard, this_error = ss[a], se[a]
        label = 'CA'+str(a+15)
        z = numpy.polyfit(this_standard, concs, 4)
        f = numpy.poly1d(z)
        x_new = f(this_standard)
        ax1.plot(concs, this_standard, 'o-')
        ax1.plot(x_new, this_standard, 'r')
        ax1.set_ylabel('Fluorescence')
        ax1.set_xlabel(r'Chitinase U $ml^{-1}$')
        plt.savefig(label+' short.pdf', bbox_inches='tight')
        plt.close()  
    return
    
def highest_means(l, dr, dg, s6, drs, dgs):
    communities = [l, dr, dg, s6, drs, dgs]
    li, dri, dgi, s6i, drsi, dgsi = read_next_generation()
    indices = [li, dri, dgi, s6i, drsi, dgsi]
    high_means, high_errors = [[], [], [], [], [], []], [[], [], [], [], [], []]
    high_all = [[], [], [], [], [], []]
    for a in range(len(indices)):
        for b in range(len(indices[a])):
            i1, i2, i3 = indices[a][b][0], indices[a][b][1], indices[a][b][2]
            community = communities[a][b]
            h1, h2, h3 = community[i1], community[i2], community[i3]
            high_all[a].append([h1, h2, h3])
            mean, error = numpy.mean([h1, h2, h3]), numpy.std([h1, h2, h3])
            high_means[a].append(mean), high_errors[a].append(error)    
        for c in range(len(high_means[a])):
            if a == 4 or a == 5:
                high_means[a][c] -= high_means[a-3][0]
            else:
                high_means[a][c] -= high_means[a][0]
    community_names = ['Light', 'Random', 'Good', 'Strain 6', 'Short random', 'Short good']
    """
    with open('High_means.csv', 'a') as f:
        writer = csv.writer(f)
        writing = ['Generation', 'Community', 'Means', 'Errors']
        writer.writerow(writing)
        for a in range(len(high_means)):
            for b in range(len(high_means[a])):
                writing = [b, community_names[a], high_means[a][b], high_errors[a][b]]
                writer.writerow(writing)
    """
    return high_means, high_errors, high_all
    
def plot_all_means(means, errors, is_high):
    #lm, sm = [means[0], means[1], means[2], means[3]], [means[4], means[5]]
    lm, sm = [means[0], means[1], means[2]], [means[4], means[5]]
    #le, se = [errors[0], errors[1], errors[2], errors[3]], [errors[4], errors[5]]
    le, se = [errors[0], errors[1], errors[2]], [errors[4], errors[5]]
    #print(is_high, 'lm = ', lm)
    #print(is_high, 'sm = ', sm)
    #print(is_high, 'le = ', le)
    #print(is_high, 'se = ', se)
    gens_l, gens_s = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20], [15, 16, 17, 18, 19, 20]
    colors_l, colors_s = ['#FFCC00', '#006699', '#CC0000', '#006600'], ['#0099FF', '#CC0099']
    labels = ['Light', 'Dark random', 'Selection (9 day)', 'Strain 6']
    labels_s = ['Dark random, \n short', 'Shortened selection \n (4 day)']
    ax1 = plt.subplot(111)
    for a in range(len(lm)):
        ax1.errorbar(gens_l, lm[a], yerr=le[a], marker='o', label=labels[a], color=colors_l[a])
    for b in range(len(sm)):
        ax1.errorbar(gens_s, sm[b], yerr=se[b], marker='o', label=labels_s[b], color=colors_s[b])
    ax1.legend(bbox_to_anchor=(1.5, 1.05))
    ax1.set_xlim([0, 21])
    ax1.set_xlabel('Generation')
    ax1.set_ylabel(r'Chitinase $\mu$M day$^{-1}$')
    plt.tight_layout()
    if is_high == 'high':
        plt.savefig('All highest no S6.pdf', bbox_inches='tight')
    else:
        plt.savefig('All means no S6.pdf', bbox_inches='tight')
    plt.close()
    colors_l[2] = '#48C9B0'
    colors_s[1] = '#F39C12'
    long_norm, short_norm = [], []
    gens_l_new = []
    plotting_here = []
    
    for c in range(len(lm[2])):
    #for c in range(8):
        long_norm.append(lm[2][c]-lm[1][c])
        plotting_here.append(c)
    for d in range(len(sm[1])):
        short_norm.append(sm[1][d]-sm[0][d])
    short_norm[0] = long_norm[15]
    ax1 = plt.subplot(111)
    ax1.plot(gens_l, long_norm, marker='o', label=labels[2], color=colors_l[2])
    #ax1.plot(plotting_here, long_norm, marker='o', label=labels[2], color=colors_l[2])
    #ax1.plot(gens_l, long_norm, marker='o', label=labels[2], color=colors_l[2])
    ax1.plot(gens_s, short_norm, marker='o', label=labels_s[1], color=colors_s[1])
    LS = [long_norm, short_norm]
    g = [gens_l, gens_s]
    ax1.legend(bbox_to_anchor=(1.52, 1.05))
    #ax1.set_xlim([0, 11])
    ax1.set_xlim([0, 21])
    #ax1.set_ylim([-0.00005, 0.0004])
    ax1.set_xticks(gens_l)
    #ax1.set_xticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    ax1.plot([0, 21], [0, 0], 'k--')
    gradient, intercept, r_value, p_value, std_err = stats.linregress(gens_l, long_norm)
    ax1.plot([0, 20], [intercept, (gradient*20 + intercept)], '-.', color=colors_l[2])
    gradient, intercept, r_value, p_value, std_err = stats.linregress(gens_s, short_norm)
    ax1.plot([15, 20], [(gradient*15 + intercept), (gradient*20 + intercept)], '-.', color=colors_s[1])
    ax1.set_xlabel('Generation')
    ax1.set_ylabel(r'Chitinase $\mu$M day$^{-1}$')
    plt.tight_layout()
    if is_high == 'high':
        plt.savefig('Normalised high.pdf', bbox_inches='tight')
    else:
        plt.savefig('Normalised means.pdf', bbox_inches='tight')
    plt.close()
    long_norm, short_norm = [], []
    gens_l_new = []
    for c in range(len(lm[2])):
        long_norm.append(lm[2][c]-lm[0][c])
    ax1 = plt.subplot(111)
    ax1.plot(gens_l, long_norm, marker='o', label=labels[2], color=colors_l[2])
    ax1.legend(bbox_to_anchor=(1.4, 1.05))
    ax1.set_xlim([0, 21])
    ax1.set_xticks(gens_l_new)
    ax1.plot([0, 21], [0, 0], 'k--')
    gradient, intercept, r_value, p_value, std_err = stats.linregress(gens_l, long_norm)
    ax1.plot([0, 20], [intercept, (gradient*20 + intercept)], '-.', color=colors_l[2])
    ax1.set_xlabel('Generation')
    ax1.set_ylabel(r'Chitinase $\mu$M day$^{-1}$')
    plt.tight_layout()
    if is_high == 'high':
        plt.savefig('Normalised high light.pdf', bbox_inches='tight')
    else:
        plt.savefig('Normalised means light.pdf', bbox_inches='tight')
    plt.close()
    return LS, g
    
    
def plot_short(means, errors):
    os.chdir(bdir+'figures/')
    gens = [15, 16, 17, 18, 19, 20]
    colors_s = ['#0099FF', '#CC0099']
    ax1 = plt.subplot(111)
    ax1.set_title('Means shortened generation time')
    ax1.set_xlabel('Generation')
    ax1.set_ylabel(r'Chitinase $\mu$M day$^{-1}$')
    ax1.errorbar(gens, means[0], yerr=errors[0], marker='o', linestyle='--', color=colors_s[0], label='Control')
    ax1.errorbar(gens, means[1], yerr=errors[1], marker='o', color=colors_s[1], label='Selection')
    ax1.legend(bbox_to_anchor=(1.4, 1.05))
    ax1.set_xlim([14.5, 20.5])
    plt.tight_layout()
    plt.savefig('Means shortened only.pdf', bbox_inches='tight')
    plt.close()
    norm_means = []
    for a in range(len(means[0])):
        norm_means.append(means[1][a]-means[0][a])
    ax1 = plt.subplot(111)
    ax1.set_title('Means shortened generation time')
    ax1.set_xlabel('Generation')
    ax1.set_ylabel(r'Chitinase $\mu$M day$^{-1}$')
    del gens[0]
    del norm_means[0]
    ax1.errorbar(gens, norm_means, marker='o', color=colors_s[1], label='Selection')
    ax1.set_xlim([15.5, 20.5])
    plt.tight_layout()
    plt.savefig('Means shortened only normalised.pdf', bbox_inches='tight')
    plt.close()
    return

def get_difference_with_errors(first, second):
    #first = good, second = random
    mean_second = []
    for a in range(len(second)):
        mean_second.append(numpy.mean(second[a]))
    difference, error = [], []
    for b in range(len(first)):
        ms = mean_second[b]
        this_gen = []
        for c in range(len(first[b])):
            this_gen.append(first[b][c]-ms)
        difference.append(numpy.mean(this_gen))
        error.append(numpy.std(this_gen))
    return difference, error
    
def get_daily():
    os.chdir('/Users/u1560915/Documents/GitHub/ChitinActivity/Experiment1/daily/')
    new_g15 = ['RW_C15_d1_1.csv', 'RW_C15_d1_2.csv', 'RW_C15_d1_3.csv', 'RW_C15_d1_4.csv',
               'RW_C15_d2_1.csv', 'RW_C15_d2_2.csv', 'RW_C15_d2_3.csv', 'RW_C15_d2_4.csv',
               'RW_C15_d4_1.csv', 'RW_C15_d4_2.csv', 'RW_C15_d4_3.csv', 'RW_C15_d4_4.csv',
               'RW_C15_d6_1.csv', 'RW_C15_d6_2.csv', 'RW_C15_d6_3.csv', 'RW_C15_d6_4.csv',
               'RW_C15_d9_1.csv', 'RW_C15_d9_2.csv', 'RW_C15_d9_3.csv', 'RW_C15_d9_4.csv',]  
    
    #g15 = ['RW_C15_d1_plate-2.csv', 'RW_C15_d1_plate-4.csv', 'RW_C15_d3_plate-2.csv', 'RW_C15_d3_plate-4.csv', 'RW_C15_d4_plate-2.csv', 'RW_C15_d4_plate-4.csv', 'RW_C15_d6_plate-2.csv', 'RW_C15_d6_plate-4.csv', 'RW_C15_d9_plate-2.csv', 'RW_C15_d9_plate-4.csv']
    g20 = ['RW_C20-4d_d1_plate-1.csv', 'RW_C20-4d_d1_plate-2.csv', 'RW_C20-4d_d2_plate-1.csv', 'RW_C20-4d_d2_plate-2.csv', 'RW_C20-4d_d3_plate-1.csv', 'RW_C20-4d_d3_plate-2.csv', 'RW_C20-4d_d4_plate-1.csv', 'RW_C20-4d_d4_plate-2.csv']
    #g20 = ['RW_C20-4d_d1_plate-1.csv', 'RW_C20-4d_d1_plate-2.csv', 'RW_C20-4d_d2_plate-1.csv', 'RW_C20-4d_d2_plate-2.csv', 'RW_C20-4d_d3_plate-1.csv', 'RW_C20-4d_d3_plate-2.csv', 'RW_C20-4d_after_plate-1.csv', 'RW_C20-4d_after_plate-2.csv']
    g15_s, g20_s = [], []
    for a in range(len(new_g15)):
        standards, community = get_this_file_txt(new_g15[a])
        new_g15[a] = community
        g15_s.append(standards)
    for a in range(len(g20)):
        standards, community = get_this_file_txt(g20[a])
        g20[a] = community
        g20_s.append(standards)
        
    r, g = 0, 3
    g15_r, g15_g = [new_g15[0+r], new_g15[4+r],new_g15[8+r],new_g15[12+r], new_g15[16+r]], [new_g15[0+g], new_g15[4+g],new_g15[8+g],new_g15[12+g], new_g15[16+g]]
    g20_r, g20_g = [g20[0], g20[2], g20[4], g20[6]], [g20[1], g20[3], g20[5], g20[7]]
    g15_s2, g15_e, g20_s2, g20_e = [], [], [], []
    for b in [0, 4, 8, 12, 16]:
        standards = [g15_s[b], g15_s[b+1], g15_s[b+2], g15_s[b+3]]
        standards, errors = get_standards_daily(standards)
        g15_s2.append(standards), g15_e.append(errors)
    for b in [0, 2, 4, 6]:
        c = b+1
        standards = [g20_s[b], g20_s[c]]
        standards, errors = get_standards_daily(standards)
        g20_s2.append(standards), g20_e.append(errors)
    for d in [0, 1, 2, 3, 4]:
        g15_r[d] = normalise_community(g15_r[d], g15_s2[d])
        g15_g[d] = normalise_community(g15_g[d], g15_s2[d])
    for d in [0, 1, 2, 3]:
        g20_r[d] = normalise_community(g20_r[d], g20_s2[d])
        g20_g[d] = normalise_community(g20_g[d], g20_s2[d])
    print('g20_g =', g20_g)
    print('g20_r =', g20_r)
    for z in range(len(g15_r)):
        for a in range(30):
            g15_r[z][a] = (g15_r[z][a])
            g15_g[z][a] = (g15_g[z][a])
    for z in range(len(g20_r)):
        for a in range(30):
            g20_r[z][a] = (g20_r[z][a])
            g20_g[z][a] = (g20_g[z][a])
    g15_rm, g15_gm, g20_rm, g20_gm = [], [], [], []
    g15_re, g15_ge, g20_re, g20_ge = [], [], [], []
    for e in g15_r:
        g15_rm.append(numpy.mean(e))
        g15_re.append(numpy.std(e))
    for e in g15_g:
        g15_gm.append(numpy.mean(e))
        g15_ge.append(numpy.std(e))
    for e in g20_r:
        g20_rm.append(numpy.mean(e))
        g20_re.append(numpy.std(e))
    for e in g20_g:
        g20_gm.append(numpy.mean(e))
        g20_ge.append(numpy.std(e))
    return [g15_rm, g15_gm, g20_rm, g20_gm], [g15_re, g15_ge, g20_re, g20_ge]
    
def get_paper_plot(LS, g):
    os.chdir(bdir+'figures/')
    col_l, col_s = '#48C9B0', '#F39C12'
    fig = plt.figure(figsize=(8, 6))
    ax1 = plt.subplot2grid((3, 20), (0, 0), colspan=13, rowspan=2)
    ax2 = plt.subplot2grid((3, 10), (0, 7), colspan=3)
    ax2.text(-3, -0.4, r'Chitinase $\mu$M day$^{-1}$', rotation=90, ha='center', va='center')
    ax3 = plt.subplot2grid((3, 10), (1, 7), colspan=3)
    [g15_r, g15_g, g20_r, g20_g], [g15_re, g15_ge, g20_re, g20_ge] = get_daily()
    os.chdir(bdir+'figures/')
    regr = make_pipeline(PolynomialFeatures(2), LinearRegression())
    with open('For regression.csv', 'w') as f:
        writer = csv.writer(f)
        writer.writerow(['X', 'Y'])
        for a in range(len(g[0])):
            writer.writerow([g[0][a], LS[0][a]])
    df = pd.read_csv('For regression.csv')
    regr.fit(df[['X']], df[['Y']])
    y_pred = regr.predict(df[['X']])
    est = sm.OLS(df[['Y']], y_pred)
    pval =  est.fit().f_pvalue
    r2 = est.fit().rsquared
    #ax1.plot(g[0], y_pred, '-.', color=col_s, linewidth=1)
    label='Selection (nine-day)'
    label2 = r'$r^2$'+'=%.2f'%r2
    label2 += ', '+r'$p$'+'=%.2f'%pval
    #ax1.plot(g[0], y_pred, '-.', color=col_l, linewidth=1)
    ax1.plot(g[0], LS[0], color=col_l, marker='o', label=label, markeredgecolor='k', linewidth=1)
    #ax1.plot(g[0], LS[0], color=col_l, label=label2, linewidth=1)
    #ax1.plot([14, 15, 16], [LS[0][14], g15_g[2]-g15_r[2], LS[0][16]], marker='o', color=col_l, markeredgecolor='k', linewidth=1)
    #ax1.annotate('C', xy=(19.8, LS[1][-1]), xytext=(18, LS[1][-1]-0.15),arrowprops={'arrowstyle':'->'})
    #ax1.text(17.5, LS[1][-1]-0.22, 'day 2', fontsize=8)
    label = 'Selection (four-day)'
    LS[1][-1] = g20_g[3]-g20_r[3]
    ax1.plot(g[1], LS[1], color=col_s, marker='^', label=label, markeredgecolor='k', linewidth=1)
    with open('For regression short day 2.csv', 'w') as f:
        writer = csv.writer(f)
        writer.writerow(['X', 'Y'])
        for a in range(len(g[1])):
            writer.writerow([g[1][a], LS[1][a]])
    df = pd.read_csv('For regression short day 2.csv')
    regr.fit(df[['X']], df[['Y']])
    y_pred = regr.predict(df[['X']])
    est = sm.OLS(df[['Y']], y_pred)
    pval_2day =  est.fit().f_pvalue
    r2_2day = est.fit().rsquared
    label2 = r'$r^2$'+'=%.2f'%r2_2day
    label2 += ', '+r'$p$'+'=%.2f'%pval_2day
    label2 += ' (day 2)'
    #ax1.plot(g[1], LS[1], color=col_s, label=label2, linewidth=1, markeredgecolor='k')
    do = '#cc5500'
    #do = '#51412d'
    #ax1.plot(g[1], y_pred, '-.', color=do, linewidth=1)
    
    day4 = LS[1]
    #diff = (g20_g[1]-g20_r[1])-(g20_g[3]-g20_r[3])
    day4[-1] = g20_g[3]-g20_r[3]
    ax1.annotate('C', xy=(19.8, day4[-1]), xytext=(19, day4[-1]-0.18),arrowprops={'arrowstyle':'->'})
    #ax1.text(18.5, day4[-1]-0.25, 'day 4', fontsize=8)
    #ax1.plot(g[1][-2:], day4[-2:], color=col_s, marker='o', markeredgecolor='k', linewidth=1)
    with open('For regression short day 4.csv', 'w') as f:
        writer = csv.writer(f)
        writer.writerow(['X', 'Y'])
        for a in range(len(g[1])):
            writer.writerow([g[1][a], day4[a]])
    df = pd.read_csv('For regression short day 4.csv')
    regr.fit(df[['X']], df[['Y']])
    y_pred = regr.predict(df[['X']])
    est = sm.OLS(df[['Y']], y_pred)
    pval_4day =  est.fit().f_pvalue
    r2_4day = est.fit().rsquared
    label3 = r'$r^2$'+'=%.2f'%r2_4day
    label3 += ', '+r'$p$'+'=%.2f'%pval_4day
    label3 += ' (day 4)'
    #ax1.plot(g[1][-2:], day4[-2:], color=col_s, label=label3, linewidth=1)  
    #ax1.plot(g[1], y_pred, '-.', color=do, linewidth=1)    
    
    ax1.annotate('B', xy=(14.8, LS[0][15]), xytext=(13, LS[0][15]-0.1),arrowprops={'arrowstyle':'->'})
    #ax1.text(12.5, LS[0][15]-0.17, 'day 9', fontsize=8)
    #ax1.annotate('B', xy=(14.8, g15_g[2]-g15_r[2]), xytext=(13, g15_g[2]-g15_r[2]-0.05),arrowprops={'arrowstyle':'->'})
    #ax1.text(12.5, g15_g[2]-g15_r[2]-0.1, 'day 4', fontsize=8)
    
    ax1.legend(loc='upper left', fontsize=8, frameon=False, numpoints=1, handlelength=0)
    ax1.plot([0, 20], [0, 0], 'k--', linewidth=1)
    ax2.plot([0, 10], [0, 0], 'k--', linewidth=1)
    ax3.plot([0, 5], [0, 0], 'k--', linewidth=1)
    ax1.set_xlim([-1, 21])
    ax1.set_ylim([-0.4, 0.3])
    ax1.set_xticks([0, 5, 10, 15, 20])
    days_20 = [1, 2, 3, 4]
    g = [g15_r, g15_g, g20_r, g20_g]
    ge = [g15_re, g15_ge, g20_re, g20_ge]
    colors = ['#66CCFF', '#FF0099']
    ax3.set_xlim([0, 5])
    ax3.errorbar(days_20, g[2], yerr=ge[2], color=colors[0], marker='s', markeredgecolor='k', capsize=2, linewidth=1)
    ax3.errorbar(days_20, g[3], yerr=ge[3], color=colors[1], marker='o', markeredgecolor='k', capsize=2, linewidth=1)
    days_15 = [1, 2, 4, 6, 9]
    ax2.errorbar(days_15, g[0], yerr=ge[0], color=colors[0], marker='s', label='Random selection', markeredgecolor='k', capsize=2, linewidth=1)
    ax2.errorbar(days_15, g[1], yerr=ge[1], color=colors[1], marker='o', label='Positive selection', markeredgecolor='k', capsize=2, linewidth=1)
    
    #print('sel_20 = ', g[3], 'sel_sd_20 = ', ge[3])
    #print('ran_20 = ', g[2], 'ran_sd_20 = ', ge[2])
    
    #print('sel_15 = ', g[1], 'sel_sd_15 = ', ge[1])
    #print('ran_15 = ', g[0], 'ran_sd_15 = ', ge[0])
    
    ax2.tick_params(axis='both', which='right', length=0)
    ax3.tick_params(axis='both', which='right', length=0)
    ax2.legend(loc='upper left', fontsize=8, frameon=False, numpoints=1, handlelength=0)
    ax1.tick_params(axis='both', which='right', length=0)
    ax3.tick_params(axis='both', which='right', length=0)
    ax1.set_ylabel(r'Chitinase $\mu$M day$^{-1}$')
    
    ax4 = ax3.twinx()
    short_good_g20 = [0.219, 0.689, 10.38, 3.646666667]
    short_good_g20_e = [0.036755952, 0.444635806, 5.191415992, 1.831402013]
    ax4.errorbar([1,2,3,4], short_good_g20, yerr=short_good_g20_e, linestyle='--', color='g', marker='^', label='DNA conc.', capsize=2, markeredgecolor='k')
    ax4.legend(loc='upper left', handlelength=0, fontsize=8)
    ax4.set_ylabel(r'DNA conc. $\mu$g mL$^{-1}$')
    
    
    #ax2.set_ylabel(r'Chitinase $\mu$M day$^{-1}$')
    #ax3.set_ylabel(r'Chitinase $\mu$M day$^{-1}$')
    ax1.set_xlabel('Generation')
    ax2.set_xlabel('Day')
    ax3.set_xlabel('Day')
    ax3.set_ylim([0, 1.5])
    ax2.set_ylim([0, 1])
    ax2.set_xlim([0, 10])
    ax3.set_yticks([0.0, 0.5, 1.0, 1.5])
    ax1.set_title('A', loc='left', fontweight='bold', fontsize=14)
    ax2.set_title('B', loc='left', fontweight='bold', fontsize=14)
    ax3.set_title('C', loc='left', fontweight='bold', fontsize=14)
    plt.subplots_adjust(wspace=20, hspace=0.8)
    plt.savefig('Paper fig R2 Jan 18.png', bbox_inches='tight', dpi=600)
    return
    
#l, dr, dg, s6, drs, dgs, ls, le, ss, se = get_all_data()
#new_means, new_errors = get_means(l, dr, dg, s6, drs, dgs)
#high_means, high_errors, high_all = highest_means(l, dr, dg, s6, drs, dgs)
#plot_standards(ls, le, ss, se)
#LS, g = plot_all_means(new_means, new_errors, 'means')
#plot_all_means(high_means, high_errors, 'high')
#long_means, long_errors = [new_means[0], new_means[1], new_means[2], new_means[3]], [new_errors[0], new_errors[1], new_errors[2], new_errors[3]]
#short_means, short_errors = [new_means[4], new_means[5]], [new_errors[4], new_errors[5]]
#plot_short(short_means, short_errors)
#get_paper_plot(LS, g)