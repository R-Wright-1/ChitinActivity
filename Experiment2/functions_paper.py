import numpy
import scipy.stats as stats
import matplotlib.pyplot as plt
import heapq
from random import *
import csv
import os
from sklearn import linear_model
from sklearn.metrics import mean_squared_error, r2_score
import pandas as pd
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
import statsmodels.api as sm

bdir='/Users/u1560915/Documents/GitHub/ChitinActivity/Experiment2/'
figures = 'New_figures/'


samples = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29]
def sorted_files(file_list):
    gens, all_gens_sorted = [], []
    for a in range(len(file_list)):
        if file_list[a][3] != '_':
            gen = int(file_list[a][2]+file_list[a][3])
        else:
            gen = int(file_list[a][2])
        adding = True
        for b in range(len(gens)):
            if gen == gens[b]:
                adding = False
        if adding == True:
            gens.append(gen)
    for c in range(len(gens)):
        this_gen = gens[c]
        this_gen_files = []
        for d in range(len(file_list)):
            if file_list[d][3] != '_':
                tg = float(file_list[d][2]+file_list[d][3])
            else:
                tg = float(file_list[d][2])
            if tg == this_gen:
                this_gen_files.append(file_list[d])
        all_gens_sorted.append(this_gen_files)
    return all_gens_sorted
    
def get_samples(fn):
    r, c = 27, 2
    rows, plate = [], []
    with open(fn, 'rU') as f:
        reader = csv.reader(f)
        for row in reader:
            rows.append(row)
        for i in range(8):
            for j in range(12):
                plate.append(float(rows[r+i][c+j]))
    standards1, standards2, selec, cont = [], [], [], []
    selection, control = [], []
    for a in range(6):
        b = a+6
        standards1.append(plate[a])
        standards2.append(plate[b])
    selec.append(plate[13:23]), selec.append(plate[25:35]), selec.append(plate[37:47])
    cont.append(plate[49:59]), cont.append(plate[61:71]), cont.append(plate[73:83])
    for a in range(3):
        for b in range(10):
            control.append(selec[a][b])
            selection.append(cont[a][b])
    return standards1, standards2, control, selection
    
def get_standards(day):
    concs = [0, 0.00016, 0.0008, 0.004, 0.02]
    for a in range(len(concs)):
        concs[a] = (concs[a])
    folder = bdir+'All_plates/'
    os.chdir(folder)
    this_day, std, standards, errors = [], [], [[], [], [], [], []], []
    for item in os.listdir(folder):
        if item.endswith('.csv') and item != 'Standards.csv':
            if item[13] != '.':
                file_day = int(item[12]+item[13])
            else:
                file_day = int(item[12])
            if file_day == day:
                this_day.append(item)
    for c in range(len(this_day)):
        s1, s2, control, selection = get_samples(this_day[c])
        std.append(s1), std.append(s2)
    for d in range(len(std)):
        for e in range(len(standards)):
            standards[e].append(std[d][e])
    for f in range(len(standards)):
        this_standard = standards[f]
        standards[f] = numpy.mean(this_standard)
        errors.append(numpy.std(this_standard))
    gradient, intercept, r_value, p_value, std_err = stats.linregress(concs, standards)
    z = numpy.polyfit(standards, concs, 4)
    f = numpy.poly1d(z)
    x_new = f(standards)
    ax1 = plt.subplot(111)
    ax1.plot(concs, standards, 'o-')
    ax1.plot(x_new, standards, 'r')
    plt.savefig('Standards '+str(day)+'.pdf', bbox_inches='tight')
    plt.close()
    return standards, errors, gradient, intercept, f

def normalise_community(comm, gradient, intercept, f):
    new_comm = []
    for a in range(len(comm)):
        num = comm[a]
        new_num = (num-intercept)/gradient
        new_comm.append(new_num)
    new_comm = f(comm)
    for a in range(len(new_comm)):
        new_comm[a] = new_comm[a]*1440
    return new_comm
    
def plot_this_plate(nc, ns, plate, gen, day):
    ax1 = plt.subplot(111)
    savestring = 'Plate '+str(plate+1)+' CC'+str(gen)+' day '+str(day)
    ax1.plot(samples, nc, color='#0099FF', marker='o', label='Control')
    ax1.plot(samples, ns, color='#FF0099', marker='o', label='Selection')
    ax1.set_xlabel('Sample')
    ax1.set_ylabel(r'Chitinase $\mu$M day$^{-1}$')
    ax1.set_title(savestring)
    ax1.legend(bbox_to_anchor=(1.4, 1))
    plt.tight_layout()
    os.chdir(bdir+'New_figures/')
    plt.savefig(savestring+'.pdf', bbox_inches='tight')
    plt.close()
    return

def plot_days(days_mean, days_std, days_control_mean, days_control_std, days_len, plate, gen, high):
    savestring = 'Plate '+str(plate+1)+' CC'+str(gen)
    ax1 = plt.subplot(111)
    if high != 8:
        high -= 1
    ax1.plot(high, days_mean[high-1], marker='*', markersize=15, color='k')
    ax1.plot(high, days_control_mean[high-1], marker='*', markersize=15, color='k')
    ax1.errorbar(days_len, days_control_mean, yerr=days_control_std, color='#0099FF', marker='o', label='Control')
    ax1.errorbar(days_len, days_mean, yerr=days_std, color='#FF0099', marker='o', label='Selection')
    ax1.set_xlabel('Day')
    ax1.set_ylabel(r'Chitinase $\mu$M day$^{-1}$')
    xmin, xmax = days_len[0]-1, days_len[-1]+1
    ax1.set_xlim([xmin, xmax])
    ax1.set_title(savestring)
    ax1.legend(bbox_to_anchor=(1.4, 1))
    plt.xticks(days_len)
    plt.tight_layout()
    os.chdir(bdir+figures)
    plt.savefig(savestring+'.pdf', bbox_inches='tight')
    plt.close()
    return
    
def get_highest(comm, samples=samples):
    highest = heapq.nlargest(3, zip(comm, samples))
    selection_numbers = [highest[0][1], highest[1][1], highest[2][1]]
    selection_highest = [highest[0][0], highest[1][0], highest[2][0]]
    index, numbers = selection_numbers, selection_highest
    return index, numbers

def get_lowest(comm, samples=samples):
    highest = heapq.nsmallest(3, zip(comm, samples))
    selection_numbers = [highest[0][1], highest[1][1], highest[2][1]]
    selection_highest = [highest[0][0], highest[1][0], highest[2][0]]
    index, numbers = selection_numbers, selection_highest
    return index, numbers

def get_random(comm, samples=samples):
    random_index = [randint(0, 29), randint(0, 29), randint(0, 29)]
    random = [comm[random_index[0]], comm[random_index[1]], comm[random_index[2]]]
    index, numbers = random_index, random
    return index, numbers
   
def next_generation(plate, gen, day, control, selection, td):
    os.chdir(bdir)
    plate += 1
    rows, new_rows = [], []
    c_index, s_index = [], []
    c_numbers, s_numbers = [], []
    with open('Next_generation_separate.csv', 'rU') as f:
        reader = csv.reader(f)
        for row in reader:
            rows.append(row)
    for a in range(len(rows)):
        if rows[a] != [] and rows[a][0] != '':
            new_rows.append(rows[a])
    rows = new_rows
    del rows[0]
    for b in range(len(new_rows)):
        if new_rows[b] != []:
            if int(new_rows[b][0]) == plate and int(new_rows[b][2]) == gen:
                td = int(new_rows[b][4])
                if new_rows[b][1] == 'Control':
                    c_index.append(int(new_rows[b][5]))
                    c_index.append(int(new_rows[b][6]))
                    c_index.append(int(new_rows[b][7]))
                elif new_rows[b][1] == 'Selection':
                    s_index.append(int(new_rows[b][5]))
                    s_index.append(int(new_rows[b][6]))
                    s_index.append(int(new_rows[b][7]))
    if len(c_index) != 0 and len(s_index) != 0:
        for c in range(len(c_index)):
            c_numbers.append(control[c_index[c]])
            s_numbers.append(selection[s_index[c]])
        return c_index, c_numbers, s_index, s_numbers, td
    if plate == 4:
        c_index, c_numbers = get_lowest(control)
    else:
        c_index, c_numbers = get_random(control)
    s_index, s_numbers = get_highest(selection)
    with open('Next_generation_separate.csv', 'a') as f:
        writer = csv.writer(f)
        writer.writerow('')
        writing = [plate, 'Control', gen, day, td, c_index[0], c_index[1], c_index[2]]
        writer.writerow(writing)
        writing = [plate, 'Selection', gen, day, td, s_index[0], s_index[1], s_index[2]]
        writer.writerow(writing)
    return c_index, c_numbers, s_index, s_numbers, td

def plot_all(all_means, all_errors, total_days, highest_means, highest_errors):
    os.chdir(bdir+figures)
    ds = [[2, 5, 8, 11, 15, 19, 26, 34, 42, 50], [2, 5, 11, 18, 26, 34, 42, 50], [2, 4, 6, 8, 10, 12, 14, 17, 20, 23, 26, 29, 32, 35, 39, 44, 50]]
    xmin, xmax = 0, 51
    colors = ['#6633FF', '#990033', '#006699']
    fig = plt.figure(figsize=(8.27, 5)) 
    ax1 = plt.subplot(221, label='A')
    ax2, ax3, ax4 = plt.subplot(222, sharex=ax1, label='B'), plt.subplot(223, sharex=ax1, label='C'), plt.subplot(224, sharex=ax1, label='D')
    fig.text(0.5, 0.00, 'Days', ha='center')
    fig.text(-0.01, 0.5, r'Chitinase $\mu$M day$^{-1}$', va='center', rotation='vertical')
    ax1.set_xlim([xmin, xmax])
    markers = ['^', 'o', '*']
    axes = [ax1, ax2, ax3]
    for d in range(len(all_means)):
        new_c = heapq.nsmallest(len(total_days[d]), zip(total_days[d], all_means[d][0]))
        new_ce = heapq.nsmallest(len(total_days[d]), zip(total_days[d], all_errors[d][0]))
        new_s = heapq.nsmallest(len(total_days[d]), zip(total_days[d], all_means[d][1]))
        new_se = heapq.nsmallest(len(total_days[d]), zip(total_days[d], all_errors[d][1]))
        new_days, new_norm = [], []
        nnew_s, nnew_c, nnew_se, nnew_ce = [], [], [], []
        days = 0
        for a in range(len(new_c)):
            days+= new_c[a][0]
            new_days.append(days)
            new_num = (new_s[a][1]-new_c[a][1])
            nnew_s.append(new_s[a][1])
            nnew_se.append(new_se[a][1])
            nnew_c.append(new_c[a][1])
            nnew_ce.append(new_ce[a][1])
            new_norm.append(new_num)
        axes[d].errorbar(ds[d], nnew_c, yerr=nnew_ce, marker=markers[d], linestyle='--', markersize=5, color=colors[d])
        axes[d].errorbar(ds[d], nnew_s, yerr=nnew_se, marker=markers[d], linestyle='-', markersize=5, color=colors[d])
        label = 'Replicate '+str(d+1)
        regr = make_pipeline(PolynomialFeatures(2), LinearRegression())
        with open('Regression '+str(d)+'.csv', 'w') as f:
            writer = csv.writer(f)
            writer.writerow(['X', 'Y'])
            for a in range(len(ds[d])):
                writer.writerow([ds[d][a], new_norm[a]])
        df = pd.read_csv('Regression '+str(d)+'.csv')
        regr.fit(df[['X']], df[['Y']])
        y_pred = regr.predict(df[['X']])
        ax4.plot(ds[d], y_pred, '-.', color=colors[d])
        ax4.errorbar(ds[d], new_norm, marker=markers[d], linestyle='-', markersize=5, color=colors[d], label=label)
    ax1.set_title('A', loc='left'), ax2.set_title('B', loc='left'), ax3.set_title('C', loc='left'), ax4.set_title('D', loc='left')
    ax1.plot([xmin, xmax], [0, 0], 'k--')
    ax2.plot([xmin, xmax], [0, 0], 'k--')
    ax3.plot([xmin, xmax], [0, 0], 'k--')
    ax4.plot([xmin, xmax], [0, 0], 'k--')
    ax4.legend(bbox_to_anchor=(0.4, 1), frameon=False, numpoints=1, handlelength=0)
    ax1.tick_params(axis='both', which='both', length=0)
    ax2.tick_params(axis='both', which='both', length=0)
    ax3.tick_params(axis='both', which='both', length=0)
    ax4.tick_params(axis='both', which='both', length=0)
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.xticks([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50])
    plt.tight_layout()
    plt.savefig('All presentation.pdf', bbox_inches='tight')
    plt.close()
    
    xmin, xmax = 0, 51
    fig = plt.figure(figsize=(8.27, 5)) 
    colors = ['#6633FF', '#990033', '#006699']
    ax1 = plt.subplot(221, label='A')
    ax2, ax3, ax4 = plt.subplot(222, sharex=ax1, label='B'), plt.subplot(223, sharex=ax1, label='C'), plt.subplot(224, sharex=ax1, label='D')
    fig.text(0.5, 0.00, 'Days', ha='center')
    fig.text(-0.01, 0.5, r'Chitinase $\mu$M day$^{-1}$', va='center', rotation='vertical')
    ax1.set_xlim([xmin, xmax])
    markers = ['^', 'o', '*']
    axes = [ax1, ax2, ax3]
    returning_days = []
    for d in range(len(all_means)):
        new_c = heapq.nsmallest(len(total_days[d]), zip(total_days[d], highest_means[d][0]))
        new_ce = heapq.nsmallest(len(total_days[d]), zip(total_days[d], highest_errors[d][0]))
        new_s = heapq.nsmallest(len(total_days[d]), zip(total_days[d], highest_means[d][1]))
        new_se = heapq.nsmallest(len(total_days[d]), zip(total_days[d], highest_errors[d][1]))
        nnew_s, nnew_c, nnew_se, nnew_ce = [], [], [], []
        new_days, new_norm = [], []
        days = 0
        for a in range(len(new_c)):
            days+= new_c[a][0]
            new_days.append(days)
            new_num = (new_s[a][1]-new_c[a][1])
            nnew_s.append(new_s[a][1])
            nnew_se.append(new_se[a][1])
            nnew_c.append(new_c[a][1])
            nnew_ce.append(new_ce[a][1])
            new_norm.append(new_num)
        returning_days.append(new_days)
        axes[d].errorbar(ds[d], nnew_c, yerr=nnew_ce, marker=markers[d], linestyle='--', markersize=5, color=colors[d])
        axes[d].errorbar(ds[d], nnew_s, yerr=nnew_se, marker=markers[d], linestyle='-', markersize=5, color=colors[d])
        regr = make_pipeline(PolynomialFeatures(2), LinearRegression())
        with open('Regression highest '+str(d)+'.csv', 'w') as f:
            writer = csv.writer(f)
            writer.writerow(['X', 'Y'])
            for a in range(len(ds[d])):
                writer.writerow([ds[d][a], new_norm[a]])
        df = pd.read_csv('Regression highest '+str(d)+'.csv')
        regr.fit(df[['X']], df[['Y']])
        y_pred = regr.predict(df[['X']])
        ax4.plot(ds[d], y_pred, '-.', color=colors[d])
        ax4.errorbar(ds[d], new_norm, marker=markers[d], linestyle='-', markersize=5, color=colors[d], label='Replicate '+str(d+1))
    ax1.set_title('A', loc='left'), ax2.set_title('B', loc='left'), ax3.set_title('C', loc='left'), ax4.set_title('D', loc='left')
    ax1.plot([xmin, xmax], [0, 0], 'k--')
    ax2.plot([xmin, xmax], [0, 0], 'k--')
    ax3.plot([xmin, xmax], [0, 0], 'k--')
    ax4.plot([xmin, xmax], [0, 0], 'k--')
    ax4.legend(bbox_to_anchor=(0.4, 1), frameon=False, numpoints=1, handlelength=0)
    ax1.tick_params(axis='both', which='both', length=0)
    ax2.tick_params(axis='both', which='both', length=0)
    ax3.tick_params(axis='both', which='both', length=0)
    ax4.tick_params(axis='both', which='both', length=0)
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.xticks([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50])
    plt.tight_layout()
    plt.savefig('All highest presentation.pdf', bbox_inches='tight')
    plt.close()
    return
    
def plot_normalised_only(all_means, all_errors, total_days, highest_means, highest_errors):
    os.chdir(bdir+figures)
    ds = [[2, 5, 8, 11, 15, 19, 26, 34, 42, 50], [2, 5, 11, 18, 26, 34, 42, 50], [2, 4, 6, 8, 10, 12, 14, 17, 20, 23, 26, 29, 32, 35, 39, 44, 50]]
    xmin, xmax = 0, 51
    colors = ['#6633FF', '#990033', '#006699']
    markers = ['^', 'o', '*']
    ax1 = plt.subplot(111)
    regr = make_pipeline(PolynomialFeatures(2), LinearRegression())
    all_norm = []
    for d in range(len(all_means)):
        new_c = heapq.nsmallest(len(total_days[d]), zip(total_days[d], all_means[d][0]))
        new_s = heapq.nsmallest(len(total_days[d]), zip(total_days[d], all_means[d][1]))
        new_days, new_norm = [], []
        days = 0
        for a in range(len(new_c)):
            days+= new_c[a][0]
            new_days.append(days)
            new_num = (new_s[a][1]-new_c[a][1])
            new_norm.append(new_num)
        all_norm.append(new_norm)
        df = pd.read_csv('Regression '+str(d)+'.csv')
        regr.fit(df[['X']], df[['Y']])
        y_pred = regr.predict(df[['X']])
        est = sm.OLS(df[['Y']], y_pred)
        pval =  est.fit().f_pvalue
        r2 = est.fit().rsquared
        ax1.plot(ds[d], y_pred, '-.', color=colors[d])
        label = 'Replicate '+str(d+1)
        label2 = r'$r^2$'+'=%.2f'%r2
        label2 += ', '+r'$p$'+'=%.2f'%pval
        ax1.errorbar(ds[d], new_norm, marker=markers[d], linestyle='-', markersize=5, color=colors[d], label=label, alpha=0.7)
        ax1.errorbar(ds[d], new_norm, linestyle='-', color=colors[d], label=label2, alpha=0.7)
    ax1.plot([xmin, xmax], [0, 0], 'k--')
    legend = ax1.legend(bbox_to_anchor=(0.28, 1), frameon=False, numpoints=1, handlelength=0)
    t1, t2, t3 = legend.get_texts()[1], legend.get_texts()[3], legend.get_texts()[5]
    for t in [t1, t2, t3]:
        t.set_fontsize(8)
    plt.xticks([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50])
    ax1.set_xlim([xmin, xmax])
    ax1.set_xlabel('Days')
    ax1.set_ylabel(r'Chitinase $\mu$M day$^{-1}$')
    ax1.set_title('A', loc='left')
    plt.tight_layout()
    plt.savefig('Presentation normalised.pdf', bbox_inches='tight')
    plt.close()
    
    ax1 = plt.subplot(111)
    for d in range(len(all_means)):
        new_c = heapq.nsmallest(len(total_days[d]), zip(total_days[d], highest_means[d][0]))
        new_s = heapq.nsmallest(len(total_days[d]), zip(total_days[d], highest_means[d][1]))
        new_days, new_norm = [], []
        days = 0
        for a in range(len(new_c)):
            days+= new_c[a][0]
            new_days.append(days)
            new_num = (new_s[a][1]-new_c[a][1])
            new_norm.append(new_num)
        df = pd.read_csv('Regression highest '+str(d)+'.csv')
        regr.fit(df[['X']], df[['Y']])
        y_pred = regr.predict(df[['X']])
        est = sm.OLS(df[['Y']], y_pred)
        pval =  est.fit().f_pvalue
        r2 = est.fit().rsquared
        ax1.plot(ds[d], y_pred, '-.', color=colors[d])
        label = 'Replicate '+str(d+1)
        label2 = r'$r^2$'+'=%.2f'%r2
        label2 += ', '+r'$p$'+'=%.2f'%pval
        ax1.errorbar(ds[d], new_norm, marker=markers[d], linestyle='-', markersize=5, color=colors[d], label=label, alpha=0.7)
        ax1.errorbar(ds[d], new_norm, linestyle='-', color=colors[d], label=label2, alpha=0.7)
        ax1.plot(ds[d], y_pred, '-.', color=colors[d])
    legend = ax1.legend(bbox_to_anchor=(0.28, 1), frameon=False, numpoints=1, handlelength=0)
    t1, t2, t3 = legend.get_texts()[1], legend.get_texts()[3], legend.get_texts()[5]
    for t in [t1, t2, t3]:
        t.set_fontsize(8)
    ax1.set_xlim([xmin, xmax])
    plt.xticks([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50])
    ax1.set_xlabel('Days')
    ax1.set_ylabel(r'Chitinase $\mu$M day$^{-1}$')
    ax1.set_title('B', loc='left')
    plt.tight_layout()
    plt.savefig('Presentation normalised highest.pdf', bbox_inches='tight')
    plt.close()
    return
    
def plot_normalised_only_broken_axes(all_means, all_errors, total_days, highest_means, highest_errors):
    os.chdir(bdir+figures)
    ds = [[2, 5, 8, 11, 15, 19, 26, 34, 42, 50], [2, 5, 11, 18, 26, 34, 42, 50], [2, 4, 6, 8, 10, 12, 14, 17, 20, 23, 26, 29, 32, 35, 39, 44, 50]]
    xmin, xmax = 0, 51
    colors = ['#6633FF', '#990033', '#006699']
    markers = ['^', 'o', '*']
    f, (ax2, ax1) = plt.subplots(2, 1, sharex=True)
    ax1 = plt.subplot2grid((2, 1), (0, 0), rowspan=1)
    ax2 = plt.subplot2grid((2, 1), (1, 0), rowspan=1)
    ax1 = plt.subplot(211)
    ax2 = plt.subplot(212, sharex=ax1)
    ax2.set_ylim(-1, 7)
    ax1.set_ylim(7, 30)
    ax1.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax1.xaxis.tick_top()
    ax1.tick_params(labeltop='off')
    ax2.xaxis.tick_bottom()
    regr = make_pipeline(PolynomialFeatures(2), LinearRegression())
    all_norm = []
    for d in range(len(all_means)):
        new_c = heapq.nsmallest(len(total_days[d]), zip(total_days[d], all_means[d][0]))
        new_s = heapq.nsmallest(len(total_days[d]), zip(total_days[d], all_means[d][1]))
        new_days, new_norm = [], []
        days = 0
        for a in range(len(new_c)):
            days+= new_c[a][0]
            new_days.append(days)
            new_num = (new_s[a][1]-new_c[a][1])
            new_norm.append(new_num)
        all_norm.append(new_norm)
        df = pd.read_csv('Regression '+str(d)+'.csv')
        regr.fit(df[['X']], df[['Y']])
        y_pred = regr.predict(df[['X']])
        est = sm.OLS(df[['Y']], y_pred)
        pval =  est.fit().f_pvalue
        r2 = est.fit().rsquared
        ax1.plot(ds[d], y_pred, '-.', color=colors[d])
        ax2.plot(ds[d], y_pred, '-.', color=colors[d])
        label = 'Replicate '+str(d+1)
        label2 = r'$r^2$'+'=%.2f'%r2
        label2 += ', '+r'$p$'+'=%.2f'%pval
        ax1.errorbar(ds[d], new_norm, marker=markers[d], linestyle='-', markersize=5, color=colors[d], label=label, alpha=0.7)
        ax1.errorbar(ds[d], new_norm, linestyle='-', color=colors[d], label=label2, alpha=0.7)
        ax2.errorbar(ds[d], new_norm, marker=markers[d], linestyle='-', markersize=5, color=colors[d], alpha=0.7)
        ax2.errorbar(ds[d], new_norm, linestyle='-', color=colors[d], alpha=0.7)
    ax2.plot([xmin, xmax], [0, 0], 'k--')
    legend = ax1.legend(bbox_to_anchor=(0.26, 1), frameon=False, numpoints=1, handlelength=0)
    t1, t2, t3 = legend.get_texts()[1], legend.get_texts()[3], legend.get_texts()[5]
    for t in [t1, t2, t3]:
        t.set_fontsize(8)
    d = 0.015
    kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
    ax1.plot((-d, +d), (-d, +d), **kwargs)
    ax1.plot((1-d, 1+d), (-d, +d), **kwargs)
    kwargs.update(transform=ax2.transAxes)
    ax2.plot((-d, +d), (1-d, 1+d), **kwargs)
    ax2.plot((1-d, 1+d), (1-d, 1+d), **kwargs)
    plt.xticks([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50])
    ax1.set_xlim([xmin, xmax])
    ax2.set_xlabel('Days')
    f.text(-0.01, 0.5, r'Chitinase $\mu$M day$^{-1}$', va='center', rotation='vertical')
    ax1.set_title('A', loc='left')
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.05, hspace=None)
    plt.tight_layout()
    plt.savefig('Presentation normalised broken axes.pdf', bbox_inches='tight')
    plt.close()
    
    f, (ax2, ax1) = plt.subplots(2, 1, sharex=True)
    ax1 = plt.subplot2grid((2, 1), (0, 0), rowspan=1)
    ax2 = plt.subplot2grid((2, 1), (1, 0), rowspan=1)
    ax1 = plt.subplot(211)
    ax2 = plt.subplot(212, sharex=ax1)
    ax2.set_ylim(-1, 15)
    ax1.set_ylim(15, 90)
    ax1.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax1.xaxis.tick_top()
    ax1.tick_params(labeltop='off')
    ax2.xaxis.tick_bottom()
    for d in range(len(all_means)):
        new_c = heapq.nsmallest(len(total_days[d]), zip(total_days[d], highest_means[d][0]))
        new_s = heapq.nsmallest(len(total_days[d]), zip(total_days[d], highest_means[d][1]))
        new_days, new_norm = [], []
        days = 0
        for a in range(len(new_c)):
            days+= new_c[a][0]
            new_days.append(days)
            new_num = (new_s[a][1]-new_c[a][1])
            new_norm.append(new_num)
        df = pd.read_csv('Regression highest '+str(d)+'.csv')
        regr.fit(df[['X']], df[['Y']])
        y_pred = regr.predict(df[['X']])
        est = sm.OLS(df[['Y']], y_pred)
        pval =  est.fit().f_pvalue
        r2 = est.fit().rsquared
        ax1.plot(ds[d], y_pred, '-.', color=colors[d])
        ax2.plot(ds[d], y_pred, '-.', color=colors[d])
        label = 'Replicate '+str(d+1)
        label2 = r'$r^2$'+'=%.2f'%r2
        label2 += ', '+r'$p$'+'=%.2f'%pval
        ax1.errorbar(ds[d], new_norm, marker=markers[d], linestyle='-', markersize=5, color=colors[d], label=label, alpha=0.7)
        ax1.errorbar(ds[d], new_norm, linestyle='-', color=colors[d], label=label2, alpha=0.7)
        ax1.plot(ds[d], y_pred, '-.', color=colors[d])
        ax2.errorbar(ds[d], new_norm, marker=markers[d], linestyle='-', markersize=5, color=colors[d], label=label, alpha=0.7)
        ax2.errorbar(ds[d], new_norm, linestyle='-', color=colors[d], label=label2, alpha=0.7)
        ax2.plot(ds[d], y_pred, '-.', color=colors[d])
    legend = ax1.legend(bbox_to_anchor=(0.26, 1), frameon=False, numpoints=1, handlelength=0)
    t1, t2, t3 = legend.get_texts()[1], legend.get_texts()[3], legend.get_texts()[5]
    for t in [t1, t2, t3]:
        t.set_fontsize(8)
    d = 0.015
    kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
    ax1.plot((-d, +d), (-d, +d), **kwargs)
    ax1.plot((1-d, 1+d), (-d, +d), **kwargs)
    kwargs.update(transform=ax2.transAxes)
    ax2.plot((-d, +d), (1-d, 1+d), **kwargs)
    ax2.plot((1-d, 1+d), (1-d, 1+d), **kwargs)
    plt.xticks([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50])
    ax1.set_xlim([xmin, xmax])
    ax2.set_xlabel('Days')
    f.text(-0.01, 0.5, r'Chitinase $\mu$M day$^{-1}$', va='center', rotation='vertical')
    plt.xticks([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50])
    ax1.set_title('B', loc='left')
    plt.tight_layout()
    plt.savefig('Presentation normalised highest broken axes.pdf', bbox_inches='tight')
    plt.close()
    return