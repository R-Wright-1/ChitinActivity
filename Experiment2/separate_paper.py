import numpy
import functions_paper as fnc
import heapq
import os


folder = '/Users/u1560915/Documents/GitHub/ChitinActivity/Experiment2/'
plates = ['Plate_1/', 'Plate_2/', 'Plate_3/']
results_files = ['Results_a.csv', 'Results_b.csv', 'Results_c.csv']
results_norm_files = ['Results_norm_a.csv', 'Results_norm_b.csv', 'Results_norm_c.csv']
communities = []

for a in range(3):
    all_files = []
    for item in os.listdir(folder+plates[a]):
        if item.endswith('.csv'):
            all_files.append(item)
    all_files = sorted(all_files)
    sorted_files = fnc.sorted_files(all_files)
    for a in range(len(sorted_files)):
        for b in range(len(sorted_files[a])):
            sorted_files[a][b] = sorted_files[a][b]
    communities.append(sorted_files)

all_means, all_errors, all_days = [[[], []], [[], []], [[], []]], [[[], []], [[], []], [[], []]], [[], [], []]
all_totals = [[], [], []]
highest_means, highest_errors = [[[], []], [[], []], [[], []]], [[[], []], [[], []], [[], []]]
differences, differences_errors = [[], [], []], [[], [], []]
for b in range(3):
    uf = folder+plates[b]
    os.chdir(uf)
    total_day = 0
    for gen in range(len(communities[b])):
        days_s, days_std, days_c, days_c_std, days_len = [], [], [], [], []
        
        this_means, this_errors, this_days = [[], []], [[], []], []
        this_nc, this_ns = [], []
        for day in range(len(communities[b][gen])):
            total_day += 1
            fn = communities[b][gen][day]
            if fn[13] != '.':
                td = int(fn[12]+fn[13])
            else:
                td = int(fn[12])
            if fn[3] != '_':
                this_gen = int(fn[2]+fn[3])
            else:
                this_gen = int(fn[2])
            sm, se, gradient, intercept, f = fnc.get_standards(td)
            control, selection = fnc.get_samples(fn)[2], fnc.get_samples(fn)[3]
            nc, ns = fnc.normalise_community(control, gradient, intercept, f), fnc.normalise_community(selection, gradient, intercept, f)
            this_nc.append(nc)
            this_ns.append(ns)
            i = [[3, 2, 3, 4, 7, 3, 8, 8, 8, 4], [3, 2, 8, 6, 8, 8, 7, 8], [3, 3, 3, 2, 3, 3, 5, 4, 2, 3, 2, 2, 6, 2, 3, 2, 2]]
            fnc.plot_this_plate(nc, ns, b, this_gen, day+1)
            days_s.append(numpy.mean(ns))
            days_std.append(numpy.std(ns))
            days_c.append(numpy.mean(nc))
            days_c_std.append(numpy.std(nc))
            days_len.append(day+1)
        fnc.plot_days(days_s, days_std, days_c, days_c_std, days_len, b, this_gen, i[b][gen])
        highest = heapq.nlargest(1, zip(days_s, days_len))
        highest_day_index = highest[0][1]
        add1 = 0
        i = [[3, 2, 3, 4, 7, 3, 8, 8, 8, 4], [3, 2, 8, 6, 8, 8, 7, 8], [3, 3, 3, 2, 3, 3, 5, 4, 2, 3, 2, 2, 6, 2, 3, 2, 2]]
        highest_day_index = i[b][gen]
        plate2gen2 = True
        if plate2gen2 == True:
            nc, ns = this_nc[highest_day_index-1], this_ns[highest_day_index-1]
            c_index, c_numbers, s_index, s_numbers, td = fnc.next_generation(b, this_gen, highest_day_index, nc, ns, td)
            all_means[b][0].append(numpy.mean(nc))
            all_means[b][1].append(numpy.mean(ns))
            all_errors[b][0].append(numpy.std(nc))
            all_errors[b][1].append(numpy.std(ns))
            all_days[b].append(highest_day_index+add1)
            all_totals[b].append(td)
            highest_means[b][0].append(numpy.mean(c_numbers))
            highest_means[b][1].append(numpy.mean(s_numbers))
            highest_errors[b][0].append(numpy.std(c_numbers))
            highest_errors[b][1].append(numpy.std(s_numbers))
fnc.plot_all(all_means, all_errors, all_days, highest_means, highest_errors)
fnc.plot_normalised_only(all_means, all_errors, all_days, highest_means, highest_errors)
fnc.plot_normalised_only_broken_axes(all_means, all_errors, all_days, highest_means, highest_errors)
fnc.plot_normalised_3_only(all_means, all_errors, all_days, highest_means, highest_errors)
