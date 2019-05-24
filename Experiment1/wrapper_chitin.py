import functions_paper as fnc

fnc.get_daily()
l, dr, dg, s6, drs, dgs, ls, le, ss, se = fnc.get_all_data()
new_means, new_errors = fnc.get_means(l, dr, dg, s6, drs, dgs)
high_means, high_errors, high_all = fnc.highest_means(l, dr, dg, s6, drs, dgs)
fnc.plot_standards(ls, le, ss, se)
LS, g = fnc.plot_all_means(new_means, new_errors, 'means')
fnc.plot_all_means(high_means, high_errors, 'high')
#print('means = ', new_means)
#print('errors = ', new_errors)
#print('high_means = ', high_means)
#print('high_errors = ', high_errors)
long_means, long_errors = [new_means[0], new_means[1], new_means[2], new_means[3]], [new_errors[0], new_errors[1], new_errors[2], new_errors[3]]
short_means, short_errors = [new_means[4], new_means[5]], [new_errors[4], new_errors[5]]
fnc.plot_short(short_means, short_errors)
g20_g1 = fnc.get_paper_plot(LS, g)

