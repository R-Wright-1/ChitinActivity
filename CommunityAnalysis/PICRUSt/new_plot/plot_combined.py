import matplotlib.pyplot as plt
import csv
import numpy

plt.tight_layout()

def get_data(f):
    order_daily = [['S20_1_4', 'S20_1_5', 'S20_1_6'], ['S20_2_4', 'S20_2_5', 'S20_2_6'], 
             ['S20_3_4', 'S20_3_5', 'S20_3_6'], ['S20_4_4', 'S20_4_5', 'S20_4_6']]
    order_daily_R = ['S20_1_A', 'S20_2_A', 'S20_3_A', 'S20_4_A']
    order_short = [['S16_4', 'S16_5', 'S16_6'], ['S17_4', 'S17_5', 'S17_6'], ['S18_4', 'S18_5', 'S18_6'],
                   ['S19_4', 'S19_5', 'S19_6'], ['S20_4_4', 'S20_4_5', 'S20_4_6']]
    order_short_R = ['S16_A', 'S17_A', 'S18_A', 'S19_A', 'S20_4_A']
    order = [['L0_7', 'L0_8', 'L0_9'], ['L1_7', 'L1_8', 'L1_9'], [], ['L3_7'], 
             ['L4_7', 'L4_8', 'L4_9'], [], ['L6_7'], [], 
             ['L8_9'], ['L9_7', 'L9_8', 'L9_9'], 
             ['L10_7', 'L10_8', 'L10_9'], ['L11_7', 'L11_8', 'L11_9'], ['L12_7', 'L12_9'], ['L13_7', 'L13_8', 'L13_9'],
             ['L14_7', 'L14_8', 'L14_9'], [], ['L16_9'], ['L17_8', 'L17_9'],
             ['L18_9'], ['L19_8'], ['L20_7', 'L20_8', 'L20_9']]
    order_R = ['L0_B', 'L1_B', 'L2_B', 'L3_B', '', 'L5_B', 'L6_B', '', 'L8_B', 'L9_B', 'L10_B',
               'L11_B', 'L12_B', 'L13_B', 'L14_B', '', 'L16_B', '', 'L18_B', 'L19_B', 'L20_B']
    data = [[], [], []]
    data_R = [[], [], []]
    data_std = [[], [], []]
    orders = [order, order_short, order_daily]
    orders_R = [order_R, order_short_R, order_daily_R]
    with open(f, 'rU') as f:
        rows = []
        for row in csv.reader(f):
            rows.append(row)
    for a in range(len(orders)):
        for b in range(len(rows)):
            if b > 0:
                taxa, taxa_std = [], []
                taxa_R = []
                for c in range(len(orders[a])):
                    this_gen = []
                    for d in range(len(orders[a][c])):
                        for e in range(len(rows[0])):
                            if rows[0][e] == orders[a][c][d]:
                                this_gen.append(float(rows[b][e])/100)
                    if this_gen != []:
                        if orders == 2:
                            print(this_gen)
                        taxa.append(numpy.mean(this_gen))
                        taxa_std.append(numpy.std(this_gen))
                    else:
                        taxa.append(0)
                        taxa_std.append(0)
                data[a].append(taxa)
                data_std.append(taxa_std)
                for f in range(len(orders_R[a])):
                    for g in range(len(rows[0])):
                        if orders_R[a][f] == rows[0][g]:
                            taxa_R.append(float(rows[b][g])/100)
                data_R.append(taxa_R)
        
    labels = []
    for h in range(len(rows)):
        if h > 0:
            labels.append(rows[h][0])
    gens, gens_text, gens_R, gens_text_R = [], [], [], []
    for a in range(len(order)):
        gens.append(a)
        if order[a] != []:
            gens_text.append(a)
        else:
            gens_text.append('*')
    for b in range(len(order_R)):
        gens_R.append(b)
        if order_R[b] != '':
            gens_text_R.append(b)
        else:
            gens_text_R.append('*')
    return labels, data, data_std, data_R, gens, gens_text, gens_R, gens_text_R

def get_distinct_colors(n):
    colors = []
    for i in numpy.arange(0., 360., 360. / n):
        h = i / 360.
        l = (50 + numpy.random.rand() * 10) / 100.
        s = (90 + numpy.random.rand() * 10) / 100.
        colors.append(hls_to_rgb(h, l, s))
    random.shuffle(colors)
    return colors

def make_plot(files):
    fig = plt.figure(figsize=(11, 10))
    ax1 = plt.subplot2grid((4,30), (0,0), colspan=21)
    ax2, ax3 = plt.subplot2grid((4,30), (0,21), colspan=5, sharey=ax1), plt.subplot2grid((4,30), (0,26), colspan=4, sharey=ax1)
    ax4 = plt.subplot2grid((4,30), (1,0), colspan=21, sharex=ax1)
    ax5, ax6 = plt.subplot2grid((4,30), (1,21), colspan=5, sharey=ax4, sharex=ax2), plt.subplot2grid((4,30), (1,26), colspan=4, sharey=ax4, sharex=ax3)
    ax7 = plt.subplot2grid((4,30), (2,0), colspan=21, sharex=ax1)
    ax8, ax9 = plt.subplot2grid((4,30), (2,21), colspan=5, sharey=ax4, sharex=ax2), plt.subplot2grid((4,30), (2,26), colspan=4, sharey=ax4, sharex=ax3)
    ax10 = plt.subplot2grid((4,30), (3,0), colspan=21, sharex=ax1)
    ax11, ax12 = plt.subplot2grid((4,30), (3,21), colspan=5, sharey=ax4, sharex=ax2), plt.subplot2grid((4,30), (3,26), colspan=4, sharey=ax4, sharex=ax3)
    ax13, ax14 = ax3.twinx(), ax9.twinx()
    chi_daily = [0.237026614506512, 0.8669166986037341, 0.80789838336188657, 0.10367476695375623]
    ax13.plot([1, 2, 3, 4], chi_daily, marker='o')
    ax14.plot([1, 2, 3, 4], chi_daily, marker='o')
    ax13.set_ylim([0, 1])
    ax14.set_ylim([0, 1])
    ax13.set_yticks([0, 0.5, 1])
    ax14.set_yticks([0, 0.5, 1])
    ax13.set_ylabel(r'Chitinase $\mu$M day$^{-1}$')
    ax14.set_ylabel(r'Chitinase $\mu$M day$^{-1}$')
    axes = [[ax1, ax2, ax3], [ax4, ax5, ax6], [ax7, ax8, ax9], [ax10, ax11, ax12]]
    for a in range(len(axes)):
        for b in range(len(axes[a])):
            if b != 0:
                plt.setp(axes[a][b].get_yticklabels(), visible=False)
            if a != 3:
                plt.setp(axes[a][b].get_xticklabels(), visible=False)
    all_data, all_data_std, all_data_R = [], [], []
    for f in range(len(files)):
        labels, data, data_std, data_R, gens, gens_text, gens_R, gens_text_R = get_data(files[f])
        all_data.append(data)
        all_data_std.append(data_std)
        all_data_R.append(data_R)
    to_plot = []
    for a in range(len(all_data[0][0])):
        adding = False
        for b in range(len(all_data)):
            for c in range(len(all_data[b])):
                if sum(all_data[b][c][a]) > 0:
                    adding = True
        if adding:
            to_plot.append(a)
    colors = get_distinct_colors(len(to_plot))
    gens = [gens, [16, 17, 18, 19, 20], [1, 2, 3, 4]]
    for a in range(len(all_data)):
        ax = axes[a]
        data = all_data[a]
        new_data = [[], [], []]
        for b in range(len(data)):
            new_labels = []
            for c in range(len(to_plot)):
                i = to_plot[c]
                new_data[b].append(data[b][i])
                new_labels.append(labels[i])
        for d in range(len(new_data)):
            this_data = numpy.array(new_data[d])
            bottom = numpy.cumsum(this_data, axis=0)
            count = [0, 0, 0, 0, 0]
            if d == 0:
                count[0] += this_data[0][16]
                count[1] += this_data[0][17]
                count[2] += this_data[0][18]
                count[3] += this_data[0][19]
                count[4] += this_data[0][20]
            elif d == 1:
                count[0] += this_data[0][0]
                count[1] += this_data[0][1]
                count[2] += this_data[0][2]
                count[3] += this_data[0][3]
                count[4] += this_data[0][4]
            ax[d].bar(gens[d], this_data[0], color=colors[0], label=new_labels[0], edgecolor='k', linewidth=1)
            for j in range(1, this_data.shape[0]):
                if d == 0:
                    count[0] += this_data[j][16]
                    count[1] += this_data[j][17]
                    count[2] += this_data[j][18]
                    count[3] += this_data[j][19]
                    count[4] += this_data[j][20]
                elif d == 1:
                    count[0] += this_data[j][0]
                    count[1] += this_data[j][1]
                    count[2] += this_data[j][2]
                    count[3] += this_data[j][3]
                    count[4] += this_data[j][4]
                ax[d].bar(gens[d], this_data[j], color=colors[j], bottom=bottom[j-1], label=new_labels[j], edgecolor='k', linewidth=0.5)
            if d == 0 or d == 1:
                m = numpy.mean(count)
                ax[d].plot([15.5, 20.5], [m, m], 'k--', markeredgecolor='k', linewidth=0.5)
                ax[d].plot([15.5, 20.5], [m, m], 'k--', markeredgecolor='k', linewidth=1)
    axes[0][2].legend(bbox_to_anchor=(2,1.1), fontsize=8)
    for a in range(len(axes)):
        plt.sca(axes[a][0])
        plt.ylabel('Gene copies\n(per bacterium)')
        plt.yticks([0, 0.5, 1, 1.5])
    plt.setp(axes[3][0], xticks=gens[0], xticklabels=gens_text)
    axes[3][0].set_xlabel('Generation')
    axes[3][0].set_xlim([-0.6, 20.6])
    axes[3][1].set_xlabel('Generation')
    axes[3][2].set_xlabel('Days')
    plt.setp(axes[3][1], xticks=gens[1])
    plt.setp(axes[3][2], xticks=gens[2])
    axes[0][0].set_title('Nine-day\nincubation')
    axes[0][1].set_title('Four-day\nincubation')
    axes[0][2].set_title('Daily\nincubation')
    axes[0][0].text(-4.5, 0.75, 'K01183', va='center', ha='center', rotation=90, fontsize=10)
    axes[1][0].text(-4.5, 0.75, 'K01452', va='center', ha='center', rotation=90, fontsize=10)
    axes[2][0].text(-4.5, 0.75, 'K01207 and K12373', va='center', ha='center', rotation=90, fontsize=10)
    axes[3][0].text(-4.5, 0.75, 'K00884, K01443, K00884,\nK18676 and K02564', va='center', ha='center', rotation=90, fontsize=10)
    ax3.set_ylim([0, 1.5])
    plt.tight_layout()
    plt.savefig('Combined PICRUSt.png', dpi=600, bbox_inches='tight')
                        
        
    return

files = ['New_file_1.csv', 'New_file_2.csv', 'New_file_3.csv', 'New_file_4.csv']
make_plot(files)