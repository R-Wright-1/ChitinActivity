import matplotlib.pyplot as plt
import scipy.stats as stats
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os

random = [0.0, 0.058375860611402985, 0.24134653330134098, 0.30966397906632548, 0.16062308163669739, 0.09440327549693904, 0.63567563474268984, 0.24953138697063479, 0.18180934172439142, 0.40243834057657946, 0.64420196547105346, 0.16333591260416894, 0.18819786067418864, 0.1463188871570934, 0.14400854585510034, 1.1438714033529274, 0.19203342058044198, 0.26308971414464594, 0.16021021625589066, 0.14637100190351432, 0.934276973755922]
random_error = [0.033071838115075067, 0.0085163916608636397, 0.10378724504401736, 0.11579195936272399, 0.061542383196122603, 0.017291368564817979, 0.3070362003620497, 0.072263295772364941, 0.062160521941519269, 0.28773701476517921, 0.26048500578777539, 0.04262759254940187, 0.054964118929279208, 0.040848253990708391, 0.028491640977595344, 0.21453232127574201, 0.062717771920727833, 0.053574764633128581, 0.029409225403674917, 0.036568014577673419, 0.19990519763619125]
good = [0.0, 0.017920536614005448, 0.10805971996112582, 0.10687069021772959, 0.25802609272303406, 0.15542142177689128, 0.66007946333179646, 0.42149028746098838, 0.2619912007764349, 0.34298475729854405, 0.29268209471849016, 0.062392440878458286, 0.070198587524897063, 0.084317671703386152, 0.095619572791297305, 0.97202012469759269, 0.10780704625230569, 0.16234494173743455, 0.16515224903390632, 0.086043385387887666, 0.96455279717118969]
good_error = [0.010957736557575124, 0.018420947112044339, 0.16445988574436288, 0.066755394726016951, 0.19373283883821379, 0.17444752070393527, 0.37381985143207536, 0.3436519053827285, 0.27309521565662265, 0.21632693042152326, 0.26007607353002499, 0.030710547903890148, 0.023369141852674614, 0.055181866234432317, 0.16451364634033055, 1.2396990954770788, 0.040952528883283407, 0.18223840611295206, 0.095448512518872333, 0.048285485972044129, 0.19664310945739275]
short = [0.97202012469759269, 0.25146295260462104, 0.29757596956213128, 0.19017258825264544, 0.13456934121016442, 0.76801885637873224] 
short_error = [1.2396990954770788, 0.086404037477184947, 0.26881617203395314, 0.064526408318092449, 0.032893954154043718, 0.22563467289662209]
short_random = [1.1438714033529274, 0.1353301929764055, 0.13663989387460368, 0.037774560738738486, 0.043134886121428247, 0.15750716402217033]
short_random_error = [0.21453232127574201, 0.045944293481570447, 0.11630722032400587, 0.015174758217427224, 0.017512119253744723, 0.075259418894834312]
del short[0]
chi_daily = [0.237026614506512, 0.8669166986037341, 0.80789838336188657, 0.10367476695375623]
chi_daily_error = [0.059736179903026806, 0.22382194238678188, 0.1198038210570813, 0.042979078330859642]
chi_daily_r = [0.060815516547624314, 0.20574092459971904, 0.28381263064226919, 0.10472441680031566]
chi_daily_r_error = [0.023864665241773816, 0.0903933660773827, 0.075946470829352422, 0.042281604809644227]

genes = [0, 17.33333333, 0, 21, 21, 0, 3, 0, 4, 12, 21.66666667, 16.66666667, 41.5, 12, 14.66666667, 0, 54, 29, 26, 16, 28.66666667]
genes_short = [194, 317, 111.6666667, 198.6666667, 584]
genes_days = [394, 584, 337.3333333, 208]
gens = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
gen_ticks = ['0', '1', '2', '3', '4', '*', '6', '*', '8', '9', '10', '11', '12', '13', '14', '*', '16', '17', '18', '19', '20']
gens_s = [16, 17, 18, 19, 20]
days = [1, 2, 3, 4]

gens_g, gens_s_g, days_g = [], [],[]
for a in range(len(gens)):
    gens_g.append(gens[a]-0.5)
for b in range(len(gens_s)):
    gens_s_g.append(gens_s[b]-0.5)
for c in range(len(days)):
    days_g.append(days[c]-0.5)

short[4] = chi_daily[3]
genes_short[4] = genes_days[3]

fig = plt.figure(figsize=(8.27, 6))
ax1 = plt.subplot2grid((3, 10), (0, 0), colspan=10)
ax2 = plt.twinx(ax1)
ax3 = plt.subplot2grid((3, 10), (1, 7), colspan=3)
pos1 = ax3.get_position() # get the original position 
pos2 = [pos1.x0 + 0.1, pos1.y0,  pos1.width / 2.0, pos1.height] 
ax3.set_position(pos2)
"""
cax.set_xticklabels(['']*10)
cax.set_yticklabels(['']*10)
cax.tick_params(top="off", right="off", bottom="off", left="off")
divider = make_axes_locatable(cax)
ax3 = divider.append_axes('right', size="35%", pad=0, axisbg='none', frameon=True)
"""
ax4 = plt.twinx(ax3)
ax5 = plt.subplot2grid((3, 10), (2, 7), colspan=3)
ax6 = plt.twinx(ax5)

ax1.bar(gens_g, genes)
ax1.xaxis.set_ticks(gens)
ax1.set_xticklabels(gen_ticks)
plt.setp(ax2.get_xticklabels(), visible=False)
ax2.plot(gens, good, 'o-', color='gray')
ax3.bar(gens_s_g, genes_short)
ax5.bar(days_g, genes_days)
ax6.plot(days, chi_daily, 'o-', color='gray')
ax1.set_xlim(-0.5, 20.75), ax3.set_xlim(15.32, 20.75), ax5.set_xlim(0.35, 4.67)
ax1.set_xlabel('Generation')
ax3.set_xticks([16, 17, 18, 19, 20])
ax3.set_xlabel('Generation')
"""
caxB = plt.subplot2grid((3, 10), (1, 0), colspan=10, frameon=False)
caxB.set_xticklabels(['']*10)
caxB.set_yticklabels(['']*10)
caxB.tick_params(top="off", right="off", bottom="off", left="off")
divider = make_axes_locatable(caxB)
ax4 = divider.append_axes('right', size="35%", pad=0, axisbg='none', frameon=True)
ax4.tick_params(top="off", right="off", bottom="off", left="off")
ax4.set_xticklabels(['']*10)
ax4.yaxis.set_label_position("right")
"""
ax4.plot(gens_s, short, 'o-', color='gray')
#ax5.set_xticks([1, 2, 3, 4])
plt.tight_layout()
#ax2.set_yticks([0, 1, 2, 3, 4])
#ax6.set_yticks([0, 0.4, 0.8, 1.2, 1.6, 2.0])
ax5.set_xlabel('Day')
ax1.set_ylabel('Gene copies'), ax3.set_ylabel('Gene copies'), ax5.set_ylabel('Gene copies')
ax2.set_ylabel('Nine-day \n')
ax4.set_ylabel('Four-day \n'+r'Chitinase $\mu$M day$^{-1}$')
ax6.set_ylabel('Daily \n')
#fig.subplots_adjust(left=0.09,right=0.96,top=0.96,bottom=0.08,wspace=0, hspace=0.3)
os.chdir('/Users/u1560915/Documents/GitHub/ChitinActivity/MiSeq/18S/Figures/')
plt.savefig('Fig supp 13.pdf', bbox_inches='tight')
plt.close()

new_genes, new_good = [], []
for a in range(len(genes)):
    if a != 2 and a != 5 and a != 7 and a != 15:
        new_genes.append(genes[a])
        new_good.append(good[a])

gradient, intercept, r_value, p_value, std_err = stats.linregress(new_genes, new_good)
r2_long = r_value**2
gradient, intercept, r_value, p_value, std_err = stats.linregress(genes_short, short)
r2_short = r_value**2
gradient, intercept, r_value, p_value, std_err = stats.linregress(genes_days, chi_daily)
r2_days = r_value**2
r2 = [r2_long, r2_short, r2_days]
#print r2