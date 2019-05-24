import matplotlib.pyplot as plt
import scipy.stats as stats
from matplotlib.lines import Line2D
import numpy

random = [0.916, 1.065333333, 0.721333333, 0.709333333, 0.692, 0.689333333, 0.706666667, 0.698666667, 1.312666667, 0.84, 0.962666667, 1.074, 1.238, 0.249, 0.261333333, 0, 0.0244, 0, 0.033733333, 0.0546, 0.048333333]
good = [1.543333333, 0.86, 0.712, 0.696, 0.702666667, 0.713333333, 0.685, 0.706666667, 0.693333333, 0.788, 1.087333333, 0.729333333, 0.7, 0.2472, 0.222, 0, 0.0252, 0.0803, 0.0204, 0.0212, 0.0432]
short_random = [0.0332, 0.108, 0.0975, 0.190666667, 0.777333333]
short_good = [0.029066667, 0.1408, 5.346666667, 1.042, 3.646666667]
short_good_g20 = [0.219, 0.689, 10.38, 3.646666667]
short_random_g20 = [0.237333333, 0.168333333, 0.470333333, 0.777333333]

random_e = [0.272939554, 0.567046147, 0.04387862, 0.028936713, 0.030199338, 0.032083225, 0.052204725, 0.025716402, 0.296784995, 0.208, 0.037166293, 0.139097088, 0.220027271, 0.152934627, 0.146677651, 0, 0, 0, 0.013724917, 0.008202439, 0.03095179]
good_e = [0.45236416, 0.242783854, 0.004, 0.010583005, 0.108024689, 0.03494758, 0.006244998, 0.012858201, 0.114146105, 0.078689262, 0.189318075, 0.032578111, 0.028, 0.319446521, 0.039344631, 0, 0, 0.080185909, 0, 0, 0.007076722]
short_random_e = [0, 0, 0.055861436, 0.087305975, 0.200013333]
short_good_e = [0.013648932, 0.133187537, 3.375519713, 0.497485678, 1.831402013]
short_random_g20_e = [0.04291076, 0.018502252, 0.11072639, 0.200013333]
short_good_g20_e = [0.036755952, 0.444635806, 5.191415992, 1.831402013]

means_lm =  [[0.0, 0.22333332919374846, 0.4168673314450024, 0.427694586822124, 0.6569287864540203, 0.29840488739094206, 0.38828377583287205, 3.829480102783878, 0.4660025383948234, 0.4442676778826085, 0.7018065951173557, 0.26936130779203654, 0.9252275008603138, 0.3697093941768678, 0.5972975631025822, 0.8037397522474726, 0.534932527325001, 1.6608330056386031, 0.26643843479880513, 0.7390817800351512, 1.5679979713684673], [0.0, 0.058375860611343436, 0.24134653330133793, 0.309663979066336, 0.16062308163667485, 0.09440327549689483, 0.6356756347426953, 0.2495313869706358, 0.1818093417243752, 0.40243834057658595, 0.6442019654710809, 0.16333591260414806, 0.18819786067417466, 0.14631888715706717, 0.14400854585507386, 1.1438714033528943, 0.19203342058042863, 0.26308971414465077, 0.16021021625586931, 0.1463710019034883, 0.934276973755944], [0.0, 0.01792053661401457, 0.10805971996112476, 0.10687069021772629, 0.2580260927230228, 0.15542142177688684, 0.6600794633317828, 0.4214902874609766, 0.26199120077642674, 0.34298475729852984, 0.29268209471847956, 0.06239244087845939, 0.0701985875248969, 0.08431767170338499, 0.09561957279129812, 0.9720201246975796, 0.10780704625230174, 0.1623449417374287, 0.16515224903389875, 0.08604338538788642, 0.9645527971711767]]
means_sm =  [[1.1438714033528943, 0.1353301929763973, 0.13663989387459582, 0.037774560738732414, 0.043134886121421995, 0.157507164022162], [0.9720201246975796, 0.25146295260461193, 0.29757596956212234, 0.19017258825263672, 0.1345693412101562, 0.7680188563787238]]
means_le =  [[0.23396465017267526, 0.0982868836541131, 0.20094772300387168, 0.16468746147983446, 0.7367847145722678, 0.062425840876510735, 0.1916156400597665, 1.0277998903888164, 0.2197206983748821, 0.47110132615972916, 0.2945742067234798, 0.07420509884536076, 0.2963598550053271, 0.12512223650588888, 0.21521684073433017, 0.32683826491220835, 0.374498430951749, 0.8010220521729765, 0.08449482668312734, 0.46763676887115374, 0.9119130320379032], [0.03307183811508742, 0.008516391660867489, 0.10378724504403757, 0.1157919593627446, 0.0615423831961383, 0.017291368564824744, 0.3070362003620273, 0.07226329577237953, 0.062160521941535984, 0.28773701476517577, 0.26048500578775563, 0.04262759254941461, 0.05496411892929351, 0.04084825399072193, 0.028491640977604503, 0.21453232127569274, 0.0627177719207433, 0.05357476463313977, 0.029409225403683736, 0.03656801457768506, 0.1999051976361676], [0.010957736557572819, 0.018420947112041182, 0.16445988574435735, 0.06675539472601204, 0.1937328388382084, 0.17444752070392783, 0.3738198514320746, 0.34365190538272494, 0.27309521565661404, 0.2163269304215169, 0.26007607353001994, 0.030710547903886297, 0.023369141852671648, 0.05518186623442764, 0.16451364634032425, 1.239699095477084, 0.040952528883279944, 0.18223840611294803, 0.09544851251886675, 0.04828548597203891, 0.1966431094573959]]
means_se =  [[0.21453232127569274, 0.0459442934815698, 0.11630722032400537, 0.01517475821742674, 0.017512119253744143, 0.07525941889483373], [1.239699095477084, 0.0864040374771847, 0.2688161720339535, 0.06452640831809206, 0.03289395415404327, 0.2256346728966229]]
high_lm =  [[0.0, 0.42942701062783706, 0.8952922113038461, 0.7845394371777487, 1.5697181290044366, 0.41866106141456494, 0.7762906498443747, 5.924010738155844, 0.6253916164345085, 1.7098709025209304, 1.3420724744033972, 0.43535499807301864, 1.5809972771886143, 0.6098375874301897, 1.0929011275628044, 1.531650815280642, 1.518243413755056, 2.9327375841122976, 0.4478109318361825, 1.548458459935044, 3.19305084723077], [0.0, 0.055787744561071904, 0.1977334385944547, 0.40289839685871714, 0.11320979376185497, 0.0834252413279753, 0.827120221352763, 0.3536946196186723, 0.1360757418467968, 0.5808161610958563, 0.6770576599488907, 0.15636653214926766, 0.14229152585042762, 0.12877847544162122, 0.13825991526190826, 0.9362907669499761, 0.18202097226892502, 0.3024113775264294, 0.15505295713126618, 0.1305904607115641, 0.922414783753914], [0.0, 0.058417920889195916, 0.49210619391890625, 0.27725964714028495, 0.3519569152443998, 0.40619853651536975, 1.0727540795423007, 1.0842478080218683, 0.8356882177016774, 0.7537390737197187, 0.8391340236810162, 0.11779409114725814, 0.10533969173918588, 0.2279655750565044, 0.46866529829455206, 3.054528435077479, 0.200108107145315, 0.56149340784675, 0.37405815803626097, 0.16648247933765556, 0.9035476321874928]]
high_sm =  [[0.9362907669499761, 0.12497142569080016, 0.07642234972106486, 0.04796358073801014, 0.050683583615966636, 0.22039840838766434], [3.054528435077479, 0.4557985537640005, 0.9092175343549779, 0.3131633845834729, 0.19110119009476803, 0.7823842740364343]]
high_le =  [[0.20883748626147158, 0.18403297642655025, 0.1696181496455233, 0.037614501228673144, 0.2855395605543355, 0.02899584471250966, 0.09702903015406612, 0.6481551184544605, 0.21829099089875398, 0.39648493515322975, 0.16576147409650593, 0.04704018065187689, 0.1966408485052826, 0.10324546388705161, 0.1594351607044061, 0.25565170681623883, 0.393774905976762, 0.10590126892208823, 0.04058402220427638, 0.5491239562784132, 0.641495161442972], [0.03971334491423845, 0.0036841630435981566, 0.08157929239084401, 0.16011234829367416, 0.011487350677031207, 0.0037005664344816874, 0.22220561355752147, 0.06239080882848145, 0.03884787016087421, 0.3205321928680813, 0.14935038057686198, 0.013627533533708231, 0.01983938631924791, 0.010397564009135163, 0.005871991076028715, 0.29556149689603595, 0.04343285993662982, 0.07342429092513868, 0.017591365382848088, 0.015378716924085625, 0.23311236040182934], [0.023247717701222068, 0.0345230389853913, 0.31149626223708493, 0.04762808137589635, 0.3587524629973905, 0.07237329311696765, 0.07043679163928816, 0.037209113542349216, 0.14547252543064154, 0.07711524144820418, 0.2434639390315709, 0.032048613393871606, 0.010630798371865124, 0.026883970599262657, 0.3146382180604768, 3.090787848259323, 0.028714725029683245, 0.36301784905569834, 0.0053450119513068755, 0.026916813082435456, 0.0987793838923987]]
high_se =  [[0.29556149689603595, 0.015219580910988552, 0.01695979225961966, 0.03562951065388353, 0.022598635670800254, 0.16093250334907508], [3.090787848259323, 0.05203781368979842, 0.51853961696074, 0.09165118994959569, 0.015322102482793405, 0.31935590747427506]]


g20_rm = [0.060815516547624696, 0.20574092459971866, 0.28381263064226686, 0.10472441680031778]
g20_gm = [0.23702661450651122, 0.8669166986037344, 0.8078983833618888, 0.10367476695375828]
g20_re = [0.023864665241773574, 0.09039336607738284, 0.07594647082935443, 0.04228160480964478]
g20_ge = [0.0597361799030266, 0.2238219423867821, 0.11980382105708018, 0.04297907833086022]


g20_g = [[0.26324135, 0.29315721, 0.19528355, 0.16024346, 0.12788805,
       0.16423899, 0.17642368, 0.24198971, 0.28022315, 0.23213431,
       0.14560475, 0.27472108, 0.13276556, 0.23797336, 0.24827147,
       0.27107701, 0.36996646, 0.30755284, 0.17912333, 0.18783982,
       0.22218193, 0.28383948, 0.24605778, 0.20633663, 0.30979065,
       0.3012303 , 0.20911751, 0.2366384 , 0.31014568, 0.29574094], [0.66061213, 0.9501739 , 0.65703477, 0.51987022, 0.43168985,
       0.82366392, 0.99493963, 0.99493963, 1.18529844, 1.06675017,
       0.73105044, 0.54859267, 0.45690861, 1.038976  , 0.95664567,
       1.17071762, 0.99797978, 0.93744314, 0.68390636, 0.76437236,
       0.93436778, 1.17558105, 1.10505223, 0.79290395, 0.9563926 ,
       0.92166171, 0.92064759, 1.19037215, 0.95981499, 0.47914161], [0.76876268, 0.65311382, 0.87348703, 0.79132574, 0.76056366,
       0.63763315, 0.87436465, 0.82473662, 0.86507347, 0.65462316,
       0.84900961, 0.66870053, 0.92543552, 0.68796403, 1.03757426,
       0.74720788, 0.82902974, 0.72721649, 0.75150961, 0.99439797,
       0.71778469, 0.83791497, 0.84770807, 0.68352696, 0.59795877,
       0.88129231, 0.7722787 , 0.94231052, 1.06572624, 0.96872065], [0.11634279, 0.2589404 , 0.1013821 , 0.10344999, 0.11146716,
       0.05329645, 0.0871048 , 0.0818458 , 0.05914638, 0.04722424,
       0.08856224, 0.12926018, 0.17060077, 0.09115218, 0.14351216,
       0.10709155, 0.06962436, 0.06077153, 0.09162303, 0.08972629,
       0.18290363, 0.12713835, 0.13691723, 0.11950547, 0.08958495,
       0.08742026, 0.07202188, 0.0952268 , 0.06902666, 0.06837338]]

g20_r = [[0.05812994, 0.03189649, 0.06409732, 0.06107695, 0.03515069,
       0.05488038, 0.05722084, 0.06900323, 0.0413578 , 0.04214756,
       0.06495181, 0.07797395, 0.060006  , 0.05607453, 0.05368982,
       0.10717691, 0.05203751, 0.0508724 , 0.06303597, 0.05356247,
       0.04015184, 0.04812561, 0.14064776, 0.12082043, 0.03736188,
       0.07266596, 0.06937343, 0.05295177, 0.04177715, 0.04624711], [0.17310858, 0.0805066 , 0.23220486, 0.11314647, 0.24008675,
       0.15527706, 0.17330204, 0.14920138, 0.32830571, 0.09478616,
       0.17402268, 0.16232798, 0.20350695, 0.16913426, 0.15808896,
       0.5385761 , 0.27755293, 0.08997329, 0.21181356, 0.21594641,
       0.25467738, 0.2202588 , 0.29163322, 0.23404781, 0.12161861,
       0.19411003, 0.37284752, 0.18247121, 0.1819482 , 0.17774621], [0.33483886, 0.54080504, 0.31047981, 0.37748398, 0.32346989,
       0.20660456, 0.28802985, 0.23303495, 0.28179785, 0.22024964,
       0.33396111, 0.31671463, 0.38230212, 0.26811343, 0.34372513,
       0.26460589, 0.25774137, 0.21925334, 0.23439546, 0.26986398,
       0.30148243, 0.3037478 , 0.40145559, 0.27739597, 0.20916136,
       0.24627142, 0.18368432, 0.2056446 , 0.19475157, 0.18331295], [0.11973646, 0.26169612, 0.1062543 , 0.10927753, 0.11188209,
       0.05549587, 0.07623004, 0.08461941, 0.06025199, 0.04857327,
       0.09158747, 0.12634533, 0.17593635, 0.09377316, 0.15285713,
       0.10964313, 0.07078308, 0.06181442, 0.09440943, 0.09222821,
       0.1506594 , 0.13023613, 0.1427817 , 0.12394591, 0.09117881,
       0.08855344, 0.07349174, 0.09804506, 0.06971456, 0.06973096]]

n = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29]
g20_d1 = sorted(zip(g20_g[1], n))
print(g20_d1)


g20_g_new, g20_r_new = [], []
g20_g_new_sd, g20_r_new_sd = [], []
for a in g20_g:
    day = []
    for b in range(len(a)):
        if b == 21 or b == 8 or b == 27:
            day.append(a[b])
    g20_g_new.append(numpy.mean(day))
    g20_g_new_sd.append(numpy.std(day))
for a in g20_r:
    day = []
    for b in range(len(a)):
        if b == 17 or b == 15 or b == 12:
            day.append(a[b])
    g20_r_new.append(numpy.mean(day))
    g20_r_new_sd.append(numpy.std(day))

gens_l, gens_s = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20], [15, 16, 17, 18, 19, 20]

plt.figure(figsize=(20,10))
ax_gm = plt.subplot2grid((2,8), (0,0), colspan=4)
ax_rm = plt.subplot2grid((2,8), (0,4), colspan=4)
ax_sgm = plt.subplot2grid((2,8), (1,0), colspan=2)
ax_srm = plt.subplot2grid((2,8), (1,4), colspan=2)
ax_gh_DNA = ax_gm.twinx()
ax_sgh_DNA = ax_sgm.twinx()
ax_rh_DNA = ax_rm.twinx()
ax_srh_DNA = ax_srm.twinx()

ax_sgm20 = plt.subplot2grid((2,8), (1,2), colspan=2)
ax_srm20 = plt.subplot2grid((2,8), (1,6), colspan=2)
ax_sgm20_DNA = ax_sgm20.twinx()
ax_srm20_DNA = ax_srm20.twinx()

all_ax = [ax_gm, ax_rm, ax_sgm, ax_srm, ax_gh_DNA, ax_sgh_DNA, ax_rh_DNA, ax_srh_DNA, ax_sgm20, ax_srm20, ax_sgm20_DNA, ax_srm20_DNA]

ax_gm.errorbar(gens_l, means_lm[2], yerr=means_le[2], marker='o', markeredgecolor='k', capsize=2)
plt.sca(ax_gm)
plt.xticks(gens_l, gens_l)
plt.xlabel('Generation')
plt.ylabel(r'Chitinase $\mu$M day$^{-1}$')
ax_rm.errorbar(gens_l, means_lm[1], yerr=means_le[1], marker='o', markeredgecolor='k', capsize=2, label='Mean chitinase activity')
plt.sca(ax_rm)
plt.xticks(gens_l, gens_l)
plt.xlabel('Generation')
plt.ylabel(r'Chitinase $\mu$M day$^{-1}$')
ax_gm.errorbar(gens_l, high_lm[2], yerr=high_le[2], marker='^', color='y', markeredgecolor='k', capsize=2)
ax_rm.errorbar(gens_l, high_lm[1], yerr=high_le[1], marker='^', color='y', markeredgecolor='k', capsize=2, label='Highest chitinase activity')
ax_sgm.errorbar(gens_s[1:], means_sm[1][1:], yerr=means_se[1][1:], marker='o', markeredgecolor='k', capsize=2)
plt.sca(ax_sgm)
plt.xticks(gens_s, gens_s)
plt.xlabel('Generation')
plt.ylabel(r'Chitinase $\mu$M day$^{-1}$')
ax_srm.errorbar(gens_s[1:], means_sm[0][1:], yerr=means_se[0][1:], marker='o', markeredgecolor='k', capsize=2)
plt.sca(ax_srm)
plt.xticks(gens_l, gens_l)
plt.xlabel('Generation')
plt.ylabel(r'Chitinase $\mu$M day$^{-1}$')
ax_sgm.errorbar(gens_s[1:], high_sm[1][1:], yerr=high_se[1][1:], marker='^', color='y', markeredgecolor='k', capsize=2)
ax_srm.errorbar(gens_s[1:], high_sm[0][1:], yerr=high_se[0][1:], marker='^', color='y', markeredgecolor='k', capsize=2)

ax_gh_DNA.errorbar(gens_l, good, yerr=good_e, marker='s', color='r', markeredgecolor='k', capsize=2)
ax_gh_DNA.set_ylabel(r'DNA concentration $\mu$g mL$^{-1}$')
ax_rh_DNA.errorbar(gens_l, random, yerr=random_e, marker='s', color='r', markeredgecolor='k', capsize=2, label='DNA concentration')
ax_rh_DNA.set_ylabel(r'DNA concentration $\mu$g mL$^{-1}$')
ax_sgh_DNA.errorbar(gens_s[1:], short_good, yerr=short_good_e, marker='s', color='r', markeredgecolor='k', capsize=2)
ax_sgh_DNA.set_ylabel(r'DNA concentration $\mu$g mL$^{-1}$')
ax_srh_DNA.errorbar(gens_s[1:], short_random, yerr=short_random_e, marker='s', color='r', markeredgecolor='k', capsize=2)
ax_srh_DNA.set_ylabel(r'DNA concentration $\mu$g mL$^{-1}$')

ax_sgm20.errorbar([1,2,3,4], g20_gm, yerr=g20_ge, marker='o', markeredgecolor='k', capsize=2)
ax_sgm20.errorbar([1,2,3,4], g20_g_new, yerr=g20_g_new_sd, marker='^', color='y', markeredgecolor='k', capsize=2)
plt.sca(ax_sgm20)
plt.xticks([1,2,3,4], [1,2,3,4])
plt.xlabel('Day')
plt.ylabel(r'Chitinase $\mu$M day$^{-1}$')
ax_srm20.errorbar([1,2,3,4], g20_rm, yerr=g20_re, marker='o', markeredgecolor='k', capsize=2)
ax_srm20.errorbar([1,2,3,4], g20_r_new, yerr=g20_r_new_sd, marker='^', color='y', markeredgecolor='k', capsize=2)
plt.sca(ax_srm20)
plt.xticks([1,2,3,4], [1,2,3,4])
plt.xlabel('Day')
plt.ylabel(r'Chitinase $\mu$M day$^{-1}$')
ax_sgm20_DNA.errorbar([1,2,3,4], short_good_g20, yerr=short_good_g20_e, marker='s', color='r', markeredgecolor='k', capsize=2)
ax_sgm20_DNA.set_ylabel(r'DNA concentration $\mu$g mL$^{-1}$')
ax_srm20_DNA.errorbar([1,2,3,4], short_random_g20, yerr=short_random_g20_e, marker='s', color='r', markeredgecolor='k', capsize=2)
ax_srm20_DNA.set_ylabel(r'DNA concentration $\mu$g mL$^{-1}$')

gradient, intercept, r_value, p_value, std_err = stats.linregress(high_lm[2], good)
r2 = '$r^{2}$ = '+str(round(r_value**2, 2))
ax_gm.set_title('Positive selection\nNine day\n'+r2)
gradient, intercept, r_value, p_value, std_err = stats.linregress(high_lm[1], random)
r2 = '$r^{2}$ = '+str(round(r_value**2, 2))
ax_rm.set_title('Random selection\nNine day\n'+r2)
gradient, intercept, r_value, p_value, std_err = stats.linregress(high_sm[1][1:], short_good)
r2 = '$r^{2}$ = '+str(round(r_value**2, 2))
ax_sgm.set_title('Positive selection\nFour day\n'+r2)
gradient, intercept, r_value, p_value, std_err = stats.linregress(high_sm[0][1:], short_random)
r2 = '$r^{2}$ = '+str(round(r_value**2, 2))
ax_srm.set_title('Random selection\nFour day\n'+r2)
gradient, intercept, r_value, p_value, std_err = stats.linregress(g20_g_new, short_good_g20)
r2 = '$r^{2}$ = '+str(round(r_value**2, 2))
ax_sgm20.set_title('Positive selection\nFour day\nWithin generation 20\n'+r2)
gradient, intercept, r_value, p_value, std_err = stats.linregress(g20_r_new, short_random_g20)
r2 = '$r^{2}$ = '+str(round(r_value**2, 2))
ax_srm20.set_title('Random selection\nFour day\nWithin generation 20\n'+r2)

for a in all_ax:
    plt.sca(a)
    plt.gca().set_ylim(bottom=0)
legend_elements = [Line2D([0], [0], marker='o', label='Mean chitinase activity', markeredgecolor='k'),
                   Line2D([0], [0], marker='^', label='Next generation chitinase activity',
                          markerfacecolor='y', markeredgecolor='k'),
                   Line2D([0], [0], marker='s', label='DNA concentration',
                          markerfacecolor='r', markeredgecolor='k')]



ax_rm.legend(handles=legend_elements, bbox_to_anchor=(1.1,1.05), handlelength=0)
plt.subplots_adjust(hspace=0.5, wspace=2)
plt.savefig('Supplementary DNA conc.png', bbox_inches='tight')
plt.close()






DNA_conc_random = [0.596, 0.112, 0.03, 0.08, 0.18, 0.303, 0.263, 0.221]
DNA_conc_high = [0.976, 0.188, 0.116, 0.5, 0.738, 0.850, 0.147, 0.63]

x = [2, 5, 11, 18, 26, 34, 42, 50]
Random_mean = [0.033810812335185014, 0.009411862351428272, 0.024827207203976177, 0.0819028786459898, 0.028648271772991504, 0.03910331010560018, 0.113527316065263, 0.18339292953767153]
Random_errors = [0.008480201348960947, 0.008264783358440934, 0.007477586805943582, 0.02413392460168274, 0.02021991763289385, 0.029098794837049597, 0.04763021288613838, 0.05196808715496919]
Positive_mean = [0.026797767397356634, 0.00971705845591092, 5.843235849902216, 0.38786048290184477, 5.010631777487455, 9.931006801653002, 16.59400957401194, 26.15099167685104]
Positive_errors = [0.004139067265074105, 0.008731004559715004, 9.663143808414658, 0.14068644212078985, 7.53818377434487, 9.28847282828771, 25.034716421550566, 33.711103814460984]
Random_high_mean = [0.04290832669015363, 0.011072043055790055, 0.02733247749830303, 0.08639947012358877, 0.026146357092380532, 0.028692197807399417, 0.12061924631722415, 0.18863191659397913]
Random_high_errors = [0.016008627617656316, 0.0013539983419268552, 0.007699669892271953, 0.013488968864775405, 0.00621108305598326, 0.008204871422026826, 0.055738634903176695, 0.060236907804587285]
Positive_high_mean = [0.024988949920592186, 0.005727933158719187, 9.3723047489977, 0.5335749650122948, 21.96099696935843, 25.379669342499543, 32.94738932617725, 88.56244684536723]
Positive_high_errors = [0.002515279415819621, 0.000815963576787886, 11.419265584183243, 0.03941440272375856, 10.516210452311718, 13.05356533666748, 13.166724900782853, 21.451089276298326]

plt.figure(figsize=(10,8))
ax1 = plt.subplot(211)
ax1_DNA = ax1.twinx()
ax2 = plt.subplot(212)
ax2_DNA = ax2.twinx()

ax1.errorbar(x, Random_mean, yerr=Random_errors, marker='o', capsize=2, label='Mean chitinase activity', markeredgecolor='k')
ax1.errorbar(x, Random_high_mean, yerr=Random_high_errors, color='y', marker='^', capsize=2, label='Next generation chitinase activity', markeredgecolor='k')
ax1_DNA.plot(x, DNA_conc_random, color='r', marker='s', label='DNA concentration', markeredgecolor='k')
ax1.legend(handles=legend_elements, bbox_to_anchor=(1.1,1.05), handlelength=0)
ax1.set_xlabel('Day')
gradient, intercept, r_value, p_value, std_err = stats.linregress(Random_high_mean, DNA_conc_random)
r2 = '$r^{2}$ = '+str(round(r_value**2, 2))
ax1.set_title('Random selection\n'+r2)
ax1.set_ylabel(r'Chitinase $\mu$M day$^{-1}$')
ax1_DNA.set_ylabel(r'DNA concentration $\mu$g mL$^{-1}$')

ax2.errorbar(x, Positive_mean, yerr=Positive_errors, marker='o', capsize=2, markeredgecolor='k')
ax2.errorbar(x, Positive_high_mean, yerr=Positive_high_errors, color='y', marker='^', capsize=2, markeredgecolor='k')
ax2_DNA.plot(x, DNA_conc_high, color='r', marker='s', markeredgecolor='k')
gradient, intercept, r_value, p_value, std_err = stats.linregress(Positive_high_mean, DNA_conc_high)
r2 = '$r^{2}$ = '+str(round(r_value**2, 2))
ax2.set_title('Positive selection\n'+r2)
ax2.set_xlabel('Day')
ax2.set_ylabel(r'Chitinase $\mu$M day$^{-1}$')
ax2_DNA.set_ylabel(r'DNA concentration $\mu$g mL$^{-1}$')
#plt.tight_layout()
plt.subplots_adjust(hspace=0.5, wspace=2)
ax1.set_ylim(bottom=0)
ax1_DNA.set_ylim(bottom=0)
ax2.set_ylim(bottom=0)
ax2_DNA.set_ylim(bottom=0)
plt.savefig('DNA experiment 2.png', dpi=600, bbox_inches='tight')