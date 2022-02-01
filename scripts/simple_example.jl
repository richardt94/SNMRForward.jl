using SpecialFunctions, Hankel, Elliptic, PyPlot

#constants in SI units
mu_0 = 4*pi*10^-7

#frequency-domain field in horizontal wavenumber space

I = 1; #current (amps)
R = 10; #radius (metres)

Hz_wavenumber(κ,z) = I * R * besselj1(κ*R) * exp(- κ * abs(z))/2;


zgrid = 0:0.01*R:10*R
##
##
max_radius = 1.5*R;
hank_0 = QDHT(0,max_radius,201);

##
Hz_radius = hank_0 \ Hz_wavenumber.(hank_0.k, R)

Hzk_mat = [Hz_wavenumber(κ,z) for κ in hank_0.k, z in zgrid]
Hzr_mat = hank_0 \ Hzk_mat

plt.figure();
plt.pcolor(hank_0.r, zgrid, Hzr_mat', vmin = -0.06, vmax = 0.06);
gca().invert_yaxis();
xlabel("r (m)")
ylabel("z (m)")
display(gcf());
plt.close("all");

##

#analytic method for spatial fields
#(without numerical Hankel transform)
function Hz_analytic(r,z)
    αsq = R^2 + r^2 + z^2 - 2*r*R;
    βsq = R^2 + r^2 + z^2 + 2*r*R;
    β = sqrt(βsq);
    ksq = 1 - αsq/βsq;
    C = 1/pi;
    
    C/(2*αsq*β)*((R^2 - r^2 - z^2)*E(ksq) + αsq*K(ksq))
end
##

Hz_analytic_mat = [Hz_analytic(r,z) for r in hank_0.r, z in zgrid]

##
plt.figure();
plt.pcolor(hank_0.r, zgrid, Hz_analytic_mat', vmin = -0.06, vmax = 0.06);
gca().invert_yaxis();
xlabel("r (m)")
ylabel("z (m)")
display(gcf());
plt.close("all");

##
zgrid = 0:0.01*R:10*R
dipole_approximation = [25/((r^2+z^2)^(3/2)) * (3*z^2/(z^2 + r^2) - 1) for r in hank_0.r, z in zgrid]

##
plt.figure();
plt.pcolor(hank_0.r, zgrid, dipole_approximation', vmin = -0.06, vmax = 0.06);
gca().invert_yaxis();
xlabel("r (m)")
ylabel("z (m)")
display(gcf());
plt.close("all");

##
#far-field comparison
figure()
plot(zgrid[end-800:end-400], dipole_approximation[1,end-800:end-400])
plot(zgrid[end-800:end-400], Hzr_mat[1,end-800:end-400])
plot(zgrid[end-800:end-400], Hz_analytic_mat[1,end-800:end-400])
legend(["dipole", "numerical hankel transform", "analytic"])
xlabel("z (m)")
ylabel("Hz field (1/m)")
display(gcf())
close("all")

##
#Digital filter Hankel inverse with "lagged convolution" method
#radius min and max in m
rmin = 0.1
rmax = 15

#abscissa (for k*r)
const Filter_base = [4.118588707535708e-6,4.662307783048404e-6,5.2778064058937756e-6,5.97456061553286e-6,6.7632974390298035e-6,7.656160041269817e-6,8.666894677627122e-6,9.811062327352187e-6,1.1106278265924788e-5,1.2572483264760042e-5,1.42322505935799e-5,1.6111133551969887e-5,1.8238058880617237e-5,2.0645772109076137e-5,2.3371341696506312e-5,2.645672972698995e-5,2.994943794568826e-5,3.390323908202406e-5,3.8379004719130766e-5,4.3445642455208216e-5,4.918115678505129e-5,5.567385003477687e-5,6.302368183899532e-5,7.134380809054521e-5,8.076232305601659e-5,9.142423147817327e-5,0.00010349368102719458,0.00011715648947091054,0.00013262300547161834,0.00015013134705348249,0.0001699510675990275,0.00019238730581535294,0.00021778548356175115,0.0002465366238651324,0.00027908337099788406,0.00031592680530155527,0.0003576341576754271,0.00040484754250000483,0.0004582938434450188,0.0005187959043609725,0.0005872851975459913,0.0006648161644249478,0.0007525824494258378,0.000851935276985483,0.0009644042546116447,0.0010917209222795108,0.001235845410722264,0.0013989966190390965,0.0015836863762264436,0.001792758112573506,0.002029430636295734,0.0022973476893787268,0.002600634045580001,0.0029439590142573526,0.0033326083277104117,0.0037725655187922052,0.004270604041657015,0.0048343915539089625,0.005472607965649308,0.006195079072871529,0.007012927832585425,0.007938745608658544,0.008986786024826483,0.010173184409377162,0.011516206210016314,0.013036528203437736,0.014757556829019875,0.01670578854762277,0.018911217773465227,0.021407798659484324,0.024233967845691123,0.027433236218606032,0.03105485879233143,0.03515459302455746,0.039795557242315927,0.0450492023935578,0.050996412085361466,0.05772874784464333,0.06534985877304568,0.07397707729864236,0.08374322559219596,0.09479866045903013,0.10731358818908403,0.12148068500391276,0.13751806344428075,0.15567263036799733,0.1762238882567611,0.19948823835583873,0.22582385189647586,0.2556361843969779,0.28938421793905067,0.32758752752368947,0.37083428029819565,0.41979029080811425,0.47520927168614446,0.5379444375946745,0.6089616410728969,0.6893542425242224,0.7803599432780343,0.8833798408827509,1.0,1.1320158709991752,1.2814599321940212,1.4506329812931589,1.6421395578187052,1.858928041846342,2.104336046415478,2.382141802457978,2.6966223273530128,3.0526192726543444,3.4556134647626755,3.911809286149797,4.428230196243525,5.012826862585463,5.6745995670177445,6.423736771429134,7.271771976378781,8.231761287547819,9.318484423780736,10.548672261378394,11.941264417849103,13.517700840802913,15.302251891207787,17.32239200287436,19.609222670922968,22.197951281441636,25.12843315425841,28.445785143962375,32.20108024599797,36.45213390178773,41.26439410861079,46.71194903811227,52.878667676447755,59.85949104702991,67.76189389517093,76.70753933829559,86.83415195624421,98.29763815922249,111.27448647797397,125.96448473034971,142.5937958969891,161.41844006140866,182.72823602144376,206.85126325595743,234.15891294197226,265.0716057862269,300.06526470124555,339.67864197737873,384.5216137578392,435.2845695160886,492.7490410932563,557.7997349371907,631.4381527880333,714.7980105104556,809.1626924564707,915.9850100811499,1036.9095690092008,1173.7980889093292,1328.7580659938624,1504.175219423221,1702.7502211507524,1927.5402746900081,2182.006182939198,2470.0656297055034,2796.1534952362003,3165.290134357197,3583.1586684094545,4056.1924809477755,4591.674264260404,5197.848141601231,5884.046591336165,6660.8341270911405,7540.169945960111,8535.592048857829,9662.425667681435,10938.019208165191,12382.011340936813,14016.633352832258,15867.051413382505,17961.754025908867,20332.99062831218,23017.268096126892,26055.912791858576,29495.706813754354,33389.60823950846,37797.56645356835,42787.44511058541,48436.06694468877,54830.396510166145,62068.8790626859,70262.95619408888,79538.78155502846,90039.16308022855,101925.76161830177,115381.57981559623,130613.77957221285,147856.8714469329,167376.3251142129,189472.6564588066,214486.05423174356,242801.6174983236]

const Filter_J0 = [0.001502009920951996,-0.010381698214761684,0.036840860097595164,-0.0899033803922747,0.1708228653683386,-0.2711574965683628,0.37649328091859574,-0.47220778569122657,0.5477821108964709,-0.598235168530351,0.6234579161218533,-0.6265043664825772,0.6119722535117335,-0.5847017386614045,0.5491105578961686,-0.5087886168486084,0.466521063454302,-0.42426478870289325,0.3833953891193336,-0.34472344791723936,0.3087689146689051,-0.2757043936836481,0.24562331000974616,-0.21839207265400126,0.19393194321380827,-0.17198162962644575,0.15242270280410272,-0.13494945901825492,0.11946763189654602,-0.10565870880444576,0.09348235554803391,-0.08261252527983624,0.07307776352617497,-0.06453648128170122,0.057096310621587334,-0.0503848596052138,0.04459955739673305,-0.039316995063444334,0.03483854475429614,-0.030664946477420088,0.027221072240291456,-0.0239015861084346,0.021281566646364738,-0.018612484595521447,0.0166553863034614,-0.014472420076504743,0.013057747649606617,-0.0112262860010603,0.01026681855793083,-0.008673802242165359,0.008110336867084816,-0.006657362814208594,0.0064551189810889585,-0.005052401511783908,0.005198901427050626,-0.0037597430078747883,0.004264054516529264,-0.0026994965688974474,0.003592800056058723,-0.001806129128998902,0.0031436470898586703,-0.0010244200001620755,0.0028888297674091743,-0.00030605070741955765,0.0028125907445050334,0.00039337862066047565,0.0029102041487017727,0.0011170884769911178,0.003187676341595336,0.0019097806429850762,0.003662102032146541,0.002820378358784995,0.004362688559939043,0.00390500249156868,0.005332492348686608,0.005230339939427049,0.006630934845839698,0.006877554068703745,0.008337168998822609,0.00894685647454101,0.010554326601376745,0.011562765604509096,0.013414536265660639,0.014879837399735491,0.01708424133651957,0.01908809265444127,0.021768513715084332,0.024416120601223806,0.027711234465350204,0.031127164234404547,0.03518414022027745,0.03949791311363894,0.04444976215720499,0.04975843326998389,0.05566745711499243,0.06194979043415182,0.06868206490557184,0.07561674395011839,0.08258498356357226,0.08919363526298704,0.09487468694200528,0.09889168186909404,0.1000429465449573,0.09701684432980286,0.08778959614991438,0.07050976759285542,0.04277885348484935,0.0035584532926218175,-0.04721045326487935,-0.10489787743225988,-0.1602095040734828,-0.19459781573132096,-0.1849077459954238,-0.1075416502002519,0.03603772748747661,0.19759013047489976,0.26132313321851336,0.11713996822458939,-0.1875877928130144,-0.3023811499746215,0.04816313568456773,0.36399529664885466,-0.14910233461562913,-0.26373490348543854,0.403626618077187,-0.3140979465010458,0.1817936940513108,-0.09073871804263177,0.04294648754516024,-0.020586135807067835,0.010392667161913182,-0.005611784806872302,0.0032402025511569896,-0.0019858724388273777,0.0012807317326135252,-0.0008625379175606825,0.0006029659078214355,-0.00043548936996943465,0.00032375891570874245,-0.0002469821224005998,0.00019279062274925971,-0.00015357911509105972,0.0001245378784936744,-0.00010255126402954421,8.555848220947627e-5,-7.217092847633435e-5,6.143686328308002e-5,-5.2693349406473615e-5,4.547125514262373e-5,-3.943328415865359e-5,3.433297116479562e-5,-2.9987212220165472e-5,2.625765743561439e-5,-2.3037978256349448e-5,2.024507133101665e-5,-1.781292564438252e-5,1.5688305607849465e-5,-1.382767949211047e-5,1.2195005442258958e-5,-1.0760110818279243e-5,9.497485767095914e-6,-8.385371134393834e-6,7.405061254071319e-6,-6.54036828602285e-6,5.7772102413267325e-6,-5.103293321879673e-6,4.5078641320869195e-6,-3.981511195181711e-6,3.515999321029234e-6,-3.1041249680128287e-6,2.7395854235553632e-6,-2.4168587823583312e-6,2.1310948134441218e-6,-1.8780185742021288e-6,1.6538488252707735e-6,-1.4552318768651144e-6,1.2791890865306748e-6,-1.1230741959148155e-6,9.845363053859617e-7,-8.61485744931908e-7,7.520624269970173e-7,-6.546081033468263e-7,5.676445159602976e-7,-4.8985884555185e-7,4.200967410687925e-7,-3.5736212641225517e-7,3.008221596993967e-7,-2.498151048163967e-7,2.0385823466866512e-7,-1.6265189071584773e-7,1.260741670061131e-7,-9.415841791345086e-8,6.704391121706343e-8,-4.4891090827293947e-8,2.7761325666544702e-8,-1.5480404355710375e-8,7.532730014109875e-9,-3.0524770418657847e-9,9.587785609683078e-10,-2.0575286298055636e-10,2.2414416956474645e-11]

const Filter_J1 = [4.782787133250618e-10,-2.9784175503440788e-9,9.772383277089722e-9,-2.238234099608581e-8,4.044677432947085e-8,-6.173481585455392e-8,8.329391218560819e-8,-1.024945350228408e-7,1.1780779749909977e-7,-1.287006146083485e-7,1.35592434383499e-7,-1.3921010821521872e-7,1.406574572276967e-7,-1.407488190837528e-7,1.4051720878600928e-7,-1.404074668777783e-7,1.4127886061686993e-7,-1.4315595655055356e-7,1.4689283208027915e-7,-1.5210916706348747e-7,1.598980155013874e-7,-1.694091840791126e-7,1.8227089415749844e-7,-1.969585687860337e-7,2.160395242710676e-7,-2.3691320619292838e-7,2.6369843208466607e-7,-2.9202021404039016e-7,3.285244508632466e-7,-3.6589094553627693e-7,4.148750103686303e-7,-4.6327136995173986e-7,5.285269736975052e-7,-5.903498371995408e-7,6.77102116115609e-7,-7.551294280790103e-7,8.70624734664098e-7,-9.678853091813045e-7,1.1222658353725904e-6,-1.2417228058919743e-6,1.449346703689914e-6,-1.593245655920808e-6,1.8747045814274419e-6,-2.0433320340041385e-6,2.4285695351967146e-6,-2.617992645652053e-6,3.151172966139178e-6,-3.349241203288127e-6,4.096418661354996e-6,-4.275827575123915e-6,5.33711975522078e-6,-5.44354237736265e-6,6.97257669216721e-6,-6.9045562968161745e-6,9.139702543697744e-6,-8.714837303363596e-6,1.2029590806160379e-5,-1.0927976968519436e-5,1.59125267194553e-5,-1.3582559661331659e-5,2.1176226828087565e-5,-1.6678205993448338e-5,2.8384979250408712e-5,-2.0132088397797457e-5,3.837204531118811e-5,-2.3702184335455945e-5,5.238530885000794e-5,-2.685437394370126e-5,7.231858155790216e-5,-2.8535361884516687e-5,0.00010108123118106823,-2.6788477540644352e-5,0.00014319185407094621,-1.8108424211338017e-5,0.00020573561552327273,3.6361648565843316e-6,0.0002999126469285961,4.8993332079278846e-5,0.00044354733854670903,0.00013589101811995494,0.0006651582352127376,0.0002945199160862439,0.0010105553806271136,0.0005753396479225405,0.0015535077418254303,0.0010621193133794828,0.0024128970258747457,0.0018929698186109245,0.0037800177772191607,0.0032937343278959356,0.005961217980039154,0.005629593532055224,0.009442252680316808,0.009481022824713792,0.014979159139973408,0.015745093424331037,0.02370896637000014,0.02574059076213687,0.0372327828431175,0.04122500890061429,0.05750710321227736,0.06404464284623569,0.08609179655185725,0.09471713980445745,0.12172497389177185,0.128535970003989,0.15450777327408322,0.1475596409096932,0.15621399202016978,0.11147620703185755,0.07748983135608338,-0.02762826685014771,-0.1019873017831784,-0.2203988997111164,-0.21185762869925318,-0.1605241508315224,0.09164902579868109,0.23792823877700942,0.26075777853738125,-0.15662188259001042,-0.28932081756330175,0.01314851911624769,0.42691302759079564,-0.4000505000648904,0.11513789407450359,0.09374824435871762,-0.16037231301955096,0.15071857939129532,-0.12120369075996129,0.09411065607998234,-0.07374223843458433,0.059038567576124905,-0.04828811752847585,0.04019705429957688,-0.03391978772064108,0.02891824715676397,-0.024845271759013743,0.021470449751150148,-0.01863582802005709,0.01622957936236386,-0.014170085406786529,0.01239608412101189,-0.010860414401084047,0.009525944424535663,-0.008362857744723338,0.0073468029551253195,-0.006457604321096636,0.005678343995599488,-0.004994694916726544,0.004394425810860881,-0.003867026401966086,0.003403418035555667,-0.0029957260668529964,0.0026370977166248776,-0.002321554011737298,0.0020438677474690805,-0.0017994616759226389,0.0015843226896713463,-0.0013949288614414647,0.001228186970888625,-0.0010813786997088304,0.0009521140746075729,-0.0008382910302044814,0.0007380601822098776,-0.0006497940667124724,0.0005720602290123041,-0.0005035976454348357,0.000443296041223001,-0.00039017773206623073,0.0003433816697409889,-0.00030214941633163506,0.00026581280850704716,-0.0002337831047958839,0.00020554143583887405,-0.00018063040107732216,0.00015864667598176302,-0.0001392345123516876,0.00012208003098625932,-0.00010690622166177854,9.346858036256814e-5,-8.155132858123432e-5,7.096417463153107e-5,-6.153959246866666e-5,5.313060911614544e-5,-4.560910598312646e-5,3.886464858418121e-5,-3.2803856352344075e-5,2.735029677587626e-5,-2.24448161508053e-5,1.804607628158424e-5,-1.4130826937491561e-5,1.0693106849359383e-5,-7.741205331453028e-6,5.29105764436983e-6,-3.3552268362550323e-6,1.928295620636745e-6,-9.725371257205876e-7,4.110080763295935e-7,-1.3553176263207053e-7,3.0748587523233524e-8,-3.5668195345476294e-9]

##
#if both r and k domains are logarithmically-spaced we can take advantage of computational efficiencies by reusing computations of
#the integration kernel at the same k values. The spacing is fixed based on predefined values in Filter_base.
filterspacing = Filter_base[2]/Filter_base[1]
kmin = Filter_base[1]/rmax
kmax = Filter_base[end]/rmin

n_k_pts = ceil(Int, (log(kmax) - log(kmin))/log(filterspacing)) + 1
k_vals = exp.(log.(kmin):log(filterspacing):log(kmin)+(n_k_pts-1)*log(filterspacing))


#near field
zgrid = 0:0.01*R:1.5*R

kcut(r) = pi/((filterspacing-1)*r)
lowpass(k) = 1/(1 + (k/kcut(R))^3)

integration_kernel = [lowpass(k) * k * Hz_wavenumber.(k, z) for k in k_vals, z in zgrid]

#radial grid - number of r vals is equal to number of k vals minus base length of the filter + 1.
#the base radius is rmax
flen = length(Filter_base)
r_vals = rmax*exp.(-(0:n_k_pts-flen)*log(filterspacing))

#we're only using Filter_J0
Hz_digital = reduce(vcat, [1/r_vals[ir] * Filter_J0' * integration_kernel[ir:ir+flen-1,:] for ir in 1:length(r_vals)])


Hz_analytic_mat = [Hz_analytic(r,z) for r in r_vals, z in zgrid]

fig, ax = subplots(1,2)
ax[1].invert_yaxis()
ax[2].invert_yaxis()
sca(ax[1])
xlabel("r (m)")
ylabel("z (m)")
title("Digital filter inverse Hz field")
pcolor(r_vals, zgrid, Hz_digital', vmin=-0.06, vmax=0.06)
sca(ax[2])
xlabel("r (m)")
ylabel("z (m)")
title("Analytic Hz field")
pcolor(r_vals, zgrid, Hz_analytic_mat', vmin=-0.06, vmax=0.06)
display(gcf())
close("all")

## plot the far-field comparison for the digital filter
# zgrid = 0:0.1*R:10*R
# integration_kernel = [k * Hz_wavenumber.(k, z) for k in k_vals, z in zgrid]
# Hz_digital = reduce(vcat, [1/r_vals[ir] * Filter_J0' * integration_kernel[ir:ir+flen-1,:] for ir in 1:length(r_vals)])


# Hz_analytic_mat = [Hz_analytic(r,z) for r in r_vals, z in zgrid]

# dipole_approximation = [25/((r^2+z^2)^(3/2)) * (3*z^2/(z^2 + r^2) - 1) for r in r_vals, z in zgrid]

# figure()
# plot(zgrid[end-800:end-400], dipole_approximation[1,end-800:end-400])
# plot(zgrid[end-800:end-400], Hz_digital[1,end-800:end-400])
# # plot(zgrid[end-800:end-400], Hz_analytic_mat[1,end-800:end-400])
# legend(["dipole", "numerical hankel transform"])
# xlabel("z (m)")
# ylabel("Hz field (1/m)")
# display(gcf())
# close("all")



## response function (computes "B-response"/impedance and alpha vals for layered earth)
RealVector = Vector{T} where T <: Real
function responses(σ::RealVector, d::RealVector, κ::RealVector, ωl::Real)
    B = zeros(ComplexF64, length(κ), length(σ));
    α = [sqrt(k^2 - im * ωl * mu_0 * σm)  for k in κ, σm in σ]
    B[:,end] = α[:,end]
    m = size(B,2) - 1
    print(m)
    while m > 0
        thamdm = tanh.(α[:,m]*d[m])
        B[:,m] = α[:,m] .* (B[:,m+1] .+ α[:,m].*thamdm)./(α[:,m] .+ B[:,m+1].*thamdm)
        m -= 1
    end
    return B, α
end

##

#Schelkunoff potential at each z
function phiz(phi0, Bresponse, α, d, zgrid)
    #schelkunoff potential at the top of each layer
    phi_tops = zeros(ComplexF64, size(Bresponse)...)
    phi_tops[:,1] = phi0
    for m=1:size(phi_tops,2)-1
        phi_tops[:,m+1] = phi_tops[:,m] .* (α[:,m] + Bresponse[:,m])./(α[:,m] + Bresponse[:,m+1]) .* exp.(-α[:,m]*d[m])
    end
    #convert depth to z-values for the bottom of each layer
    z_interface = cumsum(d)

    phiz = zeros(ComplexF64, length(phi0), length(zgrid))
    phipz = zeros(ComplexF64, length(phi0), length(zgrid))

    layer_m = 1
    coeff_minus = phi_tops[:,1] .* (1 .+ Bresponse[:,1]./α[:,1])/2
    coeff_plus = 0
    if length(d) > 1
        coeff_plus = phi_tops[:,2] .* (1 .- Bresponse[:,2]./α[:,2])/2
    end
    h_m = 0
    h_mp1 = z_interface[1]

    for (iz, z) = enumerate(zgrid)
        if z > h_mp1
            coeff_minus = phi_tops[:,layer_m] .* (1 .+ Bresponse[:,layer_m]./α[:,layer_m])/2
            coeff_plus = 0
            if length(d) > layer_m
                coeff_plus = phi_tops[:,layer_m+1] .* (1 .- Bresponse[:,layer_m+1]./α[:,layer_m+1])/2
            end
            h_m = h_mp1
            h_mp1 = z_interface[layer_m]
        end
        
        phiz[:,iz] = coeff_minus .* exp.(-α[:,layer_m]*(z-h_m)) +
                        coeff_plus .* exp.(α[:,layer_m]*(z-h_mp1))
        phipz[:,iz] = α[:,layer_m].* (coeff_plus .* exp.(α[:,layer_m]*(z-h_mp1))
                        - coeff_minus .* exp.(-α[:,layer_m]*(z-h_m)))
    end

    phiz, phipz
end

## conductive half-space, at 2 kHz
ωl = 2.0e3 #Hz, typical for Earth's field strength
d = [Inf]
σ = [5]
B_halfspace, alpha_halfspace = responses(σ, d, k_vals, ωl)

#surface potential
phi0_free = (R * besselj1.(k_vals * R) ./ (2 * k_vals.^2))
phi0 = 2* phi0_free .* k_vals ./ (k_vals .+ B_halfspace[:,1]) 

#drop the free-space layer from arguments to the phiz function
#we have already done the coupling into the first earth layer
phizk_loop, _ = phiz(phi0, B_halfspace, alpha_halfspace, d, zgrid)

#
Hzk_halfspace = k_vals.^2 .* phizk_loop

integration_kernel = k_vals .* Hzk_halfspace

Hzr_halfspace = real.(reduce(vcat, [1/r_vals[ir] * Filter_J0' * integration_kernel[ir:ir+flen-1,:] for ir in 1:length(r_vals)]))

##
fig, ax = subplots(1,2)
sca(ax[1])
pcolor(r_vals, zgrid, Hz_digital', vmin = -0.06, vmax = +0.06)
title("Free space Hz field")
xlabel("r (m)")
ylabel("z (m)")
gca().invert_yaxis()
sca(ax[2])
pcolor(r_vals, zgrid, Hzr_halfspace', vmin = -0.06, vmax = +0.06)
t = "Half space with conductivity σ = $(σ[1]) S/m"
title(t)
xlabel("r (m)")
ylabel("z (m)")
gca().invert_yaxis()
display(gcf())

#reflection coefficient
γ = (k_vals .- B_halfspace[:,1])./(k_vals .+ B_halfspace[:,1])

##
#loop separated from half-space boundary
function phi_upper(κ, z, h, σ, ωl) 
    α = sqrt(κ^2 - im*ωl*mu_0*σ)
    γ0 = (κ - α)/(κ + α)
    R*besselj1(κ * R) /(2*κ^2) * (exp(-κ*(z+h)) + γ0*exp(-κ*(h-z)))
end

function phi_lower(κ, z, h, σ, ωl)
    α = sqrt(κ^2 - im*ωl*mu_0*σ)
    γ0 = (κ - α)/(κ + α)
    phi0 = R*besselj1.(κ * R) /(2*κ^2) * (1 + γ0) * exp(-κ*h)
    phi0*exp(-α*z)
end

phi(κ,z,h,σ,ωl) = z > 0 ? phi_lower(κ,z,h,σ,ωl) : phi_upper(κ,z,h,σ,ωl)

#say the loop is 10 m above the half-space
zgrid = -2:0.2:20

phi_loop_above = [phi.(κ, z, 2, 5, 2.0e3) for κ in k_vals, z in zgrid]

integration_kernel = k_vals.^3 .* phi_loop_above

Hzr_loop_above = real.(reduce(vcat, [1/r_vals[ir] * Filter_J0' * integration_kernel[ir:ir+flen-1,:] for ir in 1:length(r_vals)]))

figure()
pcolor(r_vals, zgrid, Hzr_loop_above', vmin = -0.03, vmax = 0.03)
gca().invert_yaxis()
gcf()

##
#mixed j1/j0 filter
Hz_wavenumber_j1(κ,z,r) = I * R * besselj0(κ*r) * exp(- κ * abs(z))/2;

function Hz_offset(r,z)
    retval = 0
    if r < R
        #use j1 filter
        kgrid = Filter_base/R
        integration_kernel = kgrid .* Hz_wavenumber_j1.(kgrid, z, r)
        retval = real(1/R * Filter_J1' * integration_kernel)
    else
        #j0 filter
        kgrid = Filter_base/r
        integration_kernel = kgrid .* Hz_wavenumber.(kgrid, z)
        retval = real(1/r * Filter_J0' * integration_kernel)
    end
    retval
end

rgrid = 0:0.2:1.5*R
zgrid = 0:0.1:5*R

Hzr_mixedfilters = [Hz_offset(r,z) for r in rgrid, z in zgrid]
Hzr_analytic = [Hz_analytic(r,z) for r in rgrid, z in zgrid]

##
fig,ax = subplots(1,2,figsize=(10,5))
sca(ax[1])
pcolor(rgrid, zgrid, Hzr_mixedfilters', vmin=-0.06, vmax=0.06)
gca().invert_yaxis()
title("Mixed J0 and J1 filters")
xlabel("r (m)")
ylabel("z (m)")
sca(ax[2])
pcolor(rgrid, zgrid, Hzr_analytic', vmin=-0.06, vmax=0.06)
gca().invert_yaxis()
title("Analytic result")
xlabel("r (m)")
ylabel("z (m)")



display(gcf())
close("all")

## far field for mixed filters
dipole_approximation = [25/((r^2+z^2)^(3/2)) * (3*z^2/(z^2 + r^2) - 1) for r in rgrid, z in zgrid]
figure()
plot(zgrid, Hzr_mixedfilters[1,:],zgrid,dipole_approximation[1,:], zgrid, Hzr_analytic[1,:])
ylim([0,0.05])
display(gcf())
close("all")

## debug near-surface artifacts
Hz_surface_kernel(κ,r) = R*besselj0(κ * r) /(2) * κ 

dense_k_grid = 0:0.01:100

figure()
plot(dense_k_grid, Hz_surface_kernel.(dense_k_grid, 1))

sampling_k_grid = Filter_base/R

scatter(sampling_k_grid, Hz_surface_kernel.(sampling_k_grid,1))
xlim([0,100])
display(gcf())
close("all")

Hz_depth_kernel(κ,r,z) = R*besselj0(κ*r) / 2 * κ * exp(-κ*z)
figure()
plot(dense_k_grid, Hz_depth_kernel.(dense_k_grid, 1,5))

sampling_k_grid = Filter_base/R

scatter(sampling_k_grid, Hz_depth_kernel.(sampling_k_grid,1,5))
xlim([0,100])
display(gcf())
close("all")

## low-pass cutoff
filterspacing = Filter_base[2]/Filter_base[1]
kcut(r) = pi/(2*(filterspacing-1)*r)

function Hz_offset(r,z)
    retval = 0
    if r < R
        #use j1 filter with low-pass
        kgrid = Filter_base/R
        lowpass_k = 1 ./ (1 .+ (kgrid/kcut(r)).^3)
        integration_kernel = lowpass_k .* kgrid .* Hz_wavenumber_j1.(kgrid, z, r)
        retval = real(1/R * Filter_J1' * integration_kernel)
    else
        #j0 filter
        kgrid = Filter_base/r
        integration_kernel = kgrid .* Hz_wavenumber.(kgrid, z)
        retval = real(1/r * Filter_J0' * integration_kernel)
    end
    retval
end

rgrid = 0:0.2:1.5*R
zgrid = 0:0.1:2*R

Hzr_mixedfilters = [Hz_offset(r,z) for r in rgrid, z in zgrid]
Hzr_analytic = [Hz_analytic(r,z) for r in rgrid, z in zgrid]

##
fig,ax = subplots(1,2,figsize=(10,5))
sca(ax[1])
pcolor(rgrid, zgrid, Hzr_mixedfilters', vmin=-0.06, vmax=0.06)
gca().invert_yaxis()
title("Mixed J0 and J1 filters")
xlabel("r (m)")
ylabel("z (m)")
sca(ax[2])
pcolor(rgrid, zgrid, Hzr_analytic', vmin=-0.06, vmax=0.06)
gca().invert_yaxis()
title("Analytic result")
xlabel("r (m)")
ylabel("z (m)")

display(gcf())
close("all")

## radial field - must use J1 kernel
Hr_wavenumber(κ, z) = I * R * besselj1(κ*R) * exp(- κ * abs(z))/2;

kcut_radial(r) = pi/(2*(filterspacing-1)*r)

function Hr_offset(r,z)
    #use j1 filter with low-pass
    kgrid = Filter_base/r
    lowpass_k = 1 ./ (1 .+ (kgrid/kcut_radial(r)).^3)
    integration_kernel = lowpass_k .* kgrid .* Hr_wavenumber.(kgrid, z)
    real(1/r * Filter_J1' * integration_kernel)
end

Hr_numeric = [Hr_offset(r,z) for r in rgrid, z in zgrid]

fig,ax = subplots(1,2,figsize=(10,5))
sca(ax[1])
pcolor(rgrid, zgrid, Hzr_mixedfilters', vmin=-0.06, vmax=0.06)
gca().invert_yaxis()
title("Hz")
xlabel("r (m)")
ylabel("z (m)")
sca(ax[2])
pcolor(rgrid, zgrid, Hr_numeric', vmin=-0.06, vmax=0.06)
gca().invert_yaxis()
title("Hr")
xlabel("r (m)")
ylabel("z (m)")

display(gcf())

##

function Hr_analytic(r,z)
    αsq = R^2 + r^2 + z^2 - 2*r*R;
    βsq = R^2 + r^2 + z^2 + 2*r*R;
    β = sqrt(βsq);
    ksq = 1 - αsq/βsq;
    C = 1/pi;
    
    C * z / (2*αsq*β*r) * ((R^2 + r^2 + z^2)*E(ksq) - αsq*K(ksq))
end

Hr_analytic_result = [Hr_analytic(r,z) for r in rgrid, z in zgrid]
##
fig,ax = subplots(1,2,figsize=(10,5))
sca(ax[1])
pcolor(rgrid, zgrid, Hr_analytic_result', vmin=-0.06, vmax=0.06)
gca().invert_yaxis()
title("Hr analytic")
xlabel("r (m)")
ylabel("z (m)")
sca(ax[2])
pcolor(rgrid, zgrid, Hr_numeric', vmin=-0.06, vmax=0.06)
gca().invert_yaxis()
title("Hr numeric")
xlabel("r (m)")
ylabel("z (m)")

display(gcf())
##

# find components of mag. field perpendicular
# to earth's field, given inclination φ
# and azimuthal coordinate θ
function perpendicular_field(Bz, Br, φ, θ)
    #y is defined as magnetic west.
    #x is defined as
    #cross product of y and the unit earth field vector
    #(perpendicular to the earth field in the north-down plane)
    Bperp_x = Bz * cos(φ) + Br * cos(θ) * sin(φ)
    Bperp_y = Br * sin(θ) 
    [Bperp_x; Bperp_y]
end

# find co-rotating and counter-rotating field parameters
# need |B+|, |B-| and the phase lag ζt

function co_counter_field(Bz, Br, φ, θ)
    Bperp = perpendicular_field(Bz, Br, φ, θ)
    Bstar = conj(Bperp)
    BdotB = transpose(Bperp) * Bperp
    BdotBstar =  transpose(Bperp) * Bstar
    BcrossBstar =  real(im * (Bperp[1] * Bstar[2] - Bperp[2] * Bstar[1]))

    exp_i_ζ = sqrt(BdotB/abs(BdotB))
    αt = sqrt((BdotBstar + abs(BdotB))/2)
    βt = sqrt((BdotBstar - abs(BdotB))/2)
    
    if BcrossBstar < 0
        βt *= -1
    end
    
    Bplus = αt - βt
    Bminus = αt + βt

    [Bplus; Bminus; exp_i_ζ]
end

##

H_params_loop = co_counter_field.(Hzr_mixedfilters, Hr_numeric, π/3, 0)

H_co = first.(H_params_loop)
ζ = last.(H_params_loop)

fig,ax = subplots(1,2)
sca(ax[1])
pcolor(rgrid, zgrid, H_co',vmin=-0.1, vmax=0.1)
gca().invert_yaxis()
sca(ax[2])
pcolor(rgrid, zgrid, ζ', vmin=-1, vmax=1)
gca().invert_yaxis()
display(gcf())
