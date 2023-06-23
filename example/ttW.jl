using HS3
using BenchmarkTools
dict = file_to_dict("/ceph/groups/e4/users/cburgard/public/forRobin/ttW/models/ttW_diff-asym_DEtall_SS.json")
specs = HS3.generate_specs(dict)
@benchmark HS3.generate_specs(dict)
ll_spec = specs.likelihoods[1]
functional_specs  = merge(specs.distributions, specs.functions)
data_specs = specs.data
likelihood = HS3.make_likelihood(ll_spec, functional_specs, data_specs)
analysis = HS3.make_analyses(specs.analyses[1], specs)
val1 = rand(analysis.prior)
params = unique!(reduce(vcat, likelihood.free_parameters))
VA = (;) 
using Distributions, Random
rand(1)
for x in keys(vars)
    VA = merge(VA, (x => Uniform(0, 5),))
end
VA = merge( VA, val1)
using DensityInterface
logdensityof(analysis.likelihood, vars)
using BAT, ValueShapes, Distributions
posterior = PosteriorMeasure(analysis.likelihood, merge(analysis.prior, VA))
bat_findmode(posterior).result




vars = (Lumi = 4.267029718294233, Norm_VV_HF = 0.1493903736944362, Norm_fake_ExtConv = 3.248949912531856, Norm_fake_HF_el_TandM = 3.3441864715866174, Norm_fake_HF_mu_TandM = 1.5182253816602789, Norm_fake_IntConv = 3.766394319071051, Norm_ttZ_LF = 4.349161577379562, alpha_ATLAS_EG_RES = 2.02035103568273, alpha_ATLAS_EG_SCALE = 3.3991596767996444, alpha_ATLAS_EL_ID = 3.355536440363518, alpha_ATLAS_EL_Isol = 0.021176984362652727, alpha_ATLAS_EL_Reco = 1.6327189966495026, alpha_ATLAS_FTAG_B0 = 2.8677316414777936, alpha_ATLAS_FTAG_B1 = 3.270016669331438, alpha_ATLAS_FTAG_B10 = 3.604022623079681, alpha_ATLAS_FTAG_B11 = 1.3965600122275332, alpha_ATLAS_FTAG_B12 = 1.9978137923630415, alpha_ATLAS_FTAG_B13 = 2.8247368830914805, alpha_ATLAS_FTAG_B14 = 1.3522463910116356, alpha_ATLAS_FTAG_B15 = 0.7074752758962475, alpha_ATLAS_FTAG_B16 = 0.23435908787592663, alpha_ATLAS_FTAG_B17 = 3.0816241388544725, alpha_ATLAS_FTAG_B18 = 1.6437368742481042, alpha_ATLAS_FTAG_B19 = 1.2382757199512708, alpha_ATLAS_FTAG_B2 = 3.1244502411448796, alpha_ATLAS_FTAG_B3 = 4.970828686444538, alpha_ATLAS_FTAG_B4 = 4.599219105581895, alpha_ATLAS_FTAG_B5 = 2.9263831677971424, alpha_ATLAS_FTAG_B6 = 3.7257009841908553, alpha_ATLAS_FTAG_B7 = 4.017904054124441, alpha_ATLAS_FTAG_B8 = 0.10449650914870108, alpha_ATLAS_FTAG_B9 = 3.4786877038803206, alpha_ATLAS_FTAG_C0 = 4.280603441019295, alpha_ATLAS_FTAG_C1 = 0.8984793651813332, alpha_ATLAS_FTAG_C10 = 2.0906272507379517, alpha_ATLAS_FTAG_C11 = 4.102440552170017, alpha_ATLAS_FTAG_C12 = 1.8574305838441412, alpha_ATLAS_FTAG_C13 = 0.3341528274193218, alpha_ATLAS_FTAG_C14 = 4.298560917046819, alpha_ATLAS_FTAG_C15 = 1.387463566550799, alpha_ATLAS_FTAG_C16 = 0.05249410627099269, alpha_ATLAS_FTAG_C17 = 4.390987418742104, alpha_ATLAS_FTAG_C18 = 4.6235678688521675, alpha_ATLAS_FTAG_C19 = 1.507216927572957, alpha_ATLAS_FTAG_C2 = 3.9111142408442148, alpha_ATLAS_FTAG_C3 = 2.402707063505234, alpha_ATLAS_FTAG_C4 = 1.4596824151398788, alpha_ATLAS_FTAG_C5 = 2.1453100936267977, alpha_ATLAS_FTAG_C6 = 4.579044123595924, alpha_ATLAS_FTAG_C7 = 4.7931626609806, alpha_ATLAS_FTAG_C8 = 1.7096267659078113, alpha_ATLAS_FTAG_C9 = 3.0598023538492085, alpha_ATLAS_FTAG_L0 = 2.455431428216102, alpha_ATLAS_FTAG_L1 = 1.7911168388772158, alpha_ATLAS_FTAG_L10 = 4.201349061130297, alpha_ATLAS_FTAG_L11 = 2.9929235615905223, alpha_ATLAS_FTAG_L12 = 1.7589690124464914, alpha_ATLAS_FTAG_L13 = 0.9129073053857837, alpha_ATLAS_FTAG_L14 = 2.3517652953197787, alpha_ATLAS_FTAG_L15 = 2.6255264689808784, alpha_ATLAS_FTAG_L16 = 1.5172449733128015, alpha_ATLAS_FTAG_L17 = 3.8050974284397916, alpha_ATLAS_FTAG_L18 = 2.148793535956451, alpha_ATLAS_FTAG_L19 = 3.5226166361665476, alpha_ATLAS_FTAG_L2 = 2.5946738288762914, alpha_ATLAS_FTAG_L3 = 2.5063026266977757, alpha_ATLAS_FTAG_L4 = 4.291045723176507, alpha_ATLAS_FTAG_L5 = 2.4232748853169572, alpha_ATLAS_FTAG_L6 = 2.8487486841849634, alpha_ATLAS_FTAG_L7 = 1.460037922320199, alpha_ATLAS_FTAG_L8 = 0.9366442013500123, alpha_ATLAS_FTAG_L9 = 0.8051416675862717, alpha_ATLAS_JER_DataVsMC = 0.43390219866549395, alpha_ATLAS_JER_Eff1 = 3.33384671762535, alpha_ATLAS_JER_Eff10 = 3.4043671263340443, alpha_ATLAS_JER_Eff11 = 3.8019925954885636, alpha_ATLAS_JER_Eff12 = 4.136809643177686, alpha_ATLAS_JER_Eff2 = 1.7491676126298386, alpha_ATLAS_JER_Eff3 = 0.9891191763714282, alpha_ATLAS_JER_Eff4 = 1.9491751614984119, alpha_ATLAS_JER_Eff5 = 2.473930051391876, alpha_ATLAS_JER_Eff6 = 1.394953699846537, alpha_ATLAS_JER_Eff7 = 0.8574999477806953, alpha_ATLAS_JER_Eff8 = 2.052773799675887, alpha_ATLAS_JER_Eff9 = 0.24158681830176243, alpha_ATLAS_JES_BJES = 1.9392659000847516, alpha_ATLAS_JES_EtaInter_Model = 4.28163500896678, alpha_ATLAS_JES_EtaInter_NonClosureHighE = 4.162136383539165, alpha_ATLAS_JES_EtaInter_NonClosureNegEta = 2.3866424590061723, alpha_ATLAS_JES_EtaInter_NonClosurePosEta = 2.279205002058777, alpha_ATLAS_JES_EtaInter_NonClosure_2018data = 0.917344026979486, alpha_ATLAS_JES_EtaInter_Stat = 0.06843199169621092, alpha_ATLAS_JES_Flavor_Comp = 4.365437428381378, alpha_ATLAS_JES_Flavor_Resp = 0.33701946186978143, alpha_ATLAS_JES_NP_Det1 = 0.9007848192071337, alpha_ATLAS_JES_NP_Det2 = 3.8902774849477764, alpha_ATLAS_JES_NP_Mix1 = 2.134101349999501, alpha_ATLAS_JES_NP_Mix2 = 2.952811756642836, alpha_ATLAS_JES_NP_Mix3 = 4.982856942023799, alpha_ATLAS_JES_NP_Mod1 = 0.06569306142201264, alpha_ATLAS_JES_NP_Mod2 = 0.8865091507904654, alpha_ATLAS_JES_NP_Mod3 = 0.9423303365040026, alpha_ATLAS_JES_NP_Mod4 = 3.2874963743893204, alpha_ATLAS_JES_NP_Stat1 = 1.0059991905091947, alpha_ATLAS_JES_NP_Stat2 = 3.0822368858075153, alpha_ATLAS_JES_NP_Stat3 = 0.6894913397793415, alpha_ATLAS_JES_NP_Stat4 = 3.630540903488153, alpha_ATLAS_JES_NP_Stat5 = 3.5255069703639186, alpha_ATLAS_JES_NP_Stat6 = 3.051580330665439, alpha_ATLAS_JES_PU_OffsetMu = 0.26855975232764673, alpha_ATLAS_JES_PU_OffsetNPV = 3.0834730912208332, alpha_ATLAS_JES_PU_PtTerm = 2.7431349322029575, alpha_ATLAS_JES_PU_Rho = 3.2420836699108984, alpha_ATLAS_JES_PunchThrough = 1.364029890600915, alpha_ATLAS_JES_SinglePart = 2.5999498198629474, alpha_ATLAS_JVT = 1.1062623140649435, alpha_ATLAS_MET_Para = 0.11752869713334824, alpha_ATLAS_MET_Perp = 1.951218253276815, alpha_ATLAS_MET_Scale = 0.041234030693137864, alpha_ATLAS_MU_ID = 2.4197371413209146, alpha_ATLAS_MU_ID_STAT = 3.8015202399187724, alpha_ATLAS_MU_ID_STAT_LOWPT = 1.4426188498709451, alpha_ATLAS_MU_ID_SYST = 4.751211982623862, alpha_ATLAS_MU_ID_SYST_LOWPT = 3.8683343439386406, alpha_ATLAS_MU_Isol_STAT = 3.2812732881920827, alpha_ATLAS_MU_Isol_SYST = 3.9261407584924592, alpha_ATLAS_MU_MS = 1.1900056726819588, alpha_ATLAS_MU_RESBIAS = 2.640878974412793, alpha_ATLAS_MU_SCALE = 1.1105040718018695, alpha_ATLAS_MU_TTVA_STAT = 0.49910824266623033, alpha_ATLAS_MU_TTVA_SYST = 0.9185050857480846, alpha_ATLAS_PRW_DATASF = 3.840260597702138, alpha_ATLAS_TRIG_EL = 3.8686774746386514, alpha_ATLAS_TRIG_MU_STAT = 3.3668143121842937, alpha_ATLAS_TRIG_MU_SYST = 3.1807988662266435, alpha_ATLAS_lumi = 3.0729963043390667, alpha_ConvIntExtrap = 1.6038144288788927, alpha_ConvMatExtrap = 0.060471081473668, alpha_FakesElLead_PLIV_PtFrac = 4.970843976211812, alpha_FakesElSubLead_PLIV_PtFrac = 0.2296319537682499, alpha_FakesEl_NBjetCorr = 2.3024224705441765, alpha_FakesEl_Tight_MtoTextrap = 1.4021598394105517, alpha_FakesMuLead_PLIV_RelCaloCluster = 0.6560973105394436, alpha_FakesMuLead_PLIV_SVLongSignif = 4.305652057372955, alpha_FakesMuSubLead_PLIV_RelCaloCluster = 2.8836775503334446, alpha_FakesMuSubLead_PLIV_SVLongSignif = 1.4710932198567679, alpha_FakesMu_NBjetCorr = 1.5813336586194178, alpha_FakesMu_Tight_MtoTextrap = 3.627774814379713, alpha_PLIV_El_ID_VeryTight = 2.400864775087011, alpha_PLIV_El_Iso_VeryTight = 4.729765266990068, alpha_PLIV_El_JetModeling_Tight = 0.6444065268696766, alpha_PLIV_El_JetModeling_VeryTight = 4.791225040174218, alpha_PLIV_El_MllWindow_VeryTight = 1.7241191264189586, alpha_PLIV_El_Pileup_Tight = 3.5938546813337013, alpha_PLIV_El_Pileup_VeryTight = 0.7235929603812615, alpha_PLIV_El_Stat_Tight = 4.762943237305686, alpha_PLIV_El_Stat_VeryTight = 2.933412194675823, alpha_PLIV_El_TemplateCut_VeryTight = 4.948731365703891, alpha_PLIV_El_Trigger_Tight = 1.7107401039744752, alpha_PLIV_El_Trigger_VeryTight = 0.5449573611228025, alpha_PLIV_Mu_BkgFraction_VeryTight = 4.638553258954151, alpha_PLIV_Mu_DRMuJet_VeryTight = 1.1641858828688818, alpha_PLIV_Mu_JetModeling_Tight = 2.5065824364512586, alpha_PLIV_Mu_JetModeling_VeryTight = 1.6315243262845356, alpha_PLIV_Mu_Luminosity_VeryTight = 4.025286830983386, alpha_PLIV_Mu_MCxSec_VeryTight = 1.3289013458242895, alpha_PLIV_Mu_MllWindow_Tight = 1.3527268963559378, alpha_PLIV_Mu_MllWindow_VeryTight = 1.9111651594511652, alpha_PLIV_Mu_ProbeQuality_Tight = 1.2836891436502378, alpha_PLIV_Mu_ProbeQuality_VeryTight = 3.754237577617268, alpha_PLIV_Mu_QCDTemplate_VeryTight = 3.7151999942456233, alpha_PLIV_Mu_Stat_Tight = 0.9784491617742797, alpha_PLIV_Mu_Stat_VeryTight = 2.299223493711199, alpha_PLIV_Mu_SuppressionScale_VeryTight = 2.668842627027597, alpha_QMisIDXsec = 4.081765043437891, alpha_QMisID_MM_SYST = 1.7764300590951223, alpha_QMisID_MT_SYST = 0.5508403967425256, alpha_QMisID_TM_SYST = 3.0219520890418803, alpha_VHXsec = 1.301459532774253, alpha_VVVXsec = 2.987858208944391, alpha_VV_varRF = 4.4716324048198395, alpha_VVlightXsec = 2.2906326410384175, alpha_VVnJets = 1.0469010674426669, alpha_WtZXsec = 0.02309965295211485, alpha_fourTopXsec = 0.46102482174971093, alpha_tZXsec = 0.031040630550480012, alpha_threeTopXsec = 2.4766205837287973, alpha_ttHXsec = 2.8380093097726493, alpha_ttH_Gen_Acceptance = 4.285062616855545, alpha_ttH_PS_Acceptance = 1.9670369443928544, alpha_ttH_varRF = 4.487867753539196, alpha_ttWWXsec = 0.016191494251662872, alpha_ttW_EW_fraction = 4.211743008822156, alpha_ttW_EW_varRF = 3.0064433598167737, alpha_ttW_ME = 0.8885352379061271, alpha_ttW_ME_EW = 1.2629256508896987, alpha_ttW_PDFaS = 0.20892573905224734, alpha_ttW_PDFalternate = 1.4700704155406863, alpha_ttW_PS_EW = 1.3133043502921264, alpha_ttW_PS_QCD = 3.6905725953716746, alpha_ttW_varRF = 1.6482184721759154, alpha_ttZHFXsec = 2.8477314602901256, alpha_ttZ_Shower = 0.9200863031331925, alpha_ttZ_Var3c = 3.094854865024121, alpha_ttZ_varRF = 0.597148775311287, alpha_ttbb_XS = 1.51978323549395, alpha_ttcc_XS = 2.042459481854017, nominalLumi = 1.0002354827286402, ttW_Bin_001_Ac = 1.423871497057771, ttW_Bin_001_muInc = 2.22465459924006, ttW_Bin_002_Ac = 2.1866601486406743, ttW_Bin_002_muInc = 3.260193903642498, ttW_Bin_003_Ac = 1.8236257891259045, ttW_Bin_003_muInc = 0.22540190903798124, ttW_Bin_004_Ac = 2.985778540896263, ttW_Bin_004_muInc = 1.6555965527279959, ttW_Bin_005_Ac = -0.2920555360108111, ttW_Bin_005_muInc = 3.552473311840479, ttW_Bin_006_Ac = 2.6400244625436824, ttW_Bin_006_muInc = 2.774701233027782, ttW_Bin_007_Ac = -0.3768471020576964, ttW_Bin_007_muInc = 0.6326478772725345, ttW_Bin_008_Ac = 1.6208142503920095, ttW_Bin_008_muInc = 0.4614989047561798, ttW_Bin_009_Ac = 1.558286235976788, ttW_Bin_009_muInc = 0.7331812491119388, ttW_Bin_010_Ac = -0.027666156223725053, ttW_Bin_010_muInc = 1.5286704081315619, ttW_Bin_011_Ac = -0.21913372066832176, ttW_Bin_011_muInc = 1.761978943201553, ttW_Bin_012_Ac = -1.5927190763565586, ttW_Bin_012_muInc = 0.22275457974255422, ttW_Bin_013_Ac = 2.60680100319505, ttW_Bin_013_muInc = 0.5732610149556504, staterror_bin1 = 3.308805427637473, staterror_bin2 = 1.873677759473629, staterror_bin3 = 3.537547568068721, staterror_bin30 = 0.01790369418092533, staterror_bin14 = 4.801553630623356, shape_stat_bin1 = 3.8540313766480807, staterror_bin41 = 4.719570577254592, shape_stat_bin2 = 4.621485709389548, staterror_bin9 = 0.3515306208099395, staterror_bin22 = 0.5420515093154226, shape_stat_bin6 = 1.4781028422705529, shape_stat_bin44 = 1.5526198534245406, staterror_bin5 = 3.3513509784060127, shape_stat_bin13 = 2.4719561951670155, shape_stat_bin47 = 2.6618980188909678, shape_stat_bin10 = 2.7723133226270686, shape_stat_bin23 = 4.496750153798186, staterror_bin18 = 0.9526643069759486, shape_stat_bin18 = 3.6452130282469577, staterror_bin12 = 0.7217210323981245, staterror_bin21 = 0.6152321968665566, staterror_bin6 = 3.1184616583868516, staterror_bin46 = 1.3439375336453052, staterror_bin45 = 0.3195320134999824, shape_stat_bin48 = 2.2593501416825843, shape_stat_bin4 = 4.658286035750209, shape_stat_bin31 = 3.5874201482259145, staterror_bin32 = 1.4081927034733945, shape_stat_bin15 = 1.7207050731516202, staterror_bin27 = 3.6571888666402064, shape_stat_bin41 = 4.619320112193787, staterror_bin24 = 4.58680806159905, shape_stat_bin26 = 0.6838596568750264, staterror_bin37 = 4.038628436663305, shape_stat_bin9 = 3.306607080574696, staterror_bin4 = 2.895937536751199, staterror_bin25 = 0.18721254330615503, shape_stat_bin35 = 1.7077397832398664, shape_stat_bin43 = 2.836888747890847, shape_stat_bin3 = 2.183351064298128, shape_stat_bin38 = 2.0416722464183583, shape_stat_bin39 = 2.4114641883062538, staterror_bin38 = 3.955058214703085, shape_stat_bin16 = 1.389484366989466, staterror_bin35 = 0.8444606437341586, shape_stat_bin7 = 2.2878126412446327, staterror_bin23 = 4.859602243060119, staterror_bin26 = 0.5939730395197121, shape_stat_bin40 = 1.5504360658397782, shape_stat_bin20 = 0.566164884953962, staterror_bin42 = 0.19493955209962038, staterror_bin48 = 1.3083461393958404, staterror_bin11 = 2.9104113284217132, shape_stat_bin14 = 0.9165720706624734, shape_stat_bin36 = 1.2431753241051908, staterror_bin36 = 3.713496235180622, staterror_bin34 = 2.103180482276733, shape_stat_bin22 = 2.4948337320610148, staterror_bin40 = 1.3687317270126191, staterror_bin33 = 3.9238418111555893, shape_stat_bin46 = 0.06301981061429941, staterror_bin16 = 1.5748897069493502, staterror_bin8 = 0.41145727156903755, staterror_bin31 = 4.783062961611475, shape_stat_bin8 = 2.4003685800209347, shape_stat_bin21 = 2.701461273562679, shape_stat_bin12 = 3.3230768681638856, shape_stat_bin42 = 3.6812028820823697, staterror_bin28 = 0.18705306706494815, staterror_bin50 = 3.3070824351269863, staterror_bin10 = 3.8834442939429312, shape_stat_bin33 = 2.4089102386786094, staterror_bin19 = 3.391606146001267, shape_stat_bin49 = 1.2225331008253906, shape_stat_bin34 = 1.8477745239119352, shape_stat_bin29 = 3.8719045378541392, staterror_bin13 = 4.887748169008098, shape_stat_bin30 = 2.692023881820337, staterror_bin15 = 1.173567396852957, staterror_bin29 = 2.551742706933016, shape_stat_bin11 = 3.5924432350175044, shape_stat_bin17 = 1.545973871262736, staterror_bin17 = 0.3393702268333365, staterror_bin39 = 2.5907286676293904, staterror_bin44 = 0.9307965641171765, shape_stat_bin28 = 2.605589084966212, shape_stat_bin5 = 2.8809523211978334, staterror_bin7 = 3.3703393681540734, staterror_bin20 = 2.2919728191121878, shape_stat_bin27 = 2.9407040751202587, staterror_bin43 = 4.072766377810427, shape_stat_bin24 = 4.99049000374366, shape_stat_bin50 = 4.350483356967466, staterror_bin47 = 3.3490506844670977, shape_stat_bin37 = 3.82335005101768, staterror_bin49 = 4.074942330924208, shape_stat_bin45 = 0.03163247721187448, shape_stat_bin19 = 1.96563108222552, shape_stat_bin32 = 1.9615828548381637, shape_stat_bin25 = 0.34791013069313115)
function has_staterror_bins(list)
    pattern = r"_bin\d+$"  # Regular expression pattern to match "_bin" followed by a number
    matching_elements = []
    for element in list
        match_result = match(pattern, string(element))
        if match_result !== nothing
            name = replace(string(element), pattern => s"")  # Extract the captured group representing the name
            matching_elements = push!(matching_elements, name)
        end
    end
    return matching_elements
end

# Example usage
my_list = [
    :staterror_bin1,
    :some_other_element,
    :staterror_bin2,
    :another_element
]

has_staterror_bins(my_list)
#    println("The list contains elements that match the pattern '_binX', where X is a number.")
#else
#    println("The list does not contain elements that match the pattern '_binX', where X is a number.")
#end






@benchmark HS3.generate_specs(dict)
k = ["a", "b"]
eachindex(k)
specs.domains[1].staterror
dist = specs.distributions.regularization_ttW
x = unique!(reduce(vcat, likelihood.free_parameters))
y = keys(specs.domains[1])
println(x)
c = (;)
for e in x 
    c = merge(c, [e => .1],)
end
length(c)  
non_matching_elements = setdiff(y, x)
    
if length(non_matching_elements) > 0
       println("Non-matching elements found:")
       for element in non_matching_elements
           println(element)
       end
    end
@benchmark a = HS3.make_functional(dist, HS3.topological_sort(merge(specs[:distributions], specs[:functions])))
using DensityInterface
logdensityof(likelihood.likelihood, (Lumi =1,))

a[1]
d = a((mu_prime_2 = 2, mu_prime_1 =4 ))
d[1]
d[2]
using Distributions
k = MvNormal([0., 0], reduce(hcat, [[1, 0.], [0,1]]))
h = Normal()

Distributions.logpdf(h, 1)
#@btime HS3._val_content(Val(:a))
##@btime Symbol("a")
#example_dict = file_to_dict("./example/a.json") ####~ 50ms / 36ms Rastaban
@benchmark example_dict = file_to_dict("./example/a.json") 
#example_dict[:distributions][2].samples[1]
oldstd = stdout
redirect_stdout(Base.DevNull())
Val{:mu_prime_2} in [Val{:mu_prime_2}]
#redirect_stdout(oldstd)
#BenchmarkTools.DEFAULT_PARAMETERS.seconds = 5
#BenchmarkTools.DEFAULT_PARAMETERS.samples = 1000
#BenchmarkTools.DEFAULT_PARAMETERS.evals = 1
#HS3.generate_specs(example_dict) # ~9s
specs = HS3.generate_specs(example_dict) # ~9s, Rastaban: max. 0.5s
@benchmark HS3.generate_specs(example_dict)
fieldnames((a=2, b=3))
specs.likelihoods[2]
@benchmark HS3.make_likelihood(specs.likelihoods[2], merge(specs.distributions, specs.functions), specs.data)
a

typeof(valus)
    eh = ((-0.006452301828157481, 0.021704004705953794, -0.013718764463578559, -0.00015046859696732362, -0.004281092937398201, 0.02663982668771725, -0.0042507829075222325, -0.021982735947919663, -0.012532373928357154, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))
l = [x for x in eh]
typeof(l)
length(Tuple{Val{:mu_prime_2}, Val{:mu_prime_3}, Val{:mu_prime_4}, Val{:mu_prime_5}, Val{:mu_prime_6}, Val{:mu_prime_9}, Val{:mu_prime_10}, Val{:mu_prime_11}, Val{:mu_prime_12}, Val{:mu_prime_15}, Val{:mu_prime_16}, Val{:mu_prime_17}, Val{:mu_prime_18}, Val{:mu_prime_19}, Val{:mu_prime_22}, Val{:mu_prime_23}, Val{:mu_prime_24}, Val{:mu_prime_25}})
#specs.distributions
ana_spec = specs.analyses[1]
a = HS3.make_analyses(ana_spec, specs)
f = unique!(reduce(vcat, a.free_parameters))
using ValueShapes
using Random
mmn = (; (i => Uniform(0, 1) for i in f)...)
prior = NamedTupleDist(mmn)
valus = rand(prior)
logdensityof(a.likelihood, valus)
k = NTuple{}()
typeof(k)
map(x -> (x),  a.free_parameters[i] for i in 1:11)
e = [merge(k, a.free_parameters[i]) for i in 1:11] 
println(keys(a.prior))
println((a.parameters_of_interest))
#@benchmark HS3.make_analyses(ana_spec, specs)
println(typeof(a.likelihood))
b = a.likelihood
using BAT
using Random
using Distributions, ValueShapes
b = rand(a.prior)
b = merge(b, (staterror= 0.,))

c = NamedTupleDist( a= Uniform(0,1), b= Uniform(2, 3))
keys(b)
posterior = PosteriorMeasure(a.likelihood, a.prior)
samples = bat_sample(posterior, MCMCSampling(mcalg = MetropolisHastings(), nsteps = 10^5, nchains = 4)).result
bat_findmode(a.likelihood)
#z
#b = merge(a.distributions, a.functions)
#HS3.topological_sort(b)
#
#####Distributions einlesen 
#data_spec = AbstractDistSpec(example_dict.distributions.data_dist)
#v = (param_lambda = 2,)
#model = make_markov_kernel(data_spec)
#model(v)
#
#prior_spec = AbstractDistSpec(example_dict.distributions.prior_dist)
#prior_model = make_markov_kernel(prior_spec)
#p = prior_model(v)
#
#example_dict[:priors]
#subst_prior(example_dict[:priors], (prior_dist = p,))