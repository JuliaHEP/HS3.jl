using HS3
using BenchmarkTools
using DensityInterface

dict = file_to_dict("/net/e4-nfs-home.e4.physik.tu-dortmund.de/home/rpelkner/Documents/ttW_diff-asym_DPhill_SS.json")

@elapsed specs = HS3.generate_specs(dict)
analysis = HS3.make_analyses(specs.analyses.simPdf_obsData, specs)
using Random
analysis.likelihood
#logdensityof(analysis.likelihood, rand(NamedTupleDist(a)))
using NamedTupleTools
a = delete(a, :staterror_bin52)
logdensityof(analysis.likelihood, (a))
a = merge(h, blinded)
a = merge(a, (Lumi = 1.,))
analysis.likelihood
specs.distributions.model_CR_HFel_2lSS_MM_lepPt1_2jincl_1b.samples.FakeMuLead_Medium.modifiers.FakesMuLead_PLIV_SVLongSignif
@elapsed a = HS3.generate_domains_specs(dict[:domains])
@elapsed data = HS3.generate_data_specs(dict[:data])
@elapsed functions = HS3.generate_functions_specs(dict[:functions])
@elapsed ll = HS3.generate_likelihoods_specs(dict[:likelihoods])
@elapsed ll = HS3.generate_parameter_points_specs(dict[:parameter_points])
@elapsed ll = HS3.generate_analyses_specs(dict[:analyses])
@elapsed dists = HS3.generate_distributions_specs(dict[:distributions])
specs
specs.domains
l = [0.8632603441351486, 0.8900055786751355, 1.008779038081359, 1.2194997057067491, 0.9879044444100515, 1.012340960385593, 1.1676631638668744, 1.1537579419773842, 0.1407917501570567, 0.1455940436231445, 0.16091393579448052, 0.16177764661775934, 0.1837543018071854, 0.22861811759424167, 0.2405370122999833, 0.2736034092301505, 0.2594943949405995, 0.28258377552727654, 0.3842129124904289, 0.34713984985073565, 0.12290518896700167, 0.1571347314745957, 0.12311814083986596, 0.15596910475673892, 0.18466848023244886, 0.22011490206555062, 0.7748175122351363, 0.8644722471012506, 1.0732560548182855, 1.1156558285520908, 0.9464930766517953, 1.0144806818096435, 1.2393677133462961, 1.2094150708096663, 0.1520844075477267, 0.1789500785186177, 0.12300436521035367, 0.17317351421005342, 0.20929563812500152, 0.27136914158637176, 0.245017299540409, 0.2380796370147264, 0.2634532809094603, 0.283664733689202, 0.3316157580172222, 0.35832547829587025, 0.12511051232910475, 0.13567885916496028, 0.14355801232723608, 0.15010567030372488, 0.13765810239749346, 0.17105715414675704]
l =  [0.0, 0.0, 0.0, 0.5342207869629692, 1.1584410547939474, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0021464802742894385, 3.5981806539077272, 0.0, 0.0, 0.0, 0.0, 0.004260597542960938, 6.0172319206787845, 0.0, 0.0, 0.0, 0.0, 0.0, 3.080092189795308, 0.002643388517055127, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
length(l)
@benchmark HS3.generate_specs(dict)
ll_spec = specs.likelihoods[1]
functional_specs  = merge(specs.distributions, specs.functions)
data_specs = specs.data
likelihood = HS3.make_likelihood(ll_spec, functional_specs, data_specs)

analysis.likelihood

length(a)
a = NamedTuple(a)
a = (;)
a
a = merge(a, (shape_stat_bin28 = Normal(1, 0.2), shape_stat_bin27 = Normal(1, 0.2), shape_stat_bin26 = Normal(1, 0.2), shape_stat_bin7 = Normal(1, 0.02), shape_stat_bin6 = Normal(1, 0.02), shape_stat_bin5 = Normal(1, 0.02), shape_stat_bin4 = Normal(1, 0.02), shape_stat_bin3 = Normal(1, 0.02), shape_stat_bin2 = Normal(1, 0.02), shape_stat_bin1 = Normal(1, 0.02), Lumi = Normal(1, 0.001),))
3768.215940804496 -2430.193465218386
9338.701 - 1338.022
8000 / 2430 
a = merge(a, (Lumi = 1,ttW_Bin_014_muInc = 1.90282))
a = merge(a, (shape_stat_bin52 = 1.0, staterror_bin52 = 1.0, shape_stat_bin51 = 1.0, staterror_bin51 = 1.0))
using DensityInterface
using BAT
using ValueShapes
using Optim
import ForwardDiff
using AutoDiffOperators
2+2
posterior = PosteriorMeasure(analysis.likelihood, NamedTupleDist(a))
set_batcontext(ad = ADModule(:ForwardDiff), autodiff = :forward)
#bat_mode_result = bat_findmode(posterior, OptimAlg(optalg = Optim.ParticleSwarm()))
bat_mode_result2 = bat_findmode(posterior, OptimAlg(optalg = Optim.NelderMead()))
bat_mode = bat_mode_result2.result
sample = bat_sample(posterior)
using Plots
plot(sample.samples)

keys(analysis)
 9287.508 / 3754.843689889635 
j = analysis.parameter_domain


















































val1 = NamedTuple(analysis.prior)
ana

VA = (;) 
using Distributions, Random
rand(1)
pattern = r"ATLAS"  
nt = (;)
for (k,v) in zip(keys(val1), val1)
    if match(pattern, string(k)) !== nothing
        val1=merge(val1, (k => Uniform(-5, 5),))
    end
end
val1
for x in keys(vars)
    VA = merge(VA, (x => Uniform(0, 5),))
end
a = (HS3.make_parameterpoints(specs.parameter_points.default_values))
dump(likelihood)[1][1]
MvNormal([1, 0.1], [0.5, 0.5])
VA = merge( VA, val1)
VA = merge(VA, (Lumi = 1.,))

using DensityInterface
logdensityof(analysis.likelihood, vars)
using BAT, ValueShapes, Distributions
posterior = PosteriorMeasure(analysis.likelihood, merge(analysis.prior, VA))
BAT.bat_findmode(posterior)
rng = MersenneTwister(1234)
BAT.bat_default_withdebug(bat_findmode, Val(:algorithm), posterior)
MersenneTwister <: AbstractRNG
using Optim
optimize(f, a; inplace=false)
algorithm = LBFGSOpt()
transformed_density, trafo = BAT.transform_and_unshape(algorithm.trafo, posterior)
initalg = BAT.apply_trafo_to_init(trafo, algorithm.init)
x_init = collect(bat_initval(rng, transformed_density, initalg).result)
f = BAT.negative(logdensityof(transformed_density))
opts = Optim.Options(store_trace = true, extended_trace=true)
optimize(f, x_init, LBFGS(), opts)

optim_result = optimize(f, x_init, algorithm)
nt = merge(vars, m)
blinded
a = merge(m, blinded)
length(a)
analysis.likelihood
a= zeros(length(vars))
f(vars)
analysis.domain

vars = (ATLAS_lumi=1.0, Lumi = 1.0, Norm_VV_HF =  0.795085, Norm_fake_ExtConv =  1.21885, Norm_fake_HF_el_TandM = 0.691499, 
Norm_fake_HF_mu_TandM = 0.908619, Norm_fake_IntConv = 1.03954, Norm_ttZ_LF = 1.2362, alpha_ATLAS_EG_RES = -0.0783328, 
alpha_ATLAS_EG_SCALE = 3.3991596767996444, alpha_ATLAS_EL_ID = 3.355536440363518, alpha_ATLAS_EL_Isol = 0.021176984362652727, 
alpha_ATLAS_EL_Reco = 1.6327189966495026, alpha_ATLAS_FTAG_B0 = 2.8677316414777936, alpha_ATLAS_FTAG_B1 = 3.270016669331438, 
alpha_ATLAS_FTAG_B10 = 3.604022623079681, alpha_ATLAS_FTAG_B11 = 1.3965600122275332, alpha_ATLAS_FTAG_B12 = 1.9978137923630415,
 alpha_ATLAS_FTAG_B13 = 2.8247368830914805, alpha_ATLAS_FTAG_B14 = 1.3522463910116356, alpha_ATLAS_FTAG_B15 = 0.7074752758962475, 
 alpha_ATLAS_FTAG_B16 = 0.23435908787592663, alpha_ATLAS_FTAG_B17 = 3.0816241388544725, alpha_ATLAS_FTAG_B18 = 1.6437368742481042, 
 alpha_ATLAS_FTAG_B19 = 1.2382757199512708, alpha_ATLAS_FTAG_B2 = 3.1244502411448796, alpha_ATLAS_FTAG_B3 = 4.970828686444538, alpha_ATLAS_FTAG_B4 = 4.599219105581895, alpha_ATLAS_FTAG_B5 = 2.9263831677971424, alpha_ATLAS_FTAG_B6 = 3.7257009841908553, alpha_ATLAS_FTAG_B7 = 4.017904054124441, alpha_ATLAS_FTAG_B8 = 0.10449650914870108, alpha_ATLAS_FTAG_B9 = 3.4786877038803206, alpha_ATLAS_FTAG_C0 = 4.280603441019295, alpha_ATLAS_FTAG_C1 = 0.8984793651813332, alpha_ATLAS_FTAG_C10 = 2.0906272507379517, alpha_ATLAS_FTAG_C11 = 4.102440552170017, alpha_ATLAS_FTAG_C12 = 1.8574305838441412, alpha_ATLAS_FTAG_C13 = 0.3341528274193218, alpha_ATLAS_FTAG_C14 = 4.298560917046819, alpha_ATLAS_FTAG_C15 = 1.387463566550799, alpha_ATLAS_FTAG_C16 = 0.05249410627099269, alpha_ATLAS_FTAG_C17 = 4.390987418742104, alpha_ATLAS_FTAG_C18 = 4.6235678688521675, alpha_ATLAS_FTAG_C19 = 1.507216927572957, alpha_ATLAS_FTAG_C2 = 3.9111142408442148, alpha_ATLAS_FTAG_C3 = 2.402707063505234, alpha_ATLAS_FTAG_C4 = 1.4596824151398788, alpha_ATLAS_FTAG_C5 = 2.1453100936267977, alpha_ATLAS_FTAG_C6 = 4.579044123595924, alpha_ATLAS_FTAG_C7 = 4.7931626609806, alpha_ATLAS_FTAG_C8 = 1.7096267659078113, alpha_ATLAS_FTAG_C9 = 3.0598023538492085, alpha_ATLAS_FTAG_L0 = 2.455431428216102, alpha_ATLAS_FTAG_L1 = 1.7911168388772158, alpha_ATLAS_FTAG_L10 = 4.201349061130297, alpha_ATLAS_FTAG_L11 = 2.9929235615905223, alpha_ATLAS_FTAG_L12 = 1.7589690124464914, alpha_ATLAS_FTAG_L13 = 0.9129073053857837, alpha_ATLAS_FTAG_L14 = 2.3517652953197787, alpha_ATLAS_FTAG_L15 = 2.6255264689808784, alpha_ATLAS_FTAG_L16 = 1.5172449733128015, alpha_ATLAS_FTAG_L17 = 3.8050974284397916, alpha_ATLAS_FTAG_L18 = 2.148793535956451, alpha_ATLAS_FTAG_L19 = 3.5226166361665476, alpha_ATLAS_FTAG_L2 = 2.5946738288762914, alpha_ATLAS_FTAG_L3 = 2.5063026266977757, alpha_ATLAS_FTAG_L4 = 4.291045723176507, alpha_ATLAS_FTAG_L5 = 2.4232748853169572, alpha_ATLAS_FTAG_L6 = 2.8487486841849634, alpha_ATLAS_FTAG_L7 = 1.460037922320199, alpha_ATLAS_FTAG_L8 = 0.9366442013500123, alpha_ATLAS_FTAG_L9 = 0.8051416675862717, alpha_ATLAS_JER_DataVsMC = 0.43390219866549395, alpha_ATLAS_JER_Eff1 = 3.33384671762535, alpha_ATLAS_JER_Eff10 = 3.4043671263340443, alpha_ATLAS_JER_Eff11 = 3.8019925954885636, alpha_ATLAS_JER_Eff12 = 4.136809643177686, alpha_ATLAS_JER_Eff2 = 1.7491676126298386, alpha_ATLAS_JER_Eff3 = 0.9891191763714282, alpha_ATLAS_JER_Eff4 = 1.9491751614984119, alpha_ATLAS_JER_Eff5 = 2.473930051391876, alpha_ATLAS_JER_Eff6 = 1.394953699846537, alpha_ATLAS_JER_Eff7 = 0.8574999477806953, alpha_ATLAS_JER_Eff8 = 2.052773799675887, alpha_ATLAS_JER_Eff9 = 0.24158681830176243, alpha_ATLAS_JES_BJES = 1.9392659000847516, alpha_ATLAS_JES_EtaInter_Model = 4.28163500896678, alpha_ATLAS_JES_EtaInter_NonClosureHighE = 4.162136383539165, alpha_ATLAS_JES_EtaInter_NonClosureNegEta = 2.3866424590061723, alpha_ATLAS_JES_EtaInter_NonClosurePosEta = 2.279205002058777, alpha_ATLAS_JES_EtaInter_NonClosure_2018data = 0.917344026979486, alpha_ATLAS_JES_EtaInter_Stat = 0.06843199169621092, alpha_ATLAS_JES_Flavor_Comp = 4.365437428381378, alpha_ATLAS_JES_Flavor_Resp = 0.33701946186978143, alpha_ATLAS_JES_NP_Det1 = 0.9007848192071337, alpha_ATLAS_JES_NP_Det2 = 3.8902774849477764, alpha_ATLAS_JES_NP_Mix1 = 2.134101349999501, alpha_ATLAS_JES_NP_Mix2 = 2.952811756642836, alpha_ATLAS_JES_NP_Mix3 = 4.982856942023799, alpha_ATLAS_JES_NP_Mod1 = 0.06569306142201264, alpha_ATLAS_JES_NP_Mod2 = 0.8865091507904654, alpha_ATLAS_JES_NP_Mod3 = 0.9423303365040026, alpha_ATLAS_JES_NP_Mod4 = 3.2874963743893204, alpha_ATLAS_JES_NP_Stat1 = 1.0059991905091947, alpha_ATLAS_JES_NP_Stat2 = 3.0822368858075153, alpha_ATLAS_JES_NP_Stat3 = 0.6894913397793415, alpha_ATLAS_JES_NP_Stat4 = 3.630540903488153, alpha_ATLAS_JES_NP_Stat5 = 3.5255069703639186, alpha_ATLAS_JES_NP_Stat6 = 3.051580330665439, alpha_ATLAS_JES_PU_OffsetMu = 0.26855975232764673, alpha_ATLAS_JES_PU_OffsetNPV = 3.0834730912208332, alpha_ATLAS_JES_PU_PtTerm = 2.7431349322029575, alpha_ATLAS_JES_PU_Rho = 3.2420836699108984, alpha_ATLAS_JES_PunchThrough = 1.364029890600915, alpha_ATLAS_JES_SinglePart = 2.5999498198629474, alpha_ATLAS_JVT = 1.1062623140649435, alpha_ATLAS_MET_Para = 0.11752869713334824, alpha_ATLAS_MET_Perp = 1.951218253276815, alpha_ATLAS_MET_Scale = 0.041234030693137864, alpha_ATLAS_MU_ID = 2.4197371413209146, alpha_ATLAS_MU_ID_STAT = 3.8015202399187724, alpha_ATLAS_MU_ID_STAT_LOWPT = 1.4426188498709451, alpha_ATLAS_MU_ID_SYST = 4.751211982623862, alpha_ATLAS_MU_ID_SYST_LOWPT = 3.8683343439386406, alpha_ATLAS_MU_Isol_STAT = 3.2812732881920827, alpha_ATLAS_MU_Isol_SYST = 3.9261407584924592, alpha_ATLAS_MU_MS = 1.1900056726819588, alpha_ATLAS_MU_RESBIAS = 2.640878974412793, alpha_ATLAS_MU_SCALE = 1.1105040718018695, alpha_ATLAS_MU_TTVA_STAT = 0.49910824266623033, alpha_ATLAS_MU_TTVA_SYST = 0.9185050857480846, alpha_ATLAS_PRW_DATASF = 3.840260597702138, alpha_ATLAS_TRIG_EL = 3.8686774746386514, alpha_ATLAS_TRIG_MU_STAT = 3.3668143121842937, alpha_ATLAS_TRIG_MU_SYST = 3.1807988662266435, alpha_ATLAS_lumi = 3.0729963043390667, alpha_ConvIntExtrap = 1.6038144288788927, alpha_ConvMatExtrap = 0.060471081473668, alpha_FakesElLead_PLIV_PtFrac = 4.970843976211812, alpha_FakesElSubLead_PLIV_PtFrac = 0.2296319537682499, alpha_FakesEl_NBjetCorr = 2.3024224705441765, alpha_FakesEl_Tight_MtoTextrap = 1.4021598394105517, alpha_FakesMuLead_PLIV_RelCaloCluster = 0.6560973105394436, alpha_FakesMuLead_PLIV_SVLongSignif = 4.305652057372955, alpha_FakesMuSubLead_PLIV_RelCaloCluster = 2.8836775503334446, alpha_FakesMuSubLead_PLIV_SVLongSignif = 1.4710932198567679, alpha_FakesMu_NBjetCorr = 1.5813336586194178, alpha_FakesMu_Tight_MtoTextrap = 3.627774814379713, alpha_PLIV_El_ID_VeryTight = 2.400864775087011, alpha_PLIV_El_Iso_VeryTight = 4.729765266990068, alpha_PLIV_El_JetModeling_Tight = 0.6444065268696766, alpha_PLIV_El_JetModeling_VeryTight = 4.791225040174218, alpha_PLIV_El_MllWindow_VeryTight = 1.7241191264189586, alpha_PLIV_El_Pileup_Tight = 3.5938546813337013, alpha_PLIV_El_Pileup_VeryTight = 0.7235929603812615, alpha_PLIV_El_Stat_Tight = 4.762943237305686, alpha_PLIV_El_Stat_VeryTight = 2.933412194675823, alpha_PLIV_El_TemplateCut_VeryTight = 4.948731365703891, alpha_PLIV_El_Trigger_Tight = 1.7107401039744752, alpha_PLIV_El_Trigger_VeryTight = 0.5449573611228025, alpha_PLIV_Mu_BkgFraction_VeryTight = 4.638553258954151, alpha_PLIV_Mu_DRMuJet_VeryTight = 1.1641858828688818, alpha_PLIV_Mu_JetModeling_Tight = 2.5065824364512586, alpha_PLIV_Mu_JetModeling_VeryTight = 1.6315243262845356, alpha_PLIV_Mu_Luminosity_VeryTight = 4.025286830983386, alpha_PLIV_Mu_MCxSec_VeryTight = 1.3289013458242895, alpha_PLIV_Mu_MllWindow_Tight = 1.3527268963559378, alpha_PLIV_Mu_MllWindow_VeryTight = 1.9111651594511652, alpha_PLIV_Mu_ProbeQuality_Tight = 1.2836891436502378, alpha_PLIV_Mu_ProbeQuality_VeryTight = 3.754237577617268, alpha_PLIV_Mu_QCDTemplate_VeryTight = 3.7151999942456233, alpha_PLIV_Mu_Stat_Tight = 0.9784491617742797, alpha_PLIV_Mu_Stat_VeryTight = 2.299223493711199, alpha_PLIV_Mu_SuppressionScale_VeryTight = 2.668842627027597, alpha_QMisIDXsec = 4.081765043437891, alpha_QMisID_MM_SYST = 1.7764300590951223, alpha_QMisID_MT_SYST = 0.5508403967425256, alpha_QMisID_TM_SYST = 3.0219520890418803, alpha_VHXsec = 1.301459532774253, alpha_VVVXsec = 2.987858208944391, alpha_VV_varRF = 4.4716324048198395, alpha_VVlightXsec = 2.2906326410384175, alpha_VVnJets = 1.0469010674426669, alpha_WtZXsec = 0.02309965295211485, alpha_fourTopXsec = 0.46102482174971093, alpha_tZXsec = 0.031040630550480012, alpha_threeTopXsec = 2.4766205837287973, alpha_ttHXsec = 2.8380093097726493, alpha_ttH_Gen_Acceptance = 4.285062616855545, alpha_ttH_PS_Acceptance = 1.9670369443928544, alpha_ttH_varRF = 4.487867753539196, alpha_ttWWXsec = 0.016191494251662872, alpha_ttW_EW_fraction = 4.211743008822156, alpha_ttW_EW_varRF = 3.0064433598167737, alpha_ttW_ME = 0.8885352379061271, alpha_ttW_ME_EW = 1.2629256508896987, alpha_ttW_PDFaS = 0.20892573905224734, alpha_ttW_PDFalternate = 1.4700704155406863, alpha_ttW_PS_EW = 1.3133043502921264, alpha_ttW_PS_QCD = 3.6905725953716746, alpha_ttW_varRF = 1.6482184721759154, alpha_ttZHFXsec = 2.8477314602901256, alpha_ttZ_Shower = 0.9200863031331925, alpha_ttZ_Var3c = 3.094854865024121, alpha_ttZ_varRF = 0.597148775311287, alpha_ttbb_XS = 1.51978323549395, alpha_ttcc_XS = 2.042459481854017, nominalLumi = 1.0002354827286402, staterror_bin1 = 3.308805427637473, staterror_bin2 = 1.873677759473629, staterror_bin3 = 3.537547568068721, staterror_bin30 = 0.01790369418092533, staterror_bin14 = 4.801553630623356, shape_stat_bin1 = 3.8540313766480807, staterror_bin41 = 4.719570577254592, shape_stat_bin2 = 4.621485709389548, staterror_bin9 = 0.3515306208099395, staterror_bin22 = 0.5420515093154226, shape_stat_bin6 = 1.4781028422705529, shape_stat_bin44 = 1.5526198534245406, staterror_bin5 = 3.3513509784060127, shape_stat_bin13 = 2.4719561951670155, shape_stat_bin47 = 2.6618980188909678, shape_stat_bin10 = 2.7723133226270686, shape_stat_bin23 = 4.496750153798186, staterror_bin18 = 0.9526643069759486, shape_stat_bin18 = 3.6452130282469577, staterror_bin12 = 0.7217210323981245, staterror_bin21 = 0.6152321968665566, staterror_bin6 = 3.1184616583868516, staterror_bin46 = 1.3439375336453052, staterror_bin45 = 0.3195320134999824, shape_stat_bin48 = 2.2593501416825843, shape_stat_bin4 = 4.658286035750209, shape_stat_bin31 = 3.5874201482259145, staterror_bin32 = 1.4081927034733945, shape_stat_bin15 = 1.7207050731516202, staterror_bin27 = 3.6571888666402064, shape_stat_bin41 = 4.619320112193787, staterror_bin24 = 4.58680806159905, shape_stat_bin26 = 0.6838596568750264, staterror_bin37 = 4.038628436663305, shape_stat_bin9 = 3.306607080574696, staterror_bin4 = 2.895937536751199, staterror_bin25 = 0.18721254330615503, shape_stat_bin35 = 1.7077397832398664, shape_stat_bin43 = 2.836888747890847, shape_stat_bin3 = 2.183351064298128, shape_stat_bin38 = 2.0416722464183583, shape_stat_bin39 = 2.4114641883062538, staterror_bin38 = 3.955058214703085, shape_stat_bin16 = 1.389484366989466, staterror_bin35 = 0.8444606437341586, shape_stat_bin7 = 2.2878126412446327, staterror_bin23 = 4.859602243060119, staterror_bin26 = 0.5939730395197121, shape_stat_bin40 = 1.5504360658397782, shape_stat_bin20 = 0.566164884953962, staterror_bin42 = 0.19493955209962038, staterror_bin48 = 1.3083461393958404, staterror_bin11 = 2.9104113284217132, shape_stat_bin14 = 0.9165720706624734, shape_stat_bin36 = 1.2431753241051908, staterror_bin36 = 3.713496235180622, staterror_bin34 = 2.103180482276733, shape_stat_bin22 = 2.4948337320610148, staterror_bin40 = 1.3687317270126191, staterror_bin33 = 3.9238418111555893, shape_stat_bin46 = 0.06301981061429941, staterror_bin16 = 1.5748897069493502, staterror_bin8 = 0.41145727156903755, staterror_bin31 = 4.783062961611475, shape_stat_bin8 = 2.4003685800209347, shape_stat_bin21 = 2.701461273562679, shape_stat_bin12 = 3.3230768681638856, shape_stat_bin42 = 3.6812028820823697, staterror_bin28 = 0.18705306706494815, staterror_bin50 = 3.3070824351269863, staterror_bin10 = 3.8834442939429312, shape_stat_bin33 = 2.4089102386786094, staterror_bin19 = 3.391606146001267, shape_stat_bin49 = 1.2225331008253906, shape_stat_bin34 = 1.8477745239119352, shape_stat_bin29 = 3.8719045378541392, staterror_bin13 = 4.887748169008098, shape_stat_bin30 = 2.692023881820337, staterror_bin15 = 1.173567396852957, staterror_bin29 = 2.551742706933016, shape_stat_bin11 = 3.5924432350175044, shape_stat_bin17 = 1.545973871262736, staterror_bin17 = 0.3393702268333365, staterror_bin39 = 2.5907286676293904, staterror_bin44 = 0.9307965641171765, shape_stat_bin28 = 2.605589084966212, shape_stat_bin5 = 2.8809523211978334, staterror_bin7 = 3.3703393681540734, staterror_bin20 = 2.2919728191121878, shape_stat_bin27 = 2.9407040751202587, staterror_bin43 = 4.072766377810427, shape_stat_bin24 = 4.99049000374366, shape_stat_bin50 = 4.350483356967466, staterror_bin47 = 3.3490506844670977, shape_stat_bin37 = 3.82335005101768, staterror_bin49 = 4.074942330924208, shape_stat_bin45 = 0.03163247721187448, shape_stat_bin19 = 1.96563108222552, shape_stat_bin32 = 1.9615828548381637, shape_stat_bin25 = 0.34791013069313115)

m = (Lumi = 1.0, ttW_Bin_001_Ac = 0.453084, 
ttW_Bin_001_muInc = 1.56116, 
ttW_Bin_002_Ac = 0.399231,
ttW_Bin_002_muInc = 1.58663, 
ttW_Bin_003_Ac = 0.258398 ,
ttW_Bin_003_muInc = 1.66066,
ttW_Bin_004_Ac = 0.209334 ,
ttW_Bin_004_muInc = 1.73644,
ttW_Bin_005_Ac = 0.282084 ,
ttW_Bin_005_muInc = 1.66377,
ttW_Bin_006_Ac = 0.410482 ,
ttW_Bin_006_muInc = 1.67524,
ttW_Bin_007_Ac = 0.45381 ,
ttW_Bin_007_muInc = 1.64433, 
ttW_Bin_008_Ac = 0.535322,
ttW_Bin_008_muInc = 1.71448 ,
ttW_Bin_009_Ac = 0.121239,
ttW_Bin_009_muInc = 2.49353,
ttW_Bin_010_Ac = 0.0122867,
ttW_Bin_010_muInc = 2.3246 ,
ttW_Bin_011_Ac = -0.0939805 ,
ttW_Bin_011_muInc = 1.61365 ,
ttW_Bin_012_Ac = 0.191556 ,
ttW_Bin_012_muInc = 1.62983,
ttW_Bin_013_Ac = 0.205394 ,
ttW_Bin_013_muInc = 1.79733, 
ttW_Bin_014_Ac = 0.635568 ,
ttW_Bin_014_muInc = 1.90282 )

using NamedTupleTools
using Random
using ValueShapes
using Distributions

mutable struct Param
    name::String
    nominal::Float64
    stdev::Float64
end


function parse_line(line::String)
    parts = split(line)
    name = parts[1]
    if occursin("gamma_stat_CR_H", name)
    
    elseif occursin("gamma_stat_Z", name)

    elseif occursin("gamma_stat", name)
        bin_number = match(r"bin_(\d+)$", name).captures[1]
        name = "staterror_bin$(bin_number)"
    elseif occursin("gamma_shape_stat", name)
        println("here")
        bin_number = match(r"bin_(\d+)$", name).captures[1]
        println(bin_number)
        name = "shape_stat_bin$(bin_number)"
    end
    nominal = parse(Float64, parts[2])
    stdev = parse(Float64, parts[3])
    return Param(name, nominal, stdev)
end

function parse_file(filename::String)
    params = Param[]
    open(filename, "r") do f
        for line in eachline(f)
            if occursin("NUISANCE_PARAMETERS", line)
                continue
            elseif occursin("+", line) && occursin("-", line) # Checks if line contains a parameter
                push!(params, parse_line(line))
            elseif !occursin("+", line) && !occursin("-", line) # If line doesn't contain parameter, then it's the correlation matrix
                break
            end
        end
    end
    return params
end

function generate_nominal(named_params)
    return NamedTuple{Tuple(Symbol.(map(p->p.name, named_params)))}(map(p->p.nominal, named_params))
end

function generate_random(named_params)
    random_values = map(p->p.nominal + randn()*p.stdev, named_params)
    return NamedTuple{Tuple(Symbol.(map(p->p.name, named_params)))}(random_values)
end

function generate_distribution(named_params)
    dists = map(named_params) do p
        min_val = p.nominal - p.stdev # Prevents negative values
        max_val = p.nominal + p.stdev
        Uniform(min_val, max_val)
    end
    return NamedTupleDist(NamedTuple{Tuple(Symbol.(map(p->p.name, named_params)))}(dists))
end

# Test the script
named_params = parse_file("../ttW/results/diff-asym/ttWUnfolding_CombTruth_ttW_unfold_ttWnom_DPhill_SS_TopApproval_unblinded.txt")
blinded =  generate_nominal(named_params)
blinded = (generate_distribution(named_params))
using Distributions
using NamedTupleTools
function generate_namedtuple()
    #names_stat = [Symbol("staterror_bin$(i)") for i in 1:52]
    names_shape = [Symbol("shape_stat_bin$(i)") for i in 1:52]
    values = [Uniform(0.8, 1.2) for i in 1:52]
    return namedtuple(names_shape, values)
end
h = generate_namedtuple()

println(generate_random(named_params))
println(generate_distribution(named_params))



function parse_line(line::String)
    parts = split(line)
    name = parts[1]
    nominal = parse(Float64, parts[2])
    return (name, nominal)
end

function parse_file(filename::String)
    open(filename, "r") do f
        for line in eachline(f)
            if occursin("NUISANCE_PARAMETERS", line)
                continue
            elseif occursin("+", line) && occursin("-", line) # Checks if line contains a parameter
                name, nominal = parse_line(line)
                print("-p $(name) $(nominal) ")
            elseif !occursin("+", line) && !occursin("-", line) # If line doesn't contain parameter, then it's the correlation matrix
                break
            end
        end
    end
end

# Test the script
filename = "../ttW/results/diff-asym/ttWUnfolding_CombTruth_ttW_unfold_ttWnom_DPhill_SS_TopApproval_unblinded.txt"
parse_file(filename)


















































nt = (
Norm_VV_HF = 0.795085,
Norm_fake_ExtConv = 1.21885,
Norm_fake_HF_el_TandM = 0.691499,
Norm_fake_HF_mu_TandM = 0.908619, 
Norm_fake_IntConv = 1.03954, 
Norm_ttZ_LF = 1.2362, 
alpha_ATLAS_EG_RES = -0.0783328, 
alpha_ATLAS_EG_SCALE = -0.142793, 
alpha_ATLAS_EL_ID = 0.00805027, 
alpha_ATLAS_EL_Isol = 0.00193761, 
alpha_ATLAS_EL_Reco = 0.00726534,
alpha_ATLAS_FTAG_B0 = -0.396246, 
alpha_ATLAS_FTAG_B1 = -0.211822, 
alpha_ATLAS_FTAG_B10 = 0.00135416, 
alpha_ATLAS_FTAG_B11 = -0.00113703, 
alpha_ATLAS_FTAG_B12 = -0.00504314, 
alpha_ATLAS_FTAG_B13 = -0.00171574, 
alpha_ATLAS_FTAG_B14 = 4.51001 * 10^(-5), 
alpha_ATLAS_FTAG_B15 = -0.000628973,
alpha_ATLAS_FTAG_B16 = -0.00266549,
alpha_ATLAS_FTAG_B17 = -0.00159277,
alpha_ATLAS_FTAG_B18 = -0.00189272,
alpha_ATLAS_FTAG_B19 = 0.000653854,
alpha_ATLAS_FTAG_B2 =  0.0102706,
alpha_ATLAS_FTAG_B3 = -0.0107347,
alpha_ATLAS_FTAG_B4 = 0.00271318,
alpha_ATLAS_FTAG_B5 = 0.00473005,
alpha_ATLAS_FTAG_B6 = 0.000559886,
alpha_ATLAS_FTAG_B7 = 0.00341106,
alpha_ATLAS_FTAG_B8 = -0.00306363,
alpha_ATLAS_FTAG_B9 = -0.00195058,
alpha_ATLAS_FTAG_C0 = -0.097896,
alpha_ATLAS_FTAG_C1 = -0.0111419,
alpha_ATLAS_FTAG_C10 = -0.000656775,
alpha_ATLAS_FTAG_C11 = 0.000298785,
alpha_ATLAS_FTAG_C12 = -0.000794505,
alpha_ATLAS_FTAG_C13 = -0.000911324,
alpha_ATLAS_FTAG_C14 = -0.00324669,
alpha_ATLAS_FTAG_C15 = 0.00107841,
alpha_ATLAS_FTAG_C16 = 0.00326949,
alpha_ATLAS_FTAG_C17 = 0.00507692,
alpha_ATLAS_FTAG_C18 = -0.00191331,
alpha_ATLAS_FTAG_C19 = 0.00199327,
alpha_ATLAS_FTAG_C2 = 0.012606,  
alpha_ATLAS_FTAG_C3 = 0.0142724, 
alpha_ATLAS_FTAG_C4 = -0.00140911,
alpha_ATLAS_FTAG_C5 = -0.00067549,
alpha_ATLAS_FTAG_C6 = -0.00225739,
alpha_ATLAS_FTAG_C7 = -0.00141468,
alpha_ATLAS_FTAG_C8 = 5.78161e-05,
alpha_ATLAS_FTAG_C9 = 0.000627221,
alpha_ATLAS_FTAG_L0 = 0.0246624,
alpha_ATLAS_FTAG_L1 = 0.0232931,
alpha_ATLAS_FTAG_L10 = -0.000959045,
alpha_ATLAS_FTAG_L11 = 0.00221524,
alpha_ATLAS_FTAG_L12 = -0.00271024,
alpha_ATLAS_FTAG_L13 = -0.000630025,
alpha_ATLAS_FTAG_L14 = 0.002078330,                 
alpha_ATLAS_FTAG_L15 = 0.000368447,                 
alpha_ATLAS_FTAG_L16 = 0.004818490,                                                             
alpha_ATLAS_FTAG_L17 = 4.54547e-05,                 
alpha_ATLAS_FTAG_L18 = 5.46692e-06,                 
alpha_ATLAS_FTAG_L19 = 0.000135397,                 
alpha_ATLAS_FTAG_L2 = -0.00149133,
alpha_ATLAS_FTAG_L3 = 0.0109903,
alpha_ATLAS_FTAG_L4 = 0.0149955,
alpha_ATLAS_FTAG_L5 = -0.00107064,
alpha_ATLAS_FTAG_L6 = 0.000183499,
alpha_ATLAS_FTAG_L7 = -0.0195806,
alpha_ATLAS_FTAG_L8 = -0.00737304,
alpha_ATLAS_FTAG_L9 = 0.00241739,
ATLAS_JER_DataVsMC = 0.0710019,
alpha_ATLAS_JER_Eff1 = -0.0683332,
alpha_ATLAS_JER_Eff10 = -0.0511076, 
alpha_ATLAS_JER_Eff11 = -0.00818135, 
alpha_ATLAS_JER_Eff12 = -0.121801, 
alpha_ATLAS_JER_Eff2 = -0.0583368, 
alpha_ATLAS_JER_Eff3 = 0.103501, 
alpha_ATLAS_JER_Eff4 = 0.0372809, 
alpha_ATLAS_JER_Eff5 = 0.0183146,
alpha_ATLAS_JER_Eff6 = 0.051876, 
alpha_ATLAS_JER_Eff7 = 0.0016204,
alpha_ATLAS_JER_Eff8 = -0.113716,
alpha_ATLAS_JER_Eff9 = -0.0626828,
alpha_ATLAS_JES_BJES = -0.041382,
alpha_ATLAS_JES_EtaInter_Model = -0.279227,
alpha_ATLAS_JES_EtaInter_NonClosureHighE = -0.000154187,
alpha_ATLAS_JES_EtaInter_NonClosureNegEta =  0.00720441,
alpha_ATLAS_JES_EtaInter_NonClosurePosEta = -0.0434353,
alpha_ATLAS_JES_EtaInter_NonClosure_2018data = 0.018299,
alpha_ATLAS_JES_EtaInter_Stat  -0.00299782,
alpha_ATLAS_JES_Flavor_Comp  -0.166258,
alpha_ATLAS_JES_Flavor_Resp  0.255133,
alpha_ATLAS_JES_NP_Det1  0.00162799,
alpha_ATLAS_JES_NP_Det2  0.000989139,
alpha_ATLAS_JES_NP_Mix1  0.00478393,
alpha_ATLAS_JES_NP_Mix2  -9.99438e-05,
alpha_ATLAS_JES_NP_Mix3  0.00410553,
alpha_ATLAS_JES_NP_Mod1  -0.15331,
alpha_ATLAS_JES_NP_Mod2  -0.0102389,
alpha_ATLAS_JES_NP_Mod3  0.00657461,
alpha_ATLAS_JES_NP_Mod4  0.000922986,
alpha_ATLAS_JES_NP_Stat1  0.00193169,
alpha_ATLAS_JES_NP_Stat2  -0.00710011,
alpha_ATLAS_JES_NP_Stat3  0.000797382 ,
alpha_ATLAS_JES_NP_Stat4  -0.00610105,
alpha_ATLAS_JES_NP_Stat5  0.00154454 ,
alpha_ATLAS_JES_NP_Stat6  0.00223964,
alpha_ATLAS_JES_PU_OffsetMu  -0.0406879,
alpha_ATLAS_JES_PU_OffsetNPV  -0.143962,
alpha_ATLAS_JES_PU_PtTerm  -0.00370718,
alpha_ATLAS_JES_PU_Rho  -0.437293,
alpha_ATLAS_JES_PunchThrough  0.000786611,
alpha_ATLAS_JES_SinglePart  -0.00313452,
alpha_ATLAS_JVT  0.00232121,
alpha_ATLAS_MET_Para  -0.258547,
alpha_ATLAS_MET_Perp  0.0107446,
alpha_ATLAS_MET_Scale  0.00493141,
alpha_ATLAS_MU_ID  0.035623,
alpha_ATLAS_MU_ID_STAT  0.00375435,
alpha_ATLAS_MU_ID_STAT_LOWPT  -0.00187914,
alpha_ATLAS_MU_ID_SYST  0.0163597,
alpha_ATLAS_MU_ID_SYST_LOWPT  -0.0086716,
alpha_ATLAS_MU_Isol_STAT  0.00339818,
alpha_ATLAS_MU_Isol_SYST  0.0128649,
alpha_ATLAS_MU_MS  0.0280877,
alpha_ATLAS_MU_RESBIAS  1.59963e-05,
alpha_ATLAS_MU_SCALE  -0.0262178,
alpha_ATLAS_MU_TTVA_STAT  0.0022072,
alpha_ATLAS_MU_TTVA_SYST  0.00374264,
alpha_ATLAS_PRW_DATASF  -0.226927,
alpha_ATLAS_TRIG_EL  -0.0227215,
alpha_ATLAS_TRIG_MU_STAT  0.00463226,
alpha_ATLAS_TRIG_MU_SYST  -0.089891,
alpha_ATLAS_lumi  0.00323447,
alpha_ConvIntExtrap  -0.456217,
alpha_ConvMatExtrap  -0.0722369,
alpha_FakesElLead_PLIV_PtFrac  -0.106938,
alpha_FakesElSubLead_PLIV_PtFrac  -0.0557645,
alpha_FakesEl_NBjetCorr  -0.0105628,
alpha_FakesEl_Tight_MtoTextrap  -0.0346501,
alpha_FakesMuLead_PLIV_RelCaloCluster  0.325006,
alpha_FakesMuLead_PLIV_SVLongSignif  -0.263497,
alpha_FakesMuSubLead_PLIV_RelCaloCluster  -0.0774655,
alpha_FakesMuSubLead_PLIV_SVLongSignif  0.125008,
alpha_FakesMu_NBjetCorr  -0.01402,
alpha_FakesMu_Tight_MtoTextrap  -0.0210883,
alpha_PLIV_El_ID_VeryTight  0.0153886,
alpha_PLIV_El_Iso_VeryTight  0.0158845,
alpha_PLIV_El_JetModeling_Tight  -0.0756701 +0.986527 -0.986527
alpha_PLIV_El_JetModeling_VeryTight  -0.0381042 +0.998346 -0.998346
alpha_PLIV_El_MllWindow_VeryTight  0.015204 +0.998412 -0.998412
alpha_PLIV_El_Pileup_Tight  -0.0408464 +0.992062 -0.992062
alpha_PLIV_El_Pileup_VeryTight  0.0234225 +0.998309 -0.998309
alpha_PLIV_El_Stat_Tight  0.000543556 +0.993249 -0.993249
alpha_PLIV_El_Stat_VeryTight  0.0186413 +0.99833 -0.99833
alpha_PLIV_El_TemplateCut_VeryTight  0.0151815 +0.998422 -0.998422
alpha_PLIV_El_Trigger_Tight  -0.0375833 +0.989141 -0.989141
alpha_PLIV_El_Trigger_VeryTight  -0.108469 +0.98931 -0.98931
alpha_PLIV_Mu_BkgFraction_VeryTight  0.000712551 +0.993347 -0.993347
alpha_PLIV_Mu_DRMuJet_VeryTight  0.0147352 +0.998414 -0.998414
alpha_PLIV_Mu_JetModeling_Tight  -0.0483434 +0.987156 -0.987156
alpha_PLIV_Mu_JetModeling_VeryTight  0.156237 +0.990378 -0.990378
alpha_PLIV_Mu_Luminosity_VeryTight  0.0166909 +0.998368 -0.998368
alpha_PLIV_Mu_MCxSec_VeryTight  0.000117397 +0.993347 -0.993347
alpha_PLIV_Mu_MllWindow_Tight  -9.30215e-06 +0.993347 -0.993347
alpha_PLIV_Mu_MllWindow_VeryTight  0.00395776 +0.993342 -0.993342
alpha_PLIV_Mu_ProbeQuality_Tight  4.10466e-05 +0.993347 -0.993347
alpha_PLIV_Mu_ProbeQuality_VeryTight  0.022877 +0.998285 -0.998285
alpha_PLIV_Mu_QCDTemplate_VeryTight  0.00223977 +0.993345 -0.993345
alpha_PLIV_Mu_Stat_Tight  -0.0296234 +0.992796 -0.992796
alpha_PLIV_Mu_Stat_VeryTight  0.0035162 +0.993344 -0.993344
alpha_PLIV_Mu_SuppressionScale_VeryTight  0.000299516 +0.993347 -0.993347
alpha_QMisIDXsec  -0.00107939 +0.993242 -0.993242
alpha_QMisID_MM_SYST  1.61077e-05 +0.970848 -0.970848
alpha_QMisID_MT_SYST  -0.399013 +0.914571 -0.914571
alpha_QMisID_TM_SYST  0.0649161 +0.972 -0.972
alpha_VHXsec  0.148397 +0.957935 -0.957935
alpha_VVVXsec  -0.00631841 +0.994156 -0.994156
alpha_VV_varRF  0.0626069 +0.987098 -0.987098
alpha_VVlightXsec  -0.000279404 +0.993346 -0.993346
alpha_VVnJets  0.0362769 +0.993605 -0.993605
alpha_WtZXsec  -0.397132 +0.944428 -0.944428
alpha_fourTopXsec  0.0686055 +0.94656 -0.94656
alpha_tZXsec  0.0416696 +0.991432 -0.991432
alpha_threeTopXsec  -0.000228912 +0.993381 -0.993381
alpha_ttHXsec  -0.0513412 +0.996043 -0.996043
alpha_ttH_Gen_Acceptance  -0.011596 +0.979626 -0.979626
alpha_ttH_PS_Acceptance  0.0554918 +0.991196 -0.991196
alpha_ttH_varRF  0.0463539 +1.01452 -1.01452
alpha_ttWWXsec  -0.0194861 +0.997281 -0.997281
alpha_ttW_EW_fraction  0.0308431 +0.991211 -0.991211
alpha_ttW_EW_varRF  0.00583278 +0.991384 -0.991384
alpha_ttW_ME  0.602712 +0.614745 -0.614745
alpha_ttW_ME_EW  0.23755 +0.999273 -0.999273
alpha_ttW_PDFaS  0.000401612 +0.993101 -0.993101
alpha_ttW_PDFalternate  -0.0515659 +0.992617 -0.992617
alpha_ttW_PS_EW  0.0442807 +0.996674 -0.996674
alpha_ttW_PS_QCD  0.455915 +0.962285 -0.962285
alpha_ttW_varRF  -0.0109105 +1.03198 -1.03198
alpha_ttZHFXsec  0.271095 +0.858363 -0.858363
alpha_ttZ_Shower  -0.00804616 +0.987806 -0.987806
alpha_ttZ_Var3c  -0.983658 +1.02501 -1.02501
alpha_ttZ_varRF  0.0251098 +0.99257 -0.99257
alpha_ttbb_XS  0.058699 +0.989185 -0.989185
alpha_ttcc_XS  -0.0117024 +0.992809 -0.992809
gamma_shape_stat_ttW_CombReco_Truth_bin_10_CombReco_bin_15  1.00202 +0.0549613 -0.0549613
gamma_shape_stat_ttW_CombReco_Truth_bin_10_CombReco_bin_21  0.998524 +0.0777261 -0.0777261
gamma_shape_stat_ttW_CombReco_Truth_bin_10_CombReco_bin_9  0.993485 +0.084272 -0.084272
gamma_shape_stat_ttW_CombReco_Truth_bin_11_CombReco_bin_10  0.990771 +0.0595494 -0.0595494
gamma_shape_stat_ttW_CombReco_Truth_bin_11_CombReco_bin_16  1.00032 +0.0417086 -0.0417086
gamma_shape_stat_ttW_CombReco_Truth_bin_11_CombReco_bin_22  0.999084 +0.0699882 -0.0699882
gamma_shape_stat_ttW_CombReco_Truth_bin_12_CombReco_bin_11  0.999061 +0.0567844 -0.0567844
gamma_shape_stat_ttW_CombReco_Truth_bin_12_CombReco_bin_17  1.0012 +0.0396909 -0.0396909
gamma_shape_stat_ttW_CombReco_Truth_bin_12_CombReco_bin_23  1.00564 +0.065482 -0.065482
gamma_shape_stat_ttW_CombReco_Truth_bin_13_CombReco_bin_12  1.00246 +0.037872 -0.037872
gamma_shape_stat_ttW_CombReco_Truth_bin_13_CombReco_bin_18  1.0008 +0.0280634 -0.0280634
gamma_shape_stat_ttW_CombReco_Truth_bin_13_CombReco_bin_24  0.996698 +0.0418044 -0.0418044
gamma_shape_stat_ttW_CombReco_Truth_bin_21_CombReco_bin_32  0.991866 +0.0811324 -0.0811324
gamma_shape_stat_ttW_CombReco_Truth_bin_21_CombReco_bin_38  1.00197 +0.0614457 -0.0614457
gamma_shape_stat_ttW_CombReco_Truth_bin_21_CombReco_bin_44  1.00232 +0.0963785 -0.0963785
gamma_shape_stat_ttW_CombReco_Truth_bin_22_CombReco_bin_33  0.985215 +0.0795357 -0.0795357
gamma_shape_stat_ttW_CombReco_Truth_bin_22_CombReco_bin_39  0.999515 +0.0624903 -0.0624903
gamma_shape_stat_ttW_CombReco_Truth_bin_22_CombReco_bin_45  1.01076 +0.0941695 -0.0941695
gamma_shape_stat_ttW_CombReco_Truth_bin_23_CombReco_bin_34  1.00421 +0.0693641 -0.0693641
gamma_shape_stat_ttW_CombReco_Truth_bin_23_CombReco_bin_40  1.004 +0.0540262 -0.0540262
gamma_shape_stat_ttW_CombReco_Truth_bin_23_CombReco_bin_46  1.00288 +0.0799271 -0.0799271
gamma_shape_stat_ttW_CombReco_Truth_bin_24_CombReco_bin_35  0.986723 +0.0716223 -0.0716223
gamma_shape_stat_ttW_CombReco_Truth_bin_24_CombReco_bin_41  1.00089 +0.0569114 -0.0569114
gamma_shape_stat_ttW_CombReco_Truth_bin_24_CombReco_bin_47  1.01073 +0.0938187 -0.0938187
gamma_shape_stat_ttW_CombReco_Truth_bin_25_CombReco_bin_36  0.999996 +0.075921 -0.075921
gamma_shape_stat_ttW_CombReco_Truth_bin_25_CombReco_bin_42  0.997338 +0.0520458 -0.0520458
gamma_shape_stat_ttW_CombReco_Truth_bin_25_CombReco_bin_48  1.01026 +0.0813893 -0.0813893
gamma_shape_stat_ttW_CombReco_Truth_bin_26_CombReco_bin_37  0.998024 +0.0504316 -0.0504316
gamma_shape_stat_ttW_CombReco_Truth_bin_26_CombReco_bin_43  0.999676 +0.0373257 -0.0373257
gamma_shape_stat_ttW_CombReco_Truth_bin_26_CombReco_bin_49  1.00191 +0.057964 -0.057964
gamma_shape_stat_ttW_CombReco_Truth_bin_8_CombReco_bin_13  1.00322 +0.0437828 -0.0437828
gamma_shape_stat_ttW_CombReco_Truth_bin_8_CombReco_bin_19  1.00283 +0.0666523 -0.0666523
gamma_shape_stat_ttW_CombReco_Truth_bin_8_CombReco_bin_7  1.00385 +0.0569496 -0.0569496
gamma_shape_stat_ttW_CombReco_Truth_bin_9_CombReco_bin_14  0.996113 +0.0420375 -0.0420375
gamma_shape_stat_ttW_CombReco_Truth_bin_9_CombReco_bin_20  0.999088 +0.0655785 -0.0655785
gamma_shape_stat_ttW_CombReco_Truth_bin_9_CombReco_bin_8  1.0024 +0.0571833 -0.0571833
gamma_stat_CR_HFel_2lSS_MM_lepPt1_2jincl_1b_bin_0  1.00038 +0.116963 -0.116963
gamma_stat_CR_HFel_2lSS_MM_lepPt1_2jincl_1b_bin_1  1.00141 +0.0871377 -0.0871377
gamma_stat_CR_HFel_2lSS_MT_lepPt1_2jincl_1b_lowMTlepMET_bin_0  0.958111 +0.10088 -0.10088
gamma_stat_CR_HFel_2lSS_MT_lepPt1_2jincl_1b_lowMTlepMET_bin_1  0.983441 +0.0650199 -0.0650199
gamma_stat_CR_HFel_2lSS_TM_lepPt1_2jincl_1b_lowMTlepMET_bin_0  1.00556 +0.0535969 -0.0535969
gamma_stat_CR_HFel_2lSS_TM_lepPt1_2jincl_1b_lowMTlepMET_bin_1  1.00085 +0.0424247 -0.0424247
gamma_stat_CR_HFmu_2lSS_MM_lepPt1_2jincl_1b_bin_0  0.9966 +0.123085 -0.123085
gamma_stat_CR_HFmu_2lSS_MM_lepPt1_2jincl_1b_bin_1  0.991552 +0.0999333 -0.0999333
gamma_stat_CR_HFmu_2lSS_MT_lepPt1_2jincl_1b_lowMTlepMET_bin_0  0.992001 +0.062177 -0.062177
gamma_stat_CR_HFmu_2lSS_MT_lepPt1_2jincl_1b_lowMTlepMET_bin_1  1.00147 +0.0516016 -0.0516016
gamma_stat_CR_HFmu_2lSS_TM_lepPt1_2jincl_1b_lowMTlepMET_bin_0  1.18814 +0.140007 -0.140007
gamma_stat_CR_HFmu_2lSS_TM_lepPt1_2jincl_1b_lowMTlepMET_bin_1  1.00207 +0.0463213 -0.0463213
gamma_stat_CombReco_bin_1  1.00021 +0.0237213 -0.0237213
gamma_stat_CombReco_bin_10  0.984674 +0.0720823 -0.0720823
gamma_stat_CombReco_bin_11  0.998878 +0.0682687 -0.0682687
gamma_stat_CombReco_bin_12  1.00436 +0.0638377 -0.0638377
gamma_stat_CombReco_bin_13  1.0074 +0.0607966 -0.0607966
gamma_stat_CombReco_bin_14  0.98529 +0.0612291 -0.0612291
gamma_stat_CombReco_bin_15  1.00103 +0.0284251 -0.0284251
gamma_stat_CombReco_bin_16  1.0004 +0.0332203 -0.0332203
gamma_stat_CombReco_bin_17  1.00365 +0.0540956 -0.0540956
gamma_stat_CombReco_bin_18  1.00105 +0.031378 -0.031378
gamma_stat_CombReco_bin_19  1.00092 +0.0281731 -0.0281731
gamma_stat_CombReco_bin_2  1.00106 +0.0203207 -0.0203207
gamma_stat_CombReco_bin_20  0.999608 +0.0305933 -0.0305933
gamma_stat_CombReco_bin_21  0.999427 +0.0316146 -0.0316146
gamma_stat_CombReco_bin_22  0.999381 +0.0349312 -0.0349312
gamma_stat_CombReco_bin_23  1.01925 +0.0767348 -0.0767348
gamma_stat_CombReco_bin_24  0.995365 +0.0411627 -0.0411627
gamma_stat_CombReco_bin_25  1.00039 +0.0207922 -0.0207922
gamma_stat_CombReco_bin_26  1.00018 +0.0245784 -0.0245784
gamma_stat_CombReco_bin_28  1 +0.0249988 -0.0249988
gamma_stat_CombReco_bin_29  1.00288 +0.0297751 -0.0297751
gamma_stat_CombReco_bin_3  0.99922 +0.0228942 -0.0228942
gamma_stat_CombReco_bin_30  0.990116 +0.0384356 -0.0384356
gamma_stat_CombReco_bin_32  0.994608 +0.062872 -0.062872
gamma_stat_CombReco_bin_33  0.992165 +0.0581959 -0.0581959
gamma_stat_CombReco_bin_34  1.00199 +0.0493405 -0.0493405
gamma_stat_CombReco_bin_35  0.987002 +0.0727623 -0.0727623
gamma_stat_CombReco_bin_36  0.999997 +0.0699776 -0.0699776
gamma_stat_CombReco_bin_37  0.996318 +0.0633996 -0.0633996
gamma_stat_CombReco_bin_38  1.00457 +0.0662171 -0.0662171
gamma_stat_CombReco_bin_39  0.998727 +0.0651608 -0.0651608
gamma_stat_CombReco_bin_4  0.997019 +0.0263348 -0.0263348
gamma_stat_CombReco_bin_40  1.0023 +0.0310954 -0.0310954
gamma_stat_CombReco_bin_41  1.00064 +0.0361734 -0.0361734
gamma_stat_CombReco_bin_42  0.992231 +0.056932 -0.056932
gamma_stat_CombReco_bin_43  0.999457 +0.0324359 -0.0324359
gamma_stat_CombReco_bin_44  1.00053 +0.0294804 -0.0294804
gamma_stat_CombReco_bin_45  1.00298 +0.0328786 -0.0328786
gamma_stat_CombReco_bin_46  1.00116 +0.0336299 -0.0336299
gamma_stat_CombReco_bin_47  1.00432 +0.0371153 -0.0371153
gamma_stat_CombReco_bin_48  1.00408 +0.0319606 -0.0319606
gamma_stat_CombReco_bin_49  1.00325 +0.0458225 -0.0458225
gamma_stat_CombReco_bin_5  1.0068 +0.0362783 -0.0362783
gamma_stat_CombReco_bin_7  1.00326 +0.0625072 -0.0625072
gamma_stat_CombReco_bin_8  1.00213 +0.0579893 -0.0579893
gamma_stat_CombReco_bin_9  0.997938 +0.0489743 -0.0489743
gamma_stat_Zgamma_Ext_3l_bin_0  1.00886 +0.0936175 -0.0936175
gamma_stat_Zgamma_Int_3l_bin_0  1.00661 +0.138331 -0.138331
ttW_Bin_001_Ac  0.327371 +0.142797 -0.142797
ttW_Bin_001_muInc  1.57569 +0.30582 -0.30582
ttW_Bin_002_Ac  0.412374 +0.109402 -0.109402
ttW_Bin_002_muInc  1.51757 +0.260151 -0.260151
ttW_Bin_003_Ac  0.462804 +0.109688 -0.109688
ttW_Bin_003_muInc  1.4114 +0.247928 -0.247928
ttW_Bin_004_Ac  0.408978 +0.120397 -0.120397
ttW_Bin_004_muInc  1.24454 +0.23608 -0.23608
ttW_Bin_005_Ac  0.391153 +0.113491 -0.113491
ttW_Bin_005_muInc  1.28059 +0.243524 -0.243524
ttW_Bin_006_Ac  0.437791 +0.088965 -0.088965
ttW_Bin_006_muInc  1.5849 +0.249418 -0.249418
ttW_Bin_007_Ac  0.309153 +0.0756092 -0.0756092
ttW_Bin_007_muInc  1.75747 +0.255414 -0.255414
ttW_Bin_008_Ac  0.229588 +0.210003 -0.210003
ttW_Bin_008_muInc  2.39738 +0.65225 -0.65225
ttW_Bin_009_Ac  0.127138 +0.171567 -0.171567
ttW_Bin_009_muInc  2.0233 +0.532183 -0.532183
ttW_Bin_010_Ac  0.0214413 +0.174765 -0.174765
ttW_Bin_010_muInc  1.75371 +0.477505 -0.477505
ttW_Bin_011_Ac  0.0241357 +0.183036 -0.183036
ttW_Bin_011_muInc  1.57277 +0.440651 -0.440651
ttW_Bin_012_Ac  0.169722 +0.164302 -0.164302
ttW_Bin_012_muInc  1.52623 +0.390583 -0.390583
ttW_Bin_013_Ac  0.368124 +0.196405 -0.196405
ttW_Bin_013_muInc  1.49806 +0.388936 -0.388936




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