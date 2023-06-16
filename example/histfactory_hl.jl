using HS3
using BenchmarkTools, DensityInterface

dict = file_to_dict("./example/histfactory_hl.json")
specs = HS3.generate_specs(dict)
analysis = HS3.make_analyses(specs.analyses[1], specs)
logdensityof(analysis.likelihood, (mu=1, alpha_syst2 =0, staterror_bin2 = .1, alpha_syst1 = 0, alpha_syst3 = 0,  staterror_bin1 =.1))

using BAT, ValueShapes, Distributions
posterior = PosteriorMeasure(analysis.likelihood, NamedTupleDist(mu=Uniform(-3, 5), alpha_syst2 = Uniform(-5, 5), staterror_bin2 = .1, alpha_syst1 = Uniform(0, 1), alpha_syst3 = Uniform(-5, 5),  staterror_bin1 = Uniform(0, 1)))
bat_findmode(posterior).result

#make analysis directly from file
@benchmark HS3.make_analyses("./example/histfactory_hl.json", "primary_analysis")