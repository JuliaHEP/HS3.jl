using HS3
using BenchmarkTools, DensityInterface, Random, Distributions

dict = file_to_dict("./example/simple_signal.json")
specs = HS3.generate_specs(dict)
analysis = HS3.make_analyses(specs.analyses[1], specs)
logdensityof(analysis.likelihood, rand(analysis.prior))

using BAT, ValueShapes, Distributions
posterior = PosteriorMeasure(analysis.likelihood, analysis.prior)
bat_findmode(posterior).result

#Validator
@benchmark validate_file("./example/simple_signal.json")

