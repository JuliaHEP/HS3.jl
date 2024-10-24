# This file is a part of HS3.jl, licensed under the MIT License (MIT).
using HS3
using Test
using DensityInterface
using LiteHF

dict = load_pyhfjson("./test/jsons/rf515_hf_LiteHF.json")

pyhfmodel = build_pyhf(dict)
dump(priors(pyhfmodel))
LL = pyhf_logjointof(pyhfmodel)
LL( [1, 1., 1., 1., 1., 1.])

dict = file_to_dict("./test/jsons/rf515_hfJSON.json")
specs = HS3.generate_specs(dict)

analysis = HS3.make_analyses(specs.analyses[1], specs)

mu_1 = logdensityof(analysis.likelihood, (mu=1.,  syst3 = 0., syst1 = 0., syst2 = 0., mcstat_bin2 = 1., mcstat_bin1 =1.,))


303/19
185.1122/16.529
