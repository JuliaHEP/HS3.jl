# This file is a part of HS3.jl, licensed under the MIT License (MIT).

using HS3
using Test
using LiteHF, BAT
using DensityInterface
#@testset "LiteHF" begin
    LiteHF_pydict = load_pyhfjson("./test/sample_staterror_shapesys_LiteHF.json");
    LiteHF_pyhfmodel = build_pyhf(LiteHF_pydict);
    LiteHF_LL = pyhf_logjointof(LiteHF_pyhfmodel)
    LiteHF_likelihood(αs) = PosteriorMeasure(LL(αs))
    LiteHF_posterior = PosteriorDensity(LiteHF_LL, LiteHF_pyhfmodel.priors)
    LiteHF_LL([1., 1., 1., 1., 1.])
    HS3_dict = file_to_dict("./test/sample_staterror_shapesys_HS3.json")
    HS3_specs = HS3.generate_specs(HS3_dict)
    #likelihood = HS3.make_likelihood(HS3_specs.likelihoods[1], HS3_specs.distributions, HS3_specs.data)
    analysis = HS3.make_analyses(HS3_specs.analyses[1], HS3_specs)
    logdensityof(analysis.likelihood, (mu = 1., shape2_bin2 = 1., shape2_bin1 = 1., MCstat_bin1 = 1., MCstat_bin2 = 1.))
    @show bat_findmode(LiteHFposterior).result
    #@test HS3.hello_world() == 42
#end
