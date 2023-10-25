using HS3
using DensityInterface
using BenchmarkTools
using CSV

using Distributions
include("ll_help.jl")
##lognormal 
dict = file_to_dict("./example/test_gg/test_gg_binned_lognormal.json")
specs = HS3.generate_specs(dict)



analysis = HS3.make_analyses(specs.analyses[1], specs)
logdensityof(analysis[:likelihood], (mu =1, sigma=1,))
ex = explore_values([:mu,], (mu = LinRange(0.0:0.1:3), sigma=LinRange(0.1:0.1:1.5)), analysis[:likelihood], (mu = 1, sigma=1) )

write_to_csv(ex, "./example/test_gg/test_gg_binned_lognormal.csv")



## uniform
dict = file_to_dict("./example/test_gg/test_gg_binned_uniform.json")
specs = HS3.generate_specs(dict)



analysis = HS3.make_analyses(specs.analyses[1], specs)
logdensityof(analysis[:likelihood], (mu =1, sigma=5,))
ex = explore_values([:mu,:sigma], (mu = LinRange(0.0:0.1:3), sigma=LinRange(4:0.1:7)), analysis[:likelihood], (mu = 1, sigma=4) )

write_to_csv(ex, "./example/test_gg/test_gg_binned_uniform.csv")



### gauss:




dict = file_to_dict("./example/test_gg/test_gg_binned_gauss.json")
specs = HS3.generate_specs(dict)



analysis = HS3.make_analyses(specs.analyses[1], specs)
logdensityof(analysis[:likelihood], (mu =1, sigma=1,))
ex = explore_values([:mu,:sigma], (mu = LinRange(0.0:0.1:3), sigma=LinRange(0.1:0.1:1.5)), analysis[:likelihood], (mu = 1, sigma=1) )

write_to_csv(ex, "./example/test_gg/csv/test_gg_binned_gauss.csv")



## Test Exp 

dict = file_to_dict("./example/test_gg/test_gg_binned_exp.json")
specs = HS3.generate_specs(dict)



analysis = HS3.make_analyses(specs.analyses[1], specs)
logdensityof(analysis[:likelihood], (mu =1,))
ex = explore_values([:mu,], (mu = LinRange(0.1:0.1:3),), analysis[:likelihood], (mu = 1,) )

write_to_csv(ex, "./example/test_gg/test_gg_binned_exp.csv")



#unbinned gauss
dict = file_to_dict("./example/test_gg/test_gg_unbinned_gauss.json")
specs = HS3.generate_specs(dict)



analysis = HS3.make_analyses(specs.analyses[1], specs)
logdensityof(analysis[:likelihood], (mu =1, sigma=1,))
ex = explore_values([:mu, :sigma], (mu = LinRange(0.0:0.1:3), sigma=LinRange(0.2:0.1:1.5)), analysis[:likelihood], (mu = 1, sigma=1) )

write_to_csv(ex, "./example/test_gg/csv/test_gg_unbinned_gauss.csv")


#weighted unbinned gauss
dict = file_to_dict("./example/test_gg/test_gg_unbinned_gauss_weighted.json")
specs = HS3.generate_specs(dict)



analysis = HS3.make_analyses(specs.analyses[1], specs)
logdensityof(analysis[:likelihood], (mu =1, sigma=1,))
ex = explore_values([:mu, :sigma], (mu = LinRange(0.0:0.1:3), sigma=LinRange(0.2:0.1:1.5)), analysis[:likelihood], (mu = 1, sigma=1) )

write_to_csv(ex, "./example/test_gg/csv/test_gg_unbinned_gauss_weighted.csv")





#unbinned exponential
dict = file_to_dict("./example/test_gg/test_gg_unbinned_exponential.json")
specs = HS3.generate_specs(dict)



analysis = HS3.make_analyses(specs.analyses[1], specs)
logdensityof(analysis[:likelihood], (c =1,))
ex = explore_values([:c], (c = LinRange(0.0:0.1:10),), analysis[:likelihood], (c =1,) )

write_to_csv(ex, "./example/test_gg/csv/test_gg_unbinned_exponential.csv")





#unbinned lognormal
dict = file_to_dict("./example/test_gg/test_gg_unbinned_lognormal.json")
specs = HS3.generate_specs(dict)



analysis = HS3.make_analyses(specs.analyses[1], specs)
logdensityof(analysis[:likelihood], (mu =1, sigma=1,))
ex = explore_values([:mu, :sigma], (mu = LinRange(0.0:0.1:3), sigma=LinRange(0.2:0.1:1.5)), analysis[:likelihood], (mu = 1, sigma=1) )

write_to_csv(ex, "./example/test_gg/csv/test_gg_unbinned_lognormal.csv")

#unbinned crystalBall_dist
dict = file_to_dict("./example/test_gg/test_gg_unbinned_crystalball.json")
specs = HS3.generate_specs(dict)



analysis = HS3.make_analyses(specs.analyses[1], specs)

@benchmark logdensityof(analysis[:likelihood], (m0 =1, sigma=1, alpha=1, n=1.1))
ex = explore_values([:m0, :sigma, :alpha, :n], (m0 = LinRange(0.0:0.1:5), sigma=LinRange(0.1:0.1:1.5), alpha=LinRange(0.1:0.1:1.5),n=LinRange(1.1:0.1:20)), analysis[:likelihood], (m0 = 1, sigma=1, alpha=1, n=10) )

write_to_csv(ex, "./example/test_gg/csv/test_gg_unbinned_crystalball.csv")


#unbinned bernstein_poly_dist
dict = file_to_dict("./example/test_gg/test_gg_unbinned_bernstein.json")
specs = HS3.generate_specs(dict)

using QuadGK
quadgk(x->exp(x * 0.000000000000001), 0, 10)

analysis = HS3.make_analyses(specs.analyses[1], specs)

logdensityof(analysis[:likelihood], (c1 =0.2, c2=0, c3=0, c4=0.0, c5=0,)) 
ex = explore_values([:c1, :c2, :c3, :c4, :c5], (c1 = LinRange(0.0:0.1:10), c2 = LinRange(0.0:0.1:10),c3 = LinRange(0.0:0.1:10),c4 = LinRange(0.0:0.1:10),c5 = LinRange(0.0:0.1:10),), analysis[:likelihood], (c1 = 1, c2=1, c3=1, c4=1, c5=1) )
ex = explore_values([:c1], (c1 = LinRange(0.0:0.1:10),), analysis[:likelihood], (c1 = 1, c2 =0, c3=0, c4=0, c5=0) )
write_to_csv(ex, "./example/test_gg/csv/test_gg_unbinned_bernstein.csv")





#unbinned bifurkated gauss
dict = file_to_dict("./example/test_gg/test_gg_unbinned_split.json")
specs = HS3.generate_specs(dict)



analysis = HS3.make_analyses(specs.analyses[1], specs)
@benchmark logdensityof(analysis[:likelihood], (mu =1, sigmaL=1, sigmaR=2.,))
ex = explore_values([:mu, :sigmaL, :sigmaR], (mu = LinRange(0.0:0.1:5), sigmaL=LinRange(0.2:0.1:2), sigmaR=LinRange(0.2:0.1:2)), analysis[:likelihood], (mu = 1., sigmaL=1., sigmaR=1.5,) )

write_to_csv(ex, "./example/test_gg/csv/test_gg_unbinned_split.csv")


#unbinned product
dict = file_to_dict("./example/test_gg/test_gg_unbinned_product.json")
specs = HS3.generate_specs(dict)



analysis = HS3.make_analyses(specs.analyses[1], specs)
logdensityof(analysis[:likelihood], (m0 =1, m1=1, sigma0=1, sigma1=2.,))
ex = explore_values([:m0, :m1, :sigma0, :sigma1], (m0 = LinRange(0.0:0.1:6), m1 = LinRange(0.0:0.1:6), sigma0=LinRange(0.2:0.1:3), sigma1=LinRange(0.2:0.1:3)), analysis[:likelihood], (m0 = 1., m1=2, sigma0=1., sigma1=1.,) )

write_to_csv(ex, "./example/test_gg/csv/test_gg_unbinned_product.csv")


#unbinned mixture model
dict = file_to_dict("./example/test_gg/test_gg_unbinned_mixture_extended.json")
specs = HS3.generate_specs(dict)



analysis = HS3.make_analyses(specs.analyses[1], specs)
logdensityof(analysis[:likelihood], (m0 =1, m1=2, sigma0=1, sigma1=2., a0 = 1.0, a1 = 0.8))
ex = explore_values([:m0, :m1, :sigma0, :sigma1, :a0, :a1 ], (m0 = LinRange(0.1:0.1:6), m1 = LinRange(0.1:0.1:6), sigma0=LinRange(0.2:0.1:3), sigma1=LinRange(0.2:0.1:3), a0=LinRange(0.1:0.1:6), a1=LinRange(0.1:0.1:6)), analysis[:likelihood], (m0 = 1., m1=2, sigma0=1., sigma1=2., a0 = 1, a1 = 0.8) )

write_to_csv(ex, "./example/test_gg/csv/test_gg_unbinned_mixture_extended.csv")


#unbinned argus
dict = file_to_dict("./example/test_gg/test_gg_unbinned_argus.json")
specs = HS3.generate_specs(dict)



analysis = HS3.make_analyses(specs.analyses[1], specs)
logdensityof(analysis[:likelihood], (m0 =8, slope=1, n=0.5,))
ex = explore_values([:m0, :slope, :n], (m0 = LinRange(4:0.1:10), slope=LinRange(0.2:0.1:2), n=LinRange(0.2:0.1:2)), analysis[:likelihood], (m0 = 7, slope=1., n=1.5,) )

write_to_csv(ex, "./example/test_gg/csv/test_gg_unbinned_argus.csv")







a = truncated(Normal(1, 1), 0, 1)
logpdf(a, 1)
b = Normal(1, 1)
logpdf(b, 116)

b = Normal(1, 0.1)
pdf(b, 2.5)
a = Poisson(5.530709549844416e-49)
c = HS3.RelaxedPoisson(5.530709549844416e-49)
d = HS3.ContinousPoisson(5.530709549844416e-49)
logpdf(d, 3.2)
logpdf(c, 3)
logpdf(a, 3)
x =[0, 0.6, 1, 2, 3, 4, 5,6 ,7 ,8 , 9 , 10]
LinRange(0:0.5:10)


