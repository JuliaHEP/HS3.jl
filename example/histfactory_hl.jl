using HS3
using BenchmarkTools, DensityInterface
using JSON3


dict = file_to_dict("./example/rf515_hfJSON.json")
specs = HS3.generate_specs(dict)

analysis = HS3.make_analyses(specs.analyses[1], specs)
mu_1 = logdensityof(analysis.likelihood, (mcstat_bin1 =1, mu=1,  syst3 = 0, syst1 = 0, syst2 =0, mcstat_bin2 = 1, ))
mu_0 = logdensityof(analysis.likelihood, (mu=0, syst2 =0, mcstat_bin2 = 1, syst1 = 0, syst3 = 0,  mcstat_bin1 =1))
mu_2 = logdensityof(analysis.likelihood, (mu = 2, syst2 =0, mcstat_bin2 = 1, syst1 = 0, syst3 = 0, mcstat_bin1 = 1))
mu_n1 = logdensityof(analysis.likelihood, (mcstat_bin1 =1, mu=-1,  syst3 = 0, syst1 = 0, syst2 =0, mcstat_bin2 = 1, ))
mu_3 = logdensityof(analysis.likelihood, (mu=3, syst2 =0, mcstat_bin2 = 1, syst1 = 0, syst3 = 0,  mcstat_bin1 =1))
mu_n2 = logdensityof(analysis.likelihood, (mu = -2, syst2 =0, mcstat_bin2 = 1, syst1 = 0, syst3 = 0, mcstat_bin1 = 1))
mu_n3 = logdensityof(analysis.likelihood, (mu=-3, syst2 =0, mcstat_bin2 = 1, syst1 = 0, syst3 = 0,  mcstat_bin1 =1))

map(x -> x in [1,2,4,5,8,9,10,11] ,[1,3,5,7,9,4])
[x in [1,2,4,5,8,9,10,11]  for x = [1,3,5,7,9,4]]

mu = [-3, -2, -1, 0, 1, 2, 3]
jul = [mu_n3, mu_n2, mu_n1, mu_0, mu_1, mu_2, mu_3]
root = [81.182296, 46.760037 , 28.471124, 19.447233, 16.529263, 17.977606, 22.721993]

[jul .+ root]

HS3.BinnedDataSpec{3} <: HS3.BinnedDataSpec
mu = [0, 1, 2]
root = [16.529, 19.4472, 17.97606]

using LiteHF
using Plots
plot(mu, jul, label="Julia LL", linestyle=:solid)
plot!(mu, root, label="Root LL")
plot!(mu, -root, label="negative Root LL", linestyle=:dash)
plot!(mu, abs.(jul+root), label="Abs Diff Julia & Neg. Root")
xlabel!("mu")
ylabel!("LL")
savefig("LogLikelihoods1.pdf")
abs(a - b)

using BAT, ValueShapes, Distributions
posterior = PosteriorMeasure(analysis.likelihood, NamedTupleDist(mu=Uniform(-3, 5),   mcstat_bin1 = Uniform(0, 2), mcstat_bin2 = Uniform(0, 2), syst1 = Uniform(-5, 5), syst2 = Uniform(-5, 5), syst3 = Uniform(-5, 5)))
bat_findmode(posterior).result
using Plots
plot(samples)
#make analysis directly from file
@benchmark HS3.make_analyses("./example/histfactory_hl.json", "primary_analysis")





a = HS3.specs_to_json_dict(specs, true)
println(a)
JSON3.write(a)
nt = (data = (data1 = s,),)
open("my_new_file.json", "w") do io
    JSON3.pretty(io, a)
end