using HS3
example_dict = file_to_dict("./example/test.json")
dict = file_to_dict("./example/test.json")
specs = HS3.generate_specs(dict)
analysis = HS3.make_analyses(specs.analyses[1], specs)
a = logdensityof(analysis.likelihood, (param_c = -2, param_mean = 0, param_sigma =1))
HS3.topological_sort(merge(specs.distributions, specs.functions))
using DensityInterface
a = logdensityof(analysis.likelihood)
for expr in code
    println(expr)
    if expr isa Expr && expr.head == :call
        try 
            operation = push!(operation, expr.args[1].name)
        catch end
        for arg in expr.args
            if typeof(arg) == QuoteNode
                params_names = push!(param_names, arg.value)
            end
        end
    end
end

dump(code[1].args)
operation
param_names
param_names = []
operation = []
code = (code_lowered(a)[1].code)
dump(a)
a = a.likelihood_functions.contents[1]
(dump(a))




b = logdensityof(analysis.likelihood, (mu=1, syst2 =0, mcstat_bin2 = 1, syst1 = 0, syst3 = 0,  mcstat_bin1 =1))
c = logdensityof(analysis.likelihood, (mu = -2.4, syst2 = -0.009487064624933319, mcstat_bin2 = 0.9999950810630367, syst1 = 0.9978812987492401, syst3 = 1.0719092640875414, mcstat_bin1 = 1.0000003487684337))

using LiteHF

abs(a - b)
dump(analysis.likelihood._log_f.likelihood_functions.contents[1])
dump(analysis.likelihood)
using BAT, ValueShapes, Distributions
posterior = PosteriorMeasure(analysis.likelihood, NamedTupleDist(param_c = Uniform(-3, 1), param_mean = Uniform(0, 1), param_sigma = Uniform(0, 1)))
bat_findmode(posterior).result

using BenchmarkTools
example_dict[:functions]
example_dict[:distributions][1].samples
HS3._typed_content(HS3._typed_content(l))
a = HS3.generate_specs(example_dict)
b = merge(a.distributions, a.functions)
HS3.topological_sort(b)
typeof(b.params.summands)
f = Base.Fix1(HS3.substitute_variables, b)
f
f((sig1 = 2, background = 4, test_abc=3,))
HS3.substitute_variables(b[end], (sig1 = 2, background = 4, :test_abc => 4, gen = 3, func1 =4))
c = HS3.make_functional(b[end])
using Distributions
c((test_abc = 4, sig1 = Normal(), func1=4, background = Exponential()))
HS3.generate_specs("./example/test.json")

using BenchmarkTools

HS3.AxesSpec( min=3, nbins=5)
HS3.AxesSpec(; (max=4, min=4)...)
validate_file("./example/test.json")
as
@btime Bool <: Integer
typeof(2) <:Integer

a = HS3.make_parameterpointspecs(example_dict[:parameter_points]) 
#HS3.@store c = a[:starting_values].par_1 




#c = a[:starting_values].par_1 
b = a[:starting_values]
#typeof(b) <: NamedTuple{(), Tuple{HS3.ParameterPointSpec}}
a = HS3.make_funcspecs(example_dict[:functions])
HS3.topological_sort(merge(a, b))
NamedTuple{(:a, :b, :c), Tuple{Int64, Int64, Int64}} <: NamedTuple{(...), <:Tuple{<:Any, <:Any, <:Any}}
a[3]
b = HS3.make_functional(a[3])
b((a=1, b=2, c=3,))
HS3.generate_function(a[3], (a=1, b=2, c=3,))
k = "((a-2*b)+c)"
split(k, ["-", "+"]) 

ex = :(((((achkb) - 2))*b) + c)
ex =:(ttW_Bin_011_mu-2*ttW_Bin_012_mu+ttW_Bin_013_mu)


likelihood_specs = HS3.make_likelihood_specs(example_dict[:likelihoods])
HS3.make_likelihood(likelihood_specs[:main_likelihood], merge(a,b))

a = HS3.make_likelihood(likelihood_specs, merge(HS3.make_distspecs(example_dict[:distributions]), HS3.make_funcspecs(example_dict[:functions])), HS3.make_dataspecs(example_dict[:data]))
a[1]
a[1]
#HS3.make_functional(merge(HS3.make_distspecs(example_dict[:distributions]), HS3.make_funcspecs(example_dict[:functions])))
using DensityInterface
logdensityof(a[1], (;))
using ValueShapes
using Distributions
NamedTupleDist((a  = Normal()), )
specs = HS3.make_analysesspecs(example_dict[:analyses])
b = HS3.make_distspecs(example_dict[:distributions])
a = HS3.topological_sort(specs)
first.(last.(a))
typeof(a)
for x in a 
    print(x[2][1])
end
specs[:model1]
k = HS3.make_functional(specs[:background], a, 2)
k((;))
using Distributions
Distributions.ProductDistribution <:Distribution{<:ArrayLikeVariate{Any}}

func_specs = HS3.make_funcspecs(example_dict[:functions])
a = HS3.topological_sort(merge(dist_specs, func_specs))
k = HS3.make_functional(func_specs[:func1], a, 2)
@btime k((;))


l = HS3.make_functional(specs[:func1])
l((;))
l = Distributions.MixtureModel([Normal(0, 1), Exponential(0.6)], [0.2, 0.8] ) 
m = product_distribution(Normal(0.1, 3), Exponential(2))
product_distribution(Exponential(2), Normal(0, 0.5), Normal(0.2, 1), m)
Distributions._logpdf(product_distribution([Exponential(2), Normal(0, 0.5), l]), [2, 2, 2])
filter((x) -> x[1] !=0, a)
a = HS3.make_markov_kernel(specs[1])
a((a=1,))
l = p -> begin 
    a(p)
    print(a)
end
l((a=1, n=2, ))
Int64 <:Int64
dist = (background = make_markov_kernel(specs[:background]), sig1 = make_markov_kernel(specs[:sig1]),)
data = HS3.make_dataspecs(example_dict[:data])
lspecs = HS3.make_likelihood_specs(example_dict[:likelihoods])
ll = HS3.make_likelihood(lspecs[:bkg_likelihood], dist, data)
a = dist[:background]
isa(a.inner.x, HS3.DistributionSpec{:poisson_dist})
using Plots
using BAT
using Distributions 
logpdf(Normal(2, 1), 0)
using LiteHF
LiteHF.Shapesys(1, [1], [2])
Normsys(1.1, 2.0)
const v_bg = [30,19,9,4]
const variations = [1,2,3,3]
v_bg .+ variations
const bkgmodis =[
                 Histosys(zeros(length(v_bg)), v_bg .+ variations, v_bg .- variations),
                 Normsys(1.1, 0.9)
        ]
Normfactor()
misc = Dict()
LiteHF.build_channel(Dict(example_dict[:distributions][4]); misc)
struct a 
    name::String 
end
eachindex([])
typeof([])# <: AbstractArray{<:Any}
k = a("gsdha")
k[:name]
plot(ll)
using DensityInterface
logdensityof(ll, ((a= 12, n=0.2,)))
l((a=2,))
domain_specs = HS3.make_domainspecs(example_dict[:domains])
HS3.generate_prior_from_domain(domain_specs[:default_domain])
parameter_point_specs = HS3.make_parameterpointspecs(example_dict[:parameter_points])
axes_specs = HS3.make_axesspecs(example_dict[:data][1][:axes])
data_spec1 = HS3.generate_binneddataspec(example_dict[:data][1])
HS3.generate_hist(data_spec1)
HS3.make_dataspecs(example_dict[:data])
for i in data_spec1.contents
    print(i)
end
typeof(axes_specs) <: NamedTuple{<:Any, <:Tuple{HS3.AxesSpec}}
1 isa Integer
domain_specs.default_domain.coef_sig.max
#domain_specs = HS3.make_distspecs(example_dict[:distributions])
domains = HS3.generate_domain(domain_specs.default_domain)
domains.coef_sig
s = make_markov_kernel(HS3.make_distspecs(example_dict[:distributions]).sig1)
s((a = 1,))

a = (a = 2, b = 4.0)
typeof(a)
Bool <: Int64
Bool <: Real

NamedTuple{(:value,), Tuple{Int64}} <: NamedTuple{S, <: Tuple{Real}} where S
a = size([[1,2]  [4,5] ])
a[1] == a[2] && a[1] == 2
[1, 2, 4]

B = Array{Float64}(undef, 4)
typeof(B)
using JSON3
JSON3.Array{Int64, Vector{UInt8}, SubArray{UInt64, 1, Vector{UInt64}, Tuple{UnitRange{Int64}}, true}} <: AbstractArray{<:Number, 1}
a= 1
length(a)
size(a)[1]
typeof([[2] [2]])
size([[2,4] [3, 5]])
hcat([[2.0,4], [3, 5]]...)
dict = Dict(:type => Val(Symbol("gaussian_uncertainty")), :sigma =>[2, 3], :correlation => [[2,4], [3, 5]])
HS3.generate_uncertaintyspec(dict)


# Tuple{2, 2} #
typeof(hcat([[[1] [4]],[[2.] [9]]]...))


a = [[1, 2], [4, 3], [2., 2], [9., 1]] 
typeof(a) <: Array{Array{Float64, 1}, 1} <: AbstractArray{Float64}
a[1]
collect(reduce(hcat, a'))
typeof(reduce(vcat, a')) <: Matrix <: Array{2}
#hvcat((1,1), a)
length(a)
Matrix{Number} <: AbstractArray{Number, 2}

isa(2, Number)

using StatsBase

a = StatsBase.Histogram([1, 3, 5], [1., 2.])
StatsBase.binindex(a, 10)
fit

using JSON3
JSON3.Array{JSON3.Array, Vector{UInt8}, SubArray{UInt64, 1, Vector{UInt64}, Tuple{UnitRange{Int64}}, true}} <: AbstractArray{<:AbstractArray{<:Any}, 1}

NamedTuple{(:variable1, :variable2), Tuple{Float64, Float64}} <: NamedTuple{<:Any, <:NTuple{N, Float64}} where N
#NamedTuple{<:Any, NTuple{N, AxesSpec}} where N
Tuple{Float16, Float16} <: NTuple{N, Float16} where {N::Integer}


using Distributions

dist = MixtureModel([Normal(3, 0.4), Exponential(0.5)], [0.2, 0.8])

data = rand(dist, 1000000)
using StatsBase
hist = append!(Histogram(0:0.1:10), data)
print(hist.weights)