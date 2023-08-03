using HS3
using NamedTupleTools
using BenchmarkTools


@elapsed a = (k = HS3.ProductDomainSpec(1, 1+1), l = HS3.ProductDomainSpec(4, 1+1),)
@elapsed b = (; zip([:name], [a])...)
@elapsed c = (; zip([:name], [_test_domain()])...)
@elapsed HS3.ProductDomainSpec(1, 1+1)

function _test_domain()
    #domains = [HS3.ProductDomainSpec((factors = [i, i+1, "a"], ))  for i in range(1, 800)]
    domains = [HS3.ProductDomainSpec(i, 0.1*i+1) for i in range(1, 1600)]
    names = [Symbol(i) for i in range(1, 1600)]
    return (; zip(names, domains)...)
end
DistributionSpec{:product_dist}((factors = data.factors, )) 

HS3.DistributionSpec{:uniform_dist}((min = 1, max = 2,))
HS3.DistributionSpec(Val{:Uniform}, (a = 1, b = 2,))
namedtuple([:a], [HS3.ProductDomainSpec(1, 1)])
function _test_domain2()
    nt = (;)
    domains1 = [HS3.DistributionSpec{:product_dist}((factors = [i, i+1, "a"], ))  for i in range(1, 800)]
    domains2 = [HS3.DistributionSpec{:uniform_dist}((min = i, max = i+1,)) for i in range(1, 800)]
    names = [Symbol(i) for i in range(1, 1600)]
    domains = vcat(domains1, domains2) 
    return (; zip(names, domains)...)
    #for i in range(1, 1600)
    #   nt = merge(nt, (Symbol(i) => HS3.ProductDomainSpec(i, i+1) ,))
    #   #HS3.ProductDomainSpec(i, i+1)
    #end
    #return nt
end
names = [Symbol("a$i") for i in range(1, 2)]
domains = [HS3.ProductDomainSpec(i, i+1)  for i in range(1, 2)]
@elapsed namedtuple(names, domains)

names = [Symbol("a$i") for i in range(1, 1600)]
domains = [HS3.ProductDomainSpec(i+1, i+2)  for i in range(1, 1600)]
@elapsed namedtuple(names, domains)
@elapsed _test_domain2()

typeof([HS3.ProductDomainSpec(1, 1)])


function namedtuple2(namesforvalues::Vector{Symbol}, valuesfornames::Vector{HS3.ProductDomainSpec})
    length(namesforvalues) == length(valuesfornames) || throw(ErrorException("lengths must match"))
    return (; zip(namesforvalues, valuesfornames)...,)
end
    @elapsed @eval(domains = [HS3.ProductDomainSpec(i, i+1) for i in range(1, 1600)])

nt = (;)
@benchmark nt = (;) ; for i in range(1, 1600)
       nt = merge(nt, (Symbol(i) => HS3.ProductDomainSpec(i, i+1) ,))
       #HS3.ProductDomainSpec(i, i+1)
    end

@elapsed @eval(_test_domain2())
@elapsed @eval(_test_domain())
nt