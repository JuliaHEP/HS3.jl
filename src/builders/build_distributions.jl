# This file is a part of HS3.jl, licensed under the MIT License (MIT).
"""
generate_distribution(spec)

Generate a specific distribution based on the provided specification.

# Arguments
- `spec`: The specification for the desired distribution.

# Returns
An instance of the requested distribution.

"""
function generate_distribution(spec::DistributionSpec{:exponential_dist})
    Distributions.Exponential(-1/spec.params.c)
end

function generate_distribution(spec::DistributionSpec{:gaussian_dist})
    Distributions.Normal(spec.params.mean, spec.params.sigma)
end

function generate_distribution(spec::DistributionSpec{:lognormal_dist})
    Distributions.LogNormal(spec.params.mean, spec.params.k)
end

function generate_distribution(spec::DistributionSpec{:multinormal_dist}) 
    (Distributions.MvNormal(spec.params.mean, spec.params.covariances), spec.params.x)
end#

function generate_distribution(spec::DistributionSpec{:poisson_dist}) 
    Distributions.Poisson(spec.params.mean)
end

function generate_distribution(spec::DistributionSpec{:argus_dist})
    @error "ARGUS disitrbution is not yet implemented"
end

function generate_distribution(spec::DistributionSpec{:mixture_dist}) 
    Distributions.MixtureModel(collect(spec.params.summands), collect(spec.params.coefficients)/(sum(collect((spec.params.coefficients))))) 
end
function generate_distribution(spec::DistributionSpec{:product_dist}) 
    println(typeof(collect(spec.params.factors)))
    return Distributions.product_distribution(collect(spec.params.factors))
end

function generate_distribution(spec::DistributionSpec{:uniform_dist})
    return Distributions.Uniform(spec.params.min, spec.params.max)
end

function generate_distribution(spec::DistributionSpec{:polynomial_dist}) 
    @warn "Polynomial distribution is not yet implemented"
end

# this function is depricated an should not be used. Use make_histfact(histfact_spec::HistFactorySpec, channel_name::Symbol) directly instead
function generate_distribution(spec::DistributionSpec{:histfactory_dist}) 
    make_histfact(spec, :default_name)
end







