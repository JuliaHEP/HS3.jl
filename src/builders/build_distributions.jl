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
    Distributions.Exponential(1/spec.params.c)
end

"""
NOTE: At the moment two different functions for normalisation are used: `_normalize_logpdf` and `_norm`. 

- _normalize_logpdf: uses the truncation of a PDF if CDF is defined.
- _norm: if CDF is not defined use numerical integration here (QuadGK)
When creating PDFs as part of likelihoods, both functions are always called although always only one actually normalises the PDF. 

This will undergo changes in the future.
"""

_normalize_logpdf(d::Distributions.Exponential, lo::Real, hi::Real, x::Real) = Distributions.logpdf(Distributions.truncated(d, lo, hi), x)

_norm(d::Distributions.Exponential, lo::Real, hi::Real) = 1.0


function generate_distribution(spec::DistributionSpec{:gaussian_dist})
   Distributions.Normal(spec.params.mean, spec.params.sigma)
end

_normalize_logpdf(d::Distributions.Normal, lo::Real, hi::Real, x::Real) = Distributions.logpdf(Distributions.truncated(d, lo, hi),x)

_norm(d::Distributions.Normal, lo::Real, hi::Real) = 1.0

generate_distribution(spec::DistributionSpec{:normal_dist}) = generate_distribution(spec::DistributionSpec{:gaussian_dist})

function generate_distribution(spec::DistributionSpec{:lognormal_dist})
    Distributions.LogNormal(spec.params.mu, spec.params.sigma)
end

_normalize_logpdf(d::Distributions.LogNormal, lo::Real, hi::Real, x::Real) = Distributions.logpdf(Distributions.truncated(d, lo, hi),x)

_norm(d::Distributions.LogNormal, lo::Real, hi::Real) = 1.0

function generate_distribution(spec::DistributionSpec{:multinormal_dist}) 
    (Distributions.MvNormal(spec.params.mean, spec.params.covariances), spec.params.x)
end#

function generate_distribution(spec::DistributionSpec{:poisson_dist}) 
    Distributions.Poisson(spec.params.mean)
end

_normalize_logpdf(d::Distributions.Poisson, lo::Real, hi::Real, x::Real) = Distributions.logpdf(Distributions.truncated(d, lo, hi),x)

_norm(d::Distributions.Poisson, lo::Real, hi::Real) = 1.0

include("../definitions/Argus.jl") #Not yet supported in RooFit

function generate_distribution(spec::DistributionSpec{:argus_dist})
    Argus(spec.params.resonance, spec.params.slope, spec.params.power)
end

_normalize_logpdf(d::Argus, lo::Real, hi::Real, x::Real) = Distributions.logpdf(Distributions.truncated(d, lo, hi),x)

_norm(d::Argus, lo::Real, hi::Real) = 1.0

include("../definitions/mixture_hack.jl")

function generate_distribution(spec::DistributionSpec{:mixture_dist})
    MixtureHack(collect(spec.params.summands), [spec.params.coefficients...], spec.params.extended) 
end

function _norm(d::MixtureHack, lo::Float64, hi::Float64)
    return quadgk(x-> pdf(d, (lo+(hi-lo)*x)) * (hi-lo), lo, hi, atol=1e-15, order=100, rtol=1e-15)[1] 
end

function _normalize_logpdf(d::MixtureHack, lo::Float64, hi::Float64, y::Float64) 
    return logpdf(d,y)
end

include("../definitions/product_hack.jl")

function generate_distribution(spec::DistributionSpec{:product_dist}) 
    return ProductHack(collect(spec.params.factors))
end

function _norm(d::ProductHack, lo::Float64, hi::Float64)
    quadgk(x-> pdf(d, (lo+(hi-lo)*x)) * (hi-lo), lo, hi, atol=1e-15, order=100)[1]
end

function _normalize_logpdf(d::ProductHack, lo::Float64, hi::Float64, y::Float64) 
    return logpdf(d,y)
end

function generate_distribution(spec::DistributionSpec{:uniform_dist})
    return Distributions.Uniform(spec.params.min, spec.params.max)
end

function generate_distribution(spec::DistributionSpec{:polynomial_dist}) #is implemented, needs to be included here -> tests!
    @warn "Polynomial distribution is not yet implemented"
end

# this function is depricated an should not be used. Use make_histfact(histfact_spec::HistFactorySpec, channel_name::Symbol) directly instead
function generate_distribution(spec::DistributionSpec{:histfactory_dist}) 
    make_histfact(spec, :default_name)
end

include("../definitions/Bernstein.jl")
function generate_distribution(spec::DistributionSpec{:bernstein_poly_dist})
    return BernsteinPolynomial(collect(spec.params.coefficients))
end
#
function _normalize_logpdf(d::BernsteinPolynomial, lo::Float64, hi::Float64, y::Float64) 
    return logpdf(d,y)  
end
#
_norm(d::BernsteinPolynomial, lo::Real, hi::Real) = return quadgk(x-> pdf(d, (lo+(hi-lo)*x)) * (hi-lo), lo, hi, atol=1e-15, order=100, rtol=1e-15)[1]

function _find_val_entry(array::Vector{Any})
    for entry in array
        if typeof(entry).name.wrapper == Val
            return _val_content(entry)
        end
    end
    return nothing  # Return nothing if there are no Val entries
end

include("../definitions/Generic.jl")

function generate_distribution(spec::DistributionSpec{:generic_dist}) 
    x = _find_val_entry(collect(spec.params.var))
    Generic(spec.params.expression, collect(spec.params.var), Symbol(x))
end

function _normalize_logpdf(d::Generic, lo::Float64, hi::Float64, y::Float64) 
    return HS3.logpdf(d, y)
end
#
function _norm(d::Generic, lo::Real, hi::Real) 
    return quadgk(x-> pdf(d, (lo+(hi-lo)*x)) * (hi-lo), lo, hi, atol=1e-15, order=100, rtol=1e-15)[1]
end


include("../definitions/CrystalBall.jl")

function generate_distribution(spec::DistributionSpec{:crystalball_dist}) 
    return CrystalBall(spec.params.m0, spec.params.sigma, spec.params.alpha, spec.params.n)
end

function _norm(d::CrystalBall, lo::Float64, hi::Float64)
    quadgk(x-> pdf(d, (lo+(hi-lo)*x)) * (hi-lo), lo, hi, atol=1e-15, order=7)[1]
end

function _normalize_logpdf(d::CrystalBall, lo::Float64, hi::Float64, y::Float64) 
    return logpdf(d,y)
end
 
include("../definitions/SplitGaussian.jl")

function generate_distribution(spec::DistributionSpec{:bifurkated_gaussian_dist})
    return SplitGaussian(spec.params.mean, spec.params.sigmaL, spec.params.sigmaR)
end

function _norm(d::SplitGaussian, lo::Real, hi::Real)
    quadgk(x-> pdf(d, (lo+(hi-lo)*x)) * (hi-lo), lo, hi, atol=1e-15, order=9)[1]
end

function _normalize_logpdf(d::SplitGaussian, lo::Real, hi::Real, y::Real) 
    return logpdf(d,y)
end
