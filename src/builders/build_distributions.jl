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

_normalize_logpdf(d::Distributions.Exponential, lo::Real, hi::Real, x::Real) = Distributions.logpdf(Distributions.truncated(d, lo, hi), x)

_norm(d::Distributions.Exponential, lo::Real, hi::Real) = 1.0

#struct scuffed_exp <: Distributions.ContinuousUnivariateDistribution
#    c::Float64
#end 
#
#Distributions.pdf(d::scuffed_exp, x::Real) =  d.c * exp(-x * d.c) 
#
##Distributions.logpdf(d::scuffed_exp, x::Real) = -d.c - x 
#
#Distributions.cdf(d::scuffed_exp, x::Real) = 1. - exp(-d.c * x)

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

include("../definitions/Argus.jl")

function generate_distribution(spec::DistributionSpec{:argus_dist})
    Argus(spec.params.resonance, spec.params.slope, spec.params.power)
end

_normalize_logpdf(d::Argus, lo::Real, hi::Real, x::Real) = Distributions.logpdf(Distributions.truncated(d, lo, hi),x)

_norm(d::Argus, lo::Real, hi::Real) = 1.0

include("../definitions/mixture_hack.jl")

function generate_distribution(spec::DistributionSpec{:mixture_dist})
    #coff = spec.params.coefficients
    #if spec.params.extended == false 
    #    push!(coff, 1-sum(coff))
    #end
    MixtureHack(collect(spec.params.summands), collect(spec.params.coefficients), spec.params.extended) 
end

function _norm(d::MixtureHack, lo::Float64, hi::Float64)
    #@fastmath vegas((x,f) -> f[1] = pdf(d, (lo+(hi-lo)*x[1])) * (hi-lo))[1][1]
    #s = 0.0
    #for i in d.summands
    #    s +=  quadgk(x-> pdf(i,(lo+(hi-lo)*x)) * (hi-lo), lo, hi, atol=1e-15, order=10)[1]
    #end
    ##s = 1.0/s
    #println(s)
    return quadgk(x-> pdf(d, (lo+(hi-lo)*x)) * (hi-lo), lo, hi, atol=1e-15, order=10000000, rtol=1e-15)[1]
    #s
    #1.0
end

function _normalize_logpdf(d::MixtureHack, lo::Float64, hi::Float64, y::Float64) 
    return logpdf(d,y)
end

include("../definitions/product_hack.jl")
function generate_distribution(spec::DistributionSpec{:product_dist}) 
    #println(typeof(collect(spec.params.factors)))
    return ProductHack(collect(spec.params.factors))
end

function _norm(d::ProductHack, lo::Float64, hi::Float64)
    #@fastmath vegas((x,f) -> f[1] = pdf(d, (lo+(hi-lo)*x[1])) * (hi-lo))[1][1]
    quadgk(x-> pdf(d, (lo+(hi-lo)*x)) * (hi-lo), lo, hi, atol=1e-15, order=7)[1]
end

function _normalize_logpdf(d::ProductHack, lo::Float64, hi::Float64, y::Float64) 
    return logpdf(d,y)
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

struct BernsteinPolynomial <: Distributions.ContinuousUnivariateDistribution
    c::Vector{Float64}
end

function generate_distribution(spec::DistributionSpec{:bernstein_poly_dist})
    return BernsteinPolynomial(collect(spec.params.coefficients))
end

function _normalize_logpdf(d::BernsteinPolynomial, lo::Float64, hi::Float64, y::Float64) 
    #norm =  vegas((x,f) -> f[1] = pdf(d, (lo+(hi-lo)*x[1])) * (hi-lo))[1][1]
    return logpdf(d,y, hi, lo) #- log(norm)
end

_norm(d::BernsteinPolynomial, lo::Real, hi::Real) = 1.0

bernstein_basic_polynomials(n::Integer, x::Real) = [(binomial(n, k) * x^k * (1- x)^(n-k)) for k in 0:n]

highest_order(d::BernsteinPolynomial) = findlast(x -> x != 0, d.c) -1

lowest_order(d::BernsteinPolynomial) = findfirst(x -> x != 0, d.c) -1

function pdf(d::BernsteinPolynomial, x::Real, xmax, xmin)  
    #xmax = 1
    #xmin =0
    x = (x-xmin)/(xmax-xmin)
    #println(d.c, length(d.c)-1)
    #println(x, " here ", sum(bernstein_basic_polynomials(length(d.c)-1, x) )) 
    #sum(d.c .* bernstein_basic_polynomials(length(d.c)-1, x) ) 
    res = 0.0
    for i in 0:length(d.c)-1
        res += x^i * binomial(length(d.c)-1, i+1) * d.c[i+1] * (1-x)^(5-i)
        println(i, " ",res)
        #x *= x
    end
    #res += x * d.c[end] 
    println("y ", res)
    res *_n_val(d)
end
logpdf(d::BernsteinPolynomial, x::Real, xmax, xmin) = log(pdf(d::BernsteinPolynomial, x::Real, xmax, xmin))
#TODO: not valid for all intervals, might throws errors
#function logpdf(d::BernsteinPolynomial, x::Real) 
#    #println(bernstein_basic_polynomials(highest_order(d), x) )
#    log(sum(d.c .* bernstein_basic_polynomials(highest_order(d), x) )) -  log(_n_val(d))
#end

function _n_val(d::BernsteinPolynomial) 
    abs(length(d.c) -1 / sum(d.c))
end


#struct CrystalBall <: Distributions.ContinuousUnivariateDistribution
#    x0::Float64    
#    sigma::Float64
#    alpha::Float64
#    n::Float64
#end
#
include("../definitions/CrystalBall.jl")

function generate_distribution(spec::DistributionSpec{:crystalball_dist}) 
    return CrystalBall(spec.params.m0, spec.params.sigma, spec.params.alpha, spec.params.n)
end

function _norm(d::CrystalBall, lo::Float64, hi::Float64)
    #@fastmath vegas((x,f) -> f[1] = pdf(d, (lo+(hi-lo)*x[1])) * (hi-lo))[1][1]
    quadgk(x-> pdf(d, (lo+(hi-lo)*x)) * (hi-lo), lo, hi, atol=1e-15, order=7)[1]
end

function _normalize_logpdf(d::CrystalBall, lo::Float64, hi::Float64, y::Float64) 
    return logpdf(d,y)
end
#
#
#
#function pdf(d::CrystalBall, x::Float64)
#    (x-d.x0)/d.sigma > - d.alpha ?  exp(- (x - d.x0)^2/(2 * d.sigma^2)) * _n_val(d) : _a_val(d) * (_b_val(d) - (x - d.x0)/d.sigma)^(-d.n) *_n_val(d)
#end
#
##normalisation
#function _n_val(d::CrystalBall)
#    return 1/(d.sigma*(sqrt(pi/2) *( 1+ erf(abs(d.alpha)/sqrt(2))) + d.n /abs(d.alpha) * 1/(d.n-1) * exp(- abs(d.alpha)^2 / 2)))
#end
#
##helper functions for often used values
#function _a_val(d::CrystalBall)
#    return (d.n/abs(d.alpha))^(d.n) * exp(- abs(d.alpha)^2 / 2)
#end
#
#function _b_val(d::CrystalBall)
#    return d.n / abs(d.alpha) - abs(d.alpha)
#end
#
#function logpdf(d::CrystalBall, x::Float64)
#    (x-d.x0)/d.sigma > - d.alpha ? (- (x - d.x0)^2/(2 * d.sigma^2)) + _log_n_val(d) : _log_a_val(d) + (-d.n) * log(_b_val(d) - (x - d.x0)/d.sigma) + _log_n_val(d)
#end
#
#function _log_a_val(d::CrystalBall)
#    return (d.n) * log(d.n/abs(d.alpha)) + (- abs(d.alpha)^2 / 2)
#end
#
#function _log_n_val(d::CrystalBall)
#    return - (log(d.sigma) + log(sqrt(pi/2) *( 1+ erf(abs(d.alpha)/sqrt(2))) + d.n /abs(d.alpha) * 1/(d.n-1) * exp(- abs(d.alpha)^2 / 2)))
#end
#
#
#struct SplitGaussian <: Distributions.ContinuousUnivariateDistribution
#    μ::Float64
#    σ1::Float64 
#    σ2::Float64
#end
include("../definitions/SplitGaussian.jl")
function generate_distribution(spec::DistributionSpec{:bifurkated_gaussian_dist})
    return SplitGaussian(spec.params.mean, spec.params.sigmaL, spec.params.sigmaR)
end

#function pdf(d::SplitGaussian, x::Real) 
#    if x < d.μ
#        return exp(-(x-d.μ)^2 / (2 * d.σ1^2)) * _n_val(d)
#    else 
#        return exp(-(x-d.μ)^2 / (2 * d.σ2^2)) * _n_val(d)
#    end
#end
#
#function logpdf(d::SplitGaussian, x::Real) 
#    if x < d.μ
#        return (-(x-d.μ)^2 / (2 * d.σ1^2)) + log(_n_val(d))
#    else 
#        return (-(x-d.μ)^2 / (2 * d.σ2^2)) + log(_n_val(d))
#    end
#end
#
#_n_val(d::SplitGaussian) = sqrt(2 / π) * 1/(d.σ1 + d.σ2)



function _norm(d::SplitGaussian, lo::Real, hi::Real)
    quadgk(x-> pdf(d, (lo+(hi-lo)*x)) * (hi-lo), lo, hi, atol=1e-15, order=9)[1]
end

function _normalize_logpdf(d::SplitGaussian, lo::Real, hi::Real, y::Real) 
    return logpdf(d,y)
end
