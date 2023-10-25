# This file is a part of HEPDistributions.jl, licensed under the MIT License (MIT).
"""
CrystalBall distribution implementation in Julia.

TODO: Add PDF Math here.

The `CrystalBall` module provides functionality for working with the CrystalBall distribution. This distribution is commonly used in statistical analysis, particularly in the field of high-energy physics, to model signals with tails extending into negative or positive infinity.

The distribution is parameterized by `x0`, `sigma`, `alpha`, and `n`, representing the location, width, tail shape, and tail exponent, respectively.

## Usage

The `CrystalBall` module provides the following functions:

- `CrystalBall(x0::T, sigma::T, alpha::T, n::T)` - Constructs a `CrystalBall` distribution with the specified parameters.
- `pdf(d::CrystalBall, x::Real)` - Computes the probability density function (PDF) at the given value.
- `logpdf(d::CrystalBall, x::Real)` - Computes the logarithm of the PDF at the given value.
- `params(d::CrystalBall)` - Returns a tuple containing the distribution parameters.
- `location(d::CrystalBall)` - Returns the location parameter of the distribution.
- `width(d::CrystalBall)` - Returns the width parameter of the distribution.
- `alpha(d::CrystalBall)` - Returns the tail shape parameter of the distribution.
- `exponent(d::CrystalBall)` - Returns the tail exponent parameter of the distribution.

External Links: 
* [Wikipedia](https://en.wikipedia.org/wiki/Crystal_Ball_function)
"""

import Distributions: location, params, pdf, logpdf
using Distributions: @check_args, @distr_support, ContinuousUnivariateDistribution
using SpecialFunctions: erf

struct CrystalBall{T<:Real} <: ContinuousUnivariateDistribution
    x0::T       
    sigma::T 
    alpha::T
    n::T 
    function CrystalBall{T}(x0::T, sigma::T, alpha::T,  n::T) where {T <: Real}
        new{T}(x0, sigma, alpha, n)
    end
end

function CrystalBall(x0::T, sigma::T, alpha::T, n::T; check_args::Bool=true) where {T <: Real}
    @check_args CrystalBall (alpha, alpha != zero(alpha)) (sigma, sigma > zero(sigma)) (n, n > 1.0) 
    return CrystalBall{T}(x0, sigma, alpha, n)
end

CrystalBall(x0::Real, sigma::Real, alpha::Real, n::Real; check_args::Bool=true) = CrystalBall(promote(x0, sigma, alpha, n)...; check_args=check_args)
CrystalBall() = CrystalBall(0., 1., 1., 2.; check_args=false)

convert(::Type{CrystalBall{T}}, x0::Real, sigma::Real, alpha::Real, n::Real) where {T<:Real} = CrystalBall(T(x0),  T(sigma), T(alpha), T(n))
Base.convert(::Type{CrystalBall{T}}, d::CrystalBall) where {T<:Real} = CrystalBall{T}(T(d.x0),  T(d.sigma), T(d.alpha), T(d.n))
Base.convert(::Type{CrystalBall{T}}, d::CrystalBall{T}) where {T<:Real} = d

params(d::CrystalBall) = ((d.x0),  (d.sigma), (d.alpha), (d.alpha), (d.n))
location(d::CrystalBall) = d.x0 
width(d::CrystalBall) =  d.sigma
alpha(d::CrystalBall) = d.alpha 
exponent(d::CrystalBall) = d.n #maybe rename this one

mode(d::CrystalBall) = d.x0


function pdf(d::CrystalBall{T}, x::Real) where T<:Real 
    if (x-d.x0)/d.sigma > - d.alpha
        return exp(- (x - d.x0)^2/(2 * d.sigma^2)) * _n_val(d)
    else
        return _a_val(d) * (_b_val(d) - (x - d.x0)/d.sigma)^(-d.n) *_n_val(d)
    end
end

#normalisation
function _n_val(d::CrystalBall)
    c = d.n /abs(d.alpha) * 1/(d.n-1) * exp(- abs(d.alpha)^2 / 2)
    b = sqrt(pi/2) *( 1+ erf(abs(d.alpha)/sqrt(2)))
    return 1/(d.sigma*(b + c))
end

#helper functions for often used values
function _a_val(d::CrystalBall)
    return (d.n/abs(d.alpha))^(d.n) * exp(- abs(d.alpha)^2 / 2)
end

function _b_val(d::CrystalBall)
    return d.n / abs(d.alpha) - abs(d.alpha)
end

function logpdf(d::CrystalBall{T}, x::Real) where T<:Real
    if (x-d.x0)/d.sigma > - d.alpha
        return (- (x - d.x0)^2/(2 * d.sigma^2)) + _log_n_val(d)
    else
        return _log_a_val(d) + (-d.n) * log(_b_val(d) - (x - d.x0)/d.sigma) + _log_n_val(d)
    end
end

function _log_a_val(d::CrystalBall)
    return (d.n) * log(d.n/abs(d.alpha)) + (- abs(d.alpha)^2 / 2)
end

function _log_n_val(d::CrystalBall)
    c = d.n /abs(d.alpha) * 1/(d.n-1) * exp(- abs(d.alpha)^2 / 2)
    b = sqrt(pi/2) *( 1+ erf(abs(d.alpha)/sqrt(2)))
    return - (log(d.sigma) + log(b + c))
end

export mode, CrystalBall, location, width, params, pdf, logpdf, alpha, width, @check_args, @distr_support