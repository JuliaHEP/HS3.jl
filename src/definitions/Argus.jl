"""
ARGUS(m0, c, p)
The *Argus distribution* with cut-off `c`, curvature `χ` and power `p` has probability density function
```math
\\mathrm{Argus}(m; m0, c, p) = \\mathcal{N} \\cdot m \\cdot \\left[ 1 - \\left( \\frac{m}{m0} \\right)^2 \\right]^p
  \\cdot \\exp\\left[ c \\cdot \\left(1 - \\left(\\frac{m}{m0}\\right)^2 \\right) \\right]
```

```julia
Argus()        # Argus distribution with unit curvature, cut-off of 1 and power of 0.5 i.e. Argus(1, 1, 0.5)
Argus(c, χ)    # Argus distribution with curvature `χ`, cut-off `c` and power of 0.5 i.e. Argus(c, χ, 0.5), the so-called 'regular' Argus distribution
Argus(c, χ, p)    # Argus distribution with curvature `χ`, cut-off `c` and power of `p` i.e. Argus(c, χ, p)
params(d)        # Get the parameters, i.e. (c, χ, p)
cut_off(d)         # Get the cut-off parameter, i.e. c
curvature(d)        # Get the curvature parameter, i.e. χ 
power(d)        # Get the power parameter, i.e. p
```
External links
* [Argus distribution on Wikipedia](https://en.wikipedia.org/wiki/ARGUS_distribution)
"""

using SpecialFunctions

struct Argus{T<:Real} <: Distributions.ContinuousUnivariateDistribution
    c::T  # cut-off 
    χ::T   # curvature
    p::T   # power 
    function Argus{T}(c::T, χ::T, p::T) where {T <: Real}
        new{T}(c, χ, p)
    end
end

function Argus(c::T, χ::T, p::T; check_args::Bool=true) where {T <: Real}
    Distributions.@check_args Argus (c, c > zero(c)) (χ, χ > zero(χ)) (p, p > -1 * one(p))
    return Argus{T}(c, χ, p)
end

Argus(c::Real, χ::Real, p::Real; check_args::Bool=true) = Argus(promote(c, χ, p)...; check_args=check_args)
Argus(c::Real, χ::Real; check_args::Bool=true) = Argus(c, χ, 0.5; check_args=check_args)
Argus() = Argus(1., 1., 0.5; check_args=false)

#Distributions.@distr_support Argus 0.0 d.c

#### Conversions

convert(::Type{Argus{T}}, c::Real, χ::Real, p::Real) where {T<:Real} = Argus(T(c), T(χ), T(p))
Base.convert(::Type{Argus{T}}, d::Argus) where {T<:Real} = Argus{T}(T(d.c), T(d.χ), T(d.p))
Base.convert(::Type{Argus{T}}, d::Argus{T}) where {T<:Real} = d

##### Parameters
cut_off(d::Argus) = d.c  
curvature(d::Argus) = d.χ 
power(d::Argus) = d.p

params(d::Argus) = (d.c, d.χ, d.p)
@inline partype(::Argus{T}) where {T<:Real} = T

##### Statistics

mean(d::Argus) = d.c * d.p * sqrt(pi) * SpecialFunctions.gamma(d.p)/SpecialFunctions.gamma(5/2 + d.p)  * d.χ^(2* d.p+2) /(2^(d.p+2)) * HypergeometricFunctions.M(d.p + 1, 5/2 + d.p, -d.χ^2 / 2)/_n_val(d)

mode(d::Argus)  = d.c/(sqrt(2) * d.χ) * sqrt((d.χ^2 - 2* d.p -1) + sqrt(d.χ^2 * (d.χ^2 - 4* d.p +2) + (1 + 2* d.p)^2))

var(d::Argus) = d.c^2 * ( (0.5 * d.χ^2)^(d.p +1) * d.χ^(d.p+3) * exp(- 0.5 *d.χ^2) + (d.χ^2 - 2* (d.p+1)) * (SpecialFunctions.gamma(d.p + 2) - SpecialFunctions.gamma(d.p + 2, 0.5 * d.χ^2)))/(d.χ^2 * (d.p+1) * _n_val(d)) - mean(d)^2

##### Evaluation

function logpdf(d::Argus{T}, x::Real) where T<:Real
    if x <= 0.0 || x >= d.c
        1.0
    else
        c, χ, p = d.c, d.χ, d.p
        - p * log(2) + log(χ) * (2 * (p +1)) - log(_n_val(d)) + log(x/c^2) + p * log( 1 - x^2/c^2) + (- 0.5 * χ^2 * (1- x^2/c^2))
    end
end

function pdf(d::Argus{T}, x::Real) where T<:Real
    if x < 0.0 || x > d.c
        zero(T)
    else
        c, χ, p = params(d)
        2^(-p) * χ^(2 * (p +1)) /_n_val(d) * x/c^2 * ( 1 - x^2/c^2)^p * exp(- 0.5 * χ^2 * (1- x^2/c^2))
    end

end

function Distributions.cdf(d::Argus{T}, x::Real) where T <:Real
    if x < 0.0 || x > d.c
        0.0
    else
        (SpecialFunctions.gamma(d.p+1, 0.5 * d.χ^2 * (1 - x^2/d.c^2)) - SpecialFunctions.gamma(d.p+1, 0.5 * d.χ^2))/_n_val(d) 

    end
end

export cdf 
## often used value
function _n_val(d::Argus{T}) where T <: Real 
    (SpecialFunctions.gamma(d.p+1) - SpecialFunctions.gamma(d.p+1, 0.5 * d.χ^2))
end

