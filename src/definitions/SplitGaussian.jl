"""
SplitGaussian(mu, sigma1, sigma2) also known as bifurkated Gaussian

```math
\\mathrm{SplitGaussian}(x;\\mu,\\sigma_1,\\sigma_2)= A \\exp \\left(- \\frac {(x-\\mu)^2}{2 \\sigma_1^2}\\right) \\quad \\text{if } x< \\mu

\\mathrm{SplitGaussian}(x;\\mu,\\sigma_1,\\sigma_2)=   A \\exp \\left(- \\frac {(x-\\mu)^2}{2 \\sigma_2^2}\\right) \\quad \\text{otherwise}

A= \\sqrt{2/\\pi} (\\sigma_1+\\sigma_2)^{-1}.
```

External Links: 
* [Wikipedia](https://en.wikipedia.org/wiki/Split_normal_distribution)
* [RooFit](https://root.cern.ch/doc/master/classRooBifurGauss.htmls)
"""

struct SplitGaussian{T<:Real} <: ContinuousUnivariateDistribution
    μ::T
    σ1::T 
    σ2::T
    function SplitGaussian{T}(μ::T, σ1::T, σ2::T) where {T <: Real}
        new{T}(μ, σ1, σ2)
    end
end

function SplitGaussian(μ::T, σ1::T, σ2::T; check_args::Bool=true) where {T <: Real}
    @check_args SplitGaussian (σ1, σ1 > zero(σ1)) (σ2, σ2 > zero(σ2)) 
    if σ1 != σ2
        return SplitGaussian{T}(μ, σ1, σ2)
    else 
        return Distributions.Normal(μ, σ1)
    end
end

SplitGaussian(μ::Real, σ1::Real, σ2::Real; check_args::Bool=true) = SplitGaussian(promote(μ, σ1, σ2)...; check_args=check_args)
SplitGaussian(μ::Real, σ1::Real; check_args::Bool=true) = Normal(μ, σ1)
SplitGaussian() = SplitGaussian(0., 1., 2.; check_args=false)

convert(::Type{SplitGaussian{T}}, μ::Real, σ1::Real, σ2::Real) where {T<:Real} = SplitGaussian(T(μ),  T(σ1), T(σ2))
Base.convert(::Type{SplitGaussian{T}}, d::SplitGaussian) where {T<:Real} = SplitGaussian{T}(T(d.μ),  T(d.σ1), T(d.σ2))
Base.convert(::Type{SplitGaussian{T}}, d::SplitGaussian{T}) where {T<:Real} = d


location(d::SplitGaussian) = d.μ

params(d::SplitGaussian) = (d.μ, d.σ1, d.σ2)

partype(::SplitGaussian{T}) where {T<:Real} = T

mode(d::SplitGaussian) = d.μ

mean(d::SplitGaussian) = d.μ + sqrt(2/ π) * (d.σ2 - d.σ1)

var(d::SplitGaussian) = (1 - 2/π)*(d.σ2 - d.σ1)^2 + d.σ1 * d.σ2

skewness(d::SplitGaussian) = sqrt(2/π) * (d.σ2 - d.σ1)*((4/π - 1)*(d.σ2 - d.σ1)^2 + d.σ1 * d.σ2)

function pdf(d::SplitGaussian, x::Real) 
    if x < d.μ
        return exp(-(x-d.μ)^2 / (2 * d.σ1^2)) * _n_val(d)
    else 
        return exp(-(x-d.μ)^2 / (2 * d.σ2^2)) * _n_val(d)
    end
end

function logpdf(d::SplitGaussian, x::Real) 
    if x < d.μ
        return (-(x-d.μ)^2 / (2 * d.σ1^2)) + log(_n_val(d))
    else 
        return (-(x-d.μ)^2 / (2 * d.σ2^2)) + log(_n_val(d))
    end
end

_n_val(d::SplitGaussian) = sqrt(2 / π) * 1/(d.σ1 + d.σ2)