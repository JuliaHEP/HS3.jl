"""
Simple hack to create products of PDFs. WIP
"""

struct ProductHack <: Distributions.ContinuousUnivariateDistribution
    factors::AbstractArray{Distributions.ContinuousUnivariateDistribution}
end

pdf(d::ProductHack, x::Real) = prod(pdf.(d.factors, x))
logpdf(d::ProductHack, x::Real) = sum(logpdf.(d.factors, x))