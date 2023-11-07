"""
This is a hack introduced for mixture distributions (addition of PDFs); essenetially the same results as with MixtureModel
    of Distributions.jl can be achieved, but for convenience this hack will be kept in for now.
"""

struct MixtureHack  <: Distributions.ContinuousUnivariateDistribution
    summands::AbstractArray{<:Distributions.ContinuousUnivariateDistribution}
    coefficients::AbstractArray
end

function MixtureHack(summands::AbstractArray{<:Distributions.ContinuousUnivariateDistribution}, coefficients::AbstractArray, extended::Bool)
    extended == true ? MixtureHack(summands, coefficients) : MixtureHack(summands, push!(coefficients, 1 - sum(collect(coefficients))))
end

pdf(d::MixtureHack, x::Real) =   sum(d.coefficients .* pdf.(d.summands, x))


logpdf(d::MixtureHack, x::Real) = log(pdf(d, x))