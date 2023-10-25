abstract type MixtureHack <: Distributions.ContinuousUnivariateDistribution end


struct MixtureHackExtended  <: MixtureHack
    summands::AbstractArray{<:Distributions.ContinuousUnivariateDistribution}
    coefficients::AbstractArray
end

struct MixtureHackNotExtended  <: MixtureHack
    summands::AbstractArray{<:Distributions.ContinuousUnivariateDistribution}
    coefficients::AbstractArray
end

function MixtureHack(summands::AbstractArray{<:Distributions.ContinuousUnivariateDistribution}, coefficients::AbstractArray, extended::Bool)
    extended == true ? MixtureHackExtended(summands, coefficients) : MixtureHackNotExtended(summands, coefficients)
end


pdf(d::MixtureHackExtended, x::Real) = 1/sum(d.coefficients) *  sum(d.coefficients .* pdf.(d.summands, x))


logpdf(d::MixtureHackExtended, x::Real) = log(pdf(d::MixtureHackExtended, x::Real))