# This file is a part of HS3.jl, licensed under the MIT License (MIT).

generate_distributionspec(d::Distributions.Exponential, x) = DistributionSpec{:exponential_dist}((c = Distributions.scale(d), x = x,))

generate_distributionspec(d::Distributions.Normal, x) = DistributionSpec{:gaussian_dist}((mean = Distributions.mean(d), sigma = Distributions.std(d), x = x,))

#generate_distributionspec(::Val{:lognormal_dist}, data::NamedTuple) = DistributionSpec{:lognormal_dist}((mean = data.mean, k = data.k, x = data.x,))
#
#generate_distributionspec(::Val{:multinormal_dist}, data::NamedTuple) = DistributionSpec{:multinormal_dist}((covariances = reduce(hcat, data.covariances), mean = data.mean, x = data.x,) )
#
#generate_distributionspec(::Val{:poisson_dist}, data::NamedTuple) = DistributionSpec{:poisson_dist}((mean = data.mean, x = data.x,))
#
#generate_distributionspec(::Val{:argus_dist}, data::NamedTuple) = DistributionSpec{:argus_dist}((resonance = data.resonance, slope = data.slope, power = data.power, x = data.mass,))
#
#generate_distributionspec(::Val{:mixture_dist}, data::NamedTuple) = DistributionSpec{:mixture_dist}((summands = data.summands, coefficients = data.coefficients,)) 
#
#generate_distributionspec(::Val{:product_dist}, data::NamedTuple) = DistributionSpec{:product_dist}((factors = data.factors, )) 
#
#generate_distributionspec(::Val{:uniform_dist}, data::NamedTuple) = DistributionSpec{:uniform_dist}((min = data.min, max = data.max,))
#
#generate_distributionspec(::Val{:polynomial_dist}, data::NamedTuple) = DistributionSpec{:polynomial_dist}((coefficients = data.coefficients, x= data.x, ))
#
#generate_distributionspec(::Val{:CrystalBall_dist}, data::NamedTuple) = DistributionSpec{:CrystalBall_dist}((alpha = data.alpha, x= data.x, m0 = data.m0, n = data.n, sigma = data.sigma, ))
#
#function generate_distributionspec(::Val{:histfactory_dist}, data::NamedTuple)
#    HistFactorySpec(; (samples = generate_HistFactorySampleSpecs(data.samples), axes = generate_axes_specs(data.axes),)...)
#end
