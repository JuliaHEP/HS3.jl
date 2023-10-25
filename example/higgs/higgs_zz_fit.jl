import Pkg; 
Pkg.activate("/net/e4-nfs-home.e4.physik.tu-dortmund.de/home/rpelkner/Documents/HS3")
using HS3
using DensityInterface
using BenchmarkTools
using NamedTupleTools
using ValueShapes
using BAT

#init 
dict = file_to_dict("/ceph/groups/e4/users/cburgard/public/forRobin/higgs/ouput/noES_ind/zz/125.json")
specs = HS3.generate_specs(dict)
analysis = HS3.make_analyses(specs.analyses[5], specs)

prior = HS3.generate_uniform_prior_from_domain(analysis[:parameter_domain])
constant_points = HS3.make_const_parameterpoints(specs.parameter_points.default_values)
#all_points = HS3.make_parameterpoints(specs.parameter_points.default_values)

#remove unnecassary variables
prior = delete(NamedTuple(prior), (:obs_x_ATLAS_H_2e2mu_channel_2011_1250_1 , :obs_x_ATLAS_H_2e2mu_channel_2012_1250_1, :obs_x_ATLAS_H_2mu2e_channel_2011_1250_1, :obs_x_ATLAS_H_2mu2e_channel_2012_1250_1, :obs_x_ATLAS_H_4e_channel_2011_1250_1, :obs_x_ATLAS_H_4e_channel_2012_1250_1, :obs_x_ATLAS_H_4mu_channel_2011_1250_1, :obs_x_ATLAS_H_4mu_channel_2012_1250_1))
prior = merge(NamedTupleDist(prior), constant_points)

posterior = PosteriorMeasure(analysis[:likelihood], (prior))

import BAT
import BAT: bat_findmode
using Optim
using BAT: AnyMeasureOrDensity, AbstractMeasureOrDensity
using BAT: get_context, get_adselector, _NoADSelected
using BAT: bat_initval, transform_and_unshape, apply_trafo_to_init
using BAT: negative, inverse, fchain
using DensityInterface, ChangesOfVariables, InverseFunctions, FunctionChains, Distributions

function BAT.bat_findmode_impl(target::AnyMeasureOrDensity, algorithm::OptimAlg, context::BATContext)
    @info "Using my own bat_findmode_impl"
     transformed_density, trafo = transform_and_unshape(algorithm.trafo, target, context)
     inv_trafo = inverse(trafo)
 
     initalg = apply_trafo_to_init(trafo, algorithm.init)
     x_init = collect(bat_initval(transformed_density, initalg, context).result)
     # Maximize density of original target, but run in transformed space, don't apply LADJ:
     f = fchain(inv_trafo, logdensityof(target), -)
     optim_result = my_optim_minimize(f, x_init, Optim.LBFGS())
     r_optim = Optim.MaximizationWrapper(optim_result)
     transformed_mode = Optim.minimizer(r_optim.res)
     result_mode = inv_trafo(transformed_mode)
     (result = result_mode, result_trafo = transformed_mode, trafo = trafo, #=trace_trafo = trace_trafo,=# info = r_optim)
end

function my_optim_minimize(f, x_init, algorithm)
    @info "using local version of _optim_minimize for zeroth order optimizers!"
    algorithm = Optim.LBFGS(; m=100)
    opts = Optim.Options(store_trace = true, extended_trace=true, iterations = 4000, allow_f_increases=true, g_tol=1e-08 )
    Optim.optimize(f, x_init, algorithm, opts)
end
init = rand(prior)
optimizer1 = OptimAlg(init = ExplicitInit([init]), trafo=DoNotTransform())
println("DoNotTransform")
println(@benchmark bat_findmode(posterior, optimizer1))
optimizer2 = OptimAlg(init = ExplicitInit([init]), trafo=PriorToGaussian())
println("PriorToGaussian")
println(@benchmark bat_findmode(posterior, optimizer2))
optimizer3 = OptimAlg(init = ExplicitInit([init]), trafo=PriorToUniform())
println("PriorToUniform: ")
println(@benchmark bat_findmode(posterior, optimizer3))

#optimizer1 = OptimAlg(init = ExplicitInit([init]), trafo=DoNotTransform())
#mode1 = bat_findmode(posterior, optimizer1)
#println("Do not transform: ")
#println(mode1.result)
#optimizer2 = OptimAlg(init = ExplicitInit([init]), trafo=PriorToGaussian())
#mode2 = bat_findmode(posterior, optimizer2)
#println("PriorToGaussian: ")
#println(mode2.result)
#optimizer3 = OptimAlg(init = ExplicitInit([init]), trafo=PriorToUniform())
#println("PriorToUniform: ")
#mode3 = bat_findmode(posterior, optimizer3)
#println(mode3.result)

#open("/ceph/groups/e4/users/rpelkner/public/HS3/higgs_zz/out/optimized_not_transformed.txt", "w") do file
#    write(file, mode1)
#end
#
#open("/ceph/groups/e4/users/rpelkner/public/HS3/higgs_zz/out/optimized_priortogaussian.txt", "w") do file
#    write(file, mode2)
#end
#
#open("/ceph/groups/e4/users/rpelkner/public/HS3/higgs_zz/out/optimized_priortouniform.txt", "w") do file
#    write(file, mode3)
#end
