using HS3
using DensityInterface
using BenchmarkTools
using NamedTupleTools
using ValueShapes

#init 
dict = file_to_dict("/ceph/groups/e4/users/cburgard/public/forRobin/higgs/ouput/noES_ind/hh/125.json")
typeof(Dict)
Dict 
specs = HS3.generate_specs(dict)
analysis = HS3.make_analyses(specs.analyses[5], specs)

prior = HS3.generate_uniform_prior_from_domain(analysis[:parameter_domain])
constant_points = HS3.make_const_parameterpoints(specs.parameter_points.default_values)
all_points = HS3.make_parameterpoints(specs.parameter_points.default_values)

prior
#remove unnecassary variables
prior = delete(NamedTuple(prior), (:obs_x_ATLAS_H_2e2mu_channel_2011_1250_1 , :obs_x_ATLAS_H_2e2mu_channel_2012_1250_1, :obs_x_ATLAS_H_2mu2e_channel_2011_1250_1, :obs_x_ATLAS_H_2mu2e_channel_2012_1250_1, :obs_x_ATLAS_H_4e_channel_2011_1250_1, :obs_x_ATLAS_H_4e_channel_2012_1250_1, :obs_x_ATLAS_H_4mu_channel_2011_1250_1, :obs_x_ATLAS_H_4mu_channel_2012_1250_1))
prior = merge((prior), constant_points)

using BAT
all_points
posterior = PosteriorMeasure(analysis[:likelihood], (prior))
@btime logdensityof(analysis[:likelihood], rand(prior))

function unique_fields_values(nt1::NamedTuple, nt2::NamedTuple)
    unique_keys = setdiff(fieldnames(typeof(nt1)), fieldnames(typeof(nt2)))
    return NamedTuple{Tuple(unique_keys)}(tuple(map(k -> getproperty(nt1, k), unique_keys)...))
end
un = unique_fields_values(analysis[:parameter_domain], (constant_points))
value_nt = explore_values(un, analysis[:parameter_domain], analysis[:likelihood], all_points)


function transform_keys(nt::NamedTuple)
    # Use a generator to create the new key-value pairs
    new_pairs = (
        Symbol(replace(String(k), r"bin_(\d+)_hh$" => s"hh_sep_edit_bin_\1")) => v
        for (k, v) in pairs(nt)
    )
    return NamedTuple(new_pairs)
end

# Example usage
#nt = (gamma_stat_HHChannel_bin_0_hh = 1, gamma_stat_HHChannel_bin_1_hh = 2, some_other_key = 3)
new_nt = transform_keys(value_nt)
#value_nt = explore_values(new_nt, analysis[:parameter_domain], analysis[:likelihood], all_points)
write_to_csv(new_nt, "example/higgs/output_2h.csv")
