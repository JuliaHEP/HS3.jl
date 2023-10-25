# This file is a part of HS3.jl, licensed under the MIT License (MIT).
"""
    make_data(spec::AbstractBinnedDataSpec)

Create a histogram from the given binned data specification.

# Arguments
- `spec::AbstractBinnedDataSpec`: The binned data specification.

# Returns
A `StatsBase.Histogram` object representing the histogram of the binned data.

# Example
```julia
spec = BinnedDataSpec(...)
histogram = make_data(spec)

"""
function make_data(spec::BinnedDataSpec, domain::NamedTuple=(;))
    number_of_bins =  length(spec.contents) 
    min = 0
    max = number_of_bins 
    println(spec.axes)
    if spec.axes[1].max !== nothing 
        max = spec.axes[1].max 
    end
    if spec.axes[1].min !== nothing 
        min = spec.axes[1].min 
    end
    step = (max - min)/number_of_bins
    bins = range(min, stop = max, step=step)
    v = collect(spec.contents)
    StatsBase.Histogram(bins, v)
end


"""
    make_data(spec::PointDataSpec)

Create data from the given point data specification. 
*Note: the optional error is of the data point is not considered*

# Arguments
- `spec::PointDataSpec`: The point data specification that defines the data value.

# Returns
- The data value point specified by the specification.

"""
function make_data(spec::PointDataSpec, domain::NamedTuple)
    spec.value
end

@with_kw struct UnbinnedData 
    entries::AbstractArray
    entries_errors::Union{AbstractArray, Nothing} = nothing 
    axes::NamedTuple
    weights::Union{AbstractArray, Nothing} = nothing 
end

function make_data(spec::UnbinnedDataSpec, domain::NamedTuple)
    println(spec)
    if spec.weights === nothing
        return UnbinnedData(spec.entries, spec.entries_errors, make_axes(spec.axes, domain), ones(length(spec.entries)))
    else
        return UnbinnedData(spec.entries, spec.entries_errors, make_axes(spec.axes, domain), spec.weights)
    end
end

function _filter_nothing(axes_spec::AxesSpec)
    fields = fieldnames(typeof(axes_spec))
    filtered_fields = [(f, getfield(axes_spec, f)) for f in fields if getfield(axes_spec, f) !== nothing]
    return NamedTuple(filtered_fields)
end


function make_axes(axes_spec::NamedTuple, domain::NamedTuple)
    fields1 = fieldnames(typeof(axes_spec))
    fields2 = fieldnames(typeof(domain))

    merged_fields = Dict{Symbol, Any}()
    for f in fields1
        axes_spec_filtered = _filter_nothing(axes_spec[f])
        if f in fields2
            #println(domain[f], f, axes_spec_filtered)
            merged_fields[f] = merge(axes_spec_filtered, (range = domain[f], ))
        else
            merged_fields[f] = axes_spec_filtered
        end
    end

    return NamedTuple(merged_fields)
end

"""
    make_datasets(data_specs::NamedTuple)

Create datasets from the given data specifications.

# Arguments
- `data_specs::NamedTuple`: Named tuple containing the data specifications.

# Returns
- Named tuple containing the datasets, where each key corresponds to a data specification and the value is the corresponding dataset created using the `make_data` function.

"""
function make_datas(data_specs::NamedTuple, domain::NamedTuple)
    names = [Symbol(name) for name in keys(data_specs)]
    sets = [make_data(set, domain) for set in data_specs]
    (; zip(names, sets)...)
end