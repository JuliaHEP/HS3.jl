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
function make_data(spec::BinnedDataSpec)
    number_of_bins =  length(spec.contents) 
    min = 0
    max = number_of_bins 
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
function make_data(spec::PointDataSpec)
    spec.value
end

function make_data(spec::UnbinnedDataSpec)
    @warn "Unbinned data not yet implemented" 
end


"""
    make_datasets(data_specs::NamedTuple)

Create datasets from the given data specifications.

# Arguments
- `data_specs::NamedTuple`: Named tuple containing the data specifications.

# Returns
- Named tuple containing the datasets, where each key corresponds to a data specification and the value is the corresponding dataset created using the `make_data` function.

"""
function make_datas(data_specs::NamedTuple)
    datasets = NamedTuple()
    for (k, v) in zip(keys(data_specs), data_specs)
        datasets = merge(datasets, (Symbol(k) => make_data(v),))
    end
    datasets
end