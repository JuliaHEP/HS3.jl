# INTERNAL helper function for developement.
using CSV 
function write_to_csv(nt::NamedTuple, filename::String)
    # Create an array to collect rows for the CSV
    rows = []
    
    # Process each key-value pair in the NamedTuple
    for (key, values) in pairs(nt)
        for (param_value, func_output) in values
            push!(rows, (param_name=key, param_value=param_value, julia_output=func_output))
        end
    end
    
    # Write the data to a CSV file
    CSV.write(filename, rows)
end

#write_to_csv(value_nt, "output_modified_zgg.csv")

function explore_values(nt::NamedTuple, ranges::NamedTuple, ll, default::NamedTuple)
    # Initialize a Dict for collecting the results
    results_dict = Dict{Symbol, Vector{Tuple{Float64, Float64}}}()
    
    # For each key in the namedtuple
    for key in keys(nt)
        # If the key has a specified range
        if haskey(ranges, key)
            # A list to store tuples of (param_value, func_output) for this key
            key_results = []
            
            # Loop over the range
            for val in ranges[key]
                # Modify the namedtuple
                modified_nt = merge(default, (key => val,))
                
                # Compute the result
                result = logdensityof(ll, modified_nt) 
                
                # Store the parameter value and the result
                key_results = push!(key_results, (val, result))
            end
            
            # Save the results for this key in the results_dict
            results_dict[key] = key_results
        end
    end
    
    # Convert the Dict back to a NamedTuple for the final result
    return NamedTuple(results_dict)
end

function explore_values(nt::Vector{Symbol}, ranges::NamedTuple, ll, default::NamedTuple)
    # Initialize a Dict for collecting the results
    results_dict = Dict{Symbol, Vector{Tuple{Float64, Float64}}}()
    
    # For each key in the namedtuple
    for key in (nt)
        # If the key has a specified range
        if haskey(ranges, key)
            # A list to store tuples of (param_value, func_output) for this key
            key_results = []
            
            # Loop over the range
            for val in ranges[key]
                # Modify the namedtuple
                modified_nt = merge(default, (key => val,))
                
                # Compute the result
                result = logdensityof(ll, modified_nt) 
                
                # Store the parameter value and the result
                key_results = push!(key_results, (val, result))
            end
            
            # Save the results for this key in the results_dict
            results_dict[key] = key_results
        end
    end
    
    # Convert the Dict back to a NamedTuple for the final result
    return NamedTuple(results_dict)
end