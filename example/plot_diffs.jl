using CSV
function write_to_csv(nt::NamedTuple, filename::String)
    # Create an array to collect rows for the CSV
    rows = []
    
    # Process each key-value pair in the NamedTuple
    for (key, values) in pairs(nt)
        for (param_value, func_output) in values
            push!(rows, (param_name=key, param_value=param_value, func_output=func_output))
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

using CSV, Plots

# Read the CSV file directly into a struct
data = CSV.File("output.csv")

# Extract the unique param_names
param_names = unique([row.param_name for row in data])

# Loop through each unique param_name to create a plot
for param in param_names
    # Extract rows corresponding to the current param_name
    relevant_rows = filter(row -> row.param_name == param, data)
    
    param_values = [row.param_value for row in relevant_rows]
    func_outputs = [row.func_output for row in relevant_rows]
    new_values = [row.new_value for row in relevant_rows]

    # Calculate relative deviations
    relative_deviations = (new_values .- func_outputs) ./ func_outputs

    # Plot
    p1 = plot(param_values, [func_outputs, new_values], label=["Original Func Output" "New Value"], color=[:blue :red], marker=[:circle :x], xlabel="Param Value", ylabel="Function Output", title="Values for param_name = $param")
    p2 = plot(param_values, relative_deviations, label="Relative Deviation", color=:green, marker=:square, xlabel="Param Value", ylabel="Relative Deviation")

    plot(p1, p2, layout=(2,1), size=(600,400))
    display(current())
end
