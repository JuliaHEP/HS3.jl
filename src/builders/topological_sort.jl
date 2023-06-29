## This file is a part of HS3.jl, licensed under the MIT License (MIT).
"""
    node_deps(node_name, nodes)

Find the dependencies of a given node in a named tuple of nodes.

# Arguments
- `node_name`: The name of the node to find dependencies for.
- `nodes`: A named tuple of nodes, where the keys are node names and the values are tuples of dependencies.

# Returns
An array of dependencies for the given node.
"""
function node_deps(node_name, nodes)
    deps = []
    for element in keys(nodes) 
        if Val(node_name) in (nodes[element])
            deps = push!(deps, element)
        end
    end
    return deps
end


"""
    mutable struct NodeInfo

Node information for a node in the topological sort.

# Fields
- `state`: The state of the node (Visited = 2, Visiting = 1, Unvisited = 0).
- `layer`: The layer of the node in the topological sort.

"""
mutable struct NodeInfo
    state::Int
    layer::Int
end
#
"""
    _flatten_parameters(t::Tuple)

Flatten a nested tuple of parameters.

# Arguments
- `t`: The nested tuple of parameters.

# Returns
A flattened array of parameters.
"""
function _flatten_parameters(t::Tuple)
    result = []
    for x in t
        if isa(x, Tuple)
            result = (result..., _flatten_parameters(x)...)
        else
            result =  (result..., x)
        end
    end
    return result
end

"""
    compare_layers(a, b)

Comparison function for sorting based on layer values.

# Arguments
- `a`: The first element to compare.
- `b`: The second element to compare.

# Returns
`true` if the layer of `a` is less than the layer of `b`, `false` otherwise.
"""
function compare_layers(a, b)
    return a[2] < b[2] 
end

"""
    _create_nodes(specs::NamedTuple)

Create nodes, key array, and node info dictionary from functional specifications.

# Arguments
- `specs`: A named tuple of functional specifications.

# Returns
A named tuple of nodes, key array, and a node info dictionary for all elements in specs.
"""
function _create_nodes(specs::NamedTuple)
    node_array = []
    key_array = []
    node_info_dict = Dict{Symbol, NodeInfo}()
    for key in keys(specs)
        if typeof(specs[key]) <: AbstractDistributionSpec
            node_array = push!(node_array, _flatten_parameters(values(specs[key].params)))
        elseif typeof(specs[key]) <: AbstractFunctionSpec
            node_array = push!(node_array, _flatten_parameters((values(specs[key].params))))
        elseif typeof(specs[key]) <: HistFactorySpec
            node_array = push!(node_array, [])
        else
            @error "Spec Type $(typeof(spec)) not supported"
        end
        key_array = push!(key_array, key)
        node_info_dict[key] = NodeInfo(0, 0)
    end
    nodes = (; zip(keys(specs), node_array)...)
    return nodes, key_array, node_info_dict
end

"""
    _sort_specs(node_info_dict, specs)

Sort the specifications based on layer values.

# Arguments
- `node_info_dict`: A dictionary of node information.
- `specs`: A named tuple of specifications.

# Returns
A sorted array of specifications.
"""
function _sort_specs(node_info_dict, specs)
    to_sort = []
    for key in keys(node_info_dict)
        push!(to_sort, Pair(key, node_info_dict[key].layer))
    end
    sorted = sort(to_sort, lt=compare_layers)
    sorted_specs = []
    for element in sorted 
        push!(sorted_specs, element[2] => (Symbol(element[1]), specs[Symbol(element[1])]))
    end
    return sorted_specs
end

"""
    topological_sort(specs::NamedTuple)

Perform topological sorting on the specifications.

# Arguments
- `specs`: A named tuple of specifications.

# Returns
A sorted array of specifications.
"""
function topological_sort(specs::NamedTuple)
    nodes, key_array, node_info_dict = _create_nodes(specs)
    ### Visited = 2, Visiting = 1, Unvisited = 0
    while !isempty(key_array)
        node_name = key_array[end]
        deps = node_deps(node_name, nodes)
        if node_info_dict[node_name].state == 2
            ##already visited 
            pop!(key_array)
        elseif isempty(deps) 
            ## Top level node
            node_info_dict[node_name].layer = 0
            node_info_dict[node_name].state = 2
            pop!(key_array)
        else 
            allDepsVisited = true
            for dep in deps 
                if node_info_dict[dep].state == 1
                    @error("Not a DAG")
                elseif node_info_dict[dep].state == 0
                    allDepsVisited = false
                    key_array = push!(key_array, dep)
                end    
            end
            if allDepsVisited == true
                node_info_dict[node_name].layer = 1 + maximum([node_info_dict[x].layer for x in deps])
                node_info_dict[node_name].state = 2
                pop!(key_array)
            elseif node_info_dict[node_name].state != 1
                node_info_dict[node_name].state = 1
            else
                @error("Internal error in topoligical sort")
            end
        end
    end
    return _sort_specs(node_info_dict, specs)
end
