## This file is a part of HS3.jl, licensed under the MIT License (MIT).

function node_deps(node_name, nodes)
    deps = []
    for element in keys(nodes) 
        if node_name in (nodes[element])
            push!(deps, element)
        end
    end
    return deps
end

mutable struct Node_info
    state::Int
    layer::Int
end
#

function flatten(t::Tuple)
    result = ()
    for x in t
        if isa(x, Tuple)
            result = (result..., flatten(x)...)
        else
            result = (result..., x)
        end
    end
    return result
end

function compare_layers(a, b)
    return a[2] < b[2] 
end


function topological_sort(specs::NamedTuple)
    nodes =(;)  
    node_array = []
    key_array = []
    node_info_array = NamedTuple()
    for key in keys(specs)
        if typeof(specs[key]) <: AbstractDistributionsSpec
            node_array = push!(node_array, flatten((values(specs[key].params), values(specs[key].variables))))
        elseif typeof(specs[key]) <: AbstractFunctionsSpec
            node_array = push!(node_array, flatten((values(specs[key].params))))
        elseif typeof(specs[key]) <: HistFactorySpec
            node_array = push!(node_array, (;))
        else
            @error "Spec Type not supported"
        end
        key_array = push!(key_array, Val(key))
        node_info_array = merge(node_info_array, (key => Node_info(0, 0),))
    end
    nodes = (; zip(keys(specs), node_array)...)
    ### Visited = 2, Visiting = 1, Unvisited = 0
    while !isempty(key_array)
        node_name = key_array[end]
        deps = node_deps(node_name, nodes)
        if node_info_array[_val_content(node_name)].state == 2
            ##already visited 
            pop!(key_array)
        elseif isempty(deps) 
            ## Top level node
            node_info_array[_val_content(node_name)].layer = 0
            node_info_array[_val_content(node_name)].state = 2
            pop!(key_array)
        else 
            allDepsVisited = true
            for dep in deps 
                if node_info_array[dep].state == 1
                    @error("Not a DAG")
                elseif node_info_array[dep].state == 0
                    allDepsVisited = false
                    key_array = push!(key_array, Val(dep))
                end    
            end

            if allDepsVisited == true
                node_info_array[_val_content(node_name)].layer = 1 + maximum([node_info_array[x].layer for x in deps])
                node_info_array[_val_content(node_name)].state = 2
                pop!(key_array)
            elseif node_info_array[_val_content(node_name)].state != 1
                node_info_array[_val_content(node_name)].state = 1
            else
                @error("Internal error in topoligical sort")
            end
        end
    end
    #nodes
    #nodes[node_name]
    to_sort = []
    for key in keys(node_info_array)
        push!(to_sort, Pair(key, node_info_array[key].layer))
    end
    sorted = sort(to_sort, lt=compare_layers)
    sorted_specs = []
    for element in sorted 
        push!(sorted_specs, element[2] => (Symbol(element[1]), specs[Symbol(element[1])]))
    end
    return sorted_specs
end
export topological_sort
