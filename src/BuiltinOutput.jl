mutable struct Storage
    dt::Dict
end
export Storage

Storage() = Storage(Dict())

function store(st::Storage, x::Any, n::String)
    d = make_dict(x, n)
    st.dt[n] = d[n]
    return nothing
end

function store(st::Storage, x::Symbol, n::String)
    s = ":"*string(x)
    st.dt[n] = s
    return nothing
end 

function store(st::Storage, x::Expr, n::String)
    s = ":("*string(x)*")"
    st.dt[n] = s
    return nothing
end 

export store


function make_dict(var, name::String)
    return Dict(name => to_dict(var, name))
end

#function make_dict(var, m)
#    @capture(var, a_ = b_)
#    return Dict(String(a) => to_dict(m.eval(b), String(a)))
#end

function to_dict(x, above)
    dc = Dict()
    fns = fieldnames(typeof(x))

    if length(fns) == 0
        if isa_TOML_value(x)
            return x
        elseif isa(x, Union{Symbol, Expr})
            return repr(x)
        else
            return string(Base.typename(typeof(x)).wrapper)
        end
    else 
        tn = Base.typename(typeof(x)).wrapper
        #dc[string(above)] = string(tn)
        dc["⟪Type⟫"] = string(tn)
        for f in fns
            dc[string(f)]=to_dict(getfield(x, f), string(f)) #recursive call
        end
    end

    return dc
end

function default_storage()
    throw(
        ArgumentError(
        "No default LAMA.Storage object provided. Overload the function \'default_storage\' by doing:
        import LAMA.default_storage    
        st = Storage()
        LAMA.default_storage() = st ")
    )
end


macro store(var)
    parent_module = __module__
    d = make_dict(var, parent_module)

    st = default_storage()
    merge!(st.dt, d)

    return esc(:($var))    
end

#macro store(var::ParameterPointSpec)
#    #d = make_dict()
#    d = make_dict(var, fieldnames(typeof(var)))
#    print(d)
#    #st = default_storage()
#    #merge!(st.dt, d)
#
#    return esc(:($var))    
#end

export @store
