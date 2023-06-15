function file_to_dict(path::String)
    try
        f = open(path, "r")
        data = read(f)
        json_obj = JSON3.read(data)
        Dict(json_obj)
    catch e
        @error "Invalid JSON File"
        throw(e)
    end
end

export file_to_dict

