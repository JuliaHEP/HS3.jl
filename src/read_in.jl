function file_to_dict(path::String)
    f = open(path, "r")
    data = read(f)
    json_obj = JSON3.read(data)
    PropDicts.PropDict(Dict(json_obj))
end

export file_to_dict

#function make_all_specs(dict::AbstractDict)
#    
#end
#
#export make_all_specs