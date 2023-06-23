# This file is a part of HS3.jl, licensed under the MIT License (MIT).

### For now this is just a dump of code
function HS3_output()
    open("my_new_file.json", "w") do io
        JSON3.pretty(io, a)
    end
end