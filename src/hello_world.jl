# This file is a part of HS3.jl, licensed under the MIT License (MIT).


"""
    HS3.hello_world()

Prints "Hello, World!" and returns 42.

```jldoctest
using HS3

HS3.hello_world()

# output

Hello, World!
42
```
"""
function hello_world()
    println("Hello, World!")
    return 42
end
