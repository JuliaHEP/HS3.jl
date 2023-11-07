abstract type AbstractModifier end
abstract type AbstractInterp end

struct Histosys{T<:AbstractInterp} <: AbstractModifier
    interp::T
    function Histosys(interp::T) where T
        @assert T <: PiecewiseInterpolation #more possible codes in the future
        new{T}(interp)
    end
end

Histosys(nominal, down, up) = Histosys(PiecewiseInterpolation(nominal, down, up))

struct PiecewiseInterpolation <: AbstractInterp
    nominal::Vector{Float64} 
    low::Vector{Float64} 
    high::Vector{Float64}
end


#Piece-wise interpolation (Code 4)
@inline function (i::PiecewiseInterpolation)(α::Float64)::Vector{Float64}
    val = @. ( abs(α) <= 1.0 ? i.nominal + _interpolate_sixt_degree(α, i.low, i.high, i.nominal) :
          α < -1.0 ?  i.nominal + (i.nominal - i.low) * α :
          i.nominal + (i.high - i.nominal) * α )
    return @. max(0, val) - i.nominal
end

struct Normsys{T<:AbstractInterp} <: AbstractModifier
    interp::T
end

struct InterpCode4 <: AbstractInterp
    I0::Vector{Float64}
    f_ups::Float64
    f_downs::Float64
    α0::Float64
end

Normsys(nominal, ups, downs) = Normsys(InterpCode4(nominal, ups, downs)) #other codes possible in the future, but unecessary(?)

function InterpCode4(I0, I_up, I_down; α0=1.0)
    if I_up == 0 && I_down == 0
        @error ("error ", I_up, "  ", I_down)
    end
    InterpCode4(I0, Float64(I_up), Float64(I_down), α0)
end

@inline function (i::InterpCode4)(α::Float64)::Float64
    if abs(α) <= i.α0
        return exp(_interpolate_sixt_degree(α, log(i.f_downs), log(i.f_ups), 0.0, i.α0))
    elseif α > i.α0
        return i.f_ups != 0 ? (i.f_ups) ^ α : 0.0
    else
        return i.f_downs != 0 ? (i.f_downs) ^ (-α) : 0.0
    end
end


struct Normfactor <: AbstractModifier # is unconstrained
    interp::typeof(identity)
    Normfactor() = new(identity)
end


struct Shapefactor{T} <: AbstractModifier # is unconstrained, questionable whether this works
    interp::T
    function Shapefactor(nbins, nthbin)
        f = binidentity(nbins, nthbin, 1.)
        new{typeof(f)}(f)
    end
end

struct Shapesys{T} <: AbstractModifier
    σn2::Float64
    moddata::Float64
    interp::T
    function Shapesys(σneg2, moddata, nbins, nthbin)
        f = binidentity(nbins, nthbin, moddata)
        new{typeof(f)}(σneg2, moddata, f)
    end
end

struct Staterror{T} <: AbstractModifier
    σ::Float64
    interp::T
    function Staterror(σ, nbins, nthbin)
        f = binidentity(nbins, nthbin, 1.)
        new{typeof(f)}(σ, f)
    end
end

struct MultOneHot{T} <: AbstractVector{T}
    nbins::Int
    nthbin::Int
    α::T
end

Base.length(b::MultOneHot) = b.nbins
Base.getindex(b::MultOneHot{T}, n::Integer) where T = Base.ifelse(n==b.nthbin, b.α, one(T))
Base.size(b::MultOneHot) = (b.nbins, )

function binidentity(nbins, nthbin, moddata)
    α -> MultOneHot(nbins, nthbin, α)
end


struct CustomInterpolation <: AbstractInterp 
    funct::Function 
    variables::Vector{Symbol}
end

struct CustomMod <: AbstractModifier
    interp::CustomInterpolation
end

CustomMod(expression::Function, vars::Vector{Symbol}) = CustomMod(CustomInterpolation(expression, vars))


function (i::CustomInterpolation)(α::Vector{Float64})
    nt = NamedTuple{Tuple(i.variables)}((α))
    return i.funct(nt)
end

function _interpolate_sixt_degree(x::Float64, low::Float64, high::Float64, nominal::Float64, boundary::Float64=1.0)::Float64
    t = x / boundary
    eps_plus =   high - nominal
    eps_minus =   nominal - low
    S =   0.5 * (eps_plus + eps_minus)
    A =   0.0625 * (eps_plus - eps_minus)
    return x * (S + t * A * (15 + t * t * (-10 + t * t* 3)))
end

function _interpolate_sixt_degree(x::Float64, low::Vector{Float64}, high::Vector{Float64}, nominal::Vector{Float64}, boundary::Float64=1.0)::Vector{Float64}
    t = x / boundary
    eps_plus =   high .- nominal
    eps_minus =   nominal .- low
    S =   0.5 * (eps_plus .+ eps_minus)
    A =   0.0625 * (eps_plus .- eps_minus)
    return @. x * (S + t * A * (15 + t * t * (-10 + t * t* 3)))
end