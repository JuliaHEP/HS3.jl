using Distributions
@with_kw mutable struct BernsteinPolynomial{T<:Real} <: ContinuousUnivariateDistribution
    coefficients::Vector{T}
    a::Float64 = 0.0 # lower bound of the interval
    b::Float64 = 1.0 # upper bound of the interval
end

# Constructor that checks the input coefficients and interval
function BernsteinPolynomial(coefficients::Vector{T}, a::T=0.0, b::T=1.0) where T<:Real
    if isempty(coefficients)
        throw(ArgumentError("Coefficient array must be non-empty"))
    end
    if a >= b
        throw(ArgumentError("Lower bound a must be less than upper bound b"))
    end
    return BernsteinPolynomial{T}(coefficients, a, b)
end

# Transform x from [a, b] to [0, 1]
function normalize(d::BernsteinPolynomial, x::Real)
    return (x - d.a) / (d.b - d.a)
end

# Inverse transform from [0, 1] to [a, b]
function unnormalize(d::BernsteinPolynomial, x::Real)
    return x * (d.b - d.a) + d.a
end

# Implement the PDF function for the BernsteinPolynomial distribution
function pdf(d::BernsteinPolynomial, x::Real)
    if x < d.a || x > d.b
        return 0.0
    end
    x_normalized = normalize(d, x)
    n = length(d.coefficients) - 1
    p = 0.0
    for (i, coeff) in enumerate(d.coefficients)
        # Bernstein basis polynomial
        basis_poly = binomial(n, i-1) * x_normalized^(i-1) * (1 - x_normalized)^(n - (i-1))
        p += coeff * basis_poly
    end
    # Adjust the PDF by the width of the interval
    p /= (d.b - d.a)
    return p
end

# Implement the logpdf function for the BernsteinPolynomial distribution
function logpdf(d::BernsteinPolynomial, x::Real)
    pdf_val = pdf(d, x)
    return pdf_val > 0 ? log(pdf_val) : -Inf
end