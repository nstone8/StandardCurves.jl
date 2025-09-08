module StandardCurves

using LinearAlgebra, RecipesBase

export StandardCurve, LinearStandard, PolynomialStandard

"""
Abstract type representing a standard curve. The curve can be applied to
a value by calling an instance, i.e. `y = sc(x)`. Plot recipes have been
defined to allow visualization of standard curves by calling `plot(sc,X,Y)`
"""
abstract type StandardCurve end

(sc::StandardCurve)(x) = transform(sc,x)

"""
```julia
PolynomialStandard{k}(X,Y)
```
curve of the form ``y_i = sum_k(b_k*x_i^k)``
"""
struct PolynomialStandard{k} <: StandardCurve
    b_vec #weights corresponding of powers of x from 0 to k
    function PolynomialStandard{k}(X::Vector{<:Number},Y::Vector{<:Number}) where k
        @assert (k isa Int) "Degree of polynomial standard must be an Int"
        #set up problem a la: https://en.wikipedia.org/wiki/Polynomial_regression#Matrix_form_and_calculation_of_estimates
        A = [x^j for x in X, j in 0:k]
        b = pinv(A) * Y
        new(b)
    end
end

function transform(ps::PolynomialStandard{k},x) where k
    sum(enumerate(ps.b_vec)) do (i,bi)
        bi*x^(i-1) #the first term is x^0, not x^1
    end
end

"""
```julia
LinearStandard(X,Y)
```
curve of the form ``y_i = m*x_i + b``
"""
const LinearStandard = PolynomialStandard{1}

@recipe function f(ls::PolynomialStandard,X,Y)
    @series begin
        seriestype := :path
        label := "fit"
        xrange = range(start=minimum(X), length=100, stop=maximum(X))
        (xrange,ls.(xrange))
    end
    seriestype := :scatter
    label := "data"
    (X,Y)
end

end # module StandardCurves
