module StandardCurves

using LinearAlgebra, RecipesBase

export StandardCurve, LinearStandard

"""
Abstract type representing a standard curve. The curve can be applied to
a value by calling an instance, i.e. `y = sc(x)`. Plot recipes have been
defined to allow visualization of standard curves by calling `plot(sc,X,Y)`
"""
abstract type StandardCurve end

(sc::StandardCurve)(x) = transform(sc,x)

"""
```julia
LinearStandard(X,Y)
```
curve of the form ``y_i = m*x_i + b``
"""
struct LinearStandard <: StandardCurve
    m
    b
    function LinearStandard(X::Vector{<:Number},Y::Vector{<:Number})
        xmat = hcat(X,ones(eltype(X),length(X)))
        (m,b) = pinv(xmat) * Y
        new(m,b)
    end
end

transform(ls::LinearStandard,x) = ls.m*x + ls.b

@recipe function f(ls::LinearStandard,X,Y)
    @series begin
        seriestype := :path
        label := "fit"
        xrange = [maximum(X),minimum(X)]
        (xrange,ls.(xrange))
    end
    seriestype := :scatter
    label := "data"
    (X,Y)
end

end # module StandardCurves
