""" Geometrical shapes. """

export Point, Shape, Cylinder, Sphere, Plane, Empty, shapediff, inside

const Point = SVector{3, Float64}

abstract type Shape; end

""" A Sphere centered at (`xc`, `yc`, `zc`) with radius `R`. """
struct Sphere <: Shape
    xc::Float64
    yc::Float64
    zc::Float64
    R::Float64
end


""" 
A Cylinder between `bottom` < z < `top`, centered at (`xc`, `yc`) and with
radius `R`.
"""
struct Cylinder <: Shape
    bottom::Float64
    top::Float64
    xc::Float64
    yc::Float64
    R::Float64
end


""" 
A plane at `z` that defines a shape containing points above 
(if `up` is true) or below otherwise. 
"""
struct Plane <: Shape
    z::Float64
    up::Bool
end


""" Union of several shapes. """
struct ShapeUnion{T<:Tuple} <: Shape
    shapes::T
end

Base.union(s::Shape...) = ShapeUnion(s)

Empty() = ShapeUnion(())

""" Intersection of several shapes. """
struct ShapeIntersect{T<:Tuple} <: Shape
    shapes::T
end


Base.intersect(s::Shape...) = ShapeIntersect(s)

""" Shape1 - Shape2. """
struct ShapeSubstract{S1,S2} <: Shape
    shape1::S1
    shape2::S2
end

shapediff(s1, s2) = ShapeSubstract(s1, s2)

Base.:-(s1::Shape, s2::Shape) = shapediff(s1, s2)
Base.:+(s1::Shape, rest::Shape...) = union(s1, rest...)

# For some reason, dot for SVectors of size 3 is much slower than this
# explicit computation.
@inline @fastmath fdot(a::Point, b::Point) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]


# When we look for interface crossings sometimes we want a list of them,
# sometimes just the closest one.  To reuse the same function without
# creating too many closures we use accumulators for each possible use case
abstract type Accumulator; end

struct ListAccumulator{T <: AbstractVector} <: Accumulator
    list::T
end

@inline accum!(a::ListAccumulator, x) = push!(a.list, x)

mutable struct MinAccumulator{T <: Real} <: Accumulator
    value::T
end

@inline accum!(a::MinAccumulator, x) = (x < a.value) && (a.value = x)


""" 
Compute the intersection points between a shape `c` and the line joining
points `a` and `b`.  Each intersection point s is used in the accumulator vwhere the intersection point is
p = (1-s)a + s b
""" 
@fastmath function interfaces!(v::Accumulator, s::Sphere, a::Point, b::Point)
    (xc, yc, zc) = (s.xc, s.yc, s.zc)
    (xa, ya, za) = a
    (xb, yb, zb) = b

    # These are traditinal a, b, c in the quadratic eq. but a, b are already
    # in use.
    α = (xa - xb)^2 + (ya - yb)^2 + (za - zb)^2
    β = 2 * ((xb - xa) * (xa - xc) + (yb - ya) * (ya - yc) +
             (zb - za) * (za - zc))
    γ = (xa - xc)^2 + (ya - yc)^2 + (za - zc)^2
    Δ = β^2 - 4 * α * γ

    (Δ >= 0) || return nothing

    s1 = (-β + sqrt(Δ)) / 2α
    s2 = (-β - sqrt(Δ)) / 2α

    (0 < s1 < 1) && accum!(v, s1)
    (0 < s2 < 1) && accum!(v, s2)

    nothing
end


@fastmath function interfaces!(v::Accumulator, c::Cylinder, a::Point, b::Point)
    # Aliases to simplify the formulas
    (xa, ya, za) = a
    (xb, yb, zb) = b
    (xc, yc, R) = (c.xc, c.yc, c.R)

    st = (c.top - za) / (zb - za)
    sb = (c.bottom - za) / (zb - za)

    h2 = (xa - xb)^2 + (ya - yb)^2
    
    # Vertical lines.
    if (h2 == 0)
        
        # Check if we are inside or out of the cylinder.
        ((xa - xc)^2 + (ya - yc)^2 > R^2) && return nothing

        # If we are inside, check if the caps are between a and b
        (0 < st < 1) && accum!(v, st)
        (0 < sb < 1) && accum!(v, sb)

        return nothing
    end        
    
    Δ = R^2 * h2 - (xb * ya - xc * ya - xa * yb + xc * yb + xa * yc - xb * yc)^2
    m = xa^2 + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc)

    Δ < 0 && return nothing
    
    s1 = (m - sqrt(Δ)) / h2
    s2 = (m + sqrt(Δ)) / h2 

    z1 = (1 - s1) * a[3] + s1 * b[3]
    z2 = (1 - s2) * a[3] + s2 * b[3]

    (0 < st < 1) && (s1 < st < s2) && accum!(v, st)
    (0 < sb < 1) && (s1 < sb < s2) && accum!(v, sb)
    
    (0 < s1 < 1 && c.bottom < z1 < c.top) && accum!(v, s1)
    (0 < s2 < 1 && c.bottom < z2 < c.top) && accum!(v, s2)

    nothing
end


function interfaces!(v::Accumulator, p::Plane, a::Point, b::Point)
    # Aliases to simplify the formulas
    (xa, ya, za) = a
    (xb, yb, zb) = b
    zp = p.z
    
    s = (zp - za) / (zb - za)
    (0 < s < 1) && accum!(v, s)
    nothing
end

# Since a ShapeUnion or ShapeIntersect may be composed of shapes of different
# types, we use a generated function to unroll the loop and make the function
# type-stable.  This is cool.
@generated function interfaces!(v::Accumulator,
                                u::Union{ShapeUnion{T}, ShapeIntersect{T}},
                                a::Point, b::Point) where {T}
    L = fieldcount(T)
    out = quote end

    for i in 1:L
        push!(out.args,
              quote
              interfaces!(v, u.shapes[$i], a, b)
              end
              )
    end
    push!(out.args, :(return nothing))
    
    out
end


function interfaces!(v::Accumulator, s::ShapeSubstract, a::Point, b::Point)
    interfaces!(v, s.shape1, a, b)
    interfaces!(v, s.shape2, a, b)
    nothing
end



"""
Checks whether point `a` is inside a given shape.
"""
function inside(s::Sphere, a::Point)
    (xa, ya, za) = a
    (xc, yc, zc, R) = (s.xc, s.yc, s.zc, s.R)
    ((xa - xc)^2 + (ya - yc)^2 + (za - zc)^2 < R^2)
end


function inside(c::Cylinder, a::Point)
    (xa, ya, za) = a
    (xc, yc, R, z1, z2) = (c.xc, c.yc, c.R, c.bottom, c.top)

    ((xa - xc)^2 + (ya - yc)^2 > R^2) && return false
    (z1 < za < z2)
end


inside(p::Plane, a::Point) = (p.up == (a[3] > p.z))


# See above for why we use a generated function
@generated function inside(u::ShapeUnion{T}, a::Point) where {T}
    L = fieldcount(T)
    out = quote end

    for i in 1:L
        push!(out.args,
              quote
              inside(u.shapes[$i], a) && return true
              end
              )
    end
    push!(out.args, :(return false))

    out
end


@generated function inside(u::ShapeIntersect{T}, a::Point) where {T}
    L = fieldcount(T)
    out = quote end

    for i in 1:L
        push!(out.args,
              quote
              inside(u.shapes[$i], a) || return false
              end
              )
    end
    push!(out.args, :(return true))

    out
end


function inside(u::ShapeIntersect, a::Point)
    for item in u.shapes
        inside(item, a) || return false
    end
    true
end

inside(s::ShapeSubstract, a::Point) = inside(s.shape1, a) && !inside(s.shape2, a)
