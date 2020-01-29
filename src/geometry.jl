""" Geometrical shapes. """

using LinearAlgebra, CoordinateTransformations

export Shape, Cylinder, Sphere, Cone, Plane, Empty, shapediff, inside,
    TransformedShape, mathematica


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
A Cone with a vertical axis and `bottom` < z < `top`, with vertex 
at (`xv`, `yv`, `zv`) and slope m with 
sqrt((x - xv)^2 + (y - yv)^2) = m (z - zv).
"""
struct Cone <: Shape
    bottom::Float64
    top::Float64

    xv::Float64
    yv::Float64
    zv::Float64
    m::Float64
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


# Shortcuts for unions and substractiong
Base.:-(s1::Shape, s2::Shape) = shapediff(s1, s2)
Base.:+(s1::Shape, rest::Shape...) = union(s1, rest...)


""" Transformed shape after a transformation or shift. """
struct TransformedShape{T<:AbstractAffineMap, S<:Shape} <: Shape
    trans::T
    inv::T
    shape::S
    function TransformedShape(trans::T, shape::S) where
        {T<:AbstractAffineMap, S<:Shape}

        inv_ = inv(trans)
        new{T, S}(trans, inv_, shape)
    end
end


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


# Just for debugging
struct ShowAccumulator <: Accumulator
end

accum!(a::ShowAccumulator, x) = @show x


""" 
Compute the intersection points between a shape `c` and the line joining
points `a` and `b`.  Each intersection point s is used in the accumulator vwhere the intersection point is
p = (1-s)a + s b
""" 
@fastmath function interfaces!(v::Accumulator, s::Sphere, a::Point, b::Point)
    (xc, yc, zc) = (s.xc, s.yc, s.zc)
    (xa, ya, za) = a .- @SVector [s.xc, s.yc, s.zc]
    (xb, yb, zb) = b .- @SVector [s.xc, s.yc, s.zc]

    a2 = xa^2 + ya^2 + za^2
    b2 = xb^2 + yb^2 + zb^2
    ab = xa * xb + ya * yb + za * zb
    
    # These are traditional a, b, c in the quadratic eq. but a, b are already
    # in use.
    α = a2 + b2 - 2 * ab
    β = -2*a2 + 2 * ab
    γ = a2 - s.R^2

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


@fastmath function interfaces!(v::Accumulator, c::Cone, a::Point, b::Point)
    # Shift everything wrt the cone vertex
    (xa, ya, za) = a .- @SVector [c.xv, c.yv, c.zv]
    (xb, yb, zb) = b .- @SVector [c.xv, c.yv, c.zv]
    m = c.m
    m2 = m^2
    
    at2 = xa * xa + ya * ya
    bt2 = xb * xb + yb * yb
    adotb = xa * xb + ya * yb

    # with this we define a quadratic equation on s with coeffs
    α = at2 + bt2 - 2 * adotb - m2 * (za^2 + zb^2 - 2 * za * zb)
    β = -2 * at2 + 2 * adotb - m2 * (-2 * za^2 + 2 * za * zb)
    γ = at2 - m2 * za^2

    if α != 0
        Δ = sqrt(β^2 - 4 * α * γ)
        s1, s2 = (-β - Δ) / (2 * α), (-β + Δ) / (2 * α)
    else
        # Unfortunately we must allow the possibility that the line connecting
        # a and b passes though the vertex of the cone.  In that case
        # the line is either almost always inside or almost always outside
        # it.
        s1, s2 = -Inf, za / (za - zb)
    end
    
    # Make sure that s1 < s2
    (s1 < s2) || ((s1, s2) = (s2, s1))
    
    st = (c.top - za - c.zv) / (zb - za)    
    sb = (c.bottom - za - c.zv) / (zb - za)

    (0 < st < 1) && (s1 < st < s2) && accum!(v, st)
    (0 < sb < 1) && (s1 < sb < s2) && accum!(v, sb)
    
    z1 = (1 - s1) * za + s1 * zb + c.zv
    z2 = (1 - s2) * za + s2 * zb + c.zv

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

function interfaces!(v::Accumulator, s::TransformedShape,
                     a::Point, b::Point)
    at = s.inv(a)
    bt = s.inv(b)
    
    interfaces!(v, s.shape, at, bt)
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

function inside(c::Cone, a::Point)
    (c.bottom <= a[3] <= c.top) || return false
    
    (xa, ya, za) = a .- @SVector [c.xv, c.yv, c.zv]
    m = c.m

    (xa^2 + ya^2 > m^2 * za^2) && return false
    true
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

function inside(s::TransformedShape, a::Point)
    inside(s.shape, s.inv(a))
end


# To facilitate visualization we provide methods that allow the construction of
# a Mathematica expression that can be viewed with Mathematica versions above
# 11.2
mathlist(v) = "{" * join(v, ", ") * "}"
mathcall(func, v...) = func * "[" * join(v, ", ") * "]"

function mathematica(fig::Sphere)
    mathcall("Ball", mathlist([fig.xc, fig.yc, fig.zc]), fig.R)
end

function mathematica(fig::Cylinder)
    p1 = [fig.xc, fig.yc, fig.bottom]
    p2 = [fig.xc, fig.yc, fig.top]
    
    mathcall("Cylinder", mathlist([mathlist(p1), mathlist(p2)]),
             fig.R)
end

function mathematica(fig::Plane)
    mathcall("HalfSpace", mathlist([0, 0, fig.up ? -1 : 1]),
             mathlist([0, 0, fig.z]))
end

# Note that we are not currently supporting "double cones"
function mathematica(fig::Cone)
    rbottom = fig.m * abs(fig.bottom - fig.zv)
    rtop = fig.m * abs(fig.top - fig.zv)
    vertex = [fig.xv, fig.yv, fig.zv]

    if rbottom > rtop
        base = [fig.xv, fig.yv, fig.bottom]
        plane = Plane(fig.top, false)
    else
        base = [fig.xv, fig.yv, fig.top]
        plane = Plane(fig.top, true)
    end
    mathcall("RegionIntersection", 
             mathcall("Cone", mathlist([mathlist(base), mathlist(vertex)]),
                      rbottom),
             mathematica(plane))
end

function mathematica(fig::ShapeUnion)
    mathcall("RegionUnion",
             [mathematica(x) for x in fig.shapes]...)
end

function mathematica(fig::ShapeIntersect)
    mathcall("RegionIntersection",
             [mathematica(x) for x in fig.shapes]...)
end

function mathematica(fig::ShapeSubstract)
    mathcall("RegionDifference",
             mathematica(fig.shape1), mathematica(fig.shape2))
end

function mathematica(fig::TransformedShape)
    mathcall(
        mathcall("AffineTransform",
                 mathlist(
                     [mathlist([mathlist(r) for r in eachrow(fig.trans.linear)]),
                      mathlist(fig.trans.translation)])),
        mathematica(fig.shape))
end

