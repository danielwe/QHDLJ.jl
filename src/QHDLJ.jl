module QHDLJ

using RK4

type NLComponent
    m::Int # number internal degrees of freedom
    n::Int # number of external optical inputs
    q::Int # number of additional internal independent noises
    A::AbstractArray{Complex128,2}
    B::AbstractArray{Complex128,2}
    C::AbstractArray{Complex128,2}
    D::AbstractArray{Complex128,2}
    a::AbstractVector{Complex128}
    c::AbstractVector{Complex128}
    ANL_F!::Function # combined nonlinear deterministic part and stochastic
    JANL!::Function
end


ANL_F_trivial = ((t, z, w, zdot) -> nothing)
JANL_trivial = ((t, z, Jxx, Jxxc) -> nothing)

# NLComponent(m, n, q, A, B, C, D, a, c, ANL_F!=nothing, subcomponents=nothing) = NLComponent(m, n, q, A, B, C, D, a, c, ANL_F!, subcomponents)
NLComponent(m, n, q, A, B, C, D, a, c) = NLComponent(m, n, q, A, B, C, D, a, c, ANL_F_trivial, JANL_trivial)


function jacobian(C::NLComponent, t, z)
    ret = zeros(Complex128, 2C.m, 2C.m)
    C.JANL!(t, z, sub(ret, 1:C.m, 1:C.m), sub(ret, 1:C.m, C.m+1 : 2C.m))
    ret[1:C.m,1:C.m] += C.A
    # double up
    ret[C.m+1:end, 1:C.m] = conj(ret[1:C.m, C.m+1 : 2C.m])
    ret[C.m+1:end, C.m+1:end] = conj(ret[1:C.m, 1:C.m])
    ret
end


include("composition.jl")
include("components.jl")
include("dynamics.jl")


end # module
