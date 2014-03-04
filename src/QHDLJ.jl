module QHDLJ

using RK4

type NLComponent
    m::Int # number internal degrees of freedom
    n::Int # number of external optical inputs
    q::Int # number of additional internal independent noises
    A
    B
    C
    D
    a
    c
    ANL_F! # combined nonlinear deterministic part and stochastic
end

# NLComponent(m, n, q, A, B, C, D, a, c, ANL_F!=nothing, subcomponents=nothing) = NLComponent(m, n, q, A, B, C, D, a, c, ANL_F!, subcomponents)
NLComponent(m, n, q, A, B, C, D, a, c, ANL_F =nothing) = NLComponent(m, n, q, A, B, C, D, a, c, ANL_F)

include("composition.jl")
include("components.jl")
include("dynamics.jl")


end # module
