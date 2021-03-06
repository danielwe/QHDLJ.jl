module QHDLJ

export NLComponent, NLCircuit, NLSystem, nlcircuit,
       ninputs, jacobian, ode, feedback,
       linear_passive_SLH, linear_passive_static, 
       beamsplitter, phase, displace,
       single_mode_cavity, single_mode_kerr_cavity, nd_opo, two_mode_kerr_cavity,
       solve_nlsystem, make_nlode_sde, d_opo, 
       inputs, outputs, modes, internal, make_inputs,
       find_fixpoint, transferfunction_ee, transferfunction_ie, unwrap!, unwrap,
       make_ANL_F, make_JANL


using RK4

# ArrType = AbstractArray{Complex128,2}
ArrType = SparseMatrixCSC{Complex128, Int64}
VType = Vector{Complex128}

type NLComponent
    m::Int # number internal degrees of freedom
    n::Int # number of external optical inputs
    q::Int # number of additional internal independent noises
    A::ArrType
    B::ArrType
    C::ArrType
    D::ArrType
    a::VType
    c::VType
    Ci::ArrType
    Di::ArrType
    ci::VType
    ANL_F!::Function # combined nonlinear deterministic part and stochastic
    ANL_FS::Expr
    JANL!::Function
    modes::AbstractArray{AbstractString,1}
    input_ports::AbstractArray{AbstractString,1}
    output_ports::AbstractArray{AbstractString,1}
    internal::AbstractArray{AbstractString,1}
end

type NLCircuit
    components::Dict{AbstractString,NLComponent}
    connections::AbstractArray{AbstractString, 2}
    input_map::AbstractArray{AbstractString, 2}
    output_map::AbstractArray{AbstractString, 2}
end

NLSystem = Union{NLCircuit, NLComponent}

function ninputs(c::NLSystem)
    if typeof(c) == NLComponent
        c.n
    else
        sum([cc.n for cc in values(c.components)]) - size(c.connections, 1)
    end
end

ANL_F_trivial = ((t, z, w, zdot) -> nothing)
ANL_FS_trivial = quote end
JANL_trivial = ((t, z, Jxx, Jxxc) -> nothing)

# NLComponent(m, n, q, A, B, C, D, a, c, ANL_F!=nothing, subcomponents=nothing) = NLComponent(m, n, q, A, B, C, D, a, c, ANL_F!, subcomponents)
function _default_ports(prefix::AbstractString, n::Int)
    AbstractString[prefix*string(kk) for kk=1:n]
end
function _default_modes(m::Int)
    AbstractString["m"*string(kk) for kk=1:m]
end


NLComponent(m, n, q, A, B, C, D, a, c) = NLComponent(m, n, q, A, B, C, D, a, c,
                                                    spzeros(Complex128, 0, m),
                                                    spzeros(Complex128, 0, n),
                                                    zeros(Complex128, 0),
                                                    ANL_F_trivial, ANL_FS_trivial, JANL_trivial,
                                                    _default_modes(m),
                                                    _default_ports("In", n),
                                                    _default_ports("Out",n),
                                                    Array(AbstractString, (0,))
                                                    )




include("composition.jl")
include("components.jl")
include("dynamics.jl")
include("show.jl")
include("analysis.jl")
include("symbolic.jl")



end # module
