module QHDLJ

export NLComponent, NLCircuit, NLSystem, nlcircuit,
	   ninputs, jacobian, feedback,
	   linear_passive_SLH, linear_passive_static, 
	   beamsplitter, phase, displace,
	   single_mode_cavity, single_mode_kerr_cavity, nd_opo, 
	   solve_nlsystem, make_nlode_sde,
	   inputs, outputs, modes, internal


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
    JANL!::Function
	modes::AbstractArray{ASCIIString,1}
	input_ports::AbstractArray{ASCIIString,1}
	output_ports::AbstractArray{ASCIIString,1}
	internal::AbstractArray{ASCIIString,1}
end

type NLCircuit
	components::Dict{ASCIIString,NLComponent}
	connections::AbstractArray{ASCIIString, 2}
	input_map::AbstractArray{ASCIIString, 2}
	output_map::AbstractArray{ASCIIString, 2}
end

NLSystem = Union(NLCircuit, NLComponent)

function ninputs(c::NLSystem)
	if typeof(c) == NLComponent
		c.n
	else
		sum([cc.n for cc in values(c.components)]) - size(c.connections, 1)
	end
end




ANL_F_trivial = ((t, z, w, zdot) -> nothing)
JANL_trivial = ((t, z, Jxx, Jxxc) -> nothing)

# NLComponent(m, n, q, A, B, C, D, a, c, ANL_F!=nothing, subcomponents=nothing) = NLComponent(m, n, q, A, B, C, D, a, c, ANL_F!, subcomponents)
function _default_ports(prefix::ASCIIString, n::Int)
	ASCIIString[prefix*string(kk) for kk=1:n]
end
function _default_modes(m::Int)
	ASCIIString["m"*string(kk) for kk=1:m]
end


NLComponent(m, n, q, A, B, C, D, a, c) = NLComponent(m, n, q, A, B, C, D, a, c,
													spzeros(Complex128, 0, m),
													spzeros(Complex128, 0, n),
													zeros(Complex128, 0),
													ANL_F_trivial, JANL_trivial,
													_default_modes(m),
													_default_ports("In", n),
													_default_ports("Out",n),
													Array(ASCIIString, (0,))
													)




include("composition.jl")
include("components.jl")
include("dynamics.jl")
include("show.jl")
include("analysis.jl")



end # module
