function linear_passive_SLH(S, L, H, l = nothing, h = nothing)
    n = size(S, 1)
    if L == nothing
        L = zeros(Complex128, n, 0)
        @assert H == nothing
        H = zeros(Complex128, 0, 0)
        m = 0
    else
        m = size(H, 1)
        @assert size(S, 2) == n == size(L, 1)
        @assert size(H, 2) == m == size(L, 2)
    end
    if l == nothing
        l = zeros(Complex{Float64}, n, 1)
    else
        @assert size(l, 1) == n
    end
    if h == nothing
        h = zeros(Complex{Float64}, m, 1)
    else
        @assert size(h, 1) == m
    end
    A = sparse(-1.0im * H - .5 * L' * L)
    B = sparse(-(1.0+0im) * L' * S)
    C = sparse((1.0+0im) * L)
    D = sparse((1.0+0im) * S)
    a = (-1im * h - .5 * L' * l)[:]
    c = ((1.0+0im) * l)[:]
    
    @assert size(a,1) == m
    @assert size(c,1) == n
    
    
    NLComponent(m, n, 0, A, B, C, D, a, c)
end


function linear_passive_static(S, l=nothing)
    linear_passive_SLH(S, nothing, nothing, l, nothing)
end



function beamsplitter(theta, real=true)
    if real
        S = [cos(theta) -sin(theta); sin(theta) cos(theta)] * (1+0im)
    else
        S = [cos(theta) 1im * sin(theta); 1im * sin(theta)  cos(theta)]
    end
    linear_passive_static(S)
end


function phase(phi)
	linear_passive_static(speye(1) * exp(1im*phi))
end


function displace(amplitudes)
	if ndims(amplitudes) == 0
		amplitudes=[amplitudes]
	end
    n = size(amplitudes, 1)
    l = reshape(amplitudes, (n, 1))
    S = eye(n)
    linear_passive_static(S, l)
end


function single_mode_cavity(kappas, Delta)
    n = length(kappas)
    linear_passive_SLH(speye(n), sqrt(reshape(kappas, (n, 1))), speye(1) * Delta)
end

function single_mode_kerr_cavity(kappas, Delta, chi, wigner_correction = true)
    E = single_mode_cavity(kappas, Delta)
    function ANL_F!(t::Float64, z::AbstractVector{Complex128}, w, out::AbstractVector{Complex128})
        out[1] += chi/1im * (conj(z[1]) * z[1] - 1.) * z[1]
    end
    function JANL!(t::Float64, z::AbstractVector{Complex128}, out1::AbstractArray{Complex128, 2}, out2::AbstractArray{Complex128, 2})
        out1[1,1] = 2chi/1im * (conj(z[1]) * z[1] - .5)
        out2[1,1] = chi/1im * z[1]^2
    end
    
    
    E.ANL_F! = ANL_F!
    E.JANL! = JANL!
    E
end

function nd_opo(kappas, chi, Deltas=zeros(2))
    E = linear_passive_SLH(speye(3), sparse(sqrt(eye(3) .* kappas)), [Deltas[1] 0 0; 0 Deltas[2] 0; 0 0 0])
    function ANL_F!(t::Float64, z::AbstractVector{Complex128}, w, out::AbstractVector{Complex128})
        out[1] += chi * z[3] * conj(z[2])
        out[2] += chi * z[3] * conj(z[1])
        out[3] += -chi * z[1] * z[2]
    end
    function JANL!(t::Float64, z::AbstractVector{Complex128}, out1::Array{Complex128, 2}, out2::Array{Complex128, 2})
        out1[:,:] = 0
        out1[1,3] = chi * conj(z[2])
        out1[2,3] = chi * conj(z[1])
        out1[3,1] = -chi * z[2]
        out1[3,2] = -chi * z[1]

        out2[:,:] = 0
        out2[1,2] = chi * z[3]
        out2[2,1] = chi * z[3]
    end
    
	names = String["signal", "idler", "pump"]
    
    E.ANL_F! = ANL_F!
    E.JANL! = JANL!
	E.modes = names
	E.input_ports = names
	E.output_ports = names
    E
end

# TODO add free carrier models (including Kerr and thermal effects), two-mode cavities, two-way components
