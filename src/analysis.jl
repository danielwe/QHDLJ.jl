function jacobian(C::NLComponent, t, z)
    ret = zeros(Complex128, 2C.m, 2C.m)
    C.JANL!(t, z, sub(ret, 1:C.m, 1:C.m), sub(ret, 1:C.m, C.m+1 : 2C.m))
    ret[1:C.m,1:C.m] += C.A
    # double up
    ret[C.m+1:end, 1:C.m] = conj(ret[1:C.m, C.m+1 : 2C.m])
    ret[C.m+1:end, C.m+1:end] = conj(ret[1:C.m, 1:C.m])
    ret
end

function ode(C::NLComponent, t, z; u_t=nothing)
	ode = make_nlode_sde(C; u_t=u_t)
	zdot = zeros(Complex128, C.m)
	@assert size(z) == (C.m,)
	ode(t, z, zdot)
	zdot
end


function find_fixpoint(C::NLComponent, z0=nothing, t::Float64=0.; u_t=nothing, tol::Float64=1e-5, maxiter::Int=100, verbose=false)
	if z0 == nothing
		z0 = zeros(Complex128, C.m)
	end
    ode = make_nlode_sde(C;u_t=u_t)
    # print("Test")

    @assert size(z0) == (C.m,)
    # need to double up the vector
    zdu = [convert(Vector{Complex128},z0);conj(z0)]
    zdot = zeros(Complex128, C.m)
    ode(t, sub(zdu,1:C.m), zdot)

    iter::Int = 0
    while norm(zdot, 2) > tol && iter < maxiter
        # TODO check if we can exploit the special structure better than this
        zdu = zdu - (jacobian(C, t, sub(zdu, 1:C.m)) \ [zdot; conj(zdot)])
        # enforce double-up-ness
        for kk=1:C.m
            zdu[C.m+kk] = conj(zdu[kk])
        end
        ode(t, sub(zdu,1:C.m), zdot)
        iter += 1
    end
    if iter == maxiter
        println("Warning: Did not converge to within $tol in $iter iterations. Residual norm ||zdot||_2= $(norm(zdot,2)).")
    elseif verbose
        println("Converged in $iter steps")
    end
    return zdu[1:C.m]
end

function double_up(A, B)
	[A B; conj(B) conj(A)]
end


function transferfunction(C::NLComponent, t, z, omega=0)
	Atilde = jacobian(C, t, z)
	Ctilde = double_up(C.C, zeros(size(C.C)))
	Ctilde_conj = double_up(C.C', zeros(size(C.C)))
	Dtilde = double_up(C.D, zeros(size(C.D)))
	ident = eye(2*C.n)
	(ident - Ctilde * ((-1im * omega * ident - Atilde) \ Ctilde_conj)) * Dtilde
end

function quadrature_rep_transformation(n)
	[eye(n) 1im * eye(n); eye(n) (-1im * eye(n))]
end

function unwrap!(p)
    length(p) < 2 && return p
    for i = 2:length(p)
        d = p[i] - p[i-1]
        if abs(d) > pi
            p[i] -= floor((d+pi) / (2pi)) * 2pi
        end
    end
    return p
end

function unwrap(p)
	return unwrap!(copy(p))
end

