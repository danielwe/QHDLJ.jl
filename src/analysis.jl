function jacobian(C::NLComponent, t, z)
    ret = zeros(Complex128, 2C.m, 2C.m)
    C.JANL!(t, z, sub(ret, 1:C.m, 1:C.m), sub(ret, 1:C.m, C.m+1 : 2C.m))
    ret[1:C.m, 1:C.m] += full(C.A)
    # double up
    ret[C.m+1:end, 1:C.m] = conj(ret[1:C.m, C.m+1 : 2C.m])
    ret[C.m+1:end, C.m+1:end] = conj(ret[1:C.m, 1:C.m])
    ret
end

function ode(C::NLComponent, t, z; u_t=nothing)
    odefn = make_nlode_sde(C; u_t=u_t)
    zdot = zeros(Complex128, C.m)
    @assert size(z) == (C.m,)
    odefn(t, z, zdot)
    zdot
end

function find_fixpoint(C::NLComponent, z0=nothing, t::Float64=0.; u_t=nothing, tol::Float64=1e-5, maxiter::Int=100, alpha_min=.1, verbose=false)
    if z0 == nothing
        z0 = zeros(Complex128, C.m)
    end
    ode = make_nlode_sde(C;u_t=u_t)

    # print("Test")
    @assert size(z0) == (C.m,)

    # need to double up the vector
    zdu = [convert(Vector{Complex128},z0);conj(z0)]
    zdu2 = zero(zdu)
    zdot = zeros(Complex128, C.m)
    ode(t, sub(zdu,1:C.m), zdot)

    iter::Int = 0
    nzdot = norm(zdot, 2)
    alpha = 1.
    while nzdot > tol && iter < maxiter
        grad = (jacobian(C, t, sub(zdu, 1:C.m)) \ [zdot; conj(zdot)])
        nzdot2 = nzdot
        alpha = 1.
        
        while nzdot2 >= nzdot && alpha > alpha_min
            
            # perform line search by iterating alpha->1., 1./2, 1./4, ... 
            # until norm of gradient at target point decreases
            
            zdu2 = zdu - alpha * grad
            
            # enforce double-up-ness
            for kk=1:C.m
                zdu2[C.m+kk] = conj(zdu2[kk])
            end
            ode(t, sub(zdu2,1:C.m), zdot)
            nzdot2 = norm(zdot, 2)
            alpha /= 2
        end
        nzdot = nzdot2
        zdu[:] = zdu2[:]
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
    A = full(A)
    B = full(B)
    [A B; conj(B) conj(A)]
end


function transferfunction_ee(C::NLComponent, z, omega=0; t::Float64=0., quadrep=false)
    Atilde = jacobian(C, t, z)
    Ctilde = double_up(C.C, zeros(size(C.C)))
    Ctilde_conj = double_up(C.C', zeros(size(C.C')))
    Dtilde = double_up(C.D, zeros(size(C.D)))
    identn = eye(2*C.n)
    identm = eye(2*C.m)
    trafo = quadrature_rep_transformation(C.n)
    trafo_inv=inv(trafo)
    if ndims(omega)>0
        @assert ndims(omega) == 1
        res = zeros(Complex128, 2C.n, 2C.n, size(omega,1))
        for kk=1:size(omega,1)
            res[:,:,kk] = (identn - Ctilde * ((-1im * omega[kk] *identm - Atilde) \ Ctilde_conj)) * Dtilde
            if quadrep
                res[:,:,kk] = trafo_inv * res[:,:,kk] * trafo
            end
        end
        # if quadrep
        #   res = real(res)
        # end
        return res
    else
        res = (identn - Ctilde * ((-1im * omega * identm - Atilde) \ Ctilde_conj)) * Dtilde
        if quadrep
            res = trafo_inv * res * trafo
        end
        return res
    end
end

function transferfunction_ie(C::NLComponent, z, omega=0; t::Float64=0., quadrep=false)
    Atilde = jacobian(C, t, z)
    Ctilde = double_up(C.Ci, zeros(size(C.Ci)))
    Ctilde_conj = double_up(C.Ci', zeros(size(C.Ci')))
    Dtilde = double_up(C.Di, zeros(size(C.Di)))
    identn = eye(2*C.n)
    identm = eye(2*C.m)
    trafo = quadrature_rep_transformation(C.n)
    trafo_inv=inv(quadrature_rep_transformation(size(C.Di, 1)))
    if ndims(omega)>0
        @assert ndims(omega) == 1
        res = zeros(Complex128, 2*size(C.Di, 1), 2C.n, size(omega,1))
        for kk=1:size(omega,1)
            res[:,:,kk] = (identn - Ctilde * ((-1im * omega[kk] * identm - Atilde) \ Ctilde_conj)) * Dtilde
            if quadrep
                res[:,:,kk] = trafo_inv * res[:,:,kk] * trafo
            end
        end
        return res
    else
        res = (identn - Ctilde * ((-1im * omega * identm - Atilde) \ Ctilde_conj)) * Dtilde
        if quadrep
            res = trafo_inv * res * trafo
        end
        return res
    end
end

function transferfunction_ee(C::NLComponent, z, from, to, omega=0.; t::Float64=0., quadrep=false)
    
    f_indices = indexin(from, C.input_ports)
    @assert all(f_indices .> 0)
    
    t_indices = indexin(to, C.output_ports)
    @assert all(t_indices .> 0)
    
    flatten_omega = false
    
    if ndims(omega) == 0
        omega = [omega]
        flatten_omega = true
    end
    
    tfs = transferfunction_ee(C, z, omega; t=t, quadrep=quadrep)
    ret = zeros(Complex128, 2, 2, length(to), length(from), length(omega))
    
    for jj=1:length(to)
        for kk=1:length(from)
            for ll=1:length(omega)
                ret[1,1, jj, kk, ll] = tfs[t_indices[jj], f_indices[kk], ll]
                ret[1,2, jj, kk, ll] = tfs[t_indices[jj], f_indices[kk]+C.n, ll]
                ret[2,1, jj, kk, ll] = tfs[t_indices[jj]+C.n, f_indices[kk], ll]
                ret[2,2, jj, kk, ll] = tfs[t_indices[jj]+C.n, f_indices[kk]+C.n, ll]
            end
        end
    end
    
    if flatten_omega
        ret = reshape(ret, 2, 2, length(to), length(from))
    end
    
    return ret
end
        
function transferfunction_ie(C::NLComponent, z, from, to, omega=0; t::Float64=0., quadrep=false)
    
    f_indices = indexin(f, C.input_ports)
    @assert all(f_indices .> 0)
    
    t_indices = indexin(t, C.internal)
    @assert all(t_indices .> 0)
    
    flatten_omega = false
    
    if ndims(omega) == 0
        omega = [omega]
        flatten_omega = true
    end
    
    tfs = transferfunction_ie(C, z, omega; t=t, quadrep=quadrep)
    ret = zeros(Complex128, 2, 2, length(to), length(from), length(omega))
    
    for jj=1:length(to)
        for kk=1:length(from)
            for ll=1:length(omega)
                ret[1,1, jj, kk, ll] = tfs[t_indices[jj], f_indices[kk], ll]
                ret[1,2, jj, kk, ll] = tfs[t_indices[jj], f_indices[kk]+C.n, ll]
                ret[2,1, jj, kk, ll] = tfs[t_indices[jj]+C.n, f_indices[kk], ll]
                ret[2,2, jj, kk, ll] = tfs[t_indices[jj]+C.n, f_indices[kk]+C.n, ll]
            end
        end
    end
    
    if flatten_omega
        ret = reshape(ret, 2, 2, length(to), length(from))
    end
    
    return ret
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

function find_omega(tlist, data)
#     data -= mean(data)
    phis = unwrap!(angle(data))
    tbar = mean(tlist)
    phibar = mean(phis)
    
    omega = sum(tlist .* (phis-phibar))/sum(tlist .* (tlist - tbar))
    phi0 = -omega * tbar + phibar
    return omega, phi0
end


function power_spectral_density(C::NLComponent, z, omega, theta=0.; t::Float64=0.)
    tfs = transferfunction_ee(C, z, omega; t=t)
    n = C.n
    Smp = tfs[1:n,1:n,:]
    Spp = tfs[1:n,n+1:2n,:]
    Spm = tfs[n+1:2n,1:n,:]
    Smm = tfs[n+1:2n,n+1:2n,:]
    
    flatten_theta = false
    if ndims(theta)==0
        theta = [theta]
        flatten_theta = true
    end
    
    P = zeros(Float64, n, n, length(omega), length(theta))
    for kk=1:length(omega)
        
        Np = Spp[:,:,kk] * Spp[:,:,kk]'
        Nm = Spm[:,:,kk] * Spm[:,:,kk]'
        
        M1 = Smp[:,:,kk] * Spm[:,:,kk]'
        M2 = Spm[:,:,kk] * Smp[:,:,kk]'
        
        for ll=1:length(theta)
            P[:,:,kk, ll] = real(eye(n) + Np + Nm + exp(2im*theta[ll]) * M2 + exp(-2im*theta[ll]) * M1)
        end
    end
    if flatten_theta
        P = P[:,:,:,1]
    end
    P   
end

function estimate_spectrum(tlist, output_data)
    nT = length(tlist)-1
    T = tlist[end]-tlist[1]
    dt = tlist[2]-tlist[1]
    @assert abs(nT*dt - T) < dt/2
    freqs = collect(0 : nT)/T
    omegas = [freqs[int(floor(nT/2))+1:end]-freqs[end]; freqs[1:int(floor(nT/2))]]*2pi
    spec = fftshift(fft(output_data))*sqrt(dt/nT)
    omegas, spec
end

# function power_spectral_density(C::NLComponent, z, omega, outputs::AbstractArray{ASCIIString,1}, theta=0.; t::Float64=0.)
#   P = power_spectral_density(C, z, omega, theta; t=t)
#   indices = indexin(outputs, C.output_ports)
#   ret = 
# end
    
    
function linearize(C::NLComponent, z, t=0.; u_t=nothing)
    J = jacobian(C, t, z)
    zdot = ode(C, t, z; u_t=u_t)

    @assert norm(zdot)<1e-5

    A = sparse(J[1:C.m, 1:C.m])
    a = C.a * 0
    c = C.c + C.C * z
    ci = C.ci + C.Ci * z

    if u_t != nothing
        c = c + C.D * u_t(t)
        ci = ci + C.Di * u_t(theta)
    end

    Jnl = sparse(J[1:C.m, C.m+1:end])
    
    
    ANL_F! = (t, zs, w, out) -> begin 
        A_mul_B!(1., Jnl, conj(zs), 1., out)
    end

    JANL! = (t, zs, J1, J2) -> begin
        J2[1:C.m,1:C.m] = Jnl[:,:]
    end

    NLComponent(C.m, C.n, C.q, A, C.B, C.C, C.D, a, c, C.Ci, C.Di, ci, ANL_F!, JANL!, C.modes, C.input_ports, C.output_ports, C.internal)
end
    #     m::Int # number internal degrees of freedom
    #     n::Int # number of external optical inputs
    #     q::Int # number of additional internal independent noises
    #     A::ArrType
    #     B::ArrType
    #     C::ArrType
    #     D::ArrType
    #     a::VType
    #     c::VType
    # Ci::ArrType
    # Di::ArrType
    # ci::VType
    #     ANL_F!::Function # combined nonlinear deterministic part and stochastic
    #     JANL!::Function
    # modes::AbstractArray{ASCIIString,1}
    # input_ports::AbstractArray{ASCIIString,1}
    # output_ports::AbstractArray{ASCIIString,1}
    # internal::AbstractArray{ASCIIString,1}

# function io_rel(C::NLComponent;)


