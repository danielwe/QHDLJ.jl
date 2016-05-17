function make_nlode_sde(C::NLSystem; u_t=nothing, sde=false, symbolic_ANL=true)
    # Create an ode/sde function for a given component and input function u_t. 
    # pass sde=true to get an SDE function
    # pass _eval_fnq=false to get an unevaluated 
    #     expression object containing the source code for the function back.
    if typeof(C) != NLComponent
        C = NLComponent(C)
    end
    
    if symbolic_ANL
        ANL_F! = make_ANL_F(C)
    else
        ANL_F! = (C.ANL_F!) != nothing ? (C.ANL_F!) : ((t, z, w, zdot) -> nothing)
    end

    if u_t == nothing
        const_z = zeros(Complex128, C.n)
        u_t = t -> const_z
    end
    if sde
        w1ra, w1rb = 1, C.n
        w1ia, w1ib = C.n+1, 2C.n
        w2a, w2b = (2C.n+1),(2C.n + C.q)
        
        
        function nlsde!(t::Float64, z::AbstractVector{Complex128}, w::AbstractVector{Float64}, zdot::AbstractVector{Complex128}, sdeparams=nothing)
            zdot[:] = C.a[:]
            A_mul_B!(1., C.A, z, 1., zdot)
            A_mul_B!(.5, C.B, sub(w, w1ra : w1rb), 1., zdot)
            A_mul_B!(.5im, C.B, sub(w, w1ia : w1ib) , 1., zdot)
            A_mul_B!(1., C.B, u_t(t), 1., zdot)
            ANL_F!(t, z, sub(w, w2a : w2b), zdot)
            nothing
        end
        
    else
        function nlode!(t::Float64, z::AbstractVector{Complex128}, zdot::AbstractVector{Complex128}, odeparams=nothing)
            zdot[:] = C.a[:]
            A_mul_B!(1., C.A, z, 1., zdot)
            A_mul_B!(1., C.B, u_t(t), 1., zdot)
            ANL_F!(t, z, nothing, zdot)
            nothing
        end
    end
end

type NLEvolution
    nlcomponent
    nlde
    tlist
    zts
    inputs
    outputs
    internal
    u_t
    dAs
    dAouts
    dWs
    seed
end


function solve_nlsystem(C::NLSystem, z0::AbstractVector, tlist::AbstractVector{Float64}, hmax::Float64; sde=true, u_t=nothing, verbose=true, seed=0, symbolic_ANL=true)
    if typeof(C) != NLComponent
        C = NLComponent(C)
    end

    z0 = (1+0im) * z0

    nlde! = make_nlode_sde(C; u_t=u_t, sde=sde, symbolic_ANL=symbolic_ANL)

    inputs = zeros(Complex128, C.n, length(tlist))
    
    if u_t != nothing
        for kk=1:length(tlist)
            inputs[:,kk] = u_t(tlist[kk])
        end
    end
    ni = size(C.internal, 1)
    
    if sde
        zts, wts = rk4solve_stochastic(nlde!, z0, tlist, hmax, 2C.n + C.q; verbose=verbose, seed=seed)
        dAs = [wts[1:C.n,:] + 1im * wts[C.n+1:2C.n,:]   zeros(C.n)]/2
        dWs = [wts[2C.n+1:end,:]    zeros(C.q)]
        dAouts = C.D * dAs
        # print(size(C.C * zts), " ", size(C.c), " ", size(C.D * inputs), "\n")
        # print(typeof(C.C * zts), " ", typeof(C.c), " ", typeof(C.D * inputs), "\n")
        outputs = (C.C * zts .+ reshape(C.c, (C.n, 1))) + C.D * inputs
        internal = (C.Ci * zts .+ reshape(C.ci, (ni, 1))) + C.Di * inputs
        return NLEvolution(C, nlde!, tlist, zts, inputs, outputs, internal, u_t, dAs, dAouts, dWs, seed)
    else
        zts = rk4solve(nlde!, z0, tlist, hmax; verbose=verbose)
        outputs = (C.C * zts .+ reshape(C.c, (C.n, 1))) + C.D * inputs
        internal = (C.Ci * zts .+ reshape(C.ci, (ni, 1))) + C.Di * inputs
        return NLEvolution(C, nlde!, tlist, zts, inputs, outputs, internal, u_t, nothing, nothing, nothing, 0)
    end
end

function inputs(e::NLEvolution, names; transp=true, add_noise=false)
    if (typeof(names) <: AbstractString)
        names = [names]
    end 

    indices = indexin(names, e.nlcomponent.input_ports)
    if any(indices .== 0)
        error("Cannot find some names $(names[find(x->x==0, indices)])")
    end
    
    ret = e.inputs[indices,:]
    
    if add_noise
        ret += e.dAs[indices,:]
    end

    if transp
        transpose(ret)
    else
        ret
    end
end

function outputs(e::NLEvolution, names; transp=true, add_noise=false)
    if (typeof(names) <: AbstractString)
        names = [names]
    end 

    indices = indexin(names, e.nlcomponent.output_ports)
    if any(indices .== 0)
        error("Cannot find some names $(names[find(x->x==0, indices)])")
    end
    
    ret = e.outputs[indices,:]
    
    if add_noise
        ret += e.dAouts[indices,:]
    end
    
    if transp
        transpose(ret)
    else
        ret
    end

end

function modes(e::NLEvolution, names; transp=true)
    if (typeof(names) <: AbstractString)
        names = [names]
    end 

    indices = indexin(names, e.nlcomponent.modes)
    if any(indices .== 0)
        error("Cannot find some names $(names[find(x->x==0, indices)])")
    end

    if transp
        transpose(e.zts[indices,:])
    else
        e.zts[indices,:]
    end
end

function internal(e::NLEvolution, names; transp=true, add_noise=false)
    if (typeof(names) <: AbstractString)
        names = [names]
    end 

    indices = indexin(names, e.nlcomponent.internal)
    if any(indices .== 0)
        error("Cannot find some names $(names[find(x->x==0, indices)])")
    end
    ret = e.internal[indices,:]
    
    
    if add_noise
        ret += (e.nlcomponent.Di * e.dAs)[indices,:]
    end
    
    if transp
        ret = transpose(ret)
    end
    
    ret
end


function make_inputs(system, names, input_fn)
    if (typeof(names) <: AbstractString)
        names = [names]
        wrapinput = true
    else
        wrapinput = false
    end

    inputs = system.input_ports
    names = convert(Vector{AbstractString}, names)
    m = length(inputs)
    n = length(names)
    indices = indexin(names, inputs)
    @assert all(indices .> 0)
    mask = sparse(indices, 1:n, ones(Complex128, n), m, n)
    if wrapinput
        return t -> mask * [input_fn(t)]
    else
        return t-> mask * input_fn(t)
    end
end
