function make_nlode_sde(C::NLComponentl; u_t=nothing, sde=false, _eval_fnq=true)
    # Create an ode/sde function for a given component and input function u_t. 
    # pass sde=true to get an SDE function
    # pass _eval_fnq=false to get an unevaluated 
    #     expression object containing the source code for the function back.
    
    
    ANL_cmd = :(nothing)
    
    if u_t != nothing
        u_t_cmd = begin
                    A_mul_B!(1., C.B, u_t(t), 1., zdot)
        end
    else
        u_t_cmd = :(nothing)
    end
    
    if sde
        
        # pre-generate the ranges for how the noise vector decomposes
        # into the real and imaginare quadrature parts of the input noises
        rng1 = (1, C.n)
        rng2 = (C.n+1, 2C.n)
        rng3 = ((2C.n+1),(2C.n + C.q))
        w1r = :(sub(w, colon($(rng1[1]), $(rng1[2]))))
        w1i = :(sub(w, colon($(rng2[1]), $(rng2[2]))))
        w2 = :(sub(w, colon($(rng3[1]), $(rng3[2]))))
        
        sde_cmd = quote
            # extra factor of .5
            A_mul_B!(.5, C.B, $w1r, 1., zdot)
            A_mul_B!(.5im, C.B, $w1i, 1., zdot)
        end
        
        if C.ANL_F! != nothing
            ANL_cmd = quote
                # Call non-linear part
                C.ANL_F!(t, z, $w2, zdot)
            end
        end
        fnq = quote
            function nlsde!(t::Float64, z::AbstractVector{Complex128}, w::AbstractVector{Float64}, zdot::AbstractVector{Complex128}, sdeparams=nothing)
                zdotd[:] = C.a
                A_mul_B!(1., C.A, z, 1., zdot)
                $sde_cmd
                $u_t_dcmd
                $ANL_cmd
                nothing
            end
        end
    else
        if C.ANL_F! != nothing
            ANL_cmd = quote
                C.ANL_F!(t, z, nothing, zdot)
            end
        end
        fnq = quote
            function nlode!(t::Float64, z::AbstractVector{Complex128}, zdot::AbstractVector{Complex128}, sdeparams=nothing)
                zdot[:] = C.a
                A_mul_B!(1., C.A, z, 1., zdot)
                $u_t_cmd
                $ANL_cmd
                nothing
            end
        end
    end
    if _eval_fnq
        return eval(fnq) 
    else
        return fnq
    end
end


type NLEvolution
    tlist
    zts
    inputs
    outputs
    u_t
    dAs
    dAouts
    dWs
end


function solve_nlcircuit(C::NLComponent, z0::AbstractVector{Complex128}, tlist::AbstractVector{Float64}, hmax::Float64; sde=true, u_t=nothing)
    nlde! = make_nlode_sde(C; u_t=u_t, sde=sde)
    inputs = zeros(C.n, length(tlist))
    
    if u_t != nothing
        for kk=1:length(tlist)
            inputs[:,kk] = u_t(tlist[kk])
        end
    end
    outputs = C.C * zts + C.c + C.D * inputs
    
    if sde
        zts, wts = rk4solve_stochastic(nlde!, z0, tlist, hmax, 2C.n + C.q)
        dAs = [wts[1:C.n,:] + 1im * wts[C.n+1:2C.n,:]   zeros(C.n)]
        dWs = [wts[2C.n+1:end,:]    zeros(C.q)]
        dAouts = C.D * dAs
        
        return NLEvolution(tlist, zts, inputs, outputs, u_t, dAs, dAouts, dWs)
    else
        zts = rk4solve(nlde!, z0, tlist, hmax)
        return NLEvolution(tlist, zts, inputs, outputs, u_t, nothing, nothing)
    end
end