function spblocks(A,B)
    if A == nothing
        return B
    end
    if B == nothing
        return A
    end
    m, n = size(A)
    p, q = size(B)
    [A spzeros(m,q); spzeros(p,n) B]
end

function spblocks(A, B...)
    spblocks(A, apply(spblocks, B))
end

function spcat(a, b)
    if a == nothing
        return b
    end
    if b == nothing
        return a
    end
    [a; b]
end

function spcat(a, b...)
    spcat(a, apply(spcat, b))
end

function concatenate(E1::NLComponent, E2::NLComponent)
    m = E1.m + E2.m
    n = E1.n + E2.n
    q = E1.q + E2.q
    
    # Handle cases where one component is static, 
    # i.e. has no internal degrees of freedom.
    if E1.m == 0
        @assert E1.q == 0
        if E2.m > 0
            ANL_F! = E2.ANL_F!
        else 
            # both static
            @assert E2.q == 0
            ANL_F! = nothing
        end
    else
        if E2.m == 0
            @assert E2.q == 0
            ANL_F! = E1.ANL_F!
        else
            # combine nonlinear functions into 
            ANL_F! = (t, zs, w, out) -> begin 
                E1.ANL_F!(t, sub(zs,1:E1.m), sub(w,1:E1.q), sub(out,1:E1.m))
                E2.ANL_F!(t, sub(zs,E1.m:m), sub(w,E1.q:q), sub(out,E1.m:m))
                nothing
            end
        end
    end
    
    
    NLComponent(
        m,
        n,
        q,
        spblocks(E1.A, E2.A),
        spblocks(E1.B, E2.B),
        spblocks(E1.C, E2.C),
        spblocks(E1.D, E2.D),
        spcat(E1.a, E2.a),
        spcat(E1.c, E2.c),
        ANL_F!
        )
end


# helper function to find unspecified ports
function _get_missing(sorted_indices::Vector{Int}, n::Int)
    ret = zeros(Int, n-length(sorted_indices))
    jj = 1
    ll = 1
    for kk = 1:n
        if jj <= length(sorted_indices)
            while kk > sorted_indices[jj]
                jj+=1
            end
            if sorted_indices[jj] == kk
                jj += 1
                continue
            end
        end
        ret[ll] = kk
        ll += 1
    end
    ret
end

# function div_AB{T, S}(A::Factorization{T}, B::AbstractArray{S,2})
#     m = A.m
#     n = size(B, 2)
#     ret = zeros(promote_type(T, S), m, n)
#     for kk = 1:n
#         ret[:,kk] = A \ dense(B[:,kk])[:]
#     end
#     ret
# end


function portmapping(f, t, n)
    t = t[:]; f = f[:] 
    @assert size(f) == size(t)
    fbar = _get_missing(sort(f), n)
    tbar = _get_missing(sort(t), n)
    sparse([fbar; f], [tbar; t], ones(Int64,n))
end


function feedback(E::NLComponent, nfb::Int)
    # feedback of a components last nfb ports to themselves
    m=E.m; n=E.n
    @assert 0 < nfb < n
    # split up matrices into external/internal blocks.
    Be = E.B[1:m, 1:n-nfb]
    Bi = E.B[1:m, n-nfb+1:n]
    
    Ce = E.C[1:n-nfb, 1:m]
    Ci = E.C[n-nfb+1:n, 1:m]
    
    ce = E.c[1:n-nfb]
    ci = sub(E.c, n-nfb+1:n)
    
    Dee = E.D[1:n-nfb, 1:n-nfb]; Dei = E.D[1:n-nfb, n-nfb+1:n]
    Die = E.D[n-nfb+1:n, 1:n-nfb]; Dii = E.D[n-nfb+1:n, n-nfb+1:n]
    
        
    # compute the necessary inverses via an LU-factorization.
    # It turns out this can almost always (except when nfb approaches 10^4)
    # be done faster with dense representations.
    K = lufact(eye(Complex128, nfb) - dense(Dii)) 
    KCi=sparse(K \ dense(Ci))
    KDie=sparse(K \ dense(Die))
    Kc = (K \ ci)
    
    # update model matrices
    A = sparse(E.A + Bi * KCi)
    B = sparse(Be + Bi * KDie)
    
    C = sparse(Ce + Dei * KCi)
    D = sparse(Dee + Dei * KDie)

    a = E.a + Bi * Kc
    c = ce + Dei * Kc
        
    NLComponent(m, n-nfb, E.q, A, B, C, D, a, c, E.ANL_F!)
end


function feedback(E::NLComponent, f::Vector{Int}, t::Vector{Int})
    # handle general feedback case
    nfb=size(f,1)
    @assert size(t,1)==nfb
    lastnfbports = [j for j=(E.n-nfb+1:E.n)]
    p_in = portmapping(t, lastnfbports, E.n)
    p_out = portmapping(lastnfbports, f, E.n)
    
    # reduce to special case by permuting the internal channels to the end
    feedback(
        NLComponent(E.m, E.n, E.q, E.A, E.B * p_in, p_out * E.C, 
               p_out * E.D * p_in, E.a, p_out * E.c, E.ANL_F!),
        nfb)
end


function connect_components(E1::NLComponent, E2::NLComponent, f1, t2, f2, t1)
    E12 = concatenate(E1, E2)
    t2p = t2 + E1.n
    f2p = f2 + E1.n
    t = [t1[:]; t2p[:]]
    f = [f2p[:]; f1[:]]
    feedback(E12, f, t)
end

function connect_components(E1::NLComponent, E2::NLComponent, f1, t2)
    connect_components(E1, E2, f1, t2, [], [])
end

function feedforward(E1::NLComponent, E2::NLComponent)
    "E1 >> E2"
    r =[j for j=1:E1.n]
    connect_components(E1, E2, r, r)
end