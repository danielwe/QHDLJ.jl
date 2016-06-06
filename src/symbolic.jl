using Calculus: differentiate

function substitute(ex::Expr, substitutions::Dict)
    exs = get(substitutions, ex, nothing)
    if exs != nothing
        return exs
    end
    new_args = [substitute(arg, substitutions) for arg in ex.args]
    Expr(ex.head, new_args...)
end

function substitute(ex::Symbol, substitutions::Dict)
    get(substitutions, ex, ex)
end

function substitute(ex, substitutions::Dict)
    ex
end


function make_ANL_F(ANL_FS::Expr)
    eval(quote
        function ANL_F!(t::Float64, z::AbstractVector{Complex128}, w, out::AbstractVector{Complex128})
            $ANL_FS
        end
    end)
end


function make_ANL_F(C::NLComponent)
    make_ANL_F(C.ANL_FS)
end


function make_JANL(JANLS::Expr)
    eval(quote
        function JANL!(t::Float64, z::AbstractVector{Complex128}, out1::AbstractArray{Complex128, 2}, out2::AbstractArray{Complex128, 2})
            $JANLS
        end
    end)
end

function make_JANL(C::NLComponent)
    make_JANL(generate_JANLS(C.ANL_FS, C.m))
end


function compose_ANL_FS(ANL_FSs, ms, qs)
    N = size(ANL_FSs, 1)
    @assert N == length(ms) == length(qs)
    moffset = ms[1]
    qoffset = qs[1]
    @assert ANL_FSs[1].head == :block
    assignments = collect(filter(x -> isa(x, Expr), ANL_FSs[1].args))
    for jj=2:N
        substitutions = Dict([:(z[$kk]) => :(z[$(kk+moffset)]) for kk = 1:ms[jj]])
        merge!(substitutions, Dict([
            :(out[$kk]) => :(out[$(kk+moffset)]) for kk = 1:ms[jj]
            ]))
        merge!(substitutions, Dict([
            :(w[$kk]) => :(w[$(kk+qoffset)]) for kk = 1:qs[jj]
            ]))
        @assert ANL_FSs[jj].head == :block
        append!(assignments, 
                       collect(filter(x -> isa(x, Expr), substitute(ANL_FSs[jj], substitutions).args)))
        moffset += ms[jj]
        qoffset += qs[jj]
    end
    Expr(:block, assignments...)
end


function compose_ANL_FS(Cs::AbstractVector{NLComponent})
    ANL_FSs = [C.ANL_FS for C in Cs]
    ms = [C.m for C in Cs]
    qs = [C.q for C in Cs]
    compose_ANL_FS(ANL_FSs, ms, qs)
end


function generate_JANLS(ANL_FS::Expr, m::Int)
    subs = Dict([
        :(z[$kk]) => symbol("z$kk") for kk=1:m
    ])
    merge!(subs, Dict([
        :(conj(z[$kk])) => symbol("z$(kk)__conj") for kk=1:m
        ]))

    inverse_subs = Dict([subs[kkk] => kkk for kkk in keys(subs)])

    @assert ANL_FS.head == :block

    zdots = Array(Any, m)
    seen = Set{Int64}()
    for kk = 1:length(ANL_FS.args)
        if !isa(ANL_FS.args[kk], Expr)
          if !isa(ANL_FS.args[kk], LineNumberNode) 
            println("Cannot match $(ANL_FS.args[kk])")
          end
          continue
        end
        if ANL_FS.args[kk].head == :(+=)
            lhs, rhs = ANL_FS.args[kk].args
            @assert lhs.head == :ref
            @assert lhs.args[1] == :out
            ll = lhs.args[2]
            if ll in seen
                error("$ll appears to be in $(string(ANL_FS)) twice")
            end
            union!(seen, [ll])
            zdots[ll] = substitute(rhs, subs)
        end
    end
    unseen = setdiff!(Set(1:m), seen)
    for ll in unseen
        zdots[ll] = 0
    end

    J1 = [substitute(differentiate(zdots[kk], symbol("z$ll")), inverse_subs) for kk=1:m, ll=1:m]
    J2 = [substitute(differentiate(zdots[kk], symbol("z$(ll)__conj")), inverse_subs) for kk=1:m, ll=1:m]
    
    J1S = [(J1[kk,ll]!= 0 ? :(out1[$kk,$ll]=$(J1[kk,ll])) : 0) for kk=1:m, ll=1:m]
    J2S = [(J2[kk,ll]!= 0 ? :(out2[$kk,$ll]=$(J2[kk,ll])) : 0) for kk=1:m, ll=1:m]
    Expr(:block, filter(x -> x !=0 , [J1S[:]; J2S[:]])...)
end
