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
	if length(B) >= 1
		spblocks(A, apply(spblocks, B))
	else
		A
	end
end

# function spcat(a, b)
#     if a == nothing
#         return b
#     end
#     if b == nothing
#         return a
#     end
#     [a; b]
# end
# 
# function spcat(a, b...)
#     spcat(a, apply(spcat, b))
# end

# function concatenate(E1::NLComponent, E2::NLComponent)
#     m = E1.m + E2.m
#     n = E1.n + E2.n
#     q = E1.q + E2.q
#     
#     # Handle cases where one component is static, 
#     # i.e. has no internal degrees of freedom.
#             # combine nonlinear functions into 
#     ANL_F! = (t, zs, w, out) -> begin 
#         E1.ANL_F!(t, sub(zs,1:E1.m), E1.q > 0 ? sub(w,1:E1.q) : nothing, sub(out,1:E1.m))
#         E2.ANL_F!(t, sub(zs,E1.m+1:m), E2.q > 0 ? sub(w,E1.q+1:q) : nothing, sub(out,E1.m+1:m))
#         nothing
#     end
#     
#     JANL! = (t, zs, J1, J2) -> begin
#         E1.JANL!(t, sub(zs, 1:E1.m), sub(J1,1:E1.m, 1:E1.m), sub(J2, 1:E1.m, 1:E1.m))
#         E2.JANL!(t, sub(zs, E1.m+1:m), sub(J1, E1.m+1:m, E1.m+1:m), sub(J2, E1.m+1:m, E1.m+1:m))
#         nothing
#     end
#     
#     NLComponent(
#         m,
#         n,
#         q,
#         spblocks(E1.A, E2.A),
#         spblocks(E1.B, E2.B),
#         spblocks(E1.C, E2.C),
#         spblocks(E1.D, E2.D),
#         spcat(E1.a, E2.a),
#         spcat(E1.c, E2.c),
#         ANL_F!,
#         JANL!
#         )
# end


# helper function to find unspecified ports
function _get_missing(sorted_indices::Vector{Int}, n::Int)
	# dump(sorted_indices)
	# dump(n)
	@assert (length(sorted_indices) == 0) || (maximum(sorted_indices) <= n)
	@assert length(sorted_indices) <= n
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


function portmapping(f::Vector{Int}, t::Vector{Int}, n::Int)
    t = t[:]; f = f[:] 
    @assert size(f) == size(t)
    fbar = _get_missing(sort(f), n)
    tbar = _get_missing(sort(t), n)
    sparse([fbar; f], [tbar; t], ones(Complex128, n))
end

# function portmapping_component(f, t, n)
# 	linear_passive_static(portmapping(f, t, n))
# end


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

    a = (E.a + Bi * Kc)[:]
    c = (ce + Dei * Kc)[:]
	
	ni = size(E.internal, 1)
	# update internal modes
	Dip_e = E.Di[1:ni, 1:n-nfb]
	Dip_i = E.Di[1:ni, n-nfb+1:n]
	
	Cip = [E.Ci + Dip_i * KCi; KCi]
	Dip = [Dip_e + Dip_i * KDie; KDie]
	cip = [E.ci + Dip_i * Kc; Kc]
	internalp = [E.internal; E.output_ports[n-nfb+1:end]]
	
    NLComponent(m, n-nfb, E.q, A, B, C, D, a, c, 
		Cip, Dip, cip,
		E.ANL_F!, E.JANL!, E.modes, 
		E.input_ports[1:n-nfb], E.output_ports[1:n-nfb], 
		internalp)
end


function feedback(E::NLComponent, f::Vector{Int}, t::Vector{Int})
    # handle general feedback case
    nfb=size(f,1)
    @assert size(t,1)==nfb
    lastnfbports = [j for j=(E.n-nfb+1:E.n)]
    p_in = portmapping(t, lastnfbports, E.n)
    p_out = portmapping(lastnfbports, f, E.n)
    
    B = E.B * p_in
    C = p_out * E.C
    D = p_out * E.D * p_in
	Di = E.Di * p_in
    # D::SparseMatrixCSC = E.D * p_in
	# dump(p_out)
	# dump(D)
	# dump(typeof(p_out))
	# dump(typeof(D))
	# # dump(methods((*), (typeof(p_out), typeof(D))))
	# dump(p_out.n)
	# dump(size(D))
	# dump(methods(size, (typeof(D),)))
	# D = p_out * D
	# D = invoke((*), (SparseMatrixCSC, SparseMatrixCSC), p_out, D)
    c = (p_out * E.c)[:]
	
	pinind = convert(Vector{Int}, p_in' * collect(1:E.n))
	input_ports = E.input_ports[pinind]

	poutind = convert(Vector{Int},p_out * collect(1:E.n))
	output_ports = E.output_ports[poutind]
	
    # reduce to special case by permuting the internal channels to the end
    feedback(
        NLComponent(E.m, E.n, E.q, E.A, B, C, D, E.a, c, E.Ci, Di, E.ci, E.ANL_F!, E.JANL!, E.modes, input_ports, output_ports, E.internal),
        nfb)
end


# function connect_components(E1::NLComponent, E2::NLComponent, f1, t2, f2, t1)
#     E12 = concatenate(E1, E2)
#     t2p = t2 + E1.n
#     f2p = f2 + E1.n
#     t = [t1[:]; t2p[:]]
#     f = [f2p[:]; f1[:]]
#     feedback(E12, f, t)
# end
# 
# function connect_components(E1::NLComponent, E2::NLComponent, f1, t2)
#     connect_components(E1, E2, f1, t2, [], [])
# end
# 
# function feedforward(E1::NLComponent, E2::NLComponent)
#     "E1 >> E2"
#     r =[j for j=1:E1.n]
#     connect_components(E1, E2, r, r)
# end
# 
# function series(E1::NLComponent, E2::NLComponent)
#     "E1 << E2"
#     r =[j for j=1:E1.n]
#     connect_components(E1, E2, [], [], r, r)
# end


function nlcircuit(components::Dict, connections::AbstractArray{ASCIIString, 2}, input_map::AbstractArray{ASCIIString, 2}, output_map::AbstractArray{ASCIIString, 2}; flatten=true)
	all_components = Dict{ASCIIString, NLComponent}()
	all_connections = copy(connections)
	new_internal_connections=Array(ASCIIString, (0, 4))
	all_inputs = copy(input_map)
	all_outputs = copy(output_map)
	
	
	for (cname, comp) in components
		if typeof(comp) == NLCircuit
			if flatten
				
				input_dict = {comp.input_map[kk,1]=> comp.input_map[kk,2:3][:] for kk=1:size(comp.input_map,1)}
				output_dict = {comp.output_map[kk,1]=> comp.output_map[kk,2:3][:] for kk=1:size(comp.output_map,1)}
				
				for (cname_int, comp_int) in comp.components
					all_components[cname*"."*cname_int] = comp_int
				end
				for kk = 1:size(all_connections,1)
					fn, fp, tn, tp = all_connections[kk,:]
					dosth = false
					
					if fn == cname
						dosth = true
						if !haskey(output_dict, fp)
							# fn = cname*"."*fn
							error("You can only connect outputs of an NLCircuit that are declared in the output_map: $fp isn't. Either convert $fn to an NLComponent first or change its output_map.")
						else
							subfn, subfp = output_dict[fp]
							fn = cname*"."*subfn
							fp = subfp
						end
					end
					if tn == cname
						dosth = true
						if !haskey(input_dict, tp)
							error("You can only connect inputs of an NLCircuit that are declared in the input_map: $tp isn't. Either convert $tn to an NLComponent first or change its input_map")
							# tn = cname*"."*tn
						else
							subtn, subtp = input_dict[tp]
							tn = cname*"."*subtn
							tp = subtp
						end
					end
					if dosth
						all_connections[kk,:] = [fn fp tn tp]
					end
				end
				new_connections = copy(comp.connections)
				
				for kk = 1:size(new_connections, 1)
					fn, fp, tn, tp = comp.connections[kk,:]
					new_connections[kk,:] = [cname*"."*fn fp cname*"."*tn tp]
				end
				new_internal_connections = [new_internal_connections; new_connections]
				
				# handle input output map translation
				for kk=1:size(all_inputs, 1)
					ip, cn, cp = all_inputs[kk,:]
					if cn == cname
						if haskey(input_dict, cp)
							subn, subp = input_dict[cp]
							cn = cname*"."*subn
							cp = subp
						else
							cn = cname*"."*cn
						end
						all_inputs[kk,:] = [ip cn cp]
					end
				end
				
				for kk=1:size(all_outputs, 1)
					op, cn, cp = all_outputs[kk,:]
					if cn == cname
						if haskey(output_dict, cp)
							subn, subp = output_dict[cp]
							cn = cname*"."*subn
							cp = subp
						else
							cn = cname*"."*cn
						end
						all_outputs[kk,:] = [op cn cp]
					end
				end
						
			else
				all_components[cname] = NLComponent(comp)
			end
		else
			all_components[cname] = comp
		end
		
	end
	all_connections = [all_connections; new_internal_connections]
	NLCircuit(all_components, all_connections, all_inputs, all_outputs)
end

nlcircuit(components::Dict; flatten=false) = nlcircuit(components, Array(ASCIIString, (0, 4)); flatten=flatten)
nlcircuit(components::Dict, connections::AbstractArray{ASCIIString, 2}; flatten=false) = nlcircuit(components, connections, Array(ASCIIString, (0, 3)), Array(ASCIIString, (0, 3)); flatten=flatten)


function NLComponent(circuit::NLCircuit)
	names, components = zip(collect(circuit.components)...)
	
	nfb = size(circuit.connections, 1)
	
	nc = length(components)
	ms = [c.m for c in components]
	ns = [c.n for c in components]
	nis = [size(c.internal,1) for c in components]
	
	qs = [c.q for c in components]
	moffsets = [0; cumsum(ms)]
	noffsets = [0; cumsum(ns)]
	nioffsets = [0; cumsum(nis)]
	qoffsets = [0; cumsum(qs)]
	
	noffsets_by_name = {name => n for (name,n) in zip(names, noffsets)}
	
	m = sum(ms)
	n = sum(ns)
	q = sum(qs)

	A = spblocks([c.A for c in components]...)
	B = spblocks([c.B for c in components]...)
	C = spblocks([c.C for c in components]...)
	D = spblocks([c.D for c in components]...)
	a = vcat([c.a for c in components]...)
	c = vcat([c.c for c in components]...)
	Ci = spblocks([c.Ci for c in components]...)
	Di = spblocks([c.Di for c in components]...)
	ci = vcat([c.ci for c in components]...)
	
	f = zeros(Int, nfb)
	t = zeros(Int, nfb)
	for kk=1:nfb
		fn, fp, tn, tp = circuit.connections[kk,:]
		fp_index = indexin([fp], circuit.components[fn].output_ports)[1]
		tp_index = indexin([tp], circuit.components[tn].input_ports)[1]
		f[kk] = noffsets_by_name[fn] + fp_index
		t[kk] = noffsets_by_name[tn] + tp_index
	end
	
	
    ANL_F! = (t, zs, w, out) -> begin 
		local moff = 0
		local qoff = 0
		for kk=1:nc
			ccc = components[kk]
			components[kk].ANL_F!(t, 
				sub(zs, 1+moff : moff + ccc.m), 
				ccc.q > 0 ? sub(w, 1 + qoff : qoff + ccc.q) : nothing,
				sub(out, 1+moff : moff + ccc.m))
			moff += ccc.m
			qoff += ccc.q
		end
		nothing
	end
									
    #     E1.ANL_F!(t, sub(zs,1:E1.m), E1.q > 0 ? sub(w,1:E1.q) : nothing, sub(out,1:E1.m))
    #     E2.ANL_F!(t, sub(zs,E1.m+1:m), E2.q > 0 ? sub(w,E1.q+1:q) : nothing, sub(out,E1.m+1:m))
    #     nothing
    # end
    
    JANL! = (t, zs, J1, J2) -> begin
		local moff = 0
		for kk=1:nc
			ccc = components[kk]
			rng = (moff + 1: moff + ccc.m)
			ccc.JANL(t, sub(zs,  rng), sub(J1, rng, rng), sub(J2, rng, rng))
		end
		nothing
    end
	
	ni = sum(nis)
	
	input_ports = Array(ASCIIString, n)
	output_ports = Array(ASCIIString, n)
	modes = Array(ASCIIString, m)
	internal = Array(ASCIIString, ni)
	
	for kk=1:nc
		nnn = names[kk]
		ccc = components[kk]
		for jj=1:ccc.n
			input_ports[noffsets[kk] + jj] = nnn*"."*ccc.input_ports[jj]
			output_ports[noffsets[kk] + jj] = nnn*"."*ccc.output_ports[jj]
		end
		for jj=1:ccc.m
			modes[moffsets[kk]+jj] = nnn*"."*ccc.modes[jj]
		end
		for jj=1:nis[kk]
			internal[nioffsets[kk] + jj] = nnn*"."*ccc.internal[jj]
		end
	end
	
	
	EEE = NLComponent(m,n,q,A,B,C,D,a,c,Ci, Di, ci, ANL_F!,JANL!,modes,input_ports,output_ports, internal)
	# show(EEE)
	# show(f)
	# show(t)
	# if connections, apply feedback
	if nfb>0
		EEE = feedback(EEE, f, t)
	end
	
	# nci = size(circuit.input_map, 1)
	# fff = collect(1:nci)
	# ttt = zeros(Int, nci)
	# for kk=1:nci
	# 	npn, cn, pn = circuit.input_map[kk,:]
	# 	ind = indexin([cn*"."*pn], input_ports)[1]
	# 	ttt[kk] = ind
	# 	if ind > 0
	# 		input_ports[ind] = npn
	# 	else
	# 		error("Could not find input port $cn:$pn in circuit.")
	# 	end
	# end
	# P_in = portmapping(fff, ttt, n)
	# 
	# nco = size(circuit.output_map, 1)
	# fff = zeros(Int, nco)
	# ttt = collect(1:nco)
	# for kk=1:nco
	# 	npn, cn, pn = circuit.output_map[kk,:]
	# 	ind = indexin([cn*"."*pn], output_ports)[1]
	# 	fff[kk] = ind
	# 	if ind > 0
	# 		output_ports[ind] = npn
	# 	else
	# 		error("Could not find output port $cn:$pn in circuit.")
	# 	end
	# end
	# P_out = portmapping(fff, ttt, n)
	# 
	# 
	# D = P_out * (D * P_in)
	# Di = Di * P_in
	# B = B * P_in
	# C = P_out * C
	# c = P_out * c
	# 
	# input_ports = input_ports[convert(Vector{Int}, P_in' * collect(1:n))]
	# output_ports = input_ports[convert(Vector{Int}, P_out' * collect(1:n))]
	# 
	EEE
end
		
					