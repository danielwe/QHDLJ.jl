function jacobian(C::NLComponent, t, z)
    ret = zeros(Complex128, 2C.m, 2C.m)
    C.JANL!(t, z, sub(ret, 1:C.m, 1:C.m), sub(ret, 1:C.m, C.m+1 : 2C.m))
    ret[1:C.m,1:C.m] += C.A
    # double up
    ret[C.m+1:end, 1:C.m] = conj(ret[1:C.m, C.m+1 : 2C.m])
    ret[C.m+1:end, C.m+1:end] = conj(ret[1:C.m, 1:C.m])
    ret
end

function double_up(A, B)
	[A B; conj(B) conj(A)]
end


function transferfunction(C::NLComponent, t, z, omega=0)
	Atilde = jacobian(C, t, z)
	Ctilde = double_up(C.C, zeros(size(C.C)))
	Ctilde_conj = double_up(C.C', zeros(size(C.C)))
	Dtilde = double(C.D, zeros(size(C.D)))
	ident = eye(2*C.n)
	(ident - Ctilde * ((-1im * omega * ident - Atilde) \ Ctilde_conj)) * Dtilde
end