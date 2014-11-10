import Base: show

show(io::IO, c::NLComponent) = begin
    inputs = join(c.input_ports, ", ")
    outputs = join(c.output_ports, ", ")
    modes = join(c.modes, ", ")
    internal = join(c.internal, ", ")
    print(io, "[NLComponent: $(c.m) modes, $(c.n) ports, modes: [$modes], inputs: [$inputs], outputs: [$outputs], internal: [$internal]]")
end

show(io::IO, c::NLCircuit) = begin
    print(io, "[NLCircuit:\n\t Components: {\n")
    for (name, val) in c.components
        print(io, "\t\t$name => ")
        show(io, val)
        print(io, ",\n")
    end
    print(io, "\t},\n\tGlobal inputs: {\n")
    for kk = 1:size(c.input_map, 1)
        fp, cn, cp = c.input_map[kk,:]
        print(io, "\t\t$fp => $cn:$cp,\n")
    end
    print(io, "\t},\n\tGlobal outputs: {\n")
    for kk = 1:size(c.output_map, 1)
        fp, cn, cp = c.output_map[kk,:]
        print(io, "\t\t$fp => $cn:$cp,\n")
    end
    
    print(io, "\t\t},\n\tConnections: [\n")
    for kk=1:size(c.connections,1)
        fn, fp, tn, tp = c.connections[kk,:]
        print(io, "\t\t$fn:$fp -> $tn:$tp,\n")
    end
    print(io, "\t]]")
end

show(io::IO, e::NLEvolution) = begin
    sde = (e.dAs != nothing)
    has_input = (e.u_t != nothing)
    
    tmin, tmax = e.tlist[1], e.tlist[end]
    nt = length(e.tlist)
    print(io, "[NLEvolution: \n")
    print(io, "\tstochastic=$sde,\n")
    print(io, "\tinput=$has_input,\n")
    print(io, "\t$nt timesteps from $tmin to $tmax,\n")
    print(io, "\tNLSystem: \n\t\t")
    show(io, e.nlcomponent)
    print(io,"]")
end