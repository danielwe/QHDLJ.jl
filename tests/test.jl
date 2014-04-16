using QHDLJ

function showln(o)
	show(o)
	println("")
end

function test_feedback()
	kappas = 10. * ones(4)
	Delta = 1.
	c = single_mode_cavity(kappas, Delta)
	cfb = feedback(c, [2, 1],[1, 3])
	@assert abs(cfb.A[1,1] - (-.5*(sum(sqrt(kappas[1:3]))^2 + kappas[4])-1im*Delta)) < 1e-8
	@assert cfb.D == speye(Complex128, 2)
	@assert cfb.input_ports == ["In2", "In4"]
	@assert cfb.output_ports == ["Out3", "Out4"]
end

test_feedback()

function test_circuit_instantiation1()
	kappas = 10. * ones(2)
	Delta = 1.
	c = single_mode_cavity(kappas, Delta)
	cc = nlcircuit({"c1"=>c,"c2"=>c})
	@assert size(cc.connections, 1) == 0
	@assert size(cc.input_map, 1) == 0
	@assert size(cc.output_map, 1) == 0
end
test_circuit_instantiation1()

function test_circuit_to_component1()
	kappas = 10. * ones(2)
	Delta = 1.
	c = single_mode_cavity(kappas, Delta)
	cc = nlcircuit({"c1"=>c,"c2"=>c},
			["c1" "Out1" "c2" "In2";
			 "c2" "Out1" "c1" "In2"])
	ccc = NLComponent(cc)
	# show(ccc)
	@assert ccc.input_ports==["c1.In1", "c2.In1"]
	@assert ccc.output_ports==["c1.Out2", "c2.Out2"]	
end

test_circuit_to_component1()

function test_circuit_to_component2()
	kappas = 10. * ones(2)
	Delta = 1.
	c = single_mode_cavity(kappas, Delta)
	cc2 = nlcircuit({"c1"=>c,"c2"=>c},
			["c1" "Out1" "c2" "In2";
			 "c2" "Out1" "c1" "In2"],
			 ["c1In1" "c1" "In1";
			  "c2In1" "c2" "In1"],
			 ["c1Out2" "c1" "Out2";
			  "c2Out2" "c2" "Out2"]
	)
	cc2c = NLComponent(cc2)
	@assert cc2c.input_ports==["c1In1", "c2In1"]
	@assert cc2c.output_ports==["c1Out2", "c2Out2"]
end
test_circuit_to_component2()

function test_composite_circuit()
	kappas = 10. * ones(2)
	Delta = 1.
	c = single_mode_cavity(kappas, Delta)
	b = beamsplitter(pi/4)
	
	cb = nlcircuit({"c"=>c,"b"=>b},
			["c" "Out1" "b" "In2";
			 "b" "Out1" "c" "In2"],
			 ["cIn1" "c" "In1"],
			 ["cOut2" "c" "Out2"],
	)
	cbcnf = nlcircuit({"cb"=>cb, "c"=>c},
		["cb" "cOut2" "c" "In1"],
		["cbcIn1" "cb" "cIn1"],
		["cbcOut2" "cb" "cOut2"];
		flatten=false
	)
	cbcf = nlcircuit({"cb"=>cb, "c"=>c},
		["cb" "cOut2" "c" "In1"],
		["cbcIn1" "cb" "cIn1"],
		["cbcOut2" "cb" "cOut2"];
		flatten=true
	)
	@assert cbcnf.components |> keys |> collect |> sort ==["c", "cb"]
	@assert cbcf.components |> keys |> collect |> sort ==["c", "cb.b", "cb.c"]
	# 
	# showln(cbcnf)	
	# showln(cbcf)
	
	cbcnfc = NLComponent(cbcnf)
	cbcfc = NLComponent(cbcf)

	# showln(cbcnfc.input_ports)
	# showln(cbcfc.input_ports)
	@assert (cbcnfc.input_ports |> sort) == (cbcfc.input_ports |> sort)
	@assert (cbcnfc.output_ports |> sort) == (cbcfc.output_ports |> sort)
	
end
test_composite_circuit()

