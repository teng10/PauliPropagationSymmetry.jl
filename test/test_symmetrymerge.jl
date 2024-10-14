using Random

function test_symmetry_numerical(nq, nl, W, min_abs_coeff)
    """
    Tests for the symmetry merge functions.
    """
    Random.seed!(42) # for reproducibility
    thetax, thetaz, thetay = randn(nl), randn(nl), randn(nl)
    thetas = theta_periodic_brickwork(nq, nl, thetax, thetaz, thetay);
    # Observable to compute the expectation value of:
    symbs = [:I for _ in 1:nq]
    symbs[round(Integer, nq/2)] = :Z   # as symbol. Also works but is slower.
    obsint = symboltoint(symbs);  # for performance we work with bitoperations
    # Build single layer of translational invariant brickwork circuit
    topo = bricklayertopology(nq, periodic=true)  # periodic boundary conditions
    circ_layer = hardwareefficientcircuit(nq, 1; topology=topo)
    fastcirc_layer = tofastgates(circ_layer)
    thetas_layers = reshape(thetas, :, nl);
    symbfs = symmetrypropagate(fastcirc_layer, obsint, thetas_layers, nq, nl, W, min_abs_coeff)
    # Compare with numerical propagation
    full_circuit = hardwareefficientcircuit(nq, nl; topology=topo)
    fastfull_circuit = tofastgates(full_circuit)
    numbfs  = mergingbfs(fastfull_circuit, obsint, thetas, max_weight=W, min_abs_coeff=min_abs_coeff);    
    # Check the expectation value from symmetry propagation and numerical propagation are the same
    evsym = evalagainstzero(symbfs) # expectation
    evnum = evalagainstzero(numbfs) # expectation
    @assert isapprox(evnum, evsym; atol=1e-8) "Numerical value $evnum is not equal to symmetry value $evsym for nq = $nq"
    return evsym - evnum
end
    