"""
Merge the elements of the paulinodes based on the symmetry group action.
"""
function symmetricmerge(paulinodes, nq::Int, symmetryaction)
    # Initialize a dictionary to hold the representative orbits and accumulated values
    merged_paulinodes = sizehint!(typeof(paulinodes)(), length(paulinodes))
    return symmetricmerge(paulinodes, merged_paulinodes, nq, symmetryaction)
end

# Define the partitioning function
function symmetricmerge(paulinodes, merged_paulinodes, nq::Int, symmetryaction)
    """
    Merge the elements of the paulinodes based on the symmetry group action.

    Parameters:
    paulinodes (Dict{Int, Float64}): The dictionary of Pauli operators 
        and their corresponding values
    nq (Int): The number of qubits

    Returns:
    The merged nodes with the representative pauli and accumulated coefficients.
    """
    # Iterate over each `(p, c)` pair in the input dictionary
    for (p, c) in paulinodes
        test_p = p
        is_in = false
        for _ in 1:nq
            if haskey(merged_paulinodes, test_p)
                merged_paulinodes[test_p] += c
                is_in = true
                break
            else
                test_p = symmetryaction(test_p, nq)
            end
        end
        if !is_in
            merged_paulinodes[p] = c
        end
    end
    # Clear the first dict after merging to reuse the dictionary
    empty!(paulinodes)
    return merged_paulinodes, paulinodes
end

#TODO(YT): move `W` and `min_abs_coeff` to optional kwargs
function symmetrypropagate(circ, obsint, thetas, nq, nl, W, min_abs_coeff)
    """
    Pauli propagation with symmetry-aware merging.

    Parameters:
    circ (Array{FastGate, 1}): The circuit layer
    obsint (Int): The observable to compute the expectation value of
    thetas (Array{Float64, 2}): The angles for the circuit layers
    nq (Int): The number of qubits
    nl (Int): The number of layers
    W (Int): The maximal operator weight
    min_abs_coeff (Float64): The minimum absolute coefficient to consider

    Returns:
    The expectation value of the observable.
    """
    # Initialize the dictionary with the observable
    dnum = Dict(obsint=>1.0)
    dnum_merged = sizehint!(typeof(dnum)(), length(dnum))
    for l in nl:-1:1
        # Important to start from the last layer because you reverse the circuit
        # Merge by Paulis
        dnum = mergingbfs(
            circ, dnum, thetas[:, l]; max_weight=W, min_abs_coeff=min_abs_coeff
        )
        # Merge by symmetry
        dnum, dnum_merged = symmetricmerge(dnum, dnum_merged, nq, shiftbyone)
    end
    return dnum
end
