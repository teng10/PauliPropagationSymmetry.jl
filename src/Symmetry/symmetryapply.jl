"""
Symmetry operations for Pauli operators.
"""
function shiftbyone(op, nq::Int)::typeof(op)
    """
    Translate the elements of a Pauli operator by one position to the right.

    Parameters:
    op (Int): The Pauli operator to shift
    nq (Int): The number of qubits

    Returns:
    The Pauli operator with its elements shifted by one position to the right.
    """
    @assert nq >= 2 ["The number of qubits ($nq) must be at least 2."]
    for ii in 2:nq
        tmp1 = getelement(op, 1)
        tmp2 = getelement(op, ii)
        op = setelement!(op, 1, tmp2)
        op = setelement!(op, ii, tmp1)
    end
    return op
end
