### symmetries.jl
##
# This file contains functions to merge Pauli strings by symmetries.
##
###

using PauliPropagation


"""
    Greedy symmetric shift of a pstr.
    Shifts the pstr to the right until the first non-zero Pauli is found.
"""        
function greedysymmetricshift(pstr::PauliStringType)
    if pstr == 0
        return pstr
    end

    # shift until the first non-zero Pauli
    # this is non-exhaustive, but it is fast
    while getpauli(pstr, 1) == 0
        pstr = PauliPropagation._paulishiftright(pstr)
    end

    return pstr
end

"""
    Greedy symmetric shift and merge of a PauliSum.
"""
function greedysymmetricshift(psum::PauliSum)
    shifted_psum = PauliPropagation.similar(psum)
    for (pstr, coeff) in psum
        pstr = greedysymmetricshift(pstr)
        add!(shifted_psum, pstr, coeff)
    end
    return shifted_psum
end


"""
    Merges a pstr to an existing psum in-place. 
    If the shifted Pauli is contained, add the coefficient. 
    If not, set the new term in the shifted_psum.

Args:
    shifted_psum (PauliSum): The PauliSum to merge into.
    pstr (PauliStringType): The PauliStringType to merge.
    coeff (Float64): The coefficient of the pstr.

Returns:
    None
"""
function shiftandadd!(
    shifted_psum::PauliSum,
    pstr::PT,
    coeff
) where {PT<:PauliStringType}
    nq = shifted_psum.nqubits

    if pstr == 0
        add!(shifted_psum, pstr, coeff)
        return
    end

    for _ in 1:nq
        # shift periodically by one
        first_pauli = getpauli(pstr, 1)
        pstr = PauliPropagation._paulishiftright(pstr)
        pstr = setpauli(pstr, first_pauli, nq)

        # if contained, add and break shifting loop
        if haskey(shifted_psum.terms, pstr)
            add!(shifted_psum, pstr, coeff)
            return
        end
    end

    # symmetric merge failed, append the new term
    set!(shifted_psum, pstr, coeff)

    return
end


"""
    Exhaustive symmetric shift and merge of all pstrs in a psum.
    Shifts the pstr to the right until all possible shifts are found.

Args:
    psum (PauliSum): The PauliSum to perform exhaustsive symmetric shift on.

Returns:
    shifted_psum (PauliSum): The fully merged PauliSum.
"""
function fullsymmetricshift(psum::PauliSum)
    shifted_psum = PauliPropagation.similar(psum)
    for (pstr, coeff) in psum
        shiftandadd!(shifted_psum, pstr, coeff)
    end
    return shifted_psum
end


## 2D translation
"""
    Shifts a pstr up one row in a 2D grid of Paulis.
    This function shifts the entire bitstring up one row, 
    and sets the first row of Paulis to the last row of the bitstring.
Args:
    pstr (PauliStringType): The PauliStringType to shift.
    nq (Int): The number of qubits in the PauliStringType.
    nx (Int): The number of Paulis in the x direction.
    ny (Int): The number of Paulis in the y direction.
Returns:
    PauliStringType: The shifted PauliStringType.
"""
function shiftup(pstr::PauliStringType, nq, nx, ny)

    if pstr == 0
        return pstr
    end

    # check that the shift is within max allowed ny
    if ny < 1 || ny > nq / nx
        throw(ArgumentError("ny must be between 1 and $(nq / nx), but is $ny"))
    end
    row_paulis = getpauli(pstr, 1:nx) # get the first row of Paulis
    pstr = pstr >> (2 * nx)  # shift by nx many Paulis
    pstr = setpauli(pstr, row_paulis, (nq-nx+1):nq) 
    # set the last row of Paulis to the first row of the bitstring

    return pstr
end

"""
    Shifts a pstr up one row in a 2D grid and adds to the shifted_psum in-place.

Args:
    shifted_psum (PauliSum): The PauliSum to merge into.
    pstr (PauliStringType): The PauliStringType to merge.
    coeff (Float64): The coefficient of the pstr.
    nx (Int): The number of Paulis in the x direction.
    ny (Int): The number of Paulis in the y direction.
Returns:
    None
"""
function shiftupandadd!(
    shifted_psum::PauliSum,
    pstr::PT,
    coeff,
    nx,
    ny
) where {PT<:PauliStringType}
    nq = shifted_psum.nqubits

    if pstr == 0
        add!(shifted_psum, pstr, coeff)
        return
    end

    for _ in 1:ny
        # shift periodically by one row
        pstr = shiftup(pstr, nq, nx, 1)

        # if contained, add and break shifting loop
        if haskey(shifted_psum.terms, pstr)
            add!(shifted_psum, pstr, coeff)
            return
        end
    end

    # symmetric merge failed, append the new term
    set!(shifted_psum, pstr, coeff)

    return
end

"""
    fullmergeshiftup(psum::PauliSum, transform)
    Merges a PauliSum by shifting all pstrs up in a 2D grid of Paulis.

Args:
    shifted_psum (PauliSum): The PauliSum to merge into.
    transform (Function): The function to transform the pstr.
Returns:
    shifted_psum (PauliSum): The fully merged PauliSum.    
"""
function fullmergeshiftup(psum::PauliSum, nx, ny)
    shifted_psum = PauliPropagation.similar(psum)

    # Check number of qubits
    if shifted_psum.nqubits != nx * ny
        throw(
            ArgumentError("Number of qubits $(shifted_psum.nqubits) does not \n
                match grid size $(nx) x $(ny)"
            )
        )
    end    

    for (pstr, coeff) in psum
        shiftupandadd!(shifted_psum, pstr, coeff, nx, ny)
    end

    return shifted_psum
end

## visualization tool for 2D
function PauliPropagation.inttostring(pstr::PauliStringType, nx, ny)
    str = ""

    for ii in 1:ny
        row_paulis = getpauli(pstr, (ii-1)*nx+1:ii*nx)
        for pauli in row_paulis
            str *= inttostring(pauli, nx)
        end
        str *= "\n"
    end
    return str
end
