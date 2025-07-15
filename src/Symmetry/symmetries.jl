### symmetries.jl
##
# This file contains functions to merge Pauli strings by symmetries.
##
###

using PauliPropagation


"""
    greedysymmetricshift(pstr::PauliStringType)

Shift a `pstr` to the right until the first non-I Pauli, used for a system with
translational symmetry.
"""        
function greedysymmetricshift(pstr::PauliStringType)
    if pstr == 0
        return pstr
    end

    # shift until the first non-zero Pauli
    while getpauli(pstr, 1) == 0
        pstr = PauliPropagation._paulishiftright(pstr)
    end

    return pstr
end

"""
    greedysymmetricshift(psum::PauliSum)

Shift and merge of a `psum` in a system with translational symmetry.
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
    _shiftandadd!(psum::PauliSum, pstr::PauliStringType, coeff)

Merges a `pstr` to an existing `psum` in-place. 
If the shifted Pauli is contained, add the coefficient, else set the new term in the `psum`.
"""
function _shiftandadd!(
    psum::PauliSum,
    pstr::PT,
    coeff
) where {PT<:PauliStringType}
    nq = psum.nqubits

    if pstr == 0
        add!(psum, pstr, coeff)
        return
    end

    for _ in 1:nq
        # shift periodically by one
        first_pauli = getpauli(pstr, 1)
        pstr = PauliPropagation._paulishiftright(pstr)
        pstr = setpauli(pstr, first_pauli, nq)

        # if contained, add and break shifting loop
        if haskey(psum.terms, pstr)
            add!(psum, pstr, coeff)
            return
        end
    end

    # the shifted Pauli is not contained, append the new term
    set!(psum, pstr, coeff)

    return
end


"""
    fullsymmetricshift(psum::PauliSum)

Exhaustive symmetric shift and merge of all pstrs in a `psum` in a system
with translational symmetry.
"""
function fullsymmetricshift(psum::PauliSum)
    shifted_psum = PauliPropagation.similar(psum)

    for (pstr, coeff) in psum
        _shiftandadd!(shifted_psum, pstr, coeff)
    end

    return shifted_psum
end


## 2D translation
"""
    _shiftup(pstr::PauliStringType, nq, nx, ny)

Shifts a `pstr` up one row in a (`nx`, `ny`) 2D grid of `nq` qubits.
This function shifts the entire bitstring up one row, 
and sets the first row of Paulis to the last row of the Paulis.
"""
function _shiftup(pstr::PauliStringType, nq, nx, ny)
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
    _shiftupandadd!(
        shifted_psum::PauliSum, pstr::PauliStringType, coeff, nx, ny
    )

Shifts a `pstr` up one row in a (`nx`, `ny`) 2D grid of `nq` qubits. 
Adds the shifted pstr to the `shifted_psum` in-place.
"""
function _shiftupandadd!(
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
        pstr = _shiftup(pstr, nq, nx, 1)

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
    fullshiftup2d(psum::PauliSum, nx::Int, ny::Int)

Merges a `psum` by shifting all pstrs up in a (`nx`, `ny`) 2D grid of Paulis,
where `nx` is the number of periodic shifts in the x-direction,
and `ny` is the number of periodic shifts in the y-direction.
"""
function fullshiftup2d(psum::PauliSum, nx, ny)
    shifted_psum = PauliPropagation.similar(psum)

    if shifted_psum.nqubits != nx * ny
        throw(
            ArgumentError("Number of qubits $(shifted_psum.nqubits) does not \n
                match grid size $(nx) x $(ny)"
            )
        )
    end    

    for (pstr, coeff) in psum
        _shiftupandadd!(shifted_psum, pstr, coeff, nx, ny)
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
