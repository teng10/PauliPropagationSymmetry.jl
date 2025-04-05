using Bits
using BitIntegers

"""
    PauliStringType

The integer types we use to represent Pauli strings. 
Pauli strings are objects like X ⊗ Z ⊗ I ⊗ Y, where each term is a Pauli acting on a qubit.
"""
const PauliStringType = Integer

"""
    PauliType

A union type for the integer types used to represent Paulis.
Paulis, also known as Pauli operators, are objects like I, X, Y, Z acting on a single qubit.
"""
const PauliType = PauliStringType

include("datatypes.jl")
include("bitoperations.jl")
include("paulioperations.jl")
include("utils.jl")