module PauliPropagation

using Base.Threads

include("datatypes.jl")
export PathProperties, NumericPathProperties


include("Gates/Gates.jl")
export
    Gate,
    PauliGate,
    FastPauliGate,
    tofastgates,
    apply,
    applynoncummuting,
    CliffordGate,
    default_clifford_map,
    reset_clifford_map!,
    applywithmap

include("circuits.jl")
export
    bricklayertopology,
    get2dtopology,
    get2dstaircasetopology,
    hardwareefficientcircuit,
    efficientsu2circuit,
    tfitrottercircuit,
    heisenbergtrottercircuit,
    su4ansatz,
    qcnnansatz,
    appendSU4!

include("./PauliAlgebra/PauliAlgebra.jl")
export
    inttosymbol,
    symboltoint,
    inttostring,
    getelement,
    setelement!,
    show,
    containsXorY,
    containsYorZ

include("./Symmetry/Symmetry.jl")
export
    theta_periodic_brickwork,
    shiftbyone,
    symmetricmerge,
    symmetrypropagate

include("apply.jl")
export apply, applynoncummuting  # What should I export here?

include("truncations.jl")

include("Propagation/Propagation.jl")
export mergingbfs, applygatetoall!, applygatetoone!

include("stateoverlap.jl")
export evalagainstzero, evalagainsinitial, zerofilter

include("surrogate.jl")
export operatortopathdict, PauliGateNode, gettraceevalorder, expectation, resetnodes

end
