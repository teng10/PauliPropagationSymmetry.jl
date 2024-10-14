module PauliPropagation

using Base.Threads

include("datatypes.jl")
export PathProperties, NumericPathProperties


include("gates.jl")
export
    Gate,
    PauliGate,
    FastPauliGate,
    StaticGate,
    tofastgates

include("circuits.jl")
export
    bricklayertopology,
    get2dtopology,
    get2dstaircasetopology,
    hardwareefficientcircuit,
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
export mergingbfs

include("stateoverlap.jl")
export evalagainstzero, evalagainsinitial, zerofilter

include("surrogate.jl")
export operatortopathdict, PauliGateNode, gettraceevalorder, expectation, resetnodes

end
