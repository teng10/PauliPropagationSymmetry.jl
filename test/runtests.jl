using PauliPropagation
using Test
using Random

@testset "PauliPropagation.jl" begin

    include("test_propagate.jl")

    include("test_datatypes.jl")

    include("test_paulialgebra_utils.jl")

    include("test_noisechannels.jl")

    include("test_circuits.jl")

    include("test_cliffordgates.jl")

    include("test_frozengates.jl")

    include("test_miscgates.jl")

    include("test_overlaps.jl")

    include("test_paulirotations.jl")

    include("test_paulioperations.jl")

    include("test_paulitransfermaps.jl")

    include("test_truncations.jl")

    include("test_numericalcertificates.jl")


    include("test_symmetrymerge.jl")
    @test test_symmetry_numerical(4, 2, Inf, 0.0) ≈ 0.0 atol=1e-8

end
