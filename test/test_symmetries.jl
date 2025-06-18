# Test File for symmetries.jl
using Test

function get_psum(nq)
    """
    Create a PauliSum with terms that can be merged by translational symmetries.
    """
    input_psum = PauliSum(nq)
    add!(input_psum, :Z, 5)
    add!(input_psum, [:Z], [2])
    add!(input_psum, [:Z], [1], 1.5)
    add!(input_psum, [:X], [5], 0.5)
    add!(input_psum, [:X, :Z], [2, 5], 0.6)
    add!(input_psum, [:Z, :X], [2, 6], 0.5)
    
    return input_psum
end    

@testset "Full symmetric shift 1d merging" begin
    nq = 7
    input_psum = get_psum(nq)

    expected_psum = PauliSum(nq)
    add!(expected_psum, :Z, 5, 3.5)
    add!(expected_psum, [:X, :Z], [2, 5], 1.1)
    add!(expected_psum, [:X], [5], 0.5)

    merged_psum = fullsymmetricshift(input_psum)
    @test merged_psum == expected_psum

end

@testset "Greedy symmetric shift 1d merging" begin
    nq = 7
    input_psum = get_psum(nq)

    expected_psum = PauliSum(nq)
    add!(expected_psum, :Z, 1, 3.5)
    add!(expected_psum, [:Z, :X], [1, 5], 0.5)
    add!(expected_psum, [:X, :Z], [1, 4], 0.6)
    add!(expected_psum, :X, 1, 0.5)

    merged_psum = greedysymmetricshift(input_psum)
    @test merged_psum == expected_psum

end

@testset "shiftup 2d shift argument" begin
    # Test error on invalid ny
    nx, ny = 3, 2
    nq = nx * ny
    pstr = PauliString(nq, :Z, 3)
    @test_throws ArgumentError shiftup(pstr.term, pstr.nqubits, nx, 0)
    @test_throws ArgumentError shiftup(pstr.term, pstr.nqubits, nx, 3)
end

@testset "Full shiftup 2D translation merging" begin
    nx, ny = 3, 2
    nq = nx * ny
    p_init = PauliSum(nq)
    add!(p_init, :Z, 3)
    add!(p_init, :Z, 6)
    add!(p_init, :Z, 5)
    add!(p_init, [:X, :Z], [2, 5])
    add!(p_init, [:X, :Z], [5, 2], 0.5)


    expected_psum = PauliSum(nq)
    add!(expected_psum, :Z, 6, 2)
    add!(expected_psum, :Z, 5)
    add!(expected_psum, [:X, :Z], [2, 5], 1.5)

    # test merging with shiftup
    @test fullshiftup2d(p_init, nx, ny) == expected_psum

end
