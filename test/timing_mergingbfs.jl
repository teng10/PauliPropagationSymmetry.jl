using PauliPropagation
using BenchmarkTools
using Random


function timingnumericalPP()
    nq = 8
    nl = 4
    W = Inf
    min_abs_coeff = 0

    pstr = PauliString(nq, :Z, round(Int, nq / 2))
    opsum = PauliSum(nq, pstr)

    topo = bricklayertopology(nq; periodic=false)
    # topo = get2dtopology(4, 4)
    circ = hardwareefficientcircuit(nq, nl; topology=topo)

    m = length(circ)

    Random.seed!(42)
    thetas = randn(m)

    res1 = propagate(circ, pstr, thetas; max_weight=W, min_abs_coeff=min_abs_coeff)
    res2 = propagate(circ, opsum, thetas; max_weight=W, min_abs_coeff=min_abs_coeff)
    @show overlapwithzero(res1), overlapwithzero(res2)
    @btime propagate($circ, $pstr, $thetas; max_weight=$W, min_abs_coeff=$min_abs_coeff)
    @btime propagate($circ, $opsum, $thetas; max_weight=$W, min_abs_coeff=$min_abs_coeff)

    return
end

# 55.345 ms (354 allocations: 3.79 MiB)
timingnumericalPP()