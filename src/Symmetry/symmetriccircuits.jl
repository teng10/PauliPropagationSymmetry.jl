"""
Symmetric circuits.
"""

function theta_periodic_brickwork(
    nq::Integer,
    nl::Integer,
    thetax::Vector{Float64},
    thetaz::Vector{Float64},
    thetay::Vector{Float64},
)::Vector{Float64}
    """
    Build parameters for all layers in a translational invariant brickwork circuit
    with periodic boundary conditions.

    Angles are appended according to the gate sequence:
    [(thetax[1], thetaz[1], thetax[1])] * nq + [thetay[1]] * nq
    + [(thetax[2], thetaz[2], thetax[2])] * nq + [thetay[2]] * nq +...

    Parameters:
    nq: number of qubits
    nl: number of layers
    thetax: angles for X gates
    thetaz: angles for Z gates
    thetay: angles for YY gates

    Returns:
    Concatenated thetas.
    """
    thetas = []
    for l in 1:nl
        for i in 1:nq
            push!(thetas, thetax[l])
            push!(thetas, thetaz[l])
            push!(thetas, thetax[l])
        end
        for i in 1:nq
            push!(thetas, thetay[l])
        end
    end
    return thetas
end
