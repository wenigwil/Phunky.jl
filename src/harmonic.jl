struct LatticeVibrations
    fullq_freqs::Matrix{Float64}
    eigdisplacement::Array{ComplexF64,4}
    velocities::Array{Float64,3}

    function LatticeVibrations(
        ebdata::ebInputData,
        qedata::qeIfc2Output,
        deconvolution::DeconvData,
        qpoints_cryst::Matrix{Float64},
    )
        lattvecs = ebdata.crystal_info["lattvecs"]
        basisatoms2species = ebdata.crystal_info["atomtypes"]
        species2masses = ebdata.crystal_info["masses"]
        numatoms = ebdata.allocations["numatoms"]
        numbranches = 3 * numatoms

        ifc2 = qedata.properties["ifc2"]

        weightmap = deconvolution.weightmap
        uqf = deconvolution.unitpoints_qefrac_folded
        unitpoints_cart = deconvolution.unitpoints_cart

        # elphbolt input.nml has the atom mass in units of Dalton (amu) and we need 
        # them in multiples of double electron mass (Rydberg units) nm_to_bohr
        species2masses = species2masses * dalton_to_2me

        numqpoints = size(qpoints_cryst, 1)

        enforce_acoustic_sum_rule!(ifc2)

        fullq_freqs = Matrix{Float64}(undef, (numqpoints, numbranches))
        velocities = Array{ComplexF64,3}(undef, (numqpoints, numbranches, 3))
        # eigdisplacement[iq, branch, icart, iat]
        eigdisplacement =
            Array{ComplexF64,4}(undef, (numqpoints, numbranches, 3, numatoms))

        # This does not need to recomputed so it's out of the main q-point loop
        mass_prefactor = build_mass_prefactor(basisatoms2species, species2masses)

        # Reclattvecs are calculated in [1/nm]. We convert to [1/Bohr]
        reclattvecs = calc_reciprocal_lattvecs(lattvecs) / nm_to_bohr

        qpoints_cart = qpoints_cryst * permutedims(reclattvecs)

        Threads.@threads for iq in axes(qpoints_cryst, 1)
            dynmat, ∇q_dynmat = build_dynamical_matrix(
                weightmap,
                uqf,
                unitpoints_cart,
                ifc2,
                mass_prefactor,
                qpoints_cart[iq, :],
            )

            # Forcing the dynamical matrix to be hermitian
            dynmat = LinAlg.Hermitian(0.5 .* (dynmat + dynmat'))

            # This will return the eigenvalues and the diagonalizer of dynmat which 
            # is made up of column-eigenvectors. eigvecs[:,i] will yield the i-th 
            # eigenvector made up of 3*numatoms elements. 
            eigvals, eigvecs = LinAlg.eigen(dynmat)
            freqs = copysign.(sqrt.(abs.(eigvals)), eigvals)

            # Fixing the gauge
            if ~iszero(eigvecs[1, 1])
                eigvecs ./= (eigvecs[1, 1] / abs(eigvecs[1, 1]))
            end

            for ibranch in 1:numbranches
                for icart in 1:3
                    # A sesquilinear form of a hermitian matrix will be real
                    velocities[iq, ibranch, icart] = begin
                        real(
                            LinAlg.dot(
                                eigvecs[:, ibranch],
                                ∇q_dynmat[:, :, icart] * eigvecs[:, ibranch],
                            ),
                        )
                    end
                end
                velocities[iq, ibranch, :] ./= (2 * freqs[ibranch])
            end

            # First put the branch index in front and then demux cartesian and 
            # atomindex. In the muxing the cartesian index was the faster one so it 
            # will appear in as the first index of the new array after the branch 
            # index
            eigvecs = reshape(permutedims(eigvecs), (3 * numatoms, 3, numatoms))

            # Eigenvalues are the squared eigenfrequencies of the system
            fullq_freqs[iq, :] = freqs
            # According to Togo eq (6) and (7) the eigvecs just stay normalized
            eigdisplacement[iq, :, :, :] = eigvecs

            if iszero(qpoints_cryst[iq, :])
                fullq_freqs[iq, 1:3] .= [0.0, 0.0, 0.0]
                velocities[iq, :, :] .= 0.0
            end
        end

        # Convert from the Rydberg Unit System to SI
        fullq_freqs .*= Ryd_to_turnTHz
        velocities .*= Ryd_to_km_ov_s

        new(fullq_freqs, eigdisplacement, velocities)
    end
end

struct DensityOfStates
    density::Vector{Float64}
    cont_energies::Vector{Float64}
    energies::Vector{Float64}

    function DensityOfStates(
        ebdata::ebInputData,
        qedata::qeIfc2Output,
        deconvolution::DeconvData,
        numenergies::Int64,
        sampling::Tuple{Int64,Int64,Int64};
        scalebroad::Float64 = 1.0
    )
        numatoms = ebdata.allocations["numatoms"]

        qpoints_cryst = sample_cube(sampling)
        numq = size(qpoints_cryst, 1)
        harmonic = LatticeVibrations(ebdata, qedata, deconvolution, qpoints_cryst)

        # Get the energies in eV
        energies = harmonic.fullq_freqs * 1e12 * h_Js * J_to_eV
        energies = reshape(energies, 3 * numatoms * numq)

        velocities = harmonic.velocities * 10^3 * h_Js * J_to_eV
        # velocities = reshape(velocities, ())


        cont_energies = collect(range(0.0, maximum(energies) * 1.1, numenergies))


        density = zeros(Float64, size(cont_energies, 1))

        Threads.@threads for i in axes(cont_energies, 1)
            density[i] = 0.0
            for j in axes(energies, 1)
                density[i] +=
                # TODO FIX THIS
                # δ(cont_energies[i], energies[j]; smearing = smearing_type1(scalebroad, velocities[]))
            end
            density[i] /= (numq)
        end

        new(density, cont_energies, energies)
    end
end

function force_hermiticity!(mat::Matrix{ComplexF64})
    mat .= 1 / 2 * (mat + transpose(conj(mat)))
end

"""
The dynamical tensor is a rank-5 tensor with two atomic indices τ, τ' (inclusive
range 1 to `numatoms`), two cartesian indices α, α' (inclusive range 1 to 3) and one
q-point index `iq` (inclusive range 1 to `numqpoints`). The tensor will be
represented as a rank-3 tensor by muxing the two index-pairs τ, α and τ', α' into two
indices `i` and `j`. For a specified `iq` we then get a square
`3*numatoms`x`3*numatoms`-matrix.
"""
function build_dynamical_matrix(
    weightmap::Array{Float64,3},
    uqf::Array{Int64,4},
    unitpoints_cart::Matrix{Float64},
    ifc2::Array{Float64,7},
    mass_prefactor::Matrix{Float64},
    qpoint_cart::Vector{Float64},
)
    numatoms = size(mass_prefactor, 1)
    dynmat = zeros(ComplexF64, (3 * numatoms, 3 * numatoms))
    # Element-wise gradient of the dynamical matrix
    ∇q_dynmat = zeros(ComplexF64, (3 * numatoms, 3 * numatoms, 3))

    numunitpoints = size(unitpoints_cart, 1)

    for iat in 1:numatoms
        for jat in 1:numatoms
            for icart in 1:3
                i = mux2to1(iat, icart, 3)
                for jcart in 1:3
                    j = mux2to1(jat, jcart, 3)

                    for l in 1:numunitpoints
                        if weightmap[l, jat, iat] > 0
                            dynmat[i, j] += @views begin
                                ifc2[
                                    icart,
                                    jcart,
                                    iat,
                                    jat,
                                    uqf[l, jat, iat, 1],
                                    uqf[l, jat, iat, 2],
                                    uqf[l, jat, iat, 3],
                                ] *
                                exp(
                                    -im * LinAlg.dot(
                                        qpoint_cart,
                                        unitpoints_cart[l, :],
                                    ),
                                ) *
                                weightmap[l, jat, iat]
                            end

                            # We also build the derivative for the velocity!
                            ∇q_dynmat[i, j, :] += @views begin
                                -im *
                                unitpoints_cart[l, :] *
                                ifc2[
                                    icart,
                                    jcart,
                                    iat,
                                    jat,
                                    uqf[l, jat, iat, 1],
                                    uqf[l, jat, iat, 2],
                                    uqf[l, jat, iat, 3],
                                ] *
                                exp(
                                    -im * LinAlg.dot(
                                        qpoint_cart,
                                        unitpoints_cart[l, :],
                                    ),
                                ) *
                                weightmap[l, jat, iat]
                            end
                        end
                    end

                    dynmat[i, j] /= mass_prefactor[iat, jat]
                    ∇q_dynmat[i, j, :] /= mass_prefactor[iat, jat]
                end
            end
        end
    end

    return dynmat, ∇q_dynmat
end

"""
    qpoints_cryst2cart(
        qpoints_cryst::Matrix{Float64},
        reclattvecs::Matrix{Float64})

Convert a list of reciprocal lattice vectors given in crystal coordinates to ones
that are given in cartesian coordinates with units of `reclattvecs`.

# Arguments

  - `qpoints_cryst::Matrix{Float64}`: A list of q-points given in crystal
    coordinates. `qpoints_cryst[i,:]` should yield the `i`-th q-point in the list.
  - `reclattvecs::Matrix{Float64}`: Reciprocal lattice vectors. `reclattvecs[:,i]`
    will yield the `i`-th reciprocal lattice vector.
"""
function qpoints_cryst2cart(
    reclattvecs::Matrix{Float64},
    qpoints_cryst::Matrix{Float64},
)
    numqpoints = size(qpoints_cryst, 1)
    qpoints_cart = zeros(Float64, (numqpoints, 3))
    for iq in axes(qpoints_cryst, 1)
        qpoints_cart[iq, :] = reclattvecs * qpoints_cryst[iq, :]
    end

    return qpoints_cart
end

"""
    enforce_acoustic_sum_rule!(ifc2_tensor)

Enforce the acoustic sum rule on an tensor of force constants. It can be derived
from Newtons 3rd Law for the atoms in the unitcell at the origin.

The force constants should have the shape `(3, 3, nat, nat, sc[1], sc[2], sc[3])`,
where `nat` are the number of atoms in a unitcell and `sc` is a vector of integers
that correspond to the number of unitcells that make up the supercell over which
the tensor is defined.

# References

  - G. J. Ackland et. al. 1997 "Practical methods in ab initio lattice dynamics"
"""
function enforce_acoustic_sum_rule!(ifc2_tensor::Array{Float64,7})
    # Grab the number of atoms from the shape
    nat = size(ifc2_tensor, 3)

    for i in 1:3
        for j in 1:3
            for iat in 1:nat
                full_sum = @views sum(ifc2_tensor[i, j, iat, :, :, :, :])
                ifc2_tensor[i, j, iat, iat, 1, 1, 1] = @views begin
                    ifc2_tensor[i, j, iat, iat, 1, 1, 1] - full_sum
                end
            end
        end
    end

    return
end

"""
    build_mass_prefactor(
            numbasisatoms::Int64,
            basisatom2species::Vector{Int64},
            species2mass::Vector{Float64})

Build the mass prefactor that is needed in the computation of the
dynamical matrix.

The mass_prefactor is a matrix of size `(numbasisatoms,  numbasisatoms)`. Its `i`-th
diagonal element is just the mass of `species[basisatom2species[i]]` which we will
call `mass[i]` in this comment. Every other `(i,j)`-th element is
`sqrt( mass[i] * mass[j] )`
"""
function build_mass_prefactor(
    basisatom2species::Vector{Int64},
    species2mass::Vector{Float64},
)
    numbasisatoms = size(basisatom2species, 1)
    # Build a vector that holds the mass for each basisatom
    mass = Vector{Float64}(undef, (numbasisatoms))
    for basisatom in 1:numbasisatoms
        # Grab the species of each basis atom
        basisspecies = basisatom2species[basisatom]
        # Convert the species of each basis atom to its mass
        mass[basisatom] = species2mass[basisspecies]
    end

    mass_prefactor = sqrt.(mass * transpose(mass))

    return mass_prefactor
end
