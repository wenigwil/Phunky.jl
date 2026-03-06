struct Linewidth
    linewidths::Array{Float64,3}

    function Linewidth(
        ebdata::ebInputData,
        deconvolution::DeconvData,
        sodata::qeIfc2Output,
        todata::Ifc3Output,
        q1_cryst::Matrix{Float64},
        cont_freqs::Vector{Float64},
        T::Float64;
        brillouin_sampling::Tuple{Int64,Int64,Int64} = (30, 30, 30),
    )

        # System description
        numatoms = ebdata.allocations["numatoms"]
        numbranches = 3 * numatoms
        # TODO: What do we do about the units of the atomic masses?
        type2mass = ebdata.crystal_info["masses"]
        atindex2type = ebdata.crystal_info["atomtypes"]
        # Lattice Vectors come in [nm]
        lattvecs = ebdata.crystal_info["lattvecs"]
        reclattvecs = calc_reciprocal_lattvecs(lattvecs)
        # Calculate this factor in order to avoid recalculation during the loops
        # Temperature comes into here in [K]
        kbT = kb_J_ov_K * J_to_eV * T

        # Third Order Interatomic Force Constants
        # We convert from [eV/Angstrom^3] to [eV/nm^3]
        ifc3_tensor = todata.properties["ifc3_tensor"] .* 1000
        # Extract the data mapping the triplet index to...
        # Positions in [Angstrom]. We convert to [nm]
        trip2position_j = todata.properties["trip2position_j"] ./ 10
        trip2position_k = todata.properties["trip2position_k"] ./ 10
        # Atomic Indices of the whole triplet
        trip2atomindeces = todata.properties["trip2atomindices"]

        # We will include a mass normalization into the ifc3
        for itrip in axes(ifc3_tensor, 1)
            ifc3_tensor[itrip, :, :, :] ./= begin
                sqrt(
                    type2mass[atindex2type[trip2atomindeces[itrip, 1]]] *
                    type2mass[atindex2type[trip2atomindeces[itrip, 2]]] *
                    type2mass[atindex2type[trip2atomindeces[itrip, 3]]],
                )
            end
        end

        # The phonopy ifc3-file gives us the cell coordinates of the two displaced 
        # atoms. We make sure these are EXACTLY (numerics, huh) on our grid defined 
        # by the lattice vectors.
        snap_to_lattvecs!(lattvecs, trip2position_k)
        snap_to_lattvecs!(lattvecs, trip2position_j)

        # Calculate and reshape the frequencies and eigenvectors of 3 phonons by a 
        # given sampling of the brillouin zone.
        states = HarmonicStatesData(
            ebdata,
            sodata,
            deconvolution,
            q1_cryst;
            brillouin_sampling,
        )

        # Convert q-point for all states into cartesian coordinates
        q2_cart = states.q2_cryst * permutedims(reclattvecs)
        q3_abso_cart = Array{Float64,3}(undef, size(states.q3_abso_cryst))
        q3_emit_cart = Array{Float64,3}(undef, size(states.q3_emit_cryst))
        for iq1 in axes(q1_cryst, 1)
            q3_abso_cart[:, :, iq1] .=
                states.q3_abso_cryst[:, :, iq1] * permutedims(reclattvecs)
            q3_emit_cart[:, :, iq1] .=
                states.q3_emit_cryst[:, :, iq1] * permutedims(reclattvecs)
        end

        numq1 = size(q1_cryst, 1)
        numq2 = size(states.q2_cryst, 1)
        numfreq = size(cont_freqs, 1)

        @info """
        anharmonic.jl: Calculating scattering rates for all λ and frequencies...
        """
        # Calculating the scattering rates
        linewidths = Array{Float64,3}(undef, (numfreq, numq1, numbranches))
        Threads.@threads for ifreq in 1:numfreq
            println("At ifreq=", ifreq)
            ω_cont = cont_freqs[ifreq] * turnTHz_to_eV

            for λ in axes(states.q1_evec, 1)
                s1, iq1 = demux1to2(λ, numq1)
                ω = states.q1_freqs[λ] * turnTHz_to_eV

                # At every loop-nesting depth we will check if the energy of the 
                # state at the collective index vanishes. If so, we skip that 
                # iteration because the scattering rate will diverge in this case.
                if iszero(ω)
                    linewidths[ifreq, iq1, s1] = 0.0
                    continue
                end

                linewidths[ifreq, iq1, s1] = begin
                    calc_Λ(
                        λ,
                        ω_cont,
                        kbT,
                        states,
                        q2_cart,
                        q3_abso_cart,
                        q3_emit_cart,
                        ifc3_tensor,
                        trip2atomindeces,
                        trip2position_j,
                        trip2position_k,
                        numbranches,
                        reclattvecs,
                        brillouin_sampling,
                    ) *
                    hbar_Js *
                    J_to_eV *
                    pi / (2 * numq2 * ω)
                end
            end
        end

        new(linewidths)
    end
end

function calc_Λ(
    λ::Int64,
    ω_cont::Float64,
    kbT::Float64,
    states::HarmonicStatesData,
    q2_cart::Matrix{Float64},
    q3_abso_cart::Array{Float64,3},
    q3_emit_cart::Array{Float64,3},
    ifc3_tensor::Array{Float64,4},
    trip2atomindices::Matrix{Int64},
    trip2position_j::Matrix{Float64},
    trip2position_k::Matrix{Float64},
    numbranches::Int64,
    reclattvecs::Matrix{Float64},
    sampling::Tuple{Int64,Int64,Int64},
)::Float64
    numq1 = size(states.q1_cryst, 1)
    numq2 = size(states.q2_cryst, 1)

    _, iq = demux1to2(λ, numq1)
    W_λ = states.q1_evec[λ, :, :]

    Λ = 0.0
    for λ′ in axes(states.q2_evec, 1)
        _, iq′ = demux1to2(λ′, numq2)

        ω′ = states.q2_freqs[λ′] * turnTHz_to_eV

        # At every loop-nesting depth we will check if the energy of the state at the
        # collective index vanishes. If so, we skip that iteration because the 
        # scattering rate will diverge in this case.
        if iszero(ω′)
            continue
        end

        q′ = q2_cart[iq′, :]
        W_λ′ = states.q2_evec[λ′, :, :]
        v′ = states.q2_velos[λ′, :]

        for s′′ in 1:numbranches
            λ′′ = mux2to1(s′′, iq′, numq2)

            ω′′_abso = states.q3_abso_freqs[λ′′, iq] * turnTHz_to_eV
            ω′′_emit = states.q3_emit_freqs[λ′′, iq] * turnTHz_to_eV

            if iszero(ω′′_abso) || iszero(ω′′_emit)
                continue
            end

            v′′_abso = states.q3_abso_velos[λ′, :, iq]
            v′′_emit = states.q3_emit_velos[λ′, :, iq]

            smearing_abso = smearing_type3(v′, v′′_abso, reclattvecs, sampling)
            smearing_emit = smearing_type3(v′, v′′_emit, reclattvecs, sampling)

            energy_conversed_plus, delta_plus =
                check_energy_conservation_δ(ω_cont, ω′′_abso - ω′, smearing_abso)

            energy_conversed_minus, delta_minus =
                check_energy_conservation_δ(ω_cont, ω′′_emit + ω′, smearing_emit)

            q′′_abso = q3_abso_cart[iq′, :, iq]
            q′′_emit = q3_emit_cart[iq′, :, iq]

            W_λ′′_abso = states.q3_abso_evec[λ′′, :, :, iq]
            W_λ′′_emit = states.q3_emit_evec[λ′′, :, :, iq]

            if energy_conversed_plus && energy_conversed_minus
                statistics_abso = begin
                    bose(ω′, kbT) - bose(ω′′_abso, kbT)
                end

                statistics_emit = begin
                    bose(ω′, kbT) + bose(ω′′_emit, kbT) + 1
                end

                Λplus = begin
                    statistics_abso / (ω′′_abso * ω′) *
                    calc_V2plus(
                        q′,
                        q′′_abso,
                        W_λ,
                        W_λ′,
                        W_λ′′_abso,
                        ifc3_tensor,
                        trip2atomindices,
                        trip2position_j,
                        trip2position_k,
                    ) *
                    delta_plus
                end

                Λminus = begin
                    statistics_emit / (ω′′_emit * ω′) *
                    calc_V2minus(
                        q′,
                        q′′_emit,
                        W_λ,
                        W_λ′,
                        W_λ′′_emit,
                        ifc3_tensor,
                        trip2atomindices,
                        trip2position_j,
                        trip2position_k,
                    ) *
                    delta_minus
                end
                Λ += Λplus + 0.5 * Λminus
            elseif energy_conversed_plus
                statistics_abso = begin
                    bose(ω′, kbT) - bose(ω′′_abso, kbT)
                end

                Λplus = begin
                    statistics_abso / (ω′′_abso * ω′) *
                    calc_V2plus(
                        q′,
                        q′′_abso,
                        W_λ,
                        W_λ′,
                        W_λ′′_abso,
                        ifc3_tensor,
                        trip2atomindices,
                        trip2position_j,
                        trip2position_k,
                    ) *
                    delta_plus
                end
                Λ += Λplus
            elseif energy_conversed_minus
                statistics_emit = begin
                    bose(ω′, kbT) + bose(ω′′_emit, kbT) + 1
                end

                Λminus = begin
                    statistics_emit / (ω′′_emit * ω′) *
                    calc_V2minus(
                        q′,
                        q′′_emit,
                        W_λ,
                        W_λ′,
                        W_λ′′_emit,
                        ifc3_tensor,
                        trip2atomindices,
                        trip2position_j,
                        trip2position_k,
                    ) *
                    delta_minus
                end
                Λ += 0.5 * Λminus
            else
                continue
            end
        end
    end
    return Λ
end

"""
  - `ifc3_tensor::Array{Float64,3}`: Ifc3-tensor generated by ShengBTEs thirdorder.py
    indexed by a triplet-index `itrip` and three cartesian indices `α`,`α′` and `α′′`
  - `trip2atomindices::Matrix{Int64}`: Translation of `itrip` to the atomindex of
    both displaced atoms in the calculation of the `ifc3_tensor`.
  - `trip2position_j::Matrix{Float64}`: Unitcell positions of one of the two
    displaced atoms in the calculation of `ifc3_tensor`. Produced by thirdorder.py
    but snapped to the direct lattice grid for exact matching.
  - `trip2position_k::Matrix{Float64}`: Unitcell positions of one of the two
    displaced atoms in the calculation of `ifc3_tensor`. Produced by thirdorder.py
    but snapped to the direct lattice grid for exact matching.
"""
function calc_V2plus(
    q′::Vector{Float64},
    q′′::Vector{Float64},
    W_λ::Matrix{ComplexF64},
    W_λ′::Matrix{ComplexF64},
    W_λ′′::Matrix{ComplexF64},
    ifc3_tensor::Array{Float64,4},
    trip2atomindices::Matrix{Int64},
    trip2position_j::Matrix{Float64},
    trip2position_k::Matrix{Float64},
)::Float64
    V = 0.0 + im * 0.0
    for α in axes(W_λ, 1)
        for α′ in axes(W_λ′, 1)
            for α′′ in axes(W_λ′′, 1)
                for itrip in axes(ifc3_tensor, 1)
                    V += @views begin
                        W_λ[α, trip2atomindices[itrip, 1]] *
                        W_λ′[α′, trip2atomindices[itrip, 2]] *
                        conj(W_λ′′[α′′, trip2atomindices[itrip, 3]]) *
                        ifc3_tensor[itrip, α, α′, α′′] *
                        exp(im * LinAlg.dot(q′, trip2position_j[itrip, :])) *
                        exp(-im * LinAlg.dot(q′′, trip2position_k[itrip, :]))
                    end
                end
            end
        end
    end
    V2 = V * conj(V)

    return V2
end

function calc_V2minus(
    q′::Vector{Float64},
    q′′::Vector{Float64},
    W_λ::Matrix{ComplexF64},
    W_λ′::Matrix{ComplexF64},
    W_λ′′::Matrix{ComplexF64},
    ifc3_tensor::Array{Float64,4},
    trip2atomindices::Matrix{Int64},
    trip2position_j::Matrix{Float64},
    trip2position_k::Matrix{Float64},
)::Float64
    V = 0.0 + im * 0.0
    for α in axes(W_λ, 1)
        for α′ in axes(W_λ′, 1)
            for α′′ in axes(W_λ′′, 1)
                for itrip in axes(ifc3_tensor, 1)
                    V += @views begin
                        W_λ[α, trip2atomindices[itrip, 1]] *
                        conj(W_λ′[α′, trip2atomindices[itrip, 2]]) *
                        conj(W_λ′′[α′′, trip2atomindices[itrip, 3]]) *
                        ifc3_tensor[itrip, α, α′, α′′] *
                        exp(-im * LinAlg.dot(q′, trip2position_j[itrip, :])) *
                        exp(-im * LinAlg.dot(q′′, trip2position_k[itrip, :]))
                    end
                end
            end
        end
    end
    V2 = V * conj(V)
    return V2
end

"""
Snapping positions such that their coordinates are the result of an linear
combination of integer multiples of basisvectors from given lattice vectors. It is
assumed that `positions[i,:]` will yield the `i`-th position.

Lattice Vectors have to be in the form `lattvecs = [ a1 a2 a3 ]`
"""
function snap_to_lattvecs!(lattvecs::Matrix{Float64}, positions::Matrix{Float64})
    # Solve a system of linear equations to get the coefficients that make up the 
    # positions as linear combinations

    # Getting the positions in fractional coordinates as integers
    # Rounding like fortrans `anint()`. Ties are rounded away from zero
    positions_frac =
        round.(\(lattvecs, permutedims(positions)), RoundNearestTiesAway)

    # Overriding the original positions
    positions = permutedims(lattvecs * positions_frac)

    return
end

struct PhaseSpace
    phasespace::Matrix{Float64}

    function PhaseSpace(
        ebdata::ebInputData,
        deconvolution::DeconvData,
        sodata::qeIfc2Output,
        q1_cryst::Matrix{Float64},
        cont_freqs::Vector{Float64};
        brillouin_sampling::Tuple{Int64,Int64,Int64} = (6, 6, 6),
    )
        numatoms = ebdata.allocations["numatoms"]
        # Lattvecs in Angstroem
        lattvecs = ebdata.crystal_info["lattvecs"]
        reclattvecs = calc_reciprocal_lattvecs(lattvecs)
        numbranches = 3 * numatoms

        states = HarmonicStatesData(
            ebdata,
            sodata,
            deconvolution,
            q1_cryst;
            brillouin_sampling,
        )

        numq1 = size(q1_cryst, 1)
        numq2 = size(states.q2_cryst, 1)
        numfreq = size(cont_freqs, 1)

        @info "anharmonic.jl: Calculating phasespace for every frequency..."
        phasespace = Array{Float64}(undef, numfreq, numq1)
        Threads.@threads for ifreq in 1:numfreq
            ω = cont_freqs[ifreq] * turnTHz_to_eV

            for λ in axes(states.q1_evec, 1)
                _, iq = demux1to2(λ, numq1)

                Λ = 0.0
                for λ′ in axes(states.q2_evec, 1)
                    _, iq′ = demux1to2(λ′, numq2)

                    ω′ = states.q2_freqs[λ′] * turnTHz_to_eV
                    v′ = states.q2_velos[λ′, :]

                    for s′′ in 1:numbranches
                        λ′′ = mux2to1(s′′, iq′, numq2)

                        ω′′_abso = states.q3_abso_freqs[λ′′, iq] * turnTHz_to_eV
                        ω′′_emit = states.q3_emit_freqs[λ′′, iq] * turnTHz_to_eV
                        v′′_abso = states.q3_abso_velos[λ′, :, iq]
                        v′′_emit = states.q3_emit_velos[λ′, :, iq]

                        smearing_abso = smearing_type3(
                            v′,
                            v′′_abso,
                            reclattvecs,
                            brillouin_sampling,
                        )

                        smearing_emit = smearing_type3(
                            v′,
                            v′′_emit,
                            reclattvecs,
                            brillouin_sampling,
                        )

                        Λplus = δ(ω, ω′′_abso - ω′, smearing_abso)

                        Λminus = δ(ω, ω′′_emit + ω′, smearing_emit)

                        Λ += Λplus + 0.5 * Λminus
                    end
                end
                phasespace[ifreq, iq] = (Λ / numq2)
            end
        end

        new(phasespace)
    end
end
