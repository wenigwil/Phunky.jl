struct SpectralFunction
    spectral_function::Array{Float64,3}
    ω_cont_extended::Vector{Float64}

    function SpectralFunction(
        linewidths::Array{Float64,3},
        lineshifts::Array{Float64,3},
        ω_cont_extended::Vector{Float64},
        eigenenergies::Matrix{Float64},
    )
        spectral_function = Array{Float64,3}(undef, size(linewidths))

        for s in axes(spectral_function, 3)
            for iq in axes(spectral_function, 2)
                for ifreq in axes(spectral_function, 1)
                    ω_λ = eigenenergies[iq, s]
                    Γ_λ = linewidths[ifreq, iq, s]
                    Δ_λ = lineshifts[ifreq, iq, s]
                    ω = ω_cont_extended[ifreq]

                    spectral_function[ifreq, iq, s] = begin
                        1 / pi * (4 * ω_λ^2 * Γ_λ) /
                        ((ω^2 - ω_λ^2 - 2 * ω_λ * Δ_λ)^2 + (2 * ω_λ * Γ_λ)^2)
                    end
                end
            end
        end

        new(spectral_function, ω_cont_extended)
    end
end

function SpectralFunction(
    linewidth_container::Linewidth,
    lineshift_container::Lineshift,
)
    eigenenergies = linewidth_container.eigenenergies

    lineshifts = lineshift_container.lineshifts
    linewidths = lineshift_container.linewidths_extended
    ω_cont_extended = lineshift_container.ω_cont_extended

    return SpectralFunction(linewidths, lineshifts, ω_cont_extended, eigenenergies)
end
