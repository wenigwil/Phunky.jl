struct SpectralFunction
    spectral_function::Array{Float64,3}

    function SpectralFunction(
        lineshift_container::Lineshift,
        eigenenergies::Matrix{Float64},
        ω_cont::Vector{Float64},
    )
        lineshifts = lineshift_container.lineshifts .* 10
        linewidths = lineshift_container.linewidth .* 10
        spectral_function = Array{Float64,3}(undef, size(linewidths))

        for s in axes(spectral_function, 3)
            for iq in axes(spectral_function, 2)
                for ifreq in axes(spectral_function, 1)
                    ω_λ = eigenenergies[iq, s]
                    Γ_λ = linewidths[ifreq, iq, s]
                    Δ_λ = lineshifts[ifreq, iq, s]
                    ω = ω_cont[ifreq]

                    spectral_function[ifreq, iq, s] = begin
                        1 / pi * (4 * ω_λ^2 * Γ_λ) /
                        ((ω^2 - ω_λ^2 - 2 * ω_λ * Δ_λ)^2 + (2 * ω_λ * Γ_λ)^2)
                    end
                end
            end
        end

        new(spectral_function)
    end
end
