# The Hilbert Transform used in this part of Phunky.jl was adapted from
# Dr. Jae-Mo Lihm's codebase EPSpectral.jl

struct Lineshift
    lineshifts::Array{Float64,3}
    linewidth::Array{Float64,3}
    ω_cont::Vector{Float64}

    function Lineshift(ω_cont::Vector{Float64}, linewidths::Array{Float64,3})
        ω_cont_extended = extend_to_odd(ω_cont)
        shape = size(linewidths)
        lineshifts = Array{Float64,3}(undef, (2 * shape[1], shape[2], shape[3]))
        linewidths_extended_full = Array{Float64,3}(undef, size(lineshifts))
        for s in axes(lineshifts, 3)
            for iq in axes(lineshifts, 2)
                linewidths_extended = extend_to_odd(linewidths[:, iq, s])

                # The lineshifts will have the same shape as the linewidths
                lineshifts[:, iq, s] =
                    kramers_kronig(ω_cont_extended, linewidths_extended)

                linewidths_extended_full[:, iq, s] = linewidths_extended
            end
        end

        new(lineshifts, linewidths_extended_full, ω_cont_extended)
    end
end
function extend_to_odd(x::Vector{T}) where {T<:Number}
    x_reverse = (-1) .* reverse(x)
    x_extended = vcat(x_reverse, x)
    return x_extended
end
function kramers_kronig(
    ω_cont::Vector{Float64},
    lw::Vector{Float64};
    tail::Bool = false,
)
    ls = zeros(length(ω_cont))

    # Linearly fit y in (ωs[j] - dω/2, ωs[j] + dω/2) and integrate.
    for j in axes(ω_cont, 1)
        # Linear fit z = aω + b
        # Compute slope using the left and right points if available.
        if j == 1
            a = (lw[j + 1] - lw[j]) / (ω_cont[j + 1] - ω_cont[j])
        elseif j == length(ω_cont)
            a = (lw[j] - lw[j - 1]) / (ω_cont[j] - ω_cont[j - 1])
        else
            a = (lw[j + 1] - lw[j - 1]) / (ω_cont[j + 1] - ω_cont[j - 1])
        end
        b = lw[j] - a * ω_cont[j]

        if j == 1
            ωL = ω_cont[j] - (ω_cont[j + 1] - ω_cont[j]) / 2
        else
            ωL = (ω_cont[j] + ω_cont[j - 1]) / 2
        end

        if j == length(ω_cont)
            ωR = ω_cont[j] + (ω_cont[j] - ω_cont[j - 1]) / 2
        else
            ωR = (ω_cont[j] + ω_cont[j + 1]) / 2
        end

        for i in eachindex(ω_cont)
            ls[i] += a * (ωR - ωL)
            ls[i] +=
                (a * ω_cont[i] + b) * log(abs((ωR - ω_cont[i]) / (ωL - ω_cont[i])))
        end
    end
    ls ./= π

    if tail
        dω = ω_cont[2] - ω_cont[1]
        for (i, ω) in enumerate(ω_cont)
            # 1 / √ω extrapolation
            ωL = ω_cont[1] - dω / 2
            ωR = ω_cont[end] + dω / 2
            yL = sqrt(-ω_cont[1]) * lw[1]
            yR = sqrt(ω_cont[end]) * lw[end]

            if abs(ω) < sqrt(eps(ω))
                ls[i] += -yL * 2 / sqrt(-ωL) / π
                ls[i] += yR * 2 / sqrt(ωR) / π
            elseif ω > 0
                ls[i] += yL * (2 * atan(sqrt(-ωL / ω)) - π) / sqrt(ω) / π
                ls[i] += yR * 2 * atanh(sqrt(ω / ωR)) / sqrt(ω) / π

            else  # ω < 0
                ls[i] += -yL * 2 * atanh(sqrt(ω / ωL)) / sqrt(-ω) / π
                ls[i] += yR * 2 * atan(sqrt(-ω / ωR)) / sqrt(-ω) / π
            end
        end
    end

    return ls
end
