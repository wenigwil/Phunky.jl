include("../src/Phunky.jl")
using .Phunky
import Plots
using LaTeXStrings
# Plots.pgfplotsx()

ebdata = ebInputData("examples/input.nml")
numatoms = ebdata.allocations["numatoms"]
numbranches = 3 * numatoms

deconvolution = DeconvData(ebdata)
sodata = qeIfc2Output("examples/espresso.ifc2")
todata = Ifc3Output("examples/force.fc3")

q1_file = readlines("examples/12x12x12/ph.wavevecs_ibz")
numq1 = size(q1_file, 1)
q1_cryst = Matrix{Float64}(undef, (numq1, 3))
for iq in axes(q1_cryst, 1)
    q1_cryst[iq, :] = parse.(Float64, split(q1_file[iq]))
end

# Temperature in [K]
temperature = 300.0

sampling = (12, 12, 12)

scattratt_container = Phonons(
    ebdata,
    deconvolution,
    sodata,
    todata,
    q1_cryst,
    temperature;
    brillouin_sampling = sampling,
)

scattering_rates = scattratt_container.scattering_rate
# reshape for plotting
scattering_rates_reshaped = reshape(scattering_rates, (numq1 * numbranches))

# For plotting the scattering rates against energies we also need the eigenenergies 
# for every q-point
lattvibr = LatticeVibrations(ebdata, sodata, deconvolution, q1_cryst)
# get them in meV
freqs = lattvibr.fullq_freqs * Phunky.turnTHz_to_eV * 1000
# reshape for plotting
freqs = reshape(freqs, (numq1 * numbranches))

# take out all zeros and the corresponding energies as well
freqs = freqs[broadcast(~, iszero.(scattering_rates_reshaped))]
scattering_rates_reshaped =
    scattering_rates_reshaped[broadcast(~, iszero.(scattering_rates_reshaped))]

plot_aspectratio = 1 / 2 * (1 + sqrt(5))
plot_height = 300
plot_width = plot_height * plot_aspectratio
plot_size = (plot_width, plot_height)

# Polish for the plot
# extra_dict =
#     Dict("tick style" => "thick", "xtick pos" => "left", "ytick pos" => "left")

Plots.plot(
    freqs,
    scattering_rates_reshaped;
    size = plot_size,
    color = :blue,
    yscale = :log10,
    legend_position = false,
    grid = false,
    linealpha = 0.0,
    markershape = :circle,
    markerstrokewidth = 0.0,
    markerstrokealpha = 0.0,
    markersize = 2,
    framestyle = :box,
    tickdirection = :out,
    xlabel = L"\mathrm{Energy}\;\hbar\omega\;[\mathrm{meV}]",
    ylabel = L"\Gamma_{\!\lambda}(\omega_\lambda,T)\;\; [\mathrm{eV}^{-1}]",
    xtickfontsize = 12,
    ytickfontsize = 12,
    ylabelfontsize = 12,
    # extra_kwargs = Dict(:subplot => extra_dict),
)

# COMPARISON DATA
# We will read the scattering rates for comparison and the energies
energies_file = readlines("examples/12x12x12/ph.ens_ibz")
energies_comp = Matrix{Float64}(undef, (numq1, numbranches))
for iq in axes(q1_cryst, 1)
    energies_comp[iq, :] = parse.(Float64, split(energies_file[iq]))
end
energies_comp = reshape(energies_comp, (numq1 * numbranches)) .* 1000

# We will read the scattering rates for comparison and the energies
scattratt_file = readlines("examples/12x12x12/ph.W_rta_3ph");
scattratt_comparison = Matrix{Float64}(undef, (numq1, numbranches));
for iq in axes(q1_cryst, 1)
    scattratt_comparison[iq, :] = parse.(Float64, split(scattratt_file[iq]))
end
scattratt_comparison = reshape(scattratt_comparison, (numq1 * numbranches));

energies_comp = energies_comp[broadcast(~, iszero.(scattratt_comparison))];
scattratt_comparison =
    scattratt_comparison[broadcast(~, iszero.(scattratt_comparison))];

# identical plot with comp data
Plots.plot(
    energies_comp,
    scattratt_comparison;
    size = plot_size,
    color = :blue,
    yscale = :log10,
    legend_position = false,
    grid = false,
    linealpha = 0.0,
    markershape = :circle,
    markerstrokewidth = 0.0,
    markerstrokealpha = 0.0,
    markersize = 2,
    framestyle = :box,
    tickdirection = :out,
    xlabel = L"\mathrm{Energy}\;\hbar\omega\;[\mathrm{meV}]",
    ylabel = L"\Gamma_{\!\lambda}(\omega_\lambda,T)\;\; [\mathrm{eV}^{-1}]",
    xtickfontsize = 12,
    ytickfontsize = 12,
    ylabelfontsize = 12,
    # extra_kwargs = Dict(:subplot => extra_dict),
)
