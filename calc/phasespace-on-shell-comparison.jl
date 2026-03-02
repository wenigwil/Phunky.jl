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

q1_file = readlines("examples/12x12x12/ph.wavevecs_ibz")
numq1 = size(q1_file, 1)
q1_cryst = Matrix{Float64}(undef, (numq1, 3))
for iq in axes(q1_cryst, 1)
    q1_cryst[iq, :] = parse.(Float64, split(q1_file[iq]))
end

sampling = (12, 12, 12)

# For plotting the phase space against energies we also need the eigenenergies for 
# every q-point
lattvibr = LatticeVibrations(ebdata, sodata, deconvolution, q1_cryst)
# get them in meV
freqs = lattvibr.fullq_freqs * Phunky.turnTHz_to_eV * 1000
# reshape for plotting
freqs = reshape(freqs, (numq1 * numbranches))

phasespace_container = PhaseSpace(
    ebdata,
    deconvolution,
    sodata,
    q1_cryst;
    brillouin_sampling = sampling,
)

phasespace = phasespace_container.phasespace
phasespace = reshape(phasespace, (numq1 * numbranches))

plot_aspectratio = 1 / 2 * (1 + sqrt(5))
plot_height = 300
plot_width = plot_height * plot_aspectratio
plot_size = (plot_width, plot_height)

# Polish for the plot
# extra_dict =
#     Dict("tick style" => "thick", "xtick pos" => "left", "ytick pos" => "left")

p = Plots.plot(
    freqs,
    phasespace;
    size = plot_size,
    color = :blue,
    linealpha = 0.0,
    markershape = :circle,
    markerstrokewidth = 0.0,
    markerstrokealpha = 0.0,
    markersize = 1,
    legend_position = false,
    grid = false,
    framestyle = :box,
    tickdirection = :out,
    xlabel = L"\mathrm{Frequency}\;\omega\;[\mathrm{meV}]",
    ylabel = L"P^{3}(\bar{q}_1,\omega) \;\; [\mathrm{eV}^{-1}]",
    title = "12x12x12 Phunky.jl " * L"P^{(3)}",
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
phsp_minus_file = readlines("examples/12x12x12/ph.phase_space3_minus");
phsp_plus_file = readlines("examples/12x12x12/ph.phase_space3_plus");
phsp_minus_comp = Matrix{Float64}(undef, (numq1, numbranches));
phsp_plus_comp = Matrix{Float64}(undef, (numq1, numbranches));
for iq in axes(q1_cryst, 1)
    phsp_minus_comp[iq, :] = parse.(Float64, split(phsp_minus_file[iq]))
    phsp_plus_comp[iq, :] = parse.(Float64, split(phsp_plus_file[iq]))
end
phsp_minus_comp = reshape(phsp_minus_comp, (numq1 * numbranches));
phsp_plus_comp = reshape(phsp_plus_comp, (numq1 * numbranches));

phsp_comp = phsp_plus_comp + 0.5 * phsp_minus_comp

# identical plot with comp data
Plots.plot(
    energies_comp,
    phsp_comp;
    size = plot_size,
    color = :blue,
    legend_position = false,
    grid = false,
    linealpha = 0.0,
    # ylims = (0.0, 650.0),
    markershape = :circle,
    markerstrokewidth = 0.0,
    markerstrokealpha = 0.0,
    markersize = 1,
    framestyle = :box,
    tickdirection = :out,
    xlabel = L"\mathrm{Frequency}\;\hbar\omega\;[\mathrm{meV}]",
    ylabel = L"P^{3}(\bar{q}_1,\omega) \;\; [\mathrm{eV}^{-1}]",
    title = "12x12x12 elphbolt reference " * L"P^{(3)}",
    xtickfontsize = 12,
    ytickfontsize = 12,
    ylabelfontsize = 12,
    # extra_kwargs = Dict(:subplot => extra_dict),
)
