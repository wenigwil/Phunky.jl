include("../src/Phunky.jl")
using .Phunky
import Plots
using LaTeXStrings

ebdata = ebInputData("examples/input.nml")
numatoms = ebdata.allocations["numatoms"]
numbranches = 3 * numatoms

deconvolution = DeconvData(ebdata)
sodata = qeIfc2Output("examples/espresso.ifc2")
todata = Ifc3Output("examples/force.fc3")

q1_cryst = [0.5 0.5 0.5]
#
# # We give the frequencies in THz
ω_max = 40
ω_min = 0.01
numfreqs = 200
ω_cont = collect(range(ω_min, ω_max, numfreqs))

# Temperature in [K]
temperature = 300.0

sampling = (21, 21, 21)

linewidth_container = Linewidth(
    ebdata,
    deconvolution,
    sodata,
    todata,
    q1_cryst,
    ω_cont,
    temperature;
    brillouin_sampling = sampling,
)

eigenenergies = linewidth_container.eigenenergies

lineshift_container = Lineshift(linewidth_container)

spectralfunction_container =
    SpectralFunction(linewidth_container, lineshift_container)

phunky = spectralfunction_container.spectral_function
ω_cont = spectralfunction_container.ω_cont_extended

Phunky.write_to_file(
    "data/specfunc/linewidths-extended-21x21x21-300K-L.data",
    lineshift_container.linewidths_extended,
)
Phunky.write_to_file(
    "data/specfunc/lineshifts-21x21x21-300K-L.data",
    lineshift_container.lineshifts,
)
Phunky.write_to_file(
    "data/specfunc/contfreqs-extended-21x21x21-300K-L.data",
    lineshift_container.ω_cont_extended,
)
# linewidths = Phunky.read3d_from_file(
#     "data/18x18x18-gamma-point-zeroT-linewidths.data",
#     Float64,
# )

# # ================================================================================
# # PLOTING
# Plots.pgfplotsx()
# plot_aspectratio = 1 / 2 * (1 + sqrt(5))
# plot_height = 300
# plot_width = plot_height * plot_aspectratio
# plot_size = (plot_width, plot_height)
#
# # Polish for the plot
# extra_dict =
#     Dict("tick style" => "thick", "xtick pos" => "left", "ytick pos" => "left")
#
# p = Plots.plot(
#     ω_cont,
#     phunky[:, 1, 4];
#     size = plot_size,
#     color = :blue,
#     legend_position = false,
#     grid = false,
#     xlims = (0.0, maximum(ω_cont)),
#     ylims = (-0.05, maximum(phunky) * 1.05),
#     framestyle = :box,
#     tickdirection = :out,
#     xlabel = L"\mathrm{Frequency}\;\hbar\omega\;[\mathrm{meV}]",
#     ylabel = L"A_{\lambda}(\omega)\;[\mathrm{THz}]",
#     xtickfontsize = 12,
#     ytickfontsize = 12,
#     ylabelfontsize = 12,
#     extra_kwargs = Dict(:subplot => extra_dict),
# )
# Plots.vline!([eigenenergies[1, 4]]; color = :red)

# # ================================================================================
