include("../src/Phunky.jl")
using .Phunky
import Plots
using LaTeXStrings
Plots.pgfplotsx()

ebdata = ebInputData("examples/input.nml")
numatoms = ebdata.allocations["numatoms"]
numbranches = 3 * numatoms

deconvolution = DeconvData(ebdata)
sodata = qeIfc2Output("examples/espresso.ifc2")
todata = Ifc3Output("examples/force.fc3")

q1_cryst = [0.00 0.00 0.00]

# We give the frequencies in THz
ω_max = 40
ω_min = 0.01
numfreqs = 100
ω_cont = collect(range(ω_min, ω_max, numfreqs))

# Temperature in [K]
temperature = 300.0

sampling = (18, 18, 18)

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

linewidths = linewidth_container.linewidths

Phunky.write_to_file("data/24x24x24-gamma-point-linewidths.data", linewidths)

# # ================================================================================
# # PLOTING
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
#     ω_cont * Phunky.turnTHz_to_eV * 1000,
#     linewidths[:, 1, 4] * 1000;
#     size = plot_size,
#     color = :blue,
#     legend_position = false,
#     grid = false,
#     framestyle = :box,
#     tickdirection = :out,
#     xlabel = L"\mathrm{Frequency}\;\omega\;[\mathrm{meV}]",
#     ylabel = L"\Gamma_{\lambda}(\bar{q}_1,\omega)\; [\mathrm{meV}]",
#     xtickfontsize = 12,
#     ytickfontsize = 12,
#     ylabelfontsize = 12,
#     extra_kwargs = Dict(:subplot => extra_dict),
# )
#
# # ================================================================================
