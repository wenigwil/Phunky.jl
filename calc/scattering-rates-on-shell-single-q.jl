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

# plot_aspectratio = 1 / 2 * (1 + sqrt(5))
# plot_height = 300
# plot_width = plot_height * plot_aspectratio
# plot_size = (plot_width, plot_height)

# Polish for the plot
# extra_dict =
#     Dict("tick style" => "thick", "xtick pos" => "left", "ytick pos" => "left")

# p = Plots.plot(
#     ω_cont * Phunky.turnTHz_to_eV * 1000,
#     scattering_rates[:, 1, 6];
#     size = plot_size,
#     color = :blue,
#     legend_position = false,
#     grid = false,
#     framestyle = :box,
#     tickdirection = :out,
#     xlabel = L"\mathrm{Frequency}\;\omega\;[\mathrm{meV}]",
#     ylabel = L"P^{3}(\bar{q}_1,\omega) \;\; [\mathrm{eV}^{-1}]",
#     xtickfontsize = 12,
#     ytickfontsize = 12,
#     ylabelfontsize = 12,
#     # extra_kwargs = Dict(:subplot => extra_dict),
# )
