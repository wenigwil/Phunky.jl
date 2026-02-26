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

q1_cryst = [0.0 0.0 0.25]

# We give the frequencies in THz
ω_max = 35
ω_min = 0.01
numfreqs = 300
ω_cont = collect(range(ω_min, ω_max, numfreqs))

sampling = (24, 24, 24)

phasespace_container = PhaseSpace(
    ebdata,
    deconvolution,
    sodata,
    q1_cryst,
    ω_cont;
    brillouin_sampling = sampling,
)

phasespace = phasespace_container.phasespace

plot_aspectratio = 1 / 2 * (1 + sqrt(5))
plot_height = 300
plot_width = plot_height * plot_aspectratio
plot_size = (plot_width, plot_height)

# Polish for the plot
extra_dict =
    Dict("tick style" => "thick", "xtick pos" => "left", "ytick pos" => "left")

p = Plots.plot(
    ω_cont * Phunky.turnTHz_to_eV * 1000,
    phasespace[:, 1];
    size = plot_size,
    color = :blue,
    legend_position = false,
    grid = false,
    framestyle = :box,
    tickdirection = :out,
    xlabel = L"\mathrm{Frequency}\;\omega\;[\mathrm{meV}]",
    ylabel = L"P^{3}(\bar{q}_1,\omega) \;\; [\mathrm{eV}^{-1}]",
    xtickfontsize = 12,
    ytickfontsize = 12,
    ylabelfontsize = 12,
    extra_kwargs = Dict(:subplot => extra_dict),
)
