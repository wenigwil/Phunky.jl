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

# Seekpath cF1 high symmetry points
point_labels = ["Γ", "K", "L", "U", "W", "X", raw"$\mathrm{W}_2$"]
seek_path_points = [
    0.0 0.0 0.0
    0.375 0.375 0.75
    0.5 0.5 0.5
    0.625 0.25 0.625
    0.5 0.25 0.75
    0.5 0.0 0.5
    0.75 0.25 0.5
]

# Build a "walkable" path from the seekpath symmetry points
seek_path_1 = seek_path_points[[1, 6, 4], :]
seek_path_2 = seek_path_points[[2, 1, 3, 5, 6, 7], :]

# Construct a qpoint list from a path and supply it with things to make it plottable
@info "Building path..."
sympath = Sympath(
    seek_path_points,
    point_labels,
    seek_path_1,
    seek_path_2;
    numpoints_per_section = 40,
)
print("Generated a path of length", size(sympath.qpoints), "\n")

q1_cryst = sympath.qpoints
distances = sympath.distances

# We give the frequencies in THz
ω_max = 35
ω_min = 0.01
numfreqs = 150
ω_cont = collect(range(ω_min, ω_max, numfreqs))

sampling = (12, 12, 12)

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

p = Plots.heatmap(
    distances,
    ω_cont * Phunky.turnTHz_to_eV * 1000,
    phasespace;
    size = plot_size,
    colormap = :lajolla,
    clims = (0.0, 600.0),
    legend_position = true,
    grid = false,
    framestyle = :box,
    tickdirection = :out,
    xticks = (sympath.xticks_pos, sympath.xticks_labels),
    ylabel = L"\mathrm{Frequency}\;\omega\;[\mathrm{meV}]",
    tex_output_standalone = false,
    xtickfontsize = 12,
    ytickfontsize = 12,
    colorbar_tickfontsize = 12,
    colorbar_titlefontsize = 12,
    ylabelfontsize = 12,
    colorbar_title = L"P^{3}(\bar{q},\omega) \;\; [\mathrm{eV}^{-1}]",
    extra_kwargs = Dict(:subplot => extra_dict),
)
