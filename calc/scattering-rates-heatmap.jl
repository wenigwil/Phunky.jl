include("../src/Phunky.jl")
using .Phunky
import Plots
using LaTeXStrings

# # ================================================================================ # CALCULATION
#
# ebdata = ebInputData("examples/input.nml")
# numatoms = ebdata.allocations["numatoms"]
# numbranches = 3 * numatoms
#
# deconvolution = DeconvData(ebdata)
# sodata = qeIfc2Output("examples/espresso.ifc2")
# todata = Ifc3Output("examples/force.fc3")
#
# # Create path of q-points from Seekpath cF1 high symmetry points
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
seek_path_1 = seek_path_points[[2, 1, 3], :]
sympath =
    Sympath(seek_path_points, point_labels, seek_path_1, numpoints_per_section = 25)
q1_cryst = sympath.qpoints
distances = sympath.distances
#
#
# # We give the frequencies in THz
# ω_max = 40
# ω_min = 0.001
# numfreqs = 100
# ω_cont = collect(range(ω_min, ω_max, numfreqs))
#
# # Temperature in [K]
# temperature = 300.0
#
# sampling = (18, 18, 18)
#
# scattratt_container = Phonons(
#     ebdata,
#     deconvolution,
#     sodata,
#     todata,
#     q1_cryst,
#     ω_cont,
#     temperature;
#     brillouin_sampling = sampling,
# )
#
# scattering_rates = scattratt_container.scattering_rate
#
# # Write important things for plotting to a file
# Phunky.write_to_file("data/Si_18x18x18.scattrats", scattering_rates)
# Phunky.write_to_file("data/Si_18x18x18.distances", distances)
# Phunky.write_to_file("data/Si_18x18x18.w_cont", ω_cont)
# # ================================================================================

# # ================================================================================
# READING FROM WRITTEN
ω_cont = Phunky.read1d_from_file("distributed/data/w-cont-18x18x18.data", Float64)
distances =
    Phunky.read1d_from_file("distributed/data/distance-18x18x18.data", Float64)

scattering_rates =
    Phunky.read3d_from_file("distributed/data/scattrats-18x18x18.data", Float64)

# # ================================================================================

# # ================================================================================
# # PLOTTING

# Plots.pgfplotsx()
plot_aspectratio = 1 / 2 * (1 + sqrt(5))
plot_height = 300
plot_width = plot_height * plot_aspectratio
plot_size = (plot_width, plot_height)
extra_dict =
    Dict("tick style" => "thick", "xtick pos" => "left", "ytick pos" => "left")

p = Plots.heatmap(
    distances,
    ω_cont * Phunky.turnTHz_to_eV * 1000,
    scattering_rates[:, :, 4] * 1e6;
    size = plot_size,
    colormap = :lajolla,
    clims = (0.0, 10),
    ylims = (0.0, 85),
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
    # colorbar_title = L"P^{3}(\bar{q},\omega) \;\; [\mathrm{eV}^{-1}]",
    extra_kwargs = Dict(:subplot => extra_dict),
)

# # =================================================================================
