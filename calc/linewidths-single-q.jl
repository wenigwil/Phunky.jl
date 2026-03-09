include("../src/Phunky.jl")
using .Phunky
import Plots
using LaTeXStrings
#
# ebdata = ebInputData("examples/input.nml")
# numatoms = ebdata.allocations["numatoms"]
# numbranches = 3 * numatoms
#
# deconvolution = DeconvData(ebdata)
# sodata = qeIfc2Output("examples/espresso.ifc2")
# todata = Ifc3Output("examples/force.fc3")
#
q1_cryst = [0.00 0.00 0.00]
#
# # We give the frequencies in THz
ω_max = 35
ω_min = 0.01
numfreqs = 100
ω_cont = collect(range(ω_min, ω_max, numfreqs))
#
# # Temperature in [K]
# temperature = 1e-3
#
# sampling = (18, 18, 18)
#
# linewidth_container = Linewidth(
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
# linewidths = linewidth_container.linewidths
#
# Phunky.write_to_file("data/18x18x18-gamma-point-zeroT-linewidths.data", linewidths)
linewidths = Phunky.read3d_from_file(
    "data/18x18x18-gamma-point-zeroT-linewidths.data",
    Float64,
)

# getting Deinzer's Data in here
deinzer_data_file = readlines("data/deinzer-reference-Gamma-fig2.data")
deinzer_data = Matrix{Float64}(undef, (size(deinzer_data_file, 1), 2))
for i in axes(deinzer_data_file, 1)
    deinzer_data[i, :] = parse.(Float64, split(deinzer_data_file[i], ","))
end
deinzer_ω_cont = deinzer_data[:, 1]
deinzer_linewidths = deinzer_data[:, 2]

# pulling 0-values out of my data
plt_linewidths = linewidths[:, 1, 4]
plt_ω_cont = ω_cont
plt_ω_cont = plt_ω_cont[broadcast(~, iszero.(plt_linewidths))]
plt_linewidths = plt_linewidths[broadcast(~, iszero.(plt_linewidths))]

# # ================================================================================
# # PLOTING
Plots.pgfplotsx()
plot_aspectratio = 1 / 2 * (1 + sqrt(5))
plot_height = 300
plot_width = plot_height * plot_aspectratio
plot_size = (plot_width, plot_height)

# Polish for the plot
extra_dict =
    Dict("tick style" => "thick", "xtick pos" => "left", "ytick pos" => "left")

p = Plots.plot(
    plt_ω_cont * Phunky.turnTHz_to_eV * 1000,
    plt_linewidths * 84297.93663 / (4 * pi);
    size = plot_size,
    color = :blue,
    legend_position = false,
    grid = false,
    yscale = :log,
    xlims = (20, 125),
    ylims = (10^(-3), 1),
    framestyle = :box,
    tickdirection = :out,
    xlabel = L"\mathrm{Frequency}\;\hbar\omega\;[\mathrm{meV}]",
    ylabel = L"\Gamma_{\lambda}(\bar{q}_1,\omega) [\mathrm{THz}]",
    xtickfontsize = 12,
    ytickfontsize = 12,
    ylabelfontsize = 12,
    extra_kwargs = Dict(:subplot => extra_dict),
)
Plots.plot!(deinzer_ω_cont * 1000, deinzer_linewidths; color = :red, yscale = :log)

# # ================================================================================
