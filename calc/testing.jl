using Phunky
using Plots
# Plots.pgfplotsx();

extra_dict =
    Dict("tick style" => "thick", "xtick pos" => "left", "ytick pos" => "left")

plotargs = (
    yaxis = :log,
    markersize = 2,
    markerstrokealpha = 0.0,
    markerstrokewidth = 0.0,
    markercolor = :blue,
    legend_position = false,
    grid = false,
    framestyle = :box,
    tickdirection = :out,
    xtickfontsize = 12,
    ytickfontsize = 12,
    tex_output_standalone = true,
    extra_kwargs = Dict(:subplot => extra_dict),
)

smearing = 0.8
sampling = (6, 6, 6)

q1_cryst, energies, plus, minus, total = Phunky.test_phasespace(smearing, sampling)
_, energies_eb, plus_eb, minus_eb, total_eb = Phunky.get_elphbolt_phasespace()

new_length = prod(size(energies_eb))

energies_eb = reshape(energies_eb, (new_length))
total_eb = reshape(total_eb, (new_length))

energies = reshape(energies, (new_length))
total = reshape(total, (new_length))

my_plot = scatter(
    energies .* 1000,
    total;
    xlabel = "Energy E ",
    ylabel = "Phasespace",
    title = "my 1-phonon phasespace",
    plotargs...,
)

eb_plot = scatter(
    energies_eb .* 1000,
    total_eb;
    xlabel = "Energy E [meV]",
    ylabel = "Phasespace [1/eV]",
    title = "elphbolt 1-phonon phasespace",
    plotargs...,
)
savefig(eb_plot, "eb-phasespace.pdf")
savefig(my_plot, "my-phasespace.pdf")
