using Base: DenseUInt8OrInt8
include("../src/Phunky.jl")
using .Phunky

import Plots;
using LaTeXStrings
import LinearAlgebra as LinAlg
using BenchmarkTools;
Plots.pgfplotsx();

# Read-in the main description of the system
@info "Reading input.nml..."
ebdata = ebInputData("examples/input.nml")

# Construct a ifc2-deconvolution based on the system description
@info "Calculating deconvolution..."
deconvolution = DeconvData(ebdata)

# Read-in the second order force-constants
@info "Reading ifc2..."
sodata = qeIfc2Output("examples/espresso.ifc2")

# Main computation of the harmonic phonon properties
@info "Calculating Density of States..."
dos_container = DensityOfStates(
    ebdata,
    sodata,
    deconvolution,
    300,
    (30, 30, 30);
    scalebroad = 1.5,
)

dos = dos_container.density
ω_cont = dos_container.cont_energies

plot_aspectratio = 1 / 2 * (1 + sqrt(5))
plot_height = 300
plot_width = plot_height * plot_aspectratio
plot_size = (plot_width, plot_height)

# Polish for the plot
extra_dict =
    Dict("tick style" => "thick", "xtick pos" => "left", "ytick pos" => "left")

ω_cont .*= 1000

p = Plots.plot(
    ω_cont,
    dos;
    size = plot_size,
    legend_position = false,
    grid = false,
    framestyle = :box,
    tickdirection = :out,
    color = :blue,
    lw = 0.7,
    xlims = (minimum(ω_cont), maximum(ω_cont)),
    ylims = (minimum(dos), maximum(dos) * 1.1),
    xlabel = L"\mathrm{Frequency}\; \hbar\omega\;[\mathrm{meV]}",
    ylabel = L"P^{(2)}(\omega)\;[\mathrm{eV}^{-1}]",
    xtickfontsize = 12,
    ytickfontsize = 12,
    extra_kwargs = Dict(:subplot => extra_dict),
)

# Plots.savefig("graphs/phonon-disp.pdf")
