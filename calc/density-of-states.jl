include("../src/Phunky.jl")
using .Phunky
using Plots
using NumericalIntegration

# Read the system description
ebdata = ebInputData("examples/input.nml")

# Construct a deconvolution from the system description
deconvolution = DeconvData(ebdata)

# Read a quantum espresso ifc2 file
sodata = qeIfc2Output("examples/espresso.ifc2")

dense = DensityOfStates(ebdata, sodata, deconvolution, 100, (12, 12, 12))

density = dense.density;
cont_energies = dense.cont_energies;

# Should be 6 for now
println(integrate(cont_energies, density))
