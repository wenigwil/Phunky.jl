include("../src/Phunky.jl")
using .Phunky
using Plots

# Read the system description
ebdata = ebInputData("examples/input.nml")

# Construct a deconvolution from the system description
deconvolution = DeconvData(ebdata)

# Read a quantum espresso ifc2 file
sodata = qeIfc2Output("examples/espresso.ifc2")

# points = Phunky.sample_cube((3, 3, 3))
points = [0.0, 0.0]

vibr = LatticeVibrations(ebdata, sodata, deconvolution, points)

freqs = vibr.fullq_freqs
velocities = vibr.velocities

for j in axes(velocities, 3)
    println(velocities[end, :, j])
    println("")
end
