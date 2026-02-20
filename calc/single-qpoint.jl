include("../src/Phunky.jl")
using .Phunky

ebinput = ebInputData("examples/input.nml")
deconvolution = DeconvData(ebinput)
qeinput = qeIfc2Output("examples/espresso.ifc2")

qpoints_cryst = Phunky.sample_cube((2, 2, 2))

lattvibr = LatticeVibrations(ebinput, qeinput, deconvolution, qpoints_cryst)

freqs = lattvibr.fullq_freqs
velocities = lattvibr.velocities
qpoints_cart = lattvibr.qpoints_cart
