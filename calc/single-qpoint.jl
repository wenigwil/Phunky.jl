include("../src/Phunky.jl")
using .Phunky

ebdata = ebInputData("examples/input.nml")
deconvolution = DeconvData(ebdata)
sodata = qeIfc2Output("examples/espresso.ifc2")

# qpoints_cryst = Phunky.sample_cube((2, 2, 2))
qpoints_cryst = [0.5 0.5 0.5]

# lattvibr = LatticeVibrations(ebinput, qeinput, deconvolution, qpoints_cryst)

states = HarmonicStatesData(
    ebdata,
    sodata,
    deconvolution,
    qpoints_cryst;
    brillouin_sampling = (2, 2, 2),
)

