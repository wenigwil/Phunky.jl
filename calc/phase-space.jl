include("../src/Phunky.jl")
using .Phunky
using Plots

ebdata = ebInputData("examples/input.nml")
deconvolution = DeconvData(ebdata)
sodata = qeIfc2Output("examples/espresso.ifc2")

# We will read the qpoints from the 12x12x12 reference from elphbolt
q1_cryst_file = readlines("examples/reference-12x12x12/ph.wavevecs_ibz")
numq1 = size(q1_cryst_file, 1)
q1_cryst = Matrix{Float64}(undef, (numq1, 3))
for iq in axes(q1_cryst, 1)
    q1_cryst[iq, :] = parse.(Float64, split(q1_cryst_file[iq]))
end

# ω_max = 16 * 1.2
# ω_min = 0.01
# numfreqs = 100
# ω_cont = collect(range(ω_min, ω_max, numfreqs))

ω_cont = [12.0]

sampling = (3, 3, 3)

phasespace_container = PhaseSpace(
    ebdata,
    deconvolution,
    sodata,
    q1_cryst,
    ω_cont;
    brillouin_sampling = sampling,
)

phasespace = phasespace_container.phasespace
