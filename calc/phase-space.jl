include("../src/Phunky.jl")
using .Phunky
using Plots

ebdata = ebInputData("examples/input.nml")
numatoms = ebdata.allocations["numatoms"]
numbranches = 3 * numatoms

deconvolution = DeconvData(ebdata)
sodata = qeIfc2Output("examples/espresso.ifc2")

# We will read the qpoints from the 12x12x12 reference from elphbolt
# q1_cryst_file = readlines("examples/reference-12x12x12/ph.wavevecs_ibz")
# numq1 = size(q1_cryst_file, 1)
# q1_cryst = Matrix{Float64}(undef, (numq1, 3))
# for iq in axes(q1_cryst, 1)
#     q1_cryst[iq, :] = parse.(Float64, split(q1_cryst_file[iq]))
# end

q1_cryst = [0.0 0.0 0.25]

# We give the frequencies in THz
ω_max = 35
ω_min = 0.01
numfreqs = 400
ω_cont = collect(range(ω_min, ω_max, numfreqs))

sampling = (24, 24, 24)

phasespace_container = PhaseSpace(
    ebdata,
    deconvolution,
    sodata,
    q1_cryst,
    ω_cont;
    brillouin_sampling = sampling,
)

phasespace = phasespace_container.phasespace

# read the frequencies from reference
# freqs_file = readlines("examples/reference-12x12x12/ph.ens_ibz")
# freqs_ref = Matrix{Float64}(undef, (numq1, numbranches))
# for iq in axes(freqs_ref, 1)
#     freqs_ref[iq, :] = parse.(Float64, split(freqs_file[iq]))
# end
#
# # reading the phasespace
# phasespace_file = readlines("examples/reference-12x12x12/ph.phase_space3_total")
# phasespace_ref = Matrix{Float64}(undef, (numq1, numbranches))
# for iq in axes(phasespace_ref, 1)
#     phasespace_ref[iq, :] = parse.(Float64, split(phasespace_file[iq]))
# end
# phasespace_ref = reshape(phasespace_ref, numq1 * numbranches)
