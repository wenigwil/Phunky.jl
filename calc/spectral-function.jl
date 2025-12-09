using Phunky
import Plots as plt
# Read the system description
ebdata = ebInputData("examples/input.nml")

# Construct a deconvolution from the system description
deconvolution = DeconvData(ebdata)

# Read a quantum espresso ifc2 file
sodata = qeIfc2Output("examples/espresso.ifc2")
# Build a symmetry path along which the phonon lifetime will be computed
# Seekpath cF1 high symmetry points
# point_labels = ["Γ", "K", "L", "U", "W", "X", raw"$\mathrm{W}_2$"]
# seek_path_points = [
#     0.0 0.0 0.0
#     0.375 0.375 0.75
#     0.5 0.5 0.5
#     0.625 0.25 0.625
#     0.5 0.25 0.75
#     0.5 0.0 0.5
#     0.75 0.25 0.5
# ]
# # Walk the walk. Only small for now
# route = seek_path_points[[2, 1], :]
# sympath = Sympath(seek_path_points, point_labels, route; numpoints_per_section = 3)

# Read a phonopy ifc3 file
todata = Ifc3Output("examples/force.fc3")

starting = 0
ending = 20
numpoints = 2000
smearing = (ending - starting) / (16 * numpoints)

cont_freqs = collect(range(starting, ending, numpoints))
T = 300.0

# samplings = [4, 6, 8, 12, 14]
# phonon_container = Vector{Phunky.ThreePhononDensity}(undef, length(samplings))

# for i in axes(samplings, 1)
phonon_container = Phunky.ThreePhononDensity(
    ebdata,
    deconvolution,
    sodata,
    [0.0 0.0 0.0],
    cont_freqs,
    smearing;
    # brillouin_sampling = (samplings[i], samplings[i], samplings[i]),
    brillouin_sampling = (14, 14, 14),
)

density = phonon_container.density[:, 1, :]
# end

# for i in axes(samplings, 1)
#     phonon_container[i] = Phonons(
#         ebdata,
#         deconvolution,
#         sodata,
#         todata,
#         [0.0 0.0 0.0],
#         cont_freqs,
#         T,
#         smearing;
#         brillouin_sampling = (samplings[i], samplings[i], samplings[i]),
#     )
# end
