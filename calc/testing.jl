using Phunky
using Plots

smearing = 0.03
sampling = (6, 6, 6)

q1_cryst, energies, plus, minus, total = Phunky.test_phasespace(smearing, sampling)
_, energies_eb, plus_eb, minus_eb, total_eb = Phunky.get_elphbolt_phasespace()
