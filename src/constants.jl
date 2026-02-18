# ========================
# Constants from reference
# ========================

# Planksches Wirkungsquantum [J*s] https://physics.nist.gov/cgi-bin/cuu/Value?h
const h_Js = 6.62607015e-34

# Elementary charge [C] https://physics.nist.gov/cgi-bin/cuu/Value?e
const e_C = 1.602176634e-19

# Vacuum light speed [m/s] https://physics.nist.gov/cgi-bin/cuu/Value?c
const c_m_ov_s = 299792458

# Vacuum electric permittivity [F/m] https://physics.nist.gov/cgi-bin/cuu/Value?ep0
const ε0_F_ov_m = 8.8541878188e-12

# Electron Mass [kg] https://physics.nist.gov/cgi-bin/cuu/Value?me
const me_kg = 9.1093837139e-31

# Unified atomic mass unit [kg] https://physics.nist.gov/cgi-bin/cuu/Value?ukg
const mu_kg = 1.66053906892e-27

# Boltzmann Constant [J/K] https://physics.nist.gov/cgi-bin/cuu/Value?k
const kb_J_ov_K = 1.308649e-23

# =============================
# First Stage Derived Constants
# =============================

# Reduced Plank constant [Js]
const hbar_Js = h_Js / (2 * pi)

# Fine-structure constant [1]
const α = e_C^2 / (2 * ε0_F_ov_m * h_Js * c_m_ov_s)

# ==============================
# Second Stage Derived Constants
# ==============================

# Bohr radius [m]
const a0_m = hbar_Js / (me_kg * c_m_ov_s * α)

# Rydberg constant [1/m]
const Rinf_1_ov_m = me_kg * c_m_ov_s * α^2 / (2 * h_Js)

# =======================
# UNIT CONVERSION factors
# =======================

# MULTIPLY with a frequency in [10^12*tr/s] to convert into  
const turnTHz_to_eV = 10^(12) * h_Js / e_C

# MULTIPLY with an Energy in [J] to convert into [eV]
# DIVIDE with an Energy in [eV] to convert into [J]
const J_to_eV = 1 / e_C

# MULTIPLY with an Energy in [Ryd] to convert into [10^12*tr/s] or equally [THz]
# ALTERNATIVE: const Ryd_to_turnTHz = c_m_ov_s * Rinf_1_ov_m * 10^(-12)
const Ryd_to_turnTHz = 10^(-12) * c_m_ov_s * α / (4 * pi * a0_m)

# MULTIPLY with an Energy in [Ryd] to convert into [10^12*rad/s]
# ALTERNATIVE: const RydtoTHz = 2 * pi * c_m_ov_s * Rinf_1_ov_m * 10^(-12)
const Ryd_to_radTHz = 10^(-12) * c_m_ov_s * α / (2 * a0_m)

# MULTIPLY with a Length in [nm] to convert into [Bohr]
# DIVIDE with a reciprocal Length in [1/nm] to convert into [1/Bohr]
const nm_to_bohr = 1 / (10^9 * a0_m)

# MULTIPLY with a Mass in [Dalton] to convert into [DEM]
# [DEM] is a placeholder unit of mass double the electron mass. This convention comes 
# from the Rydberg Unit System. Info at the following document.
# http://ilan.schnell-web.net/physics/rydberg.pdf
const dalton_to_2me = mu_kg / (2 * me_kg)

# MULTIPLY with a Velocity in [1/Bohr] to convert into [rad * km / s]
const Ryd_to_km_ov_s = c_m_ov_s * α / 2 * 10^(-3)
