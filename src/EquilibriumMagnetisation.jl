export mag_factor
"""
    m0 = mag_factor(temp)

Returns the ratio between the magnetisation of completely saturated water (100% water content) and the ambient magnetic field at thermal equilibrium.

Parameters:
- `temp`: temperature in Kelvin.

Returns:
- `m0`: magnetisation ratio, in Amperes per (metre Tesla).
"""
mag_factor(temp) = γh^2 * ħ^2/(2*kb*temp) * nh2o
