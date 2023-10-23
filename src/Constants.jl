# -*- coding: utf-8 -*-
module Constants

export GAS_CONSTANT
export ZERO_CELSIUS
export ONE_ATM

"Ideal gas constant [J/(mol.K)]."
const GAS_CONSTANT::Float64 = 8.314_462_618_153_24

"Zero degrees Celsius in Kelvin for conversion [K]."
const ZERO_CELSIUS::Float64 = 273.15

"Atmospheric pressure at sea level [Pa]."
const ONE_ATM::Float64 = 101325.0

end # module Constants