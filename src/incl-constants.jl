# -*- coding: utf-8 -*-
export GAS_CONSTANT
export ZERO_CELSIUS
export ONE_ATM
export STEFAN_BOLTZMANN

"Ideal gas constant [$(GAS_CONSTANT) ``J mol^{-1} K^{-1}``]."
const GAS_CONSTANT::Float64 = 8.314_462_618_153_24

"Zero degrees Celsius in Kelvin [$(ZERO_CELSIUS) ``K``]."
const ZERO_CELSIUS::Float64 = 273.15

"Atmospheric pressure at sea level [$(ONE_ATM) ``Pa``]."
const ONE_ATM::Float64 = 101325.0

"Stefan-Boltzmann constant [$(STEFAN_BOLTZMANN) ``W m^{-2} K^{-4}``]"
const STEFAN_BOLTZMANN = 5.670374419e-08
