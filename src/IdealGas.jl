# -*- coding: utf-8 -*-
module IdealGas

using ModelingToolkit: Num
using DryTooling.General: makestepwise1d
using DryTooling.Transport: AbstractTransportModel
import DryTooling.Elements as Elements
import DryTooling.Transport as Transport

export AbstractGasThermo
export IdealGasThermo
export IdealGasSpecies
export IdealGasMixture
export specificheatmass
export specificheatmole
export enthalpymass
export enthalpymole
export meanmolecularmass
export massfraction2molefraction
export molefraction2massfraction
export densitymass
export massfractions
export molefractions

""" Ideal gas constant [J/(mol.K)]. """
const GAS_CONSTANT = 8.314_462_618_153_24

""" Zero degrees Celsius in Kelvin for conversion [K]. """
const ZERO_CELSIUS = 273.15

""" Atmospheric pressure at sea level [Pa]. """
const ONE_ATM = 101325.0

f64 = Float64
vf64 = Vector{f64}

#############################################################################
# Types
#############################################################################

""" Base type for thermodynamic models. """
abstract type AbstractGasThermo end

""" Ideal gas phase thermodynamics model. """
struct IdealGasThermo <: AbstractGasThermo
    model::String
    temperature_ranges::vf64
    data::Vector{vf64}
    specificheat::Function
    enthalpy::Function

    function IdealGasThermo(thermo; verbose = true)
        model = lowercase(thermo["model"])
        rngs = thermo["temperature-ranges"]
        data = thermo["data"]
        func = getthermo(model, data, rngs..., verbose)
        return new(model, rngs, data, func[1], func[2])
    end
end

""" Ideal gas phase species model. """
struct IdealGasSpecies
    name::String
    composition::Dict{String, Int64}
    transport::AbstractTransportModel
    thermo::IdealGasThermo
    molecularmass::f64

    function IdealGasSpecies(species; verbose = true)
        composition = species["composition"]
        transport = species["transport"]["model"]
        model = getfield(Transport, Transport.MODELS[transport])

        new(species["name"],
            composition,
            model(species["transport"]),
            IdealGasThermo(species["thermo"], verbose = verbose),
            computemolecularmass(composition))
    end
end

""" Ideal gas phase mixture model. """
struct IdealGasMixture
    species::Vector{IdealGasSpecies}
    nspecies::Int32

    T::f64
    P::f64
    Y::vf64

    molecularmasses::vf64

    function IdealGasMixture(data, selected; phasename = "gas")
        nspecies = length(selected)
        species = Vector{IdealGasSpecies}(undef, nspecies)
        molecularmasses = zeros(nspecies)
        Y = zeros(nspecies)

        phase = getphase(data["phases"], phasename)

        if haskey(phase, "state")
            T = getkey(phase, "T", 300.0)
            P = getkey(phase, "P", ONE_ATM)
        end

        for (i, name) in enumerate(selected)
            species[i] = IdealGasSpecies(data["species"], name)
            molecularmasses[i] = mass(species[i])
        end

        return new(species, nspecies, T, P, Y, molecularmasses)
    end
end

""" Queries database and constructs species from its name. """
function IdealGasSpecies(
        speciesdata::Vector{Dict{Any, Any}},
        name::String
    )::IdealGasSpecies
    return IdealGasSpecies(getnameditem(speciesdata, name))
end

#############################################################################
# Public (IdealGasSpecies)
#############################################################################

function specificheatmass(species::IdealGasSpecies, T::f64)::f64
    return species.thermo.specificheat(T) / mass(species)
end

function enthalpymass(species::IdealGasSpecies, T::f64)::f64
    return species.thermo.enthalpyheat(T) / mass(species)
end

function specificheatmole(species::IdealGasSpecies, T::f64)::f64
    return species.thermo.specificheat(T)
end

function enthalpymole(species::IdealGasSpecies, T::f64)::f64
    return species.thermo.enthalpyheat(T)
end

#############################################################################
# Public (by definitions)
#############################################################################

""" Mixture mean molecular mass [kg/mol]. """
function meanmolecularmass(M::vf64, Y::vf64)::f64
    return  1.0 / sum(@. Y / M)
end

""" Convert mass fractions to mole fractions. """
function massfraction2molefraction(M::vf64, Y::vf64)::vf64
    return @. Y * meanmolecularmass(M, Y) / M
end

""" Convert mole fractions to mass fractions. """
function molefraction2massfraction(M::vf64, X::vf64)::vf64
    return @. X * M / sum(X * M)
end

#############################################################################
# Public (IdealGasMixture)
#############################################################################

""" Mixture mean molecular mass [kg/mol]. """
function meanmolecularmass(mix::IdealGasMixture)::f64
    return meanmolecularmass(mix.molecularmasses, mix.Y)
end

""" Convert mass fractions to mole fractions. """
function massfraction2molefraction(mix::IdealGasMixture)::vf64
    return massfraction2molefraction(mix.molecularmasses, mix.Y)
end

""" Mixture specific mass [kg/m³]. """
function densitymass(mix::IdealGasMixture)::f64
    return mix.P * meanmolecularmass(mix) / (GAS_CONSTANT * mix.T)
end

""" Mixture composition in mole fractions. """
function massfractions(mix::IdealGasMixture)::vf64
    return mix.Y
end

""" Mixture composition in mole fractions. """
function molefractions(mix::IdealGasMixture)::vf64
    return massfraction2molefraction(mix)
end

""" Mixture mass-averaged specific heat [J/(kg.K)] """
function specificheatmass(mix::IdealGasMixture)::f64
    contrib(s, y) = specificheatmass(s, mix.T) * y
    return sum(contrib(s, y) for (s, y) ∈ zip(mix.species, mix.Y))
end

#     function [h] = enthalpy_mass(self, T, Y)
#         % Mixture mass-averaged enthalpy [J/kg].
#         h = sum((Y .* self.enthalpies_mass(T))')';
#     endfunction

#     function [hs] = enthalpies_mass(self, T)
#         % Matrix of species enthalpies [J/kg].
#         hs = [];
#         for k=1:self.n_species
#             hs = horzcat(hs, self.species{k}.enthalpy_mole(T) ./ self.mw(k));
#         endfor
#     endfunction

#     function hdot = heat_release_rate(self, h, mdotk)
#         % Heat release rate [W/m³].
#         hdot = sum((mdotk .* h)')';
#     endfunction

#############################################################################
# Private
#############################################################################

mass(s::IdealGasSpecies) = s.molecularmass / 1000
getnameditem(data, name) = first(filter(s -> s["name"] == name, data))

function computemolecularmass(composition)
    return sum(n * Elements.mass(s) for (s, n) in composition)
end

""" Molar specific heat from NASA7 polynomial [J/(mol.K)]. """
function nasa7specificheat(T, c)
    p = c[1]+T*(c[2]+T*(c[3]+T*(c[4]+c[5]*T)))
    return GAS_CONSTANT * p
end

""" Molar enthalpy from NASA7 polynomial [J/mol]. """
function nasa7enthapy(T, c)
    d = c[1:5] / collect(1:5)
    p = d[1]+T*(d[2]+T*(d[3]+T*(d[4]+d[5]*T)))+c[6]/T
    return GAS_CONSTANT * T * p
end

""" Create specific heat and enthalpy functions for species. """
function getthermo(model, data, xl, xc, xh, verbose)
    cpname = string(model, "specificheat")
    hmname = string(model, "enthapy")

    cpfun = getfield(IdealGas, Symbol(cpname))
    hmfun = getfield(IdealGas, Symbol(hmname))

    cp = makestepwise1d(T -> cpfun(T, data[1]),
                        T -> cpfun(T, data[2]),
                        xc, differentiable = true)

    hm = makestepwise1d(T -> hmfun(T, data[1]),
                        T -> hmfun(T, data[2]),
                        xc, differentiable = true)

    function prewarning(T, f)
        if !(T isa Num) && (T < xl || T > xh)
            @warn "Temperature out of range = $(T)K"
        end
        return f(T)
    end

    specificheat = verbose ? (T -> prewarning(T, cp)) : cp
    enthalpy = verbose ? (T -> prewarning(T, hm)) : hm
    return specificheat, enthalpy
end

function getphase(data, name)
    try
        return getnameditem(data, name)
    catch
        return data[1]
    end
end

end # (module IdealGas)