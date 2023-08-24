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

""" Ideal gas constant [J/(mol.K)]. """
const GAS_CONSTANT = 8.314_462_618_153_24

""" Base type for thermodynamic models. """
abstract type AbstractGasThermo end

""" Ideal gas phase thermodynamics model. """
struct IdealGasThermo <: AbstractGasThermo
    model::String
    temperature_ranges::Vector{Float64}
    data::Vector{Vector{Float64}}
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
    molecularmass::Float64

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
    species::Array{IdealGasSpecies,1}
end

""" Queries database and constructs species from its name. """
function IdealGasSpecies(
        speciesdata::Vector{Dict{Any, Any}},
        name::String
    )::IdealGasSpecies
    return IdealGasSpecies(getspeciesbyname(speciesdata, name))
end

function getspeciesbyname(speciesdata, name)
    return first(filter(s -> s["name"] == name, speciesdata))
end

function computemolecularmass(composition)
    return sum(n * Elements.mass(s) for (s, n) in composition)
end

mass(s::IdealGasSpecies) = s.molecularmass / 1000

""" Molar specific heat from NASA7 polynomial [J/(mol.K)]. """
function nasa7specificheat(T, c)
    p = c[1]+T.*(c[2]+T.*(c[3]+T.*(c[4]+c[5].*T)))
    return GAS_CONSTANT .* p
end

""" Molar enthalpy from NASA7 polynomial [J/mol]. """
function nasa7enthapy(T, c)
    d = c[1:5] ./ collect(1:5)
    p = d[1]+T.*(d[2]+T.*(d[3]+T.*(d[4]+d[5].*T)))+c[6]./T
    return GAS_CONSTANT .* T .* p
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

function specificheatmass(species, T)
    return species.thermo.specificheat(T) / mass(species)
end

function enthalpymass(species, T)
    return species.thermo.enthalpyheat(T) / mass(species)
end

function specificheatmole(species, T)
    return species.thermo.specificheat(T)
end

function enthalpymole(species, T)
    return species.thermo.enthalpyheat(T)
end

#     function [result] = n_species(self)
#         % Access to number of species in mechanism [-].
#         result = numel(self.mw);
#     endfunction

#     function [result] = molecular_masses(self)
#         % Access to array of molecular masses [kg/mol].
#         result = self.mw;
#     endfunction

#     function [mw] = mean_molecular_mass_xfrac(self, X)
#         % Mixture mean molecular mass from mole fractions [kg/mol].
#         mw = X * self.mw';
#     endfunction

#     function [mw] = mean_molecular_mass_yfrac(self, Y)
#         % Mixture mean molecular mass from mass fractions [kg/mol].
#         mw = (sum(Y' ./ self.mw').^(-1))';
#     endfunction

#     function [Y] = mole_to_mass_fraction(self, X)
#         % Convert mole to mass fractions.
#         Y = X .* self.mw ./ self.mean_molecular_mass_xfrac(X);
#     endfunction

#     function [X] = mass_to_mole_fraction(self, Y)
#         % Convert mass to mole fractions.
#         X = ((Y .* self.mean_molecular_mass_yfrac(Y))' ./ self.mw')';
#     endfunction

#     function [rho] = density_mass(self, T, P, Y)
#         % Mixture specific mass [kg/m³].
#         m = self.mean_molecular_mass_yfrac(Y);
#         rho = P .* m ./ (Thermodata.GAS_CONSTANT .* T);
#     endfunction

#     function [cp] = specific_heat_mass(self, T, Y)
#         % Mixture mass-averaged specific heat [J/(kg.K)].
#         cp = 0.0;
#         for k=1:self.n_species
#             cp = cp + self.species{k}.cp_mole(T) .* Y(:, k) ./ self.mw(k);
#         endfor
#     endfunction

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

end # (module IdealGas)