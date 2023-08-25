# -*- coding: utf-8 -*-

""" Base type for transport models. """
abstract type AbstractTransportModel end

""" Base type for thermodynamic models. """
abstract type AbstractGasThermo end

""" Named access to transport models. """
const TRANSPORT_MODELS = Dict("gas" => :LennardJonesTransport)

""" Lennard-Jones ideal gas transport model. """
struct LennardJonesTransport <: AbstractTransportModel
    geometry::String
    welldepth::Float64
    diameter::Float64
    polarizability::Float64
    rotationalrelaxation::Float64

    function LennardJonesTransport(transport)
        return new(
            transport["geometry"],
            get(transport, "well-depth", 0.0),
            get(transport, "diameter", 0.0),
            get(transport, "polarizability", 0.0),
            get(transport, "rotational-relaxation", 0.0)
        )
    end
end

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
        model = getfield(DryTooling, TRANSPORT_MODELS[transport])

        new(species["name"],
            composition,
            model(species["transport"]),
            IdealGasThermo(species["thermo"], verbose = verbose),
            sum(n * elementmass(s) for (s, n) in composition))
    end
end

""" Species specific heat in mass units [J/(kg.K)]. """
function specificheatmass(species::IdealGasSpecies, T)
    return species.thermo.specificheat(T) / mass(species)
end

""" Species specific heat in mole units [J/(mol.K)]. """
function specificheatmole(species::IdealGasSpecies, T)
    return species.thermo.specificheat(T)
end

""" Species enthalpy in mass units [J/kg]. """
function enthalpymass(species::IdealGasSpecies, T)
    return species.thermo.enthalpy(T) / mass(species)
end

""" Species enthalpy in mole units [J/mol]. """
function enthalpymole(species::IdealGasSpecies, T)
    return species.thermo.enthalpy(T)
end

""" Ideal gas phase mixture model. """
struct IdealGasMixture
    species::Vector{IdealGasSpecies}
    molecularmasses::Vector{Float64}
    nspecies::Int32

    function IdealGasMixture(
            data::Dict{Any, Any},
            selected::Vector{String}
        )
        nspecies = length(selected)
        species = Vector{IdealGasSpecies}(undef, nspecies)
        molecularmasses = zeros(nspecies)

        for (i, name) in enumerate(selected)
            thisone = getnameditem(data["species"], name)
            species[i] = IdealGasSpecies(thisone)
            molecularmasses[i] = mass(species[i])
        end

        return new(species, molecularmasses, nspecies)
    end
end

""" Mixture mean molecular mass [kg/mol]. """
function meanmolecularmass(mix::IdealGasMixture, Y)
    return meanmolecularmass(mix.molecularmasses, Y)
end

""" Mixture specific mass [kg/m³]. """
function densitymass(mix::IdealGasMixture, T, P, Y)
    return P * meanmolecularmass(mix, Y) / (GAS_CONSTANT * T)
end

""" Mixture concentration [mol/m³]. """
function concentration(mix::IdealGasMixture, T, P, Y)
    return densitymass(mix, T, P, Y) * (@. Y / mix.molecularmasses)
end

""" Mixture mass-averaged specific heat [J/(kg.K)]. """
function specificheatmass(mix::IdealGasMixture, T, Y)
    contrib(s, y) = specificheatmass(s, T) * y
    return sum(contrib(s, y) for (s, y) ∈ zip(mix.species, Y))
end

#         % Mixture mass-averaged enthalpy [J/kg].
#     function [h] = enthalpy_mass(self, T, Y)
#         h = sum((Y .* self.enthalpies_mass(T))')';
#     endfunction

#         % Matrix of species enthalpies [J/kg].
#     function [hs] = enthalpies_mass(self, T)
#         hs = [];
#         for k=1:self.n_species
#             hs = horzcat(hs, self.species{k}.enthalpy_mole(T) ./ self.mw(k));
#         endfor
#     endfunction

#         % Heat release rate [W/m³].
#     function hdot = heat_release_rate(self, h, mdotk)
#         hdot = sum((mdotk .* h)')';
#     endfunction

# mutable struct IdealGasSolution
#     mix::IdealGasMixture
#     T::Num
#     P::Num
#     Y::AbstractArray

#     function IdealGasSolution(mix::IdealGasMixture)
#         new(mix, 300.0, ONE_ATM, zeros(mix.nspecies))
#     end
# end

# """ Mixture composition in mole fractions. """
# function massfractions(gas::IdealGasSolution)
#     return gas.Y
# end

# """ Mixture composition in mole fractions. """
# function molefractions(gas::IdealGasSolution)
#     return massfraction2molefraction(gas.mix.molecularmasses, gas.Y)
# end

#############################################################################
# Private
#############################################################################

""" Retrieve atomic mass of species [kg/mol]. """
mass(s::IdealGasSpecies) = s.molecularmass / 1000

""" Query first item matching name in dictionary. """
getnameditem(data, name) = first(filter(s -> s["name"] == name, data))

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

    cpfun = getfield(DryTooling, Symbol(cpname))
    hmfun = getfield(DryTooling, Symbol(hmname))

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

#         function [wdot] = wdot_mak(self, z, T, Y, L)
#             % Mass action kinetics methane combustion rate [kg/(m³.s)].
#             k0 = 1.1e+07;
#             Ea = 83680.0;

#             X = self.mass_to_mole_fraction(Y);
#             C = (X * self.PRESSURE ./ (self.GAS_CONSTANT .* T));
#             k = k0 * exp(-Ea ./ (self.GAS_CONSTANT .* T));
#             rt = k .* C(:, 1) .* C(:, 2).^0.5;

#             wdot = rt * (self.mw .* self.SPECIES_COEFS);
#         endfunction

#         function [wdot] = wdot_ebu(self, z, T, Y, L)
#             % Eddy break-up kinetics methane combustion rate [kg/(m³.s)].
#             cr = 4.000e+00;
#             bo = 4.375e+00;
#             k0 = 1.600e+10;
#             Ea = 1.081e+05;

#             k = k0 * exp(-Ea ./ (self.GAS_CONSTANT .* T));
#             rho = self.density_mass(T, self.PRESSURE, Y);
#             yf = Y(:, 1);
#             yo = Y(:, 2);

#             % TODO implement this in ProjectData
#             ke = z ./ L;

#             R_ebu = (rho.^1) .* cr .* ke .* min(yf, yo ./ bo);
#             R_arr = (rho.^2) .* yf .* yo .* k;

#             rt = min(R_ebu, R_arr) / self.mw(1);

#             wdot = rt * (self.mw .* self.SPECIES_COEFS);
#         endfunction

#         function [mu] = viscosity(self, T, Y)
#             % Gas molecular viscosity [Pa.s].
#             mu = 1.0e-05 * (0.1672 * sqrt(T) - 1.058);
#         endfunction

#         function [k] = thermal_conductivity(self, T, Y)
#             % Gas thermal conductivity [W/(m³.K)].
#             k = 1.581e-17;
#             k = T .* k - 9.463e-14;
#             k = T .* k + 2.202e-10;
#             k = T .* k - 2.377e-07;
#             k = T .* k + 1.709e-04;
#             k = T .* k - 7.494e-03;
#         endfunction

#     % Species stoichiometric coefficients.
#     SPECIES_COEFS = [-1.0, -2.0, 1.0, 2.0, 0.0, 0.0];

#     % Operating pressure [Pa].
#     PRESSURE = 101325.0;
