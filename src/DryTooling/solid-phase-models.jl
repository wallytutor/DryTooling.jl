# -*- coding: utf-8 -*-
export MaterialShomate

""" Base type for thermodynamic models. """
abstract type AbstractSolidThermo end

""" Base type for transport models. """
abstract type AbstractSolidTransport end

""" Base type for """
abstract type AbstractSolidMaterial end

"""
    MaterialShomate

**IMPORTANT:** the implementation of `makestepwise1d` used for step-wise
function evaluation takes the mean of both ranges at the `T_ch`, while the
actual Shomate uses the high range for doing so. The function enforces the
right behavior by multiplying `T_ch` the value by `1-eps()`.

$(TYPEDFIELDS)
"""
struct MaterialShomate <: AbstractSolidThermo
    """ Molar specific heat [J/(mol.K)]. """
    cₚ::Function

    """ Molar enthalpy [J/mol]. """
    h::Function

    """ Molar entropy [J/K]. """
    s::Function
    
    """ Low temperature range Shomate coefficients. """
    a_lo::Vector{Float64}

    """ High temperature range Shomate coefficients. """
    a_hi::Vector{Float64}
    
    """ Temperature of range change for evaluation. """
    T_ch::Float64

    function MaterialShomate(;
            a_lo::Vector{Float64},
            a_hi::Vector{Float64},
            T_ch::Float64
        )
        T_ch = (1.0 - eps()) * T_ch

        c = makestepwise1d(
            (T) -> shomatespecificheat(T/1000.0, a_lo),
            (T) -> shomatespecificheat(T/1000.0, a_hi),
            T_ch; differentiable = true
        )
        h = makestepwise1d(
            (T) -> shomateenthalpy(T/1000.0, a_lo),
            (T) -> shomateenthalpy(T/1000.0, a_hi),
            T_ch; differentiable = true
        )
        s = makestepwise1d(
            (T) -> shomateentropy(T/1000.0, a_lo),
            (T) -> shomateentropy(T/1000.0, a_hi),
            T_ch; differentiable = true
        )

        return new(c, h, s, a_lo, a_hi, T_ch)
    end
end

"""
    MaterialTransportPolynomial

Transport properties for a solid material.

$(TYPEDFIELDS)
"""
struct MaterialTransportProperties <: AbstractSolidTransport
    """ Thermal conductivity [W/(m.K)]. """
    k::Function
    
    """ Emissivity [-]. """
    ε::Function
    
    function MaterialTransportProperties(;
            k::Function,
            ε::Function
        )
        return new(k, ε)
    end
end

"""
    MaterialPowderBed

Description of a powder bed material for a rotary kiln.

$(TYPEDFIELDS)
"""
struct MaterialPowderBed <: AbstractSolidMaterial
    """ Density [kg/m³]. """
    ρ::Float64

    """ Repose angle [rad]. """
    γ::Float64

    """ Solid packing fraction [-]. """
    ϕ::Float64

    """ Particle mean diameter [m]. """
    d::Float64

    """ Molar mass [kg/mol]. """
    M::Float64

    """ Thermal conductivity [W/(m.K)]. """
    k::Function

    """ Emissivity [-]. """
    ε::Function

    """ Molar specific heat [J/(mol.K)]. """
    cₚ::Function

    """ Molar enthalpy [J/mol]. """
    h::Function

    """ Molar entropy [J/K]. """
    s::Function

    """ Access to thermodynamic model. """
    thermo::AbstractSolidThermo

    """ Access to transport model. """
    transport::AbstractSolidTransport

    function MaterialPowderBed(;
            ρ::Float64,
            γ::Float64,
            ϕ::Float64,
            d::Float64,
            M::Float64,
            thermo::AbstractSolidThermo,
            transport::AbstractSolidTransport
        )
        return new(
            ρ, γ, ϕ, d, M,
            transport.k,
            transport.ε,
            thermo.cₚ,
            thermo.h,
            thermo.s,
            thermo, transport
        )
    end
end

#############################################################################
# For YAML parsing
#############################################################################

function MaterialShomate(data::Dict{Any, Any})
    return MaterialShomate(; a_lo = data["coefs_low"],
                             a_hi = data["coefs_high"],
                             T_ch = data["change_temperature"])
end

function MaterialTransportProperties(data::Dict{Any, Any})
    return MaterialTransportProperties(;
        k = eval(Meta.parse(data["thermal_conductivity"])),
        ε = eval(Meta.parse(data["emissivity"])))
end

function MaterialPowderBed(data::Dict{Any, Any})
    model = getfield(DryTooling, Symbol(data["thermo"]["type"]))

    return MaterialPowderBed(;
        ρ = data["density"],
        γ = deg2rad(data["repose_angle"]),
        ϕ = data["solid_filling"],
        d = data["particle_diam"],
        M = data["molar_mass"],
        thermo = model(data["thermo"]),
        transport = MaterialTransportProperties(data["transport"])
    )
end

#############################################################################
# Private
#############################################################################

"""
    shomatespecificheat(T::Float64, c::Vector{Float64})::Float64

Molar specific heat with Shomate equation [J/(mol.K)].
"""
function shomatespecificheat(T::Float64, c::Vector{Float64})::Float64
    p0 = c[3] + T * c[4]
    p1 = c[2] + T * p0
    return T * p1 + c[5] / T^2 + c[1]
end

"""
    shomateenthalpy(T::Float64, c::Vector{Float64})::Float64

Molar enthalpy with Shomate equation [J/mol].
"""
function shomateenthalpy(T::Float64, c::Vector{Float64})::Float64
    p0 = (c[3]/3) + T * (c[4]/4)
    p1 = (c[2]/2) + T * p0
    p2 = (c[1]/1) + T * p1
    return T * p2 - (c[5]/T) + c[6] - c[8]
end

"""
    shomateentropy(T::Float64, c::Vector{Float64})::Float64

Entropy with Shomate equation [J/K].
"""
function shomateentropy(T::Float64, c::Vector{Float64})::Float64
    p0 = (c[3]/2) + T * (c[4]/3)
    p1 = (c[2]/1) + T * p0
    return c[1] * log(T) + T * p1 + c[5]/(2T^2) + c[7]
end
