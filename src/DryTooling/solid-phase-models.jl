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
    cₚ::Function
    h::Function
    s::Function
    a_lo::Vector{Float64}
    a_hi::Vector{Float64}
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

        new(c, h, s, a_lo, a_hi, T_ch)
    end
end

struct MaterialTransportPolynomial <: AbstractSolidTransport
end

struct MaterialPowderBed <: AbstractSolidMaterial
    thermo::AbstractSolidThermo
    transport::AbstractSolidTransport

    function MaterialPowderBed(;
            thermo::AbstractSolidThermo
            transport::AbstractSolidTransport
        )

        return new(thermo, transport)
    end
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
