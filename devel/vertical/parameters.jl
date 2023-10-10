# -*- coding: utf-8 -*-
using Polynomials

##############################################################################
# STRUCTURES
##############################################################################

"Geometric description of vertical reactor bounding volume."
struct VerticalReactorGeometry
    "Reactor total height [m]"
    H::Float64

    "Reactor cross-section depth [m]"
    D::Float64

    "Reactor cross-section width [m]"
    W::Float64

    "Reactor cross-section perimeter [m]"
    P::Float64

    "Reactor cross-section area [m²]"
    A::Float64

    "Reactor total volume [m³]"
    V::Float64

    function VerticalReactorGeometry(; H, D, W)
        P = 2 * (D + W)
        A = D * W
        V = A * H
        return new(H, D, W, P, A, V)
    end
end

##############################################################################
# TEMPERATURES
##############################################################################

"Reactor external environment temperature[K]"
const Tenv::Float64 = 313.0

"Reference initial gas temperature [K]"
const Tg₀::Float64 = 2073.0

"Reference initial solids temperature [K]"
const Ts₀::Float64 = 313.0

##############################################################################
# FLOW RATES
##############################################################################

"Fumes inlet gas flow rate [kg/s]"
const ṁgₛ::Float64 = 8.0

"Mass production rate of solids [kg/s]"
const ṁsₑ::Float64 = 470_000/(24*3600)

##############################################################################
# SOLIDS
##############################################################################

"Characteristic block size [m]"
const blocksize::Float64 = 0.1

"Solids specific mass [kg/m³]"
const ρₛ::Float64 = 3000.0

"Solids phase specific heat [J/(kg.K)]"
const cₚₛ = Polynomial([900.0], :T)

"Solids phase thermal conductivity [W/(m.K)]"
const kₛ = Polynomial([5.0], :T)

##############################################################################
# GEOMETRY
##############################################################################

const REACTOR = VerticalReactorGeometry(; H = 10.0, D = 2.1, W = 6.7)

"Porosity volume fraction in reactor bed [-]"
const ϕₛ::Float64 = 0.65

@info("Check `parameters.jl` for details...")