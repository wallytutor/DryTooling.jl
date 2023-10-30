# -*- coding: utf-8 -*-
import Pkg
import Revise

if Base.current_project() != Base.active_project()
    Pkg.activate(Base.current_project())
    Pkg.resolve()
    Pkg.instantiate()

    using CairoMakie
    using CommonSolve
    using CommonSolve: solve
    using DocStringExtensions: TYPEDFIELDS
    using Trapz: trapz
    using DryTooling
    using DryTooling.Simulation
    using DryTooling: Temperature1DModelStorage
    using DryTooling: interfaceconductivity1D
end



# "Thermal diffusion in a sphere represented in temperature space."
# struct Sphere1DEnthalpyModel <: dry.AbstractDiffusionModel1D
#     "Grid over which problem will be solved."
#     grid::dry.AbstractGrid1D

#     "Memory for model linear algebra problem."
#     problem::dry.TridiagonalProblem

#     "Constant part of model coefficient α."
#     α′::Vector{Float64}

#     "Constant part of model coefficient β."
#     β′::Vector{Float64}

#     "Thermal conductivity in terms of temperature."
#     κ::Function

#     "Enthalpy in terms of temperature."
#     h::Function

#     "Global heat transfer coefficient ``U=hR²``."
#     U::Float64

#     "Surface environment temperature."
#     B::Float64

#     "Time-step used in integration."
#     τ::Base.RefValue{Float64}

#     "Memory storage for solution retrieval."
#     mem::Base.RefValue{dry.Temperature1DModelStorage}

#     function Sphere1DEnthalpyModel(;
#             grid::dry.AbstractGrid1D,
#             h::Function,
#             κ::Function,
#             ρ::Float64,
#             u::Float64,
#             B::Float64
#         )
#         problem = dry.TridiagonalProblem(grid.N)

#         rₙ = dry.tail(grid.w)
#         rₛ = dry.head(grid.w)
#         α′ = @. ρ * (rₙ^3 - rₛ^3) / 3.0

#         rₙ = dry.tail(grid.r)
#         rₛ = dry.head(grid.r)
#         wⱼ = dry.body(grid.w)
#         β′ = @. wⱼ^2 / (rₙ - rₛ)

#         U = u * last(grid.r)^2
#         τ = Ref(-Inf)
#         mem = Ref(dry.Temperature1DModelStorage(0, 0))

#         return new(grid, problem, α′, β′, κ, h, U, B, τ, mem)
#     end
# end

# "Time-step dependent updater for model."
# function sphereenthalpyouter!(
#         m::Sphere1DEnthalpyModel,
#         t::Float64,
#         n::Int64
#     )::Nothing
#     # Note the factor 4π because U = r²h only and A = 4πr²!!!!
#     m.mem[].Q[n] = 4π * m.U * (m.B - last(m.problem.x))
#     m.mem[].T[n, 1:end] = m.problem.x

#     @. m.problem.b[1:end] = m.h.(m.problem.x)
#     m.problem.b[end] += m.U * m.B

#     return nothing
# end

# "Non-linear iteration updater for model."
# function sphereenthalpyinner!(
#         m::Sphere1DEnthalpyModel,
#         t::Float64,
#         n::Int64
#     )::Nothing
#     a_p = m.problem.A.d
#     a_s = m.problem.A.dl
#     a_n = m.problem.A.du
#     T_p = m.problem.x

#     κ = interfaceconductivity1D(m.κ.(T_p))
#     β = κ .* m.β′
#     α = m.α′./ m.τ

#     a_s[1:end] = -β
#     a_n[1:end] = -β
#     a_p[1:end] = α

#     a_p[2:end-1] += tail(β) + head(β)
#     a_p[1]       += first(β)
#     a_p[end]     += last(β) + m.U

#     return nothing
# end

# function sphereenthalpysolve!()
# end
