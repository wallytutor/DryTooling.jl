# -*- coding: utf-8 -*-
export SphereTemperatureModel
export solve

#############################################################################
# Helpers
#############################################################################

"Data storage for 1D temperature solution models."
struct Temperature1DModelStorage
    "Tracker for surface heat flux."
    Q::Vector{Float64}

    "Tracker for intermediate temperature solutions."
    T::Matrix{Float64}

    function Temperature1DModelStorage(N, M)
        return new(zeros(M+1), zeros(M+1, N+1))
    end
end

#############################################################################
# Models
#############################################################################

"Thermal diffusion in a cylinder represented in temperature space."
struct Cylinder1DTemperatureModel <: AbstractDiffusionModel1D
    "Grid over which problem will be solved."
    grid::AbstractGrid1D

    "Memory for model linear algebra problem."
    problem::TridiagonalProblem

    "Constant part of model coefficient α."
    α′::Vector{Float64}

    "Constant part of model coefficient β."
    β′::Vector{Float64}

    "Thermal conductivity in terms of temperature."
    κ::Function

    "Global heat transfer coefficient ``U=hR``."
    U::Float64

    "Surface environment temperature."
    B::Float64

    "Time-step used in integration."
    τ::Base.RefValue{Float64}

    "Memory storage for solution retrieval."
    mem::Base.RefValue{Temperature1DModelStorage}

    function Cylinder1DTemperatureModel(;
            grid::AbstractGrid1D,
            κ::Function,
            ρ::Float64,
            c::Float64,
            h::Float64,
            B::Float64
        )
        problem = TridiagonalProblem(grid.N)

        rₙ = tail(grid.w)
        rₛ = head(grid.w)
        α′ = @. ρ * c * (rₙ^2 - rₛ^2) / 2.0

        rₙ = tail(grid.r)
        rₛ = head(grid.r)
        wⱼ = head(tail(grid.w))
        β′ = @. wⱼ / (rₙ - rₛ)

        U = h * last(grid.r)
        τ = Ref(-Inf)
        mem = Ref(Temperature1DModelStorage(0, 0))

        return new(grid, problem, α′, β′, κ, U, B, τ, mem)
    end
end

"Thermal diffusion in a sphere represented in temperature space."
struct SphereTemperatureModel <: AbstractPhysicalModel
    problem::TridiagonalProblem
    z::Vector{Float64}
    α::Vector{Float64}
    β::Vector{Float64}
    U::Float64
    T∞::Float64
    k::Function
    t::Float64
    τ::Float64
    Q::Vector{Float64}
    T::Matrix{Float64}

    function SphereTemperatureModel(;
            N::Int64,
            R::Float64,
            ρ::Float64,
            c::Float64,
            h::Float64,
            T∞::Float64,
            T₀::Float64,
            k::Function,
            t::Union{Float64, Nothing} = nothing,
            M::Union{Int64, Nothing} = nothing
        )
        # Compute final integration time if required.
        tend = isnothing(t) ? 2*ρ*c*R^2/k(0.5(T∞+T₀)) : t
        steps = isnothing(M) ? convert(Int64, round(tend/10)) : M

        # Space discretization.
        δr = R / N
        z = collect(0.0:δr:R)
        w = collect(0.5δr:δr:R-0.5δr)
        r = vcat(0.0, w, R)

        # Increments.
        δ = z[2:end-0] - z[1:end-1]
        τ = tend / steps

        # To handle boundaries use ``r`` for computing α.
        α = @. (r[2:end-0]^3 - r[1:end-1]^3) * ρ * c / (3τ)

        # For β use only internal walls ``w``.
        β = @. w^2 / δ

        # Heat transfer coefficient multiplied by area.
        U = R^2 * h

        # Create linear problem memory.
        problem = TridiagonalProblem(N)
        problem.x[:] .= T₀

        # Memory for storing surface fluxes.
        Q = zeros(steps+2)

        # Memory for intermediate solution.
        T = zeros(steps+2, N+1)

        return new(problem, z, α, β, U, T∞, k, tend, τ, Q, T)
    end
end

#############################################################################
# Initialization
#############################################################################

"Set initial condition of thermal diffusion model."
function initialize(
        m::AbstractDiffusionModel1D,
        mem::Temperature1DModelStorage,
        T::Float64,
        τ::Float64,
    )::Nothing
    m.problem.x[:] .= T
    m.τ[] = τ
    m.mem[] = mem
    return nothing
end

#############################################################################
# Solution
#############################################################################

"Interface for solving a `Cylinder1DTemperatureModel` instance."
function CommonSolve.solve(m::Cylinder1DTemperatureModel; kwargs...)
    kwargs_dict = Dict(kwargs)

    tend = kwargs_dict[:t]
    τ = kwargs_dict[:τ]
    T = kwargs_dict[:T]

    delete!(kwargs_dict, :t)
    delete!(kwargs_dict, :τ)
    delete!(kwargs_dict, :T)

    steps = convert(Int64, round(tend / τ))
    mem = Temperature1DModelStorage(m.grid.N, steps)
    initialize(m, mem, T, τ)

    return relaxationouterloop(;
        model = m,
        updaterouter = cylindertemperatureouter!,
        updaterinner = cylindertemperatureinner!,
        tend = tend,
        steps = steps,
        kwargs_dict...
    )
end

"Interface for solving a `SphereTemperatureModel` instance."
function CommonSolve.solve(model::SphereTemperatureModel; kwargs...)
    return relaxationouterloop(;
        model = model,
        updaterouter = spheretemperatureouter!,
        updaterinner = spheretemperatureinner!,
        kwargs...
    )
end

#############################################################################
# Internals
#############################################################################

"Non-linear iteration updater for model."
function cylindertemperatureinner!(
        m::Cylinder1DTemperatureModel,
        t::Float64,
        n::Int64
    )::Nothing
    a_p = m.problem.A.d
    a_s = m.problem.A.dl
    a_n = m.problem.A.du
    T_p = m.problem.x

    κ = interfaceconductivity1D(m.κ.(T_p))
    β = κ .* m.β′
    α = m.α′./ m.τ

    a_s[1:end] = -β
    a_n[1:end] = -β
    a_p[1:end] = α

    a_p[2:end-1] += tail(β) + head(β)
    a_p[1]       += first(β)
    a_p[end]     += last(β) + m.U

    return nothing
end

"Time-step dependent updater for model."
function cylindertemperatureouter!(
        m::Cylinder1DTemperatureModel,
        t::Float64,
        n::Int64
    )::Nothing
    # Follow surface heat flux and store partial solutions.
    # XXX: note the factor 2π because U = rh only and A = 2πrl!!!!
    m.mem[].Q[n] = 2π * m.U * (m.B - last(m.problem.x))
    m.mem[].T[n, 1:end] = m.problem.x

    @. m.problem.b[1:end] = (m.α′ / m.τ) * m.problem.x
    m.problem.b[end] += m.U * m.B
    return nothing
end

"Non-linear iteration updater for model."
function spheretemperatureinner!(
        model::SphereTemperatureModel,
        t::Float64,
        n::Int64
    )::Nothing
    # Aliases to approach mathermatical formulation.
    a_p = model.problem.A.d
    a_s = model.problem.A.dl
    a_n = model.problem.A.du

    # Wall interpolated thermal conductivities.
    kp = @. model.k(model.problem.x)
    ks = kp[1:end-1]
    kn = kp[2:end+0]
    κ = @. 2 * ks * kn / (ks + kn)

    # Update temperature dependency of β.
    βₖ = @. κ * model.β

    # Main, lower and upper diagonal elements.
    a_s[1:end] = -βₖ
    a_n[1:end] = -βₖ
    a_p[1:end] = model.α

    # Update main but leave boundaries off.
    @. a_p[2:end-1] += βₖ[2:end] + βₖ[1:end-1]

    # Boundary conditions on main.
    a_p[1] += βₖ[1]
    a_p[end] += βₖ[end] + model.U

    return nothing
end

"Time-step dependent updater for model."
function spheretemperatureouter!(
        model::SphereTemperatureModel,
        t::Float64,
        n::Int64
    )::Nothing
    # Follow surface heat flux and store partial solutions.
    # XXX: note the factor 4π because U = r²h only and A = 4πr²!!!!
    model.Q[n] = 4π * model.U * (model.T∞ - model.problem.x[end])
    model.T[n, 1:end] = model.problem.x[1:end]

    # Problem right-hand side.
    @. model.problem.b[1:end] = model.α * model.problem.x

    # Apply boundary condition.
    model.problem.b[end] += model.U * model.T∞

    return nothing
end

#############################################################################
# Utilities
#############################################################################

"Dimensional analysis time-scale of diffusion process."
function diffusiontimescale(R, α)
    return 2 * R^2 / α
end

"Dimensional analysis time-scale of diffusion process."
function diffusiontimescale(R, ρ, c, κ)
    return diffusiontimescale(R, κ / (ρ * c))
end

"Interface thermal conductivity assuming equidistant centers."
function interfaceconductivity1D(κ::Vector{Float64})
    κₛ, κₙ = head(κ), tail(κ)
    return @. 2 * κₛ * κₙ / (κₛ + κₙ)
end

"Reconstruct time axis of integrated diffusion model."
function timeaxis(m::AbstractDiffusionModel1D)
    tend = (size(m.mem[].Q)[1] - 1) * m.τ[]
    return 0.0:m.τ[]:tend
end
