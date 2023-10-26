# -*- coding: utf-8 -*-
export step!
export advance!
export fouter!
export finner!
export fsolve!
export timepoints
export relaxationstep!

"""
    step!(
        m::AbstractPhysicalModel,
        t::Float64,
        n::Int64;
        fouter!::Function,
        finner!::Function,
        fsolve!::Function,
        α::Float64 = 0.1,
        iters::Int64 = 20,
        tol::Float64 = 1.0e-10
    )

Manage the integration of a model `m` from time `t` corresponding to
step call `n` using model internal time step. All the updates of
coefficients and solution are performed through user-supplied functions.
"""
function step!(
        m::AbstractPhysicalModel,
        t::Float64,
        n::Int64;
        α::Float64 = 0.1,
        ε::Float64 = 1.0e-10,
        M::Int64 = 20
    )::Nothing
    fouter!(m, t, n)
    finner!(m, t, n)

    if α <= 0.0
        solve!(m.problem)
        m.res[].innersteps[n] = 1
        return nothing
    end

    for niter in 1:M
        if fsolve!(m, t, n, α) <= ε
            m.res[].innersteps[n] = niter
            return nothing
        end
        finner!(m, t, n)
    end

    @warn "Did not converge after $(M) steps"
    m.res[].innersteps[n] = M
    return nothing
end

"""
    advance!(
        m::AbstractPhysicalModel;
        α::Float64 = 0.1,
        ε::Float64 = 1.0e-10,
        M::Int64 = 20,
        t0::Float64 = 0.0
    )

Manage execution of `step!` over the integration time interval.
"""
function advance!(
        m::AbstractPhysicalModel;
        α::Float64 = 0.1,
        ε::Float64 = 1.0e-10,
        M::Int64 = 20,
        t0::Float64 = 0.0
    )
    if α >= 1.0 || α < 0.0
        @error """\
        Relaxation factor out-of-range (α = $(α)), this is not \
        supported / does not make physical sense in most cases!
        """
    end

    times = timepoints(m)

    # XXX: `times` contain all solution points from initial condition.
    # For time-advance use only the `head` of this to avoid performing
    # and extra step. Uncomment the `@info` for checking if in doubt!
    for (n, t) in enumerate(head(times))
        tn = t + t0
        # @info "Advancing from $(tn) to $(tn+m.τ[]) s"
        step!(m, tn, n; α, ε, M)
    end

    # @info "Reached time $(last(times))"
    fouter!(m, last(times) + t0, length(times))
end

"""
    fouter!(::AbstractPhysicalModel, ::Float64, ::Int64)

Outer loop update for [`step!`](@ref).
"""
function fouter!(::AbstractPhysicalModel, ::Float64, ::Int64)
    @error "An specialization of this method is expexted!"
end

"""
    finner!(::AbstractPhysicalModel, ::Float64, ::Int64)

Inner loop update for [`step!`](@ref).
"""
function finner!(::AbstractPhysicalModel, ::Float64, ::Int64)
    @error "An specialization of this method is expexted!"
end

"""
    fsolve!(::AbstractPhysicalModel, ::Float64, ::Int64, ::Float64)

Solution update for [`step!`](@ref).
"""
function fsolve!(::AbstractPhysicalModel, ::Float64, ::Int64, ::Float64)
    @error "An specialization of this method is expexted!"
end

"""
    timepoints(::AbstractPhysicalModel)

Get array of model time-points for use in [`step!`](@ref).
"""
function timepoints(::AbstractPhysicalModel)
    @error "An specialization of this method is expexted!"
end

function relaxationstep!(
    p::AbstractMatrixProblem,
    α::Float64,
    f::Function
    )::Float64
    "Applies relaxation to solution and compute residual."
    # Compute relaxed increment based on solution change.
    Δx = (1.0 - α) * change(p)

    # Evaluate residual metric.
    ε = f(p.x, Δx)

    # Update solution with relaxation.
    p.x[:] += Δx

    # Return residual.
    return ε
end