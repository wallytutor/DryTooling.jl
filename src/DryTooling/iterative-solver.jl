# -*- coding: utf-8 -*-

"Maximum relative change in a solution array."
function maxrelativevariation(
        x::Vector{Float64},
        Δx::Vector{Float64}
    )::Float64
    return maximum(abs.(Δx ./ x))
end

"Integrates and arbitrary problem stepping over non-linear relaxations."
function relaxationouterloop(;
        model::AbstractPhysicalModel,
        updaterouter::Function,
        updaterinner::Function,
        tend::Float64,
        tau::Float64,
        iters::Int64 = 10,
        relax::Float64 = 0.5,
        tol::Float64 = 1.0e-08,
        metric::Function = maxrelativevariation
    )::ResidualsProcessed
    times = 0.0:tau:tend
    residual = ResidualsRaw(iters, length(times))

    # TODO use `ts` for time dependent *things*!
    for (nouter, ts) in enumerate(times)
        updaterouter(model, nouter, ts)

        residual.innersteps[nouter] = relaxationinnerloop(;
            model    = model,
            updater  = updaterinner,
            residual = residual,
            iters    = iters,
            relax    = relax,
            tol      = tol,
            metric   = metric,
        )
    end

    updaterouter(model, length(model.Q), model.t)

    return ResidualsProcessed(residual)
end

"Solve an arbitrary problem through successive relaxations."
function relaxationinnerloop(;
        model::AbstractPhysicalModel,
        updater::Function,
        residual::ResidualsRaw,
        iters::Int64 = 10,
        relax::Float64 = 0.5,
        tol::Float64 = 1.0e-08,
        metric::Function = maxrelativevariation
    )::Int64
    for niter in 1:iters
        updater(model)
        ε = relaxationstep(model.problem, relax, metric)
        feedinnerresidual(residual, ε)
        if ε <= tol
            return niter
        end
    end

    @warn "Did not converge after $(iters) steps"
    return iters
end

"Applies relaxation to solution and compute residual."
function relaxationstep(
        p::AbstractMatrixProblem,
        relax::Float64,
        metric::Function
    )::Float64
    Δx = (1.0 - relax) * (p.A \ p.b - p.x)
    p.x[:] += Δx
    return metric(p.x, Δx)
end
