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

    if relax >= 1.0 || relax < 0.0
        @error """
        Relaxation factor out-of-range (α = $(relax)), this is
        not supported / does not make physical sense in most cases!
        """
    end

    if relax <= 1.0e-06
        @warn """
        Relaxation factor below threshold (α = $(relax)). If you are
        this low, that probably means that the model is not nonlinear.
        In this cases better performance (solution times) would be
        achieved with a direct solver instead, although this solver
        will manage it properly.
        """
    end

    times = 0.0:tau:tend
    residual = ResidualsRaw(iters, length(times))

    for (nouter, ts) in enumerate(times)
        updaterouter(model, ts, nouter)
        residual.innersteps[nouter] = relaxationinnerloop(;
            model    = model,
            updater  = (m)->updaterinner(m, ts, nouter),
            residual = residual,
            iters    = iters,
            relax    = relax,
            tol      = tol,
            metric   = metric,
        )
    end

    updaterouter(model, model.t, length(model.Q))

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
        p = model.problem

        if relax <= 0.0
            p.x[:] = p.A \ p.b
            ε = tol
        else
            ε = relaxationstep(p, relax, metric)
        end
        
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
    ε = metric(p.x, Δx)
    p.x[:] += Δx
    return ε
end
