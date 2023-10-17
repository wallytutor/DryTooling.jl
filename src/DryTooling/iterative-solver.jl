# -*- coding: utf-8 -*-

"Maximum relative change in a solution array."
function maxrelativevariation(
        x::Vector{Float64},
        Δx::Vector{Float64}
    )::Float64
    return maximum(abs.(Δx ./ x))
end

"Maximum change in a solution array."
function maxvariation(
        x::Vector{Float64},
        Δx::Vector{Float64}
    )::Float64
    return maximum(abs.(Δx))
end

"Integrates and arbitrary problem stepping over non-linear relaxations."
function relaxationouterloop(;
        model::AbstractPhysicalModel,
        updaterouter::Function,
        updaterinner::Function,
        tend::Float64,
        steps::Int64,
        iters::Int64 = 10,
        relax::Float64 = 0.5,
        tol::Float64 = 1.0e-08,
        metric::Function = maxvariation
    )::ResidualsProcessed

    if relax >= 1.0 || relax < 0.0
        @error """\
        Relaxation factor out-of-range (α = $(relax)), this is \
        not supported / does not make physical sense in most cases!
        """
    end

    if relax <= 1.0e-06
        @warn """\
        Relaxation factor below threshold (α = $(relax)). If you are this \
        low, that probably means that the model is not nonlinear. In \
        this cases better performance in terms of solution times is \
        possible with a direct solver instead, although the present \
        will handle the problem it properly.
        """
    end

    times = range(0.0, tend, steps)
    residual = ResidualsRaw(iters, steps)

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

    updaterouter(model, tend, steps+1)

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
