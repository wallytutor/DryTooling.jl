# -*- coding: utf-8 -*-

"Manage iterative solvers residuals storage during a simulation."
mutable struct ResidualsRaw
    "Inner loop iteration counter."
    inner::Int64

    "Outer loop iteration counter."
    outer::Int64

    "Total iteration counter."
    counter::Int64

    "Number of inner steps in solution."
    innersteps::Vector{Int64}
    
    "Store residuals of each inner step."
    residuals::Vector{Float64}
    
    function ResidualsRaw(inner::Int64, outer::Int64)
        innersteps = -ones(Int64, outer)
        residuals = -ones(Float64, outer * inner)
        return new(inner, outer, 0, innersteps, residuals)
    end
end

"Post-processing of `ResidualsRaw` for simulation residuals."
struct ResidualsProcessed
    "Total iteration counter."
    counter::Int64
    
    "Number of inner steps in solution."
    innersteps::Vector{Int64}
    
    "Store residuals of each inner step."
    residuals::Vector{Float64}
    
    "Step number of inner iteration exit."
    finalsteps::Vector{Int64}

    "Residuals at exit iteration in `finalsteps`."
    finalresiduals::Vector{Float64}

    function ResidualsProcessed(r::ResidualsRaw)
        innersteps = r.innersteps[r.innersteps .> 0.0]
        residuals = r.residuals[r.residuals .> 0.0]

        finalsteps = cumsum(innersteps)
        finalresiduals = residuals[finalsteps]

        return new(r.counter, innersteps, residuals,
                   finalsteps, finalresiduals)
    end
end

"Feed residuals in an inner relaxation loop."
function feedinnerresidual(r::ResidualsRaw, ε::Float64)
    # TODO: add resizing test here!
    r.counter += 1
    r.residuals[r.counter] = ε
end

"Plot problem residuals over iterations or steps."
function plotsimulationresiduals(
        r::ResidualsProcessed;
        ε::Union{Float64, Nothing} = nothing,
        showinner::Bool = false,
        yticks::Any = nothing,
        xbase::Number = 20
    )::Figure
    function getxticks(xv)
        δi = closestpowerofx(xv/10; x = xbase)
        imax = closestpowerofx(xv; x = xbase)
        return 0:δi:imax
    end

    xs = r.finalsteps
    ys = log10.(r.finalresiduals)

    fig = CM.Figure(resolution = (720, 500))
    ax = CM.Axis(fig[1, 1], yscale = identity)
    ax.ylabel = "log10(Residual)"
    ax.title = "Maximum of $(maximum(r.innersteps)) iterations per step"

    if !isnothing(yticks)
        ax.yticks = yticks
        CM.ylims!(ax, extrema(yticks))
    end

    if showinner
        xg = 1:r.counter
        yg = log10.(r.residuals)
        xticks = getxticks(xg[end])
        ax.xlabel = "Global iteration"
        CM.lines!(ax, xg, yg, color = :black, linewidth = 0.5)
        CM.scatter!(ax, xs, ys, color = :red)
    else
        xticks = getxticks(length(xs))
        ax.xlabel = "Outer iteration"
        CM.lines!(ax, xs, ys, color = :red)
    end

    ax.xticks = xticks
    CM.xlims!(ax, extrema(xticks))

    if !isnothing(ε)
        CM.hlines!(ax, log10(ε), color = :blue, linewidth = 2)
    end

    return fig
end