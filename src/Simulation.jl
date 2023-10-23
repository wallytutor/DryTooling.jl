# -*- coding: utf-8 -*-
module Simulation

using CairoMakie
using DryTooling

export TimeSteppingSimulationResiduals
export addresidual!
export plotsimulationresiduals

struct TimeSteppingSimulationResiduals
    """
        TimeSteppingSimulationResiduals
    
    Manage iterative solvers residuals storage during a simulation.
    The memory is initialized with a given number of inner and outer
    iterations and resizing is not under the scope of this structure.
    
    $(TYPEDFIELDS)
    """

    "Number of variables being tracked in problem."
    N::Int64

    "Total iteration counter."
    counter::Base.RefValue{Int64}

    "Number of inner steps per outer loop in solution."
    innersteps::Vector{Int64}

    "Store residuals of each inner step, one variable per column."
    residuals::Matrix{Float64}
end

function TimeSteppingSimulationResiduals(N::Int64, inner::Int64, outer::Int64)
    """
        TimeSteppingSimulationResiduals(N::Int64, inner::Int64, outer::Int64)
    
    Outer constructor for starting a simulation with pre-allocated memory.
    """
    innersteps = -ones(Int64, outer)
    residuals = -ones(Float64, (outer * inner, N))
    return TimeSteppingSimulationResiduals(N, Ref(0), innersteps, residuals)
end

function TimeSteppingSimulationResiduals(r::TimeSteppingSimulationResiduals)
    """
        TimeSteppingSimulationResiduals(r::TimeSteppingSimulationResiduals)
    
    Outer constructor for post-processing an already filled object.
    """
    # XXX: the equal sign below is required! In some cases, *e.g.*
    # when you try to simulate a system that is already at constant
    # state or when using a closed boundary condition, the error
    # can be exactly zero because of the form of the equations!
    innersteps = r.innersteps[r.innersteps .>= 0.0]
    residuals = r.residuals[r.residuals .>= 0.0]

    N = last(size(r.residuals))
    residuals = reshape(residuals, r.counter[], N)

    return TimeSteppingSimulationResiduals(N, r.counter, innersteps, residuals)
end

function addresidual!(r::TimeSteppingSimulationResiduals, ε::Vector{Float64})::Nothing
    """
        step!(r::TimeSteppingSimulationResiduals, ε::Vector{Float64})::Nothing

    Utility to increment iteration counter and store residuals.
    """
    # TODO: add resizing test here!
    r.counter[] += 1
    r.residuals[r.counter[], :] = ε
    return nothing
end

function finaliterationdata(
        r::TimeSteppingSimulationResiduals
    )::Tuple{Vector{Int64}, Matrix{Float64}}
    """
        finaliterationdata(
            r::TimeSteppingSimulationResiduals
        )::Tuple{Vector{Int64}, Matrix{Float64}}

    Retrieve data at iterations closing an outer loop of solution.
    """
    steps = cumsum(r.innersteps)
    residuals = r.residuals[steps, :]
    return steps, residuals
end

function plotsimulationresiduals(
        r::TimeSteppingSimulationResiduals;
        ε::Union{Float64, Nothing} = nothing,
        showinner::Bool = false,
        scaler::Function = log10,
        resolution::Tuple{Int64, Int64} = (720, 500)
    )::Tuple{Figure, Axis, Vector}
    """
        plotsimulationresiduals(
            r::TimeSteppingSimulationResiduals;
            ε::Union{Float64, Nothing} = nothing,
            showinner::Bool = false,
            resolution::Tuple{Int64, Int64} = (720, 500)
        )::Tuple{Figure, Axis, Vector}

    Plot problem residuals over iterations or steps. It performs the basic
    figure setup, configuration of axis and details beign left to the user.
    """
    xs, ys = finaliterationdata(r)
    ys = scaler.(ys)
    xs .-= 1

    units, unitp = axesunitscaler(last(xs))

    fig = Figure(resolution = resolution)
    ax = Axis(fig[1, 1], yscale = identity)
    ax.ylabel = "$(scaler)(Residual)"
    ax.xlabel = "Global iteration $(units)"

    p = []

    if showinner
        xg = collect(0:(r.counter[]-1))
        yg = scaler.(r.residuals)

        for col in 1:last(size(ys))
            lines!(ax, xg ./ unitp, yg[:, col], color = :black, linewidth = 0.5)
            push!(p, scatter!(ax, xs ./ unitp, ys[:, col]))
        end
    else
        for col in 1:last(size(ys))
            push!(p, lines!(ax, xs ./ unitp, ys[:, col]))
        end
    end

    if !isnothing(ε)
        hlines!(ax, scaler(ε), color = :blue, linewidth = 2)
    end

    return (fig, ax, p)
end

end # module Simulation