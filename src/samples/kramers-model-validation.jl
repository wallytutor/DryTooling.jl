### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 27971906-767d-4638-b7b1-bfb6aee9e7d8
begin
    import Pkg
    Pkg.activate(Base.current_project())
    Pkg.instantiate()

    using CSV
    using DataFrames
    using Ipopt
    using JuMP
    using Plots
    using PlutoUI
    using Printf
    import DryTooling as dry

    TableOfContents()
end

# ╔═╡ bbc70420-5554-11ee-035b-f7d156a056b2
md"""
# Kramers equation validation

## Problem specification

Bed profile $h(z)$ can be evaluated from *volume* conservation through

```math
\dfrac{dh}{dz} = C₁ f\left(\frac{h}{R}\right) - C₂
```

where the following partial expressions are used

```math
\begin{align}
C₁          &= \frac{3}{4}\dfrac{Φ\tan{γ}}{π R^3 ω}\\
C₂          &= \dfrac{\tan{β}}{\cos{γ}}\\
f(r)        &= (2r - r²)^{-\frac{3}{2}}
\end{align}
```

and coordinate $z=0$ represents the discharge.
"""

# ╔═╡ ce518a29-75b1-4e2a-9b9a-8ce73e07de51
md"""
## Comparison with analytical

The final step in model validation is to compare the approximate analytical solution proposed by Kramers and the results of numerical integration. It is worth mentioning that numerical integration remains the recommended method because one does not need to verify the ranges of validity of analytical approximation for every use case.
"""

# ╔═╡ 58b08548-efb9-41f6-8b8f-b2d7e4e1bcb6
md"""
## Industrial cases
"""

# ╔═╡ 27cacb32-fb32-4a8b-a293-ea02296a0054
md"""
### Alumina kiln

The following illustrates a practical use case of the model. Next we scan a parameter space to confirm once again the model suitability as an alternative to analytical engineering estimations as per Peray's notebook.
"""

# ╔═╡ 0d184dcf-3592-40cc-9cd3-743bebcbf184
md"""
The following table confirms the expected values as per Peray.
"""

# ╔═╡ 0cafec10-ed46-49fe-9388-dc9ab84411fc
md"""
## Appendix
"""

# ╔═╡ 436f5f65-065c-4b63-b138-95d5a999b870
md"""
### Data and functions
"""

# ╔═╡ 39f66a95-0a70-43ef-957a-e3fe7f07529a
"Partial data from Kramers (1952) Table 3"
const DATA_TABLE3 = """\
ρ,γ,tan(β),n,ṁ,prod_dimless,η̄ᵣ,hold_real
1480.0,36.0,0.0094,0.059,5.15e-03,18.3,0.111,8.10
1480.0,36.0,0.0094,0.090,2.68e-03,6.25,0.054,5.00
1480.0,36.0,0.0094,0.195,1.32e-02,14.2,0.088,7.75
1480.0,36.0,0.0094,0.232,7.24e-03,6.55,0.043,3.85
1480.0,36.0,0.0100,0.040,6.38e-03,29.7,0.169,13.3
1480.0,36.0,0.0100,0.040,5.00e-03,23.2,0.144,11.2
1480.0,36.0,0.0100,0.069,9.20e-03,24.8,0.150,10.6
1480.0,36.0,0.0100,0.069,6.53e-03,17.6,0.113,8.50
1480.0,36.0,0.0100,0.106,1.50e-02,27.8,0.162,12.2
1480.0,36.0,0.0100,0.159,1.20e-02,14.0,0.092,7.49
1480.0,36.0,0.0100,0.238,1.55e-02,12.1,0.083,7.48
1480.0,36.0,0.0100,0.238,1.19e-02,9.22,0.068,6.13
"""

# ╔═╡ c5d1f88d-90ab-4d3e-b878-4c281081494f
"Kramers (1952) dimensionless group NΦ."
function dimlessNΦ(R, β, ω, Φ, γ)
    return Φ * sin(γ) / (ω * R^3 * tan(β))
end

# ╔═╡ 7a839722-9dee-4cb5-b900-042140cf4eeb
"Kramers (1952) dimensionless group Nₖ."
function dimlessNₖ(L, R, β, γ)
    return R * cos(γ) / (L * tan(β))
end

# ╔═╡ cb33d18a-e11e-4596-86a3-1d1d7389f982
"Sullivans approximation to kiln filling."
function sullivansηₘ(R, β, ω, Φ, γ)
    return 3.8 * dimlessNΦ(R, β, ω, Φ, γ) * sqrt(γ) / sin(γ)
end

# ╔═╡ 049056cf-f754-4fdc-ab20-a6b177942c11
"Nonlinear formulation of Kramers model approximate solution."
function kramersanalytical(; z, R, Φ, ω, β, γ, d)
    L = z[end]
    N = length(z)

    NΦ = dimlessNΦ(R, β, ω, Φ, γ)
    Nₖ = dimlessNₖ(L, R, β, γ)

    C₁ = R * NΦ
    C₂ = 3C₁ / (4π * 1.24)
    C₃ = C₁ / (L * NΦ * Nₖ)

    optim = JuMP.Model(Ipopt.Optimizer)
    JuMP.set_silent(optim)

    @JuMP.variable(optim, h[1:N])
    @JuMP.NLconstraint(
        optim,
        [i = 1:N],
        C₂ * log((d - C₂) / (h[i] - C₂)) - C₃ * z[i] - h[i] + d == 0,
    )
    @JuMP.NLconstraint(optim, [i = 1:N], h[i] >= 0.0)
    JuMP.optimize!(optim)

    return dry.RotaryKilnBedSolution(z, JuMP.value.(h), β, R, Φ)
end

# ╔═╡ cd18eabc-3c5f-40ec-bb84-523ef0577871
"Compute residence time from Peray's equation."
function perrayresidence(L, ω, D, β)
    return 0.19 * L / (ω * D * tan(β))
end

# ╔═╡ 75228bbe-4da9-4f70-a065-c96f47adc704
"Compares approximate analytical to numerical solution."
function solvekiln(; L, D, Φ, ω, β, γ, d, show = true)
    model = dry.solvelinearkramersmodel(;
        model = dry.SymbolicLinearKramersModel(),
        L     = L,
        R     = D / 2.0,
        Φ     = Φ / 3600.0,
        ω     = ω / 60.0,
        β     = deg2rad(β),
        γ     = deg2rad(γ),
        d     = d / 1000.0
    )

    optim = kramersanalytical(;
        z = model.z,
        R = D / 2.0,
        Φ = Φ / 3600.0,
        ω = ω / 60.0,
        β = deg2rad(β),
        γ = deg2rad(γ),
        d = d / 1000.0
    )

    p = nothing;

    if show
        p = plot()
        plot!(p, 100model.z/L, 100model.h,
              linewidth = 3, label = "Numerical")
        plot!(p, 100optim.z/L, 100optim.h,
              linewidth = 3, label = "Analytical")

        a = @sprintf("%.1f", model.ηₘ)
        b = @sprintf("%.1f", optim.ηₘ)
        title = "Loading: $(a)% (numerical) | $(b)% (analytical)"

        plot!(p,
              title = title,
              xaxis = "Coordinate [%]",
              yaxis = "Bed height [cm]",
              label = nothing,
              xlims = (0.0, 100.0),
              xticks = 0.0:20.0:100.0,
        )
    end

    return model, optim, p
end

# ╔═╡ b258541d-fca9-4c2e-b6c1-d602d7d077ed
let
    _, _, p = solvekiln(
        L = 10.0,
        D = 1.0,
        Φ = 1.0,
        ω = 1.0,
        β = 3.0,
        γ = 45.0,
        d = 0.001
    )

    plot!(p,
          ylims  = (0.0, 20.0),
          yticks = 0.0:4.0:20.0)
end

# ╔═╡ 3a855976-5150-47f6-b3ef-cf1ba4a11239
"Reference case for alumina kiln testing."
function aluminakiln(ṁ, ω; show = false)
    # Density of bed [kg/m³]
    ρ = 800.0
    L = 34.0
    D = 1.5
    β = atan(0.025)

    model, optim, p = solvekiln(
        L = L,
        D = D,
        Φ = (1000// 24) * ṁ / ρ,
        ω = ω,
        β = rad2deg(β),
        γ = 33.0,
        d = 0.050,
        show = show
    )

    τₚ = perrayresidence(L, ω, D, β)

    return model, optim, p, τₚ
end

# ╔═╡ e2d3263d-bec5-4c6a-83da-e6e6f4deaeab
let
    ṁ = 33.6
    ω = 0.85
    _, _, p, _ = aluminakiln(ṁ, ω, show = true)

    plot!(p,
          ylims  = (0.0, 30.0),
          yticks = 0.0:6.0:30.0)
end

# ╔═╡ c09fb848-1e24-43a1-be49-07deeb3e76cc
"Run `aluminakiln` against some known conditions."
function scanaluminakiln()
    ṁlist = [33.6, 43.2]
    ωlist = [0.85, 1.20]

    df = DataFrame(
        ṁ = Float64[],
        ω = Float64[],
        η̄ = Float64[],
        τᵢ = Float64[],
        τₚ = Float64[]
    )

    for ṁ ∈ ṁlist, ω ∈ ωlist
        model, _, _, τ = aluminakiln(ṁ, ω, show = false)



        η̄ = round(model.ηₘ, digits = 0)
        τᵢ = round(model.τ / 60.0, digits = 0)
        τₚ = round(τ, digits = 0)
        push!(df, [ṁ ω η̄ τᵢ τₚ])
    end

    return df
end

# ╔═╡ fe157b40-e2d7-404c-87ad-ad3f2f4a6004
scanaluminakiln()

# ╔═╡ 8d01b94b-f443-41e2-8f55-04253d737ba0
md"""
### Problems setup
"""

# ╔═╡ 9b77c70d-06d2-485d-bbc4-087ed0ffaad1
let
    println("Solution of reference case")

    in1_to_m1(v) = 0.0254 * v
    ft1_to_m1(v) = in1_to_m1(12.0) * v
    ft3_to_m3(v) = ft1_to_m1(1.0)^3 * v

    # Kiln length [m]
    L = ft1_to_m1(45.0)

    # Kiln diameter [m]
    D = 2 * ft1_to_m1(3.1)

    # Volume flow rate [m³/h]
    Φ = ft3_to_m3(6.1) * 60

    # Rotation rate (+0.0005) [rev/min]
    ω = 0.0505 * 60.0

    # Kiln slope (0.5in/ft) [°]
    β = rad2deg(atan(0.5 / 12))

    # Repose angle [°]
    γ = 45.0

    # Particle size [mm]
    d = 0.050

    # Conversions to match model inputs.
    R = D / 2.0
    Φ = Φ / 3600.0
    ω = ω / 60.0
    β = deg2rad(β)
    γ = deg2rad(γ)
    d = d / 1000.0

    # Create problem container.
    kramers = dry.solvelinearkramersmodel(;
        model = dry.SymbolicLinearKramersModel(),
        L     = L,
        R     = R,
        Φ     = Φ,
        ω     = ω,
        β     = β,
        γ     = γ,
        d     = d
    )

    optim = kramersanalytical(;
        z = kramers.z,
        R = R,
        Φ = Φ,
        ω = ω,
        β = β,
        γ = γ,
        d = d
    )

    global kramers_NΦ = dimlessNΦ(R, β, ω, Φ, γ)
    global kramers_Nₖ = dimlessNₖ(L, R, β, γ)
    global kramers_η̄ₛ = sullivansηₘ(R, β, ω, Φ, γ)
    global kramers_ref = kramers
    global optim_ref = optim
end;

# ╔═╡ e8e2a667-fc8d-47be-b8ae-9d71c7c3c036
md"""
## Sample reference case

Here we make use of the current implementation to check if it correctly approximates the last example provided in reference paper from [Kramers (1952)](https://doi.org/10.1016/0009-2509(52)87019-8). To minimize rounding errors causes by unit conversions, we provide the required functions to convert from imperial to international system in the solution process.

The next table summarizes the results. It is seen that the dimensionless numbers are well approximated. It must be emphasized that the reference estimates η̄ᵣ by a graphical method -- it was 1952 -- and the current value is considered a good enough approximation. Additionally, the equation was not integrated numerically as done here, but engineering relationships were used in the approximation. That said, the proper loading to be considered in our days is η̄ᵢ, what gives $(@sprintf("%.2f", 100 * (kramers_ref.ηₘ - 5.65) / 5.65))% more loading than the conventional estimation.

| Quantity | Reference | Computed |
| -------- | :-------: | :------: |
| NΦ       | 1.15      | $(@sprintf("%.2f", kramers_NΦ)) |
| Nₖ       | 1.17      | $(@sprintf("%.2f", kramers_Nₖ)) |
| η̄ᵣ       | 5.65      | $(@sprintf("%.2f", kramers_η̄ₛ)) |
| η̄ᵢ       | $(@sprintf("%.2f", optim_ref.ηₘ)) [^1] | $(@sprintf("%.2f", kramers_ref.ηₘ)) |

[^1]: This is not provided in Kramers' paper but computed from the approximate analytical solution provided by the authors. As we see here, it may get >20% error under some circumstances.
"""

# ╔═╡ f8e6391a-c942-4b6f-9359-e36c196a63c0
const TABLE3 = let
    println("Verification of *Table 3*")

    Dₖ = 0.197
    Lₖ = 1.780
    dₖ = 0.0012

    table3 = DataFrame(CSV.File(IOBuffer(DATA_TABLE3)))
    table3[!, "η̄ᵢ"] = zeros(length(table3[!, "η̄ᵣ"]))
    table3[!, "η̄ᵣ"] *= 100

    model = dry.SymbolicLinearKramersModel()

    for (i, row) in enumerate(eachrow(table3))
        Φ = 3600.0 * row["ṁ"] / row["ρ"]
        ω = row["n"] * 60.0
        β = rad2deg(atan(row["tan(β)"]))
        γ = row["γ"]

        kramers = dry.solvelinearkramersmodel(;
            model = model,
            L     = Lₖ,
            R     = Dₖ / 2.0,
            Φ     = Φ / 3600.0,
            ω     = ω / 60.0,
            β     = deg2rad(β),
            γ     = deg2rad(γ),
            d     = dₖ / 1000.0
        )

        table3[i, "η̄ᵢ"] = round(kramers.ηₘ, digits = 1)
    end

    exclude = ["ρ", "γ", "prod_dimless", "hold_real"]
    select(table3, Not(exclude))
end;

# ╔═╡ f9c3f3ee-fbe5-4ff0-a7c6-614f21e34c76
md"""
## Verification of *Table 3*

In the next cell we provide the kiln dimensions used by Kramers (1952) to experimentally validate the model. Some data from their Tab. 3 is then loaded and all rows are simulated with current model. Fractional hold-up seems to be well correlated at least to a few percent of the reference value.

$(TABLE3)
"""

# ╔═╡ 15143f31-e958-466a-bb1d-31a68f52351e
const DIMLESSPLOT = let
    println("Dimensionless profiles solution")

    ρ = 1480.0
    L = 20.0
    D = 0.197
    Φ = 5.15e-03 / ρ * 3600
    ω = 0.059 * 60
    β = rad2deg(atan(0.0094))
    γ = 36.0

    # Conversions to match model inputs.
    R = D / 2.0
    Φ = Φ / 3600.0
    ω = ω / 60.0
    β = deg2rad(β)
    γ = deg2rad(γ)

    # Things held constant in loop.
    NΦ = dimlessNΦ(R, β, ω, Φ, γ)
    Nₖ = dimlessNₖ(L, R, β, γ)
    model = dry.SymbolicLinearKramersModel()

    p = plot()

    for d in [0.05, 0.10, 0.15, 0.193, 0.25]
        kramers = dry.solvelinearkramersmodel(;
            model = model,
            L     = L,
            R     = R,
            Φ     = Φ,
            ω     = ω,
            β     = β,
            γ     = γ,
            d     = d * R * NΦ
        )

        # Dimensionless axes.
        z = kramers.z
        h = kramers.h / (R * NΦ)
        z = @. (L - z) / L * 1 / (NΦ * Nₖ)
        z = @. z[1] - z

        plot!(p, z, h,
              linewidth = 2,
              label = @sprintf("%.3f", d))
    end

    plot!(p,
         title = "Dimensionless loading curves",
         xaxis = "Coordinate",
         yaxis = "Bed height",
         label = nothing,
         xlims = (0.0, 0.5),
         ylims = (0.0, 0.25),
         xticks = 0.0:0.1:0.5,
         yticks = 0.0:0.05:0.25
    )

    p
end;

# ╔═╡ 0ff54742-abf4-4b28-8ce9-ee8572c85285
md"""
## Dimensionless profiles

Next step in validation is to check profiles in dimensionless format, as done by Kramers in their Fig. 3. Notice that here we used the numerical integration curves instead of the analytical approximation of profiles, so reproduction and consequences of results are not exactly the same.

$(DIMLESSPLOT)
"""

# ╔═╡ Cell order:
# ╟─bbc70420-5554-11ee-035b-f7d156a056b2
# ╟─e8e2a667-fc8d-47be-b8ae-9d71c7c3c036
# ╟─f9c3f3ee-fbe5-4ff0-a7c6-614f21e34c76
# ╟─0ff54742-abf4-4b28-8ce9-ee8572c85285
# ╟─ce518a29-75b1-4e2a-9b9a-8ce73e07de51
# ╟─b258541d-fca9-4c2e-b6c1-d602d7d077ed
# ╟─58b08548-efb9-41f6-8b8f-b2d7e4e1bcb6
# ╟─27cacb32-fb32-4a8b-a293-ea02296a0054
# ╟─e2d3263d-bec5-4c6a-83da-e6e6f4deaeab
# ╟─0d184dcf-3592-40cc-9cd3-743bebcbf184
# ╟─fe157b40-e2d7-404c-87ad-ad3f2f4a6004
# ╟─0cafec10-ed46-49fe-9388-dc9ab84411fc
# ╟─436f5f65-065c-4b63-b138-95d5a999b870
# ╟─39f66a95-0a70-43ef-957a-e3fe7f07529a
# ╟─c5d1f88d-90ab-4d3e-b878-4c281081494f
# ╟─7a839722-9dee-4cb5-b900-042140cf4eeb
# ╟─cb33d18a-e11e-4596-86a3-1d1d7389f982
# ╟─049056cf-f754-4fdc-ab20-a6b177942c11
# ╟─cd18eabc-3c5f-40ec-bb84-523ef0577871
# ╟─75228bbe-4da9-4f70-a065-c96f47adc704
# ╟─3a855976-5150-47f6-b3ef-cf1ba4a11239
# ╟─c09fb848-1e24-43a1-be49-07deeb3e76cc
# ╟─8d01b94b-f443-41e2-8f55-04253d737ba0
# ╟─9b77c70d-06d2-485d-bbc4-087ed0ffaad1
# ╟─f8e6391a-c942-4b6f-9359-e36c196a63c0
# ╟─15143f31-e958-466a-bb1d-31a68f52351e
# ╟─27971906-767d-4638-b7b1-bfb6aee9e7d8
