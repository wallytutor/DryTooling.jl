### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ e37d9d68-9052-4e5e-9540-91e9e0ce6b27
begin
    import Pkg
    Pkg.activate(Base.current_project())
    Pkg.resolve()
    Pkg.instantiate()

    using CairoMakie
    using DryTooling
    using Trapz
    using YAML
    using PlutoUI
    include("parameters.jl")
end

# ╔═╡ f6b28d70-641d-11ee-01c7-2772f3bc7473
md"""
## Plug-flow model

$(TableOfContents())
"""

# ╔═╡ 433dcf12-39c1-4bd3-a3e3-2d3b0da6b11d
data = YAML.load_file("mixtures.yaml")
# mix = GasMixture(data, order = ["fumes", "co2"])

# W = molecularmasses(mix)

# T = ZEROCELSIUS
# P = ONEATM
# Y = [1.0, 0.0]
# ρ = idealgasdensity(T, P, Y; W = W)

# K = 20
# T = collect(range(300.0, 3000.0, K))
# P = collect(range(300.0, 3000.0, K)) .+ ONEATM
# Y = [range(0.0, 1.0, K) range(1.0, 0.0, K)]
# ρ = idealgasdensity.(T, P, eachrow(Y); W = W)

# K = 20
# T = (ZEROCELSIUS + 25.0) * ones(K)
# P = ONEATM * ones(K)
# Y = [range(0.0, 1.0, K) range(1.0, 0.0, K)]
# ρ, μ, k, c = zip(mixtureproperties(T, P, Y; m = mix, W = W)...)

# K = 20
# T = collect(range(300.0, 2000.0, K))
# P = ONEATM * ones(K)
# Y = [ones(Float64, K) zeros(Float64, K)]
# ρ, μ, k, c = zip(mixtureproperties(T, P, Y; m = mix, W = W)...)

# ╔═╡ ef85aab3-7126-4c79-a2f8-dc996b639486
md"""
## Tools
"""

# ╔═╡ Cell order:
# ╟─f6b28d70-641d-11ee-01c7-2772f3bc7473
# ╠═433dcf12-39c1-4bd3-a3e3-2d3b0da6b11d
# ╟─ef85aab3-7126-4c79-a2f8-dc996b639486
# ╠═e37d9d68-9052-4e5e-9540-91e9e0ce6b27
