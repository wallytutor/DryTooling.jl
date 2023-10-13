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
    using Trapz
    using YAML
    using PlutoUI
    import DryTooling as dry
    include("parameters.jl")
end

# ╔═╡ f6b28d70-641d-11ee-01c7-2772f3bc7473
md"""
## Plug-flow model

$(TableOfContents())
"""

# ╔═╡ bb617963-b2b0-42ef-857a-cde4ae0387f5
mix = let
    data = YAML.load_file("mixtures.yaml")
    mix = dry.GasMixturePhase(data, order = ["fumes", "co2"])
end

# ╔═╡ 851467a3-8a11-43b3-9063-6c4d002407dc
W = dry.molecularmasses(mix)

# ╔═╡ 8782664f-2f03-472d-aac4-652d00581544
let
    T = dry.ZERO_CELSIUS
    P = dry.ONE_ATM
    Y = [1.0, 0.0]
    ρ = dry.idealgasdensity(T, P, Y; W = W)
end

# ╔═╡ 665a8c34-0f33-4094-888a-d02b08e2d01f
let
    K = 20
    T = collect(range(300.0, 3000.0, K))
    P = collect(range(300.0, 3000.0, K)) .+ dry.ONE_ATM
    Y = [range(0.0, 1.0, K) range(1.0, 0.0, K)]
    ρ = dry.idealgasdensity.(T, P, eachrow(Y); W = W)
end

# ╔═╡ 1a6c601e-7600-432d-8231-bf25f4621a55
let
    K = 20
    T = (dry.ZERO_CELSIUS + 25.0) * ones(K)
    P = dry.ONE_ATM * ones(K)
    Y = [range(0.0, 1.0, K) range(1.0, 0.0, K)]
    ρ, μ, k, c = zip(dry.mixtureproperties(T, P, Y; m = mix, W = W)...)
    ρ, μ, k, c
end

# ╔═╡ 433dcf12-39c1-4bd3-a3e3-2d3b0da6b11d
let
    K = 20
    T = collect(range(300.0, 2000.0, K))
    P = dry.ONE_ATM * ones(K)
    Y = [ones(Float64, K) zeros(Float64, K)]
    ρ, μ, k, c = zip(dry.mixtureproperties(T, P, Y; m = mix, W = W)...)
    ρ, μ, k, c
end

# ╔═╡ ef85aab3-7126-4c79-a2f8-dc996b639486
md"""
## Tools
"""

# ╔═╡ Cell order:
# ╟─f6b28d70-641d-11ee-01c7-2772f3bc7473
# ╠═bb617963-b2b0-42ef-857a-cde4ae0387f5
# ╠═851467a3-8a11-43b3-9063-6c4d002407dc
# ╠═8782664f-2f03-472d-aac4-652d00581544
# ╠═665a8c34-0f33-4094-888a-d02b08e2d01f
# ╠═1a6c601e-7600-432d-8231-bf25f4621a55
# ╠═433dcf12-39c1-4bd3-a3e3-2d3b0da6b11d
# ╟─ef85aab3-7126-4c79-a2f8-dc996b639486
# ╟─e37d9d68-9052-4e5e-9540-91e9e0ce6b27
