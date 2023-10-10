### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 52d557a4-8ef1-4b50-a145-7995745502c6
begin
    import Pkg
    Pkg.activate(Base.current_project())
    Pkg.resolve()
    Pkg.instantiate()

    using Symbolics
    using PlutoUI
    include("parameters.jl")
end

# ╔═╡ f175ae20-641b-11ee-2bf8-e72783c8a834
md"""
# Porosity calculations

In what follows we develop a simple approach for estimating quantities required for the computation of quantities relevant to heat and momentum transfer in the counter-current porous reactor model. Expressions are manipulated symbolically and then numerical results are obtained.

$(TableOfContents())
"""

# ╔═╡ a677c786-874c-4cc7-9871-6016064797a3
md"""
## Geometrical considerations

We start by declaring the relevant symbols:

|     |    |
|:---:|:---|
``ϕ`` | Volumetric porosity fraction
``l`` | Characteristic size of blocks
``D`` | Depth of reactor cross-section
``W`` | Width of reactor cross-section
``H`` | Height vertical reactor
"""

# ╔═╡ 5ef67a89-2d81-4863-b16f-a0f0a3c6321f
args = @variables ϕ l D W H

# ╔═╡ 13daafac-abfa-476d-b845-a78cc864559e
md"""
The surface area ``Aᵦ`` and volume ``Vᵦ`` of a single cubic block is computed as:
"""

# ╔═╡ e48a03d9-7189-4ff0-95cf-6bb942cab369
Aᵦ = 6l^2

# ╔═╡ 915f13a0-c078-47e8-ba4e-2eb0ab8d7de8
Vᵦ = l^3

# ╔═╡ 1eabfb91-5d33-4f4a-a9b4-6f4def62a1d6
md"""
The cross section of the reactor is ``Aᵥ``, total volume ``Vᵥ``, and perimeter ``Pᵥ`` are:
"""

# ╔═╡ 4dfcfd32-3b16-48b1-81e4-4371bb11a731
Aᵥ = D * W

# ╔═╡ bd7ff6f9-1578-4780-8666-6ca6cfa0ede8
Vᵥ = H * Aᵥ

# ╔═╡ 1d886b01-f1c5-4d85-ad61-970ed562c0e8
Pᵥ = 2 * (D + W)

# ╔═╡ a68a22cc-8a2f-4bee-82d7-14e1182381ca
md"""
From porosity ``ϕ`` the total volume of solids inside the reactor is ``Vₛ`` is:
"""

# ╔═╡ d0e39584-6c13-4804-8adc-5b5eb962940f
Vₛ = (1 - ϕ) * Vᵥ

# ╔═╡ 9876504d-9732-47c1-9ea5-003e6a916c8f
md"""
The number of blocks ``Nₛ`` is estimated as the ratio between ``Vₛ`` and a block volume ``Vᵦ``
"""

# ╔═╡ a87f0302-8e77-4dea-9e10-0469ecfdf6e3
Nₛ = Vₛ / Vᵦ

# ╔═╡ 19ae6a86-684d-4dca-b75f-27b9bfe71b9c
md"""
Thus the surface area of gas-blocks can be estimated as ``Aₛ``:
"""

# ╔═╡ d003c1fb-2c61-4031-b7bc-9c11d6e2dcc2
Aₛ = Nₛ * Aᵦ

# ╔═╡ 5ff1e62f-2a34-4dee-b2de-12fa670212c1
md"""
The blocks perimeter *per unit height* required in the PFR model is simply
"""

# ╔═╡ 48d89f5c-80e7-4c79-869a-c9dcd87b5c77
Aₛ / H

# ╔═╡ b6efc1d0-f5cc-42ac-b46b-920f6ce6e72e
md"""
Because reactor cross-section is constant, we simplify the expression per unit area:
"""

# ╔═╡ 2e4ce97e-2ec1-4691-a0f2-451915f9fd2f
Pₛ = Aₛ / Vᵥ

# ╔═╡ 8d78856b-0a9c-4909-bb30-fc530197b9ca
md"""
Actually this is exactly the same result as found by Gunn (1978), thus we express their expression in terms of ``Pₛ`` to ensure a proper estimation of the channel radius which is found to be:
"""

# ╔═╡ 66de140a-3a62-4b74-9b2b-97efaff2f3a5
Dᵧ = 2 * 2ϕ / Pₛ

# ╔═╡ e118a6e8-9b69-4987-801e-f2035687c41f
md"""
The factor ``2`` in the above equation reflects the use of a diameter instead of a radius.
"""

# ╔═╡ 5fcd5364-5d57-4711-9c63-8be40607aa56
md"""
## Statistical considerations
"""

# ╔═╡ 738db443-ad85-4f10-8864-9ba9b653fae6
md"""
New symbols are introduced to compute the uncertainties on parameters:
"""

# ╔═╡ 4983a5d8-1129-423e-8427-320edb206956
delta = @variables δϕ δl

# ╔═╡ 7206a2a0-4c77-4861-988c-49249966d047
md"""
The associated differentials are declared to perform a sensitivity analysis:
"""

# ╔═╡ c758b980-9aba-41dc-a696-4fa0ed8cee9e
ddϕ, ddl = let
    ddϕ = Differential(ϕ)
    ddl = Differential(l)
    ddϕ, ddl
end

# ╔═╡ d33a74c4-bb51-44ba-b61a-3bddadb30730
md"""
In what follows, this section provides functions for the evaluation of model parameters uncertainty.
"""

# ╔═╡ c1661200-1e5c-45b2-984e-56c91100c220
md"""
### Solids perimeter uncertainty

The total derivative of perimeter with respect to ``ϕ`` and ``l`` is given by:
"""

# ╔═╡ f76010d5-491a-4e38-b111-daa78f15f34a
dPₛ = expand_derivatives(ddl(Pₛ)*δl + ddϕ(Pₛ)*δϕ)

# ╔═╡ 944c1d88-edf7-4bba-8cb6-0490347b9f88
md"""
### Channel diameter uncertainty
"""

# ╔═╡ cf128a15-dae6-44fc-a28d-c861a57fef91
md"""
The same approach is now applied to channel diameter:
"""

# ╔═╡ 88769a22-d80b-4f3e-88cd-519400c47de9

dDᵧ = expand_derivatives(ddl(Dᵧ)*δl + ddϕ(Dᵧ)*δϕ)

# ╔═╡ 5cb8cc6e-1097-47a9-a0c0-ec93fbf711ef
md"""
## Numerical evaluation
"""

# ╔═╡ f7b401c6-a28d-4a3c-9dae-0cf471d7b162
md"""

|    |    |
|---:|:---|
 Voidage relative error [%]    | $(@bind δϕn PlutoUI.Slider(0.0:1.0:40.0, default = 10.0, show_value = true))
 Block size relative error [%] | $(@bind δln PlutoUI.Slider(0.0:1.0:40.0, default = 10.0, show_value = true))

"""

# ╔═╡ eb8e27b4-a2cd-422b-995e-3465478c3f70
porosityparameters(ϕₛ, blocksize, A = REACTOR.A)

# ╔═╡ 0948d221-cbac-4f6e-ae9f-c6a3ff3f0705
md"""
## Functions
"""

# ╔═╡ 078c50a8-4efa-463d-9110-61a6a38431bd
begin
    funPₛ = build_function(Pₛ, args...; expression = Val{false})

    "Evaluates the perimeter of solids per unit height of reactor."
    function solidsperim(; ϕ, l)
        return funPₛ(ϕ, l, REACTOR.D, REACTOR.W, REACTOR.H)
    end
end

# ╔═╡ 65c99f9b-0322-4229-a108-ead8525c4091
begin
    funDᵧ = build_function(Dᵧ, args...; expression = Val{false})

    "Evaluates the channel diameter estimation in cross-section."
    function channeldiam(; ϕ, l)
        return funDᵧ(ϕ, l, REACTOR.D, REACTOR.W, REACTOR.H)
    end
end

# ╔═╡ 0d81315f-9962-450e-b42a-af7525754408
begin
    fundPₛ = build_function(dPₛ, args..., delta...; expression = Val{false})

    "Evaluate blocks perimeter uncertainty from sizes and porosity."
    function Δsolidsperim(; ϕ, l, δϕ, δl)
        return fundPₛ(ϕ, l, REACTOR.D, REACTOR.W, REACTOR.H, δϕ, δl)
    end
end

# ╔═╡ c7eaa50c-0aaf-49c0-bd71-69246cf42bff
begin
    fundDᵧ = build_function(dDᵧ, args..., delta...; expression = Val{false})

    "Evaluate blocks perimeter uncertainty from sizes and porosity."
    function Δchanneldiam(; ϕ, l, δϕ, δl)
        return fundDᵧ(ϕ, l, REACTOR.D, REACTOR.W, REACTOR.H, δϕ, δl)
    end
end

# ╔═╡ 9b99bd14-8b19-4069-bac9-ebe147b6b515
μP, μD, σP, σD = let
    ϕ = ϕₛ
    l = blocksize
    δϕ = 0.5ϕ * δϕn / 100.0
    δl = 0.5l * δln / 100.0

    μP = solidsperim(; ϕ, l)
    μD = channeldiam(; ϕ, l)

    σP = abs(Δsolidsperim(; ϕ, l, δϕ, δl = δl) / 3)
    σD = abs(Δchanneldiam(; ϕ, l, δϕ, δl = δl) / 3)

    μP *= REACTOR.A
    σP *= REACTOR.A

    Prng = sort([μP, μP + 3σP, μP - 3σP])
    Drng = sort([μD, μD + 3σD, μD - 3σD])

    @info "Perimeter range ... $(Prng)"
    @info "Diameter range .... $(Drng)"

    μP, μD, σP, σD
end

# ╔═╡ 21ab6e10-9438-46b1-8ac7-3d93e4c70d89
md"""
## Tools
"""

# ╔═╡ Cell order:
# ╟─f175ae20-641b-11ee-2bf8-e72783c8a834
# ╟─a677c786-874c-4cc7-9871-6016064797a3
# ╟─5ef67a89-2d81-4863-b16f-a0f0a3c6321f
# ╟─13daafac-abfa-476d-b845-a78cc864559e
# ╟─e48a03d9-7189-4ff0-95cf-6bb942cab369
# ╟─915f13a0-c078-47e8-ba4e-2eb0ab8d7de8
# ╟─1eabfb91-5d33-4f4a-a9b4-6f4def62a1d6
# ╟─4dfcfd32-3b16-48b1-81e4-4371bb11a731
# ╟─bd7ff6f9-1578-4780-8666-6ca6cfa0ede8
# ╟─1d886b01-f1c5-4d85-ad61-970ed562c0e8
# ╟─a68a22cc-8a2f-4bee-82d7-14e1182381ca
# ╟─d0e39584-6c13-4804-8adc-5b5eb962940f
# ╟─9876504d-9732-47c1-9ea5-003e6a916c8f
# ╟─a87f0302-8e77-4dea-9e10-0469ecfdf6e3
# ╟─19ae6a86-684d-4dca-b75f-27b9bfe71b9c
# ╟─d003c1fb-2c61-4031-b7bc-9c11d6e2dcc2
# ╟─5ff1e62f-2a34-4dee-b2de-12fa670212c1
# ╟─48d89f5c-80e7-4c79-869a-c9dcd87b5c77
# ╟─b6efc1d0-f5cc-42ac-b46b-920f6ce6e72e
# ╟─2e4ce97e-2ec1-4691-a0f2-451915f9fd2f
# ╟─8d78856b-0a9c-4909-bb30-fc530197b9ca
# ╟─66de140a-3a62-4b74-9b2b-97efaff2f3a5
# ╟─e118a6e8-9b69-4987-801e-f2035687c41f
# ╟─5fcd5364-5d57-4711-9c63-8be40607aa56
# ╟─738db443-ad85-4f10-8864-9ba9b653fae6
# ╟─4983a5d8-1129-423e-8427-320edb206956
# ╟─7206a2a0-4c77-4861-988c-49249966d047
# ╟─c758b980-9aba-41dc-a696-4fa0ed8cee9e
# ╟─d33a74c4-bb51-44ba-b61a-3bddadb30730
# ╟─c1661200-1e5c-45b2-984e-56c91100c220
# ╟─f76010d5-491a-4e38-b111-daa78f15f34a
# ╟─944c1d88-edf7-4bba-8cb6-0490347b9f88
# ╟─cf128a15-dae6-44fc-a28d-c861a57fef91
# ╟─88769a22-d80b-4f3e-88cd-519400c47de9
# ╟─5cb8cc6e-1097-47a9-a0c0-ec93fbf711ef
# ╟─f7b401c6-a28d-4a3c-9dae-0cf471d7b162
# ╟─9b99bd14-8b19-4069-bac9-ebe147b6b515
# ╟─eb8e27b4-a2cd-422b-995e-3465478c3f70
# ╟─0948d221-cbac-4f6e-ae9f-c6a3ff3f0705
# ╟─078c50a8-4efa-463d-9110-61a6a38431bd
# ╟─65c99f9b-0322-4229-a108-ead8525c4091
# ╟─0d81315f-9962-450e-b42a-af7525754408
# ╟─c7eaa50c-0aaf-49c0-bd71-69246cf42bff
# ╟─21ab6e10-9438-46b1-8ac7-3d93e4c70d89
# ╟─52d557a4-8ef1-4b50-a145-7995745502c6
