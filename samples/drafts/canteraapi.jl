cantera = "C:\\Program Files\\Cantera\\bin\\cantera_shared.dll"
	
ENV["CANTERA_SHARED"] = cantera
@assert haskey(ENV, "CANTERA_SHARED")

include("..\\..\\src\\CanteraAPI.jl")

# import DryTooling.CanteraAPI as ct;
import .CanteraAPI as ct;

@assert ct.appdelete()
@assert ct.resetstorage()
@assert ct.clearstorage()
@assert ct.suppress_thermo_warnings(true)
@assert ct.use_legacy_rate_constants(false)

sol = ct.Solution("gri30.yaml", "gri30", "mixture-averaged")
gas = ct.Solution("gri30.yaml", "gri30", "mixture-averaged")

Xᵣ = zeros(sol.nspecies)
Xᵣ[1] = 1.0

Tᵣ = 3500.0
Pᵣ = 50000.0

ct.set_TPX!(sol, Tᵣ, Pᵣ, Xᵣ; norm = true)

@assert ct.gettemperature(sol) ≈ Tᵣ
@assert ct.getpressure(sol) ≈ Pᵣ
@assert all(ct.getmolefractions(sol) ≈ Xᵣ)

ct.equilibrate!(sol, "HP", print_results = true)

# ct.delete!(sol)
# ct.gettemperature(sol)