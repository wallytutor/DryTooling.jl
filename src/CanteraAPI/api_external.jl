#############################################################################
# ERRORS
#############################################################################

mutable struct CanteraError <: Exception
    name::String
end

showerror(io::IO, e::CanteraError) = print(io, "CanteraError: $(e)")

#############################################################################
# General utilities (public)
#############################################################################

function appdelete()::Bool
    return ct_appdelete() == 0 || throw(CanteraError("appdelete"))
end

function resetstorage()::Bool
    return ct_resetStorage() == 0 || throw(CanteraError("resetstorage"))

end
function clearstorage()::Bool
    return ct_clearStorage() == 0 || throw(CanteraError("clearstorage"))
end

function suppress_thermo_warnings(flag::Bool)::Bool
    return (ct_suppress_thermo_warnings(Int32(flag)) == 0||
            throw(CanteraError("suppress_thermo_warnings")))
end

function use_legacy_rate_constants(flag::Bool)::Bool
    return (ct_use_legacy_rate_constants(Int32(flag)) == 0 ||
            throw(CanteraError("use_legacy_rate_constants")))
end

#############################################################################
# SolutionIndex (private)
#############################################################################

mutable struct SolutionIndex
    solution::Int32
    thermo::Int32
    kinetics::Int32
    transport::Int32

    function SolutionIndex(
            infile::String,
            name::String,
            transport::String
        )
        obj = new()
        obj.solution = soln_newSolution(infile, name, transport)
        if obj.solution < 0
            throw(CanteraError("SolutionIndex cannot create new solution"))
        end

        obj.thermo    = soln_thermo(obj.solution)
        obj.kinetics  = soln_kinetics(obj.solution)
        obj.transport = soln_transport(obj.solution)
        return obj
    end
end

function nelements(obj::SolutionIndex)::UInt32
    value = thermo_nElements(obj.thermo)
    if value == 0
        throw(CanteraError("nelements : no elements in thermo object"))
    end
    return convert(UInt32, value)
end

function nspecies(obj::SolutionIndex)::UInt32
    value = thermo_nSpecies(obj.thermo)
    if value == 0
        throw(CanteraError("nspecies : no species in thermo object"))
    end
    return convert(UInt32, value)
end

#############################################################################
# Solution (public)
#############################################################################

mutable struct Solution
    nelements::UInt32
    nspecies::UInt32
    
    _index::SolutionIndex
    _X::Vector{Float64}
    _Y::Vector{Float64}

    function Solution()
        return nothing
    end

    function Solution(
            infile::String,
            name::String,
            transport::String
        )
        obj = new()
        
        # Create index of internals.
        obj._index = SolutionIndex(infile, name, transport)

        # Get counters of elements/species.
        obj.nelements = nelements(obj._index)
        obj.nspecies  = nspecies(obj._index)

        # Initialize memory for mass/mole fractions.
        obj._X = zeros(obj.nspecies)
        obj._Y = zeros(obj.nspecies)

        return obj
    end
end

#############################################################################
# Solution (public - main)
#############################################################################

function set_TPX!(
        obj::Solution,
        T::Float64,
        P::Float64,
        X::Vector{Float64};
        norm::Bool = true
    )::Nothing
    settemperature(obj, T)
    setpressure(obj, P)
    setmolefractions(obj, X, norm=norm)
end

function equilibrate!(
        obj::Solution,
        XY::String;
        solver::String      = "auto",
        rtol::Float64       = 1.0e-09,
        maxsteps::Int32     = Int32(1000),
        maxiter::Int32      = Int32(100),
        loglevel::Int32     = Int32(0),
        print_results::Bool = false,
        show_thermo::Bool   = true,
        threshold::Float64  = 1.0e-14
    )::Nothing

    solver_code = Dict(
        "auto"              => -1,
        "element_potential" => 0,
        "gibbs"             => 1,
        "vcs"               => 2,
    )[solver]

    status = thermo_equilibrate(obj._index.thermo, XY, solver_code,
                                rtol, maxsteps, maxiter, loglevel)

    if status < 0
        throw(CanteraError("equilibrate! : failed ($(status))"))
    end
    
    if print_results
        status = thermo_print(obj._index.thermo, Int32(show_thermo), threshold)
        if status < 0
            throw(CanteraError("equilibrate! thermo_print : ($(status))"))
        end
    end
end

function delete!(obj::Solution)::Nothing
    # TODO: this is doing nothing!
    if soln_del(obj._index.solution) != 0
        throw(CanteraError("delete(obj::Solution)"))
    end
    obj = nothing
end

#############################################################################
# Solution (public - setters)
#############################################################################

function settemperature(obj::Solution, T::Float64)::Nothing
    if thermo_setTemperature(obj._index.thermo, T) < 0
        throw(CanteraError("settemperature"))
    end
end

function setpressure(obj::Solution, P::Float64)::Nothing
    if thermo_setPressure(obj._index.thermo, P) < 0
        throw(CanteraError("setpressure"))
    end
end

function setmolefractions(obj::Solution, X::Vector{Float64};
                          norm::Bool)::Nothing
    n = obj._index.thermo
    if thermo_setMoleFractions(n, obj.nspecies, X, Int(norm)) < 0
        throw(CanteraError("setmolefractions"))
    end
end

#############################################################################
# Solution (public - getters)
#############################################################################

function gettemperature(obj::Solution)::Float64
    T = thermo_temperature(obj._index.thermo)
    if T < 0.0
        throw(CanteraError("""
        gettemperature : get a negative temperature
        The parent solution object has been freed elsewhere!
        """))
    end
    return T
end

function getpressure(obj::Solution)::Float64
    P = thermo_pressure(obj._index.thermo)
    if P < 0.0
        throw(CanteraError("""
        getpressure : get a negative pressure
        The parent solution object has been freed elsewhere!
        """))
    end
    return P
end

function getmolefractions(obj::Solution)::Vector{Float64}
    status = thermo_getMoleFractions(
        obj._index.thermo, obj.nspecies, obj._X)
    if status < 0
        throw(CanteraError("""
        getmolefractions : get a negative status
        The parent solution object has been freed elsewhere!
        """))
    end
    return obj._X
end
