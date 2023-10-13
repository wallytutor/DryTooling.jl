# -*- coding: utf-8 -*-

""" Mixture mean molecular mass [kg/mol]. """
function meanmolecularmass(
        Y::Union{Vector{Float64},SubArray},
        W::Vector{Float64}
    )::Float64
    return 1.0 / sum(@. Y / W)
end

""" Convert mass fractions to mole fractions. """
function massfraction2molefraction(Y, M::Vector{Float64})
    return meanmolecularmass(M, Y) * @. Y / M
end

""" Convert mole fractions to mass fractions. """
function molefraction2massfraction(X, M::Vector{Float64})
    return (@. X * M) / sum(@. X * M)
end
