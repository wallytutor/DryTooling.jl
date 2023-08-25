# -*- coding: utf-8 -*-

""" Mixture mean molecular mass [kg/mol]. """
function meanmolecularmass(M::Vector{Float64}, Y)
    return  1.0 / sum(@. Y / M)
end

""" Convert mass fractions to mole fractions. """
function massfraction2molefraction(M::Vector{Float64}, Y)
    return meanmolecularmass(M, Y) * @. Y / M
end

""" Convert mole fractions to mass fractions. """
function molefraction2massfraction(M::Vector{Float64}, X)
    return (@. X * M) / sum(@. X * M)
end
