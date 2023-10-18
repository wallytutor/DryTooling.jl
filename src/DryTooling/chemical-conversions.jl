# -*- coding: utf-8 -*-

""" Mixture mean molecular mass [kg/mol]. """
function meanmolecularmass(
        Y::Union{Vector{Float64},SubArray},
        W::Vector{Float64}
    )::Float64
    return 1.0 / sum(@. Y / W)
end

# TODO maybe the above is not providing any advantage here!
function meanmolecularmass(Y, W)
    return 1.0 / sum(@. Y / W)
end

""" Convert mass fractions to mole fractions. """
function massfraction2molefraction(Y, W::Vector{Float64})
    return meanmolecularmass(Y, W) * @. Y / W
end

""" Convert mole fractions to mass fractions. """
function molefraction2massfraction(X, W::Vector{Float64})
    return (@. X * W) / sum(@. X * W)
end
