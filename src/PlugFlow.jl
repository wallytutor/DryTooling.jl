# -*- coding: utf-8 -*-
module PlugFlow
    
export  RectangularReactorGeometry

struct RectangularReactorGeometry
    """
        RectangularReactorGeometry

    Geometric description of rectangular reactor bounding volume.

    $(TYPEDFIELDS)
    """

    "Reactor total height/length [m]"
    H::Float64

    "Reactor cross-section depth [m]"
    D::Float64

    "Reactor cross-section width [m]"
    W::Float64

    "Reactor cross-section perimeter [m]"
    P::Float64

    "Reactor cross-section area [m²]"
    A::Float64

    "Reactor total volume [m³]"
    V::Float64

    function RectangularReactorGeometry(; H, D, W)
        P = 2 * (D + W)
        A = D * W
        V = A * H
        return new(H, D, W, P, A, V)
    end
end

end # module PlugFlow