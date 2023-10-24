# -*- coding: utf-8 -*-
export unknowns
export params

function unknowns(k::AbstractKineticsMechanism)
    @error "An specialization of this method is expexted!"
end

function params(k::AbstractKineticsMechanism)
    @error "An specialization of this method is expexted!"
end
