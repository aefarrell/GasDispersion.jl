# subroutine editcc
# calculates the time averaged volume concentration

function editcc(vecs::SLAB_Vecs{F,A}, params::SLAB_Params{I,F,A}, mffm<:I) where {I <: Integer, F <: AbstractFloat}
    tcmx = 2.0*vecs.bbx[end]
end