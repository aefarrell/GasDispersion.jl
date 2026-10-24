struct SLAB_Plume_Integrated{T}
# the integrated state variables defined in section 2.4 of the SLAB User's Guide
    r::T    # R = ρUBh
    bbv::T  # plume half-width B
    bv::T   # plume half-width parameter b
    g::T    # G = RVg (I think)
    gw::T   # ¯\_(ツ)_/¯
    zc::T   # centerline height
    qint::T # ¯\_(ツ)_/¯
    sft::T  # Fₜ = ∫fₜdx
    sfu::T  # Fᵤ = ∫fᵤdx
    sfy::T  # Fy = ∫fydx
    sfz::T  # Fz = ∫fzdx
end

struct SLAB_Thermo_State{T,P}
# the plume variables that are calculated by the thermo subroutine
    cm::T     # mass fraction of the released species in the plume
    cv::T     # volume fraction of the released species in the plume
    cmw::T    # mass fraction of water in the plume
    cmda::T   # mass fraction of dry air in the plume
    cmwv::T   # mass fraction of water vapor in the plume
    cmwd::T   # mass fraction of water in the plume
    cmev::T   # mass fraction of evaporated species in the plume
    cmed::T   # 
    cp::T     # specific heat capacity of the plume
    t::T      # temperature of the plume
    rho::T    # density of the plume
    phase::P  # phase of the plume
    transported_energy::T # 
end

struct SLAB_Plume_State{I,T}
    integrated::I
    thermo::T
end