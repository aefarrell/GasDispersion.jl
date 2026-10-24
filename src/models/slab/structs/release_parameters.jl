struct SLAB_Release_Source{T}
    bse::T   # B_{se} effective source half-width
    qs::T   # q_s source mass emission rate
    wms::T  # M_s molar mass of released species
    cps::T  # C_{ps} specific heat capacity of released species
    cpsl::T # C_{psl} specific heat capacity of released species in liquid phase
    ts::T   # T_s source temperature
    antoine::SLAB_Antoine{T} # Antoine equation parameters for the source species
end

struct SLAB_Release_Ambient_Air{T}
    wmae::T  # M_{ae} molar mass of air
    wma::T   # M_a molar mass dry air
    wmw::T   # M_w molar mass water
    cmwa::T  # m_{wa} mass fraction water in ambient air
    cpaa::T  # C_{paa} specific heat capacity of ambient air
    cpa::T   # C_{pa} specific heat capacity of dry air
    cpwv::T  # C_{pwv} specific heat capacity of water vapor
    cpwl::T  # C_{pwl} specific heat capacity of liquid water
    ta::T    # T_a ambient temperature
    antoine::SLAB_Antoine{T} # Antoine equation parameters for water
end

struct SLAB_Release_Parameters{T}
    source::SLAB_Release_Source{T}
    ambient_air::SLAB_Release_Ambient_Air{T}
end