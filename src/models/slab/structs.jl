struct SLAB_Input{I <: Integer, F <: Number, A <: AbstractVector{F}}
    idspl::I
    ncalc::I
    wms::F
    cps::F
    tbp::F
    cmed0::F
    dhe::F
    cpsl::F
    rhosl::F
    spb::F
    spc::F
    ts::F
    qs::F
    as::F
    tsd::F
    qtis::F
    hs::F
    tav::F
    xffm::F
    zp::A
    z0::F
    za::F
    ua::F
    ta::F
    rh::F
    stab::F
    ala::F
end

SLAB_Input(;idspl,ncalc,wms,cps,tbp,cmed0,dhe,cpsl,
            rhosl,spb,spc,ts,qs,as,tsd,qtis,hs,tav,
            xffm,zp,z0,za,ua,ta,rh,stab,ala) = SLAB_Input(idspl,ncalc,wms,cps,tbp,cmed0,dhe,cpsl,
            rhosl,spb,spc,ts,qs,as,tsd,qtis,hs,tav,
            xffm,zp,z0,za,ua,ta,rh,stab,ala)

struct SLAB_Release_Gas_Props{F <: Number}
    wms::F
    cps::F
    ts::F
    rhos::F
    tbp::F
    cmed0::F
    cpsl::F
    dhe::F
    rhosl::F
    spa::F
    spb::F
    spc::F
end

struct SLAB_Spill_Chars{I <: Integer, F <: Number}
    idspl::I
    qs::F
    tsd::F
    qtcs::F
    qtis::F
    as::F
    ws::F
    bs::F
    hs::F
    us::F
end

struct SLAB_Field_Params{F <: Number, A <: AbstractVector{F}}
    tav::F
    hmx::F
    xffm::F
    tffm::F
    zp::A
end

struct SLAB_Ambient_Met_Props{F <: Number}
    wmae::F
    cpaa::F
    rhoa::F
    za::F
    pa::F
    ua::F
    ta::F
    rh::F
    uastr::F
    stab::F
    ala::F
    z0::F
    stb::F
    phimi::F
    phgam::F
    cmwa::F
    cmdaa::F
end

struct SLAB_Additional_Params{I <: Integer, F <: Number}
    ncalc::I
    nssm::I
    grav::F
    rr::F
    xk::F
end

struct SLAB_Wind_Profile{F <: Number}
    z0::F
    ala0::F
    zl::F
    hmx::F
    zt::F
    cu1::F
    cu2::F
end

struct SLAB_Other_Params{F <: Number}
    tgon::F
    bse::F
    hrf::F
    urf::F
    cf0::F
    rcf::F
    tau0::F
    at0::F
    afa::F
end

struct SLAB_Params{I <:Integer, F <: Number, A <: AbstractVector{F}}
    rgp::SLAB_Release_Gas_Props{F}
    spl::SLAB_Spill_Chars{I,F}
    fld::SLAB_Field_Params{F,A}
    met::SLAB_Ambient_Met_Props{F}
    xtra::SLAB_Additional_Params{I,F}
    wps::SLAB_Wind_Profile{F}
    othr::SLAB_Other_Params{F}
end


# for other constants that are initialized
# and passed to the integrators
struct SLAB_Loop_Init{I <: Integer, F <: Number}
    nxi::I
    msfm::I
    mnfm::I
    mffm::I
    gam::F
    ft::F
    fu::F
    fv::F
    fw::F
    fug::F
    bbv0::F
    bv0::F
    r0::F
    cp0::F
    alfg::F
    sru0::F
    htp0::F
    ubs20::F
    rmi::F
    bx::F
    bbx::F
    bbvx0::F
    bvx0::F
    xcc0::F
    bxs0::F
end


# SLAB state vectors
# instantaneous spatially averaged cloud parameters
struct SLAB_Vecs{F <: Number, A <: AbstractVector{F}}
    x::A
    zc::A
    h::A
    bb::A
    b::A
    bbx::A
    bx::A
    cv::A
    rho::A
    t::A
    u::A
    uab::A
    cm::A
    cmev::A
    cmda::A
    cmw::A
    cmwv::A
    wc::A
    vg::A
    ug::A
    w::A
    v::A
    vx::A
    tim::A
    beta::A
    qint::A
    betax::A
    xccp::A
    tccp::A
end

SLAB_Vecs(t::Type,n::Integer) = SLAB_Vecs(
    zeros(t,n),#x
    zeros(t,n),#zc
    zeros(t,n),#h
    zeros(t,n),#bb
    zeros(t,n),#b
    zeros(t,n),#bbx
    zeros(t,n),#bx
    zeros(t,n),#cv
    zeros(t,n),#rho
    zeros(t,n),#t
    zeros(t,n),#u
    zeros(t,n),#uab
    zeros(t,n),#cm
    zeros(t,n),#cmev
    zeros(t,n),#cmda
    zeros(t,n),#cmw
    zeros(t,n),#cmwv
    zeros(t,n),#wc
    zeros(t,n),#vg
    zeros(t,n),#ug
    zeros(t,n),#w
    zeros(t,n),#v
    zeros(t,n),#vx
    zeros(t,n),#tim
    zeros(t,n),#beta
    zeros(t,n),#qint
    zeros(t,n),#betax
    zeros(t,n),#xccp
    zeros(t,n)#tccp
)

struct SLAB_CC_Vecs{F <: Number, A <: AbstractVector{F}}
    x::A
    cc::A
    b::A
    betac::A
    zc::A
    sig::A
    t::A
    xc::A
    bx::A
    betax::A
    tim::A
    tcld::A
    bbc::A
end

struct SLAB_Output{I <: Integer, F <: Number, A <: AbstractVector{F}}
    p::SLAB_Params{I,F,A}
    s::SLAB_Vecs{F,A}
    cc::SLAB_CC_Vecs{F,A}
end