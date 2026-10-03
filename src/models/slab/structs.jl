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

SLAB_Input(;idspl,ncalc,wms,cps,tbp,cmed0,dhe,cpsl,rhosl,spb,spc,ts,qs,as,tsd,qtis,hs,
            tav,xffm,zp,z0,za,ua,ta,rh,stab,ala) = SLAB_Input(idspl,ncalc,wms,cps,tbp,
            cmed0,dhe,cpsl,rhosl,spb,spc,ts,qs,as,tsd,qtis,hs,tav,xffm,zp,z0,za,ua,ta,
            rh,stab,ala)

function SLAB_Input(idspl,ncalc,wms,cps,tbp,cmed0,dhe,cpsl,rhosl,spb,spc,ts,qs,as,tsd,
                    qtis,hs,tav,xffm,zp,z0,za,ua,ta,rh,stab,ala)

    return SLAB_Input(idspl,ncalc,promote(wms,cps,tbp,cmed0,dhe,cpsl,rhosl,spb,spc,ts,
                      qs,as,tsd,qtis,hs,tav,xffm)...,zp,promote(z0,za,ua,ta,rh,stab,ala)...)
end

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
    bbx::A
    betax::A
    tim::A
    tcld::A
    bbc::A
end

struct SLAB_Interpolations{C,B,BC,Z,S,XC,BX,BX2,BXX,BBXX,TCLD,O}
    cc::C
    b::B
    betac::BC
    zc::Z
    sig::S
    xc::XC
    bx::BX
    betax::BX2
    bx_x::BXX
    bbx_x::BBXX
    tcld::TCLD
    ode_solution::O
end

function _slab_legacy_interpolations(cc::SLAB_CC_Vecs)
    xperm = sortperm(cc.x)
    tperm = sortperm(cc.t)
    return SLAB_Interpolations(
        AkimaInterpolation(cc.cc[xperm], cc.x[xperm]),
        AkimaInterpolation(cc.b[xperm], cc.x[xperm]),
        AkimaInterpolation(cc.betac[xperm], cc.x[xperm]),
        AkimaInterpolation(cc.zc[xperm], cc.x[xperm]),
        AkimaInterpolation(cc.sig[xperm], cc.x[xperm]),
        AkimaInterpolation(cc.xc[tperm], cc.t[tperm]),
        AkimaInterpolation(cc.bx[tperm], cc.t[tperm]),
        AkimaInterpolation(cc.betax[tperm], cc.t[tperm]),
        AkimaInterpolation(cc.bx[xperm], cc.x[xperm]),
        AkimaInterpolation(cc.bbx[xperm], cc.x[xperm]),
        AkimaInterpolation(cc.tcld[xperm], cc.x[xperm]),nothing)
end

abstract type AbstractSLABSolution end

mutable struct SLAB_Steady_State{F <: Number}
    r::F
    bbv::F
    bv::F
    g::F
    gw::F
    sft::F
    sfu::F
    sfy::F
    sfz::F
    zc::F
    qint::F
    h::F
    b::F
    bb::F
    rho::F
    t::F
    u::F
    uab::F
    vg::F
    wc::F
    htp::F
    cm::F
    cmw::F
    cmwv::F
    cmev::F
    cp::F
    ft::F
    fu::F
    fv::F
    fw::F
    fug::F
    ubs2::F
end

mutable struct SLAB_Transient_State{F <: Number}
    x::F
    tim::F
    r::F
    bbv::F
    bv::F
    bbvx::F
    bvx::F
    g::F
    gw::F
    gx::F
    sft::F
    sfu::F
    sfx::F
    sfy::F
    sfz::F
    zc::F
    qint::F
    h::F
    b::F
    bb::F
    bx::F
    bbx::F
    rho::F
    t::F
    u::F
    uab::F
    vg::F
    ug::F
    wc::F
    cm::F
    cmw::F
    cmwv::F
    cmev::F
    cp::F
    ft::F
    fu::F
    fv::F
    fw::F
    fug::F
    ubs2::F
end

struct SLAB_Steady_Solution{I <: Integer, F <: Number, A <: AbstractVector{F}, O, G} <: AbstractSLABSolution
    params::SLAB_Params{I,F,A}
    state::SLAB_Vecs{F,A}
    cc::SLAB_CC_Vecs{F,A}
    initial::SLAB_Steady_State{F}
    ode_solution::O
    ode_segments::G
end

struct SLAB_Transient_Solution{I <: Integer, F <: Number, A <: AbstractVector{F}} <: AbstractSLABSolution
    params::SLAB_Params{I,F,A}
    state::SLAB_Vecs{F,A}
    cc::SLAB_CC_Vecs{F,A}
    initial::SLAB_Transient_State{F}
end

struct SLAB_ODE_FieldInterpolation{S,G,P,F,BX,BBX,T,Fallback}
    solution::S
    segments::G
    params::P
    field::Symbol
    x0::F
    bx_x::BX
    bbx_x::BBX
    tcld::T
    fallback::Fallback
end

function (itp::SLAB_ODE_FieldInterpolation)(x)
    index = findlast(segment -> segment.x0 <= x <= segment.x1, itp.segments)
    index === nothing && return itp.fallback(x)
    y = itp.solution(x)
    state, _ = _slab_steady_project(y, itp.segments[index].context, x)
    if itp.field === :b
        return state.b
    elseif itp.field === :zc
        return state.zc
    end
    p = itp.params
    cv = (p.met.wmae*state.cm)/(p.rgp.wms+(p.met.wmae-p.rgp.wms)*state.cm)
    point = _slab_editcc_point(p,x,itp.x0,state.zc,state.h,state.b,state.beta,
        state.uab,state.cm,cv,itp.bx_x(x),itp.bbx_x(x),itp.tcld(x))
    return getproperty(point,itp.field)
end

function _slab_output_interpolations(cc::SLAB_CC_Vecs, steady_solution, params)
    legacy = _slab_legacy_interpolations(cc)
    if steady_solution === nothing || steady_solution.ode_solution === nothing ||
       isempty(steady_solution.ode_segments)
        return legacy
    end
    makefield(field, fallback) = SLAB_ODE_FieldInterpolation(
        steady_solution.ode_solution, steady_solution.ode_segments, params, field,
        cc.x[1], legacy.bx_x, legacy.bbx_x, legacy.tcld, fallback)
    return SLAB_Interpolations(
        makefield(:cc,legacy.cc), makefield(:b,legacy.b),
        makefield(:betac,legacy.betac), makefield(:zc,legacy.zc),
        makefield(:sig,legacy.sig), legacy.xc, legacy.bx, legacy.betax,
        legacy.bx_x, legacy.bbx_x, legacy.tcld, steady_solution.ode_solution)
end

struct SLAB_Output{I <: Integer, F <: Number, A <: AbstractVector{F}, P,
                   S <: Union{Nothing,SLAB_Steady_Solution{I,F,A}},
                   T <: Union{Nothing,SLAB_Transient_Solution{I,F,A}}}
    p::SLAB_Params{I,F,A}
    s::SLAB_Vecs{F,A}
    cc::SLAB_CC_Vecs{F,A}
    interpolations::P
    steady::S
    transient::T
end

SLAB_Output(params, state, cc) = SLAB_Output(params, state, cc,
                                              _slab_legacy_interpolations(cc), nothing, nothing)

function _slab_initial_steady_state(vecs::SLAB_Vecs{F}, vars::SLAB_Loop_Init{I,F}) where {I,F}
    return SLAB_Steady_State(vars.r0, vars.bbv0, vars.bv0, zero(F), zero(F), zero(F), zero(F),
        zero(F), zero(F), vecs.zc[1], vecs.qint[1], vecs.h[1], vecs.b[1], vecs.bb[1],
        vecs.rho[1], vecs.t[1], vecs.u[1], vecs.uab[1], vecs.vg[1], vecs.wc[1],
        vars.htp0, vecs.cm[1], vecs.cmw[1], vecs.cmwv[1], vecs.cmev[1], vars.cp0,
        vars.ft, vars.fu, vars.fv, vars.fw, vars.fug, vars.ubs20)
end

function _slab_initial_transient_state(vecs::SLAB_Vecs{F}, vars::SLAB_Loop_Init{I,F}, index::Integer) where {I,F}
    return SLAB_Transient_State(vecs.x[index], vecs.tim[index], vars.r0, vars.bbv0, vars.bv0,
        vars.bbvx0, vars.bvx0, zero(F), zero(F), zero(F), zero(F), zero(F), zero(F), zero(F), zero(F),
        vecs.zc[index], vecs.qint[index], vecs.h[index], vecs.b[index], vecs.bb[index],
        vecs.bx[index], vecs.bbx[index], vecs.rho[index], vecs.t[index], vecs.u[index],
        vecs.uab[index], vecs.vg[index], vecs.ug[index], vecs.wc[index], vecs.cm[index],
        vecs.cmw[index], vecs.cmwv[index], vecs.cmev[index], vars.cp0, vars.ft, vars.fu,
        vars.fv, vars.fw, vars.fug, vars.ubs20)
end