# original SLAB subroutines

include("subroutines/slope.jl")
include("subroutines/slopepf.jl")
include("subroutines/solve.jl")
include("subroutines/solvepf.jl")
include("subroutines/thermo.jl")
include("subroutines/eval.jl")
include("subroutines/evalpf.jl")
include("subroutines/entran.jl")
include("subroutines/store.jl")

# integration routines
include("subroutines/integrate_steadystate.jl")
include("subroutines/integrate_transient.jl")

# intialize parameters

function _slab_init_rgps(wms::F,cps::F,tbp::F,cmed0::F,dhe::F,cpsl::F,rhosl::F,spb::F,
                         spc::F,ts::F) where F <: AbstractFloat

    # c  saturation pressure constant default
    if spb < zero(F)
        spb = dhe*wms/rr
        spc = 0.0
    end

    #c  boiling point temperature check
    spa = spb/(tbp+spc)
    if ts < tbp
        ts = tbp
    end
    if cmed0 > zero(F)
        ts = tbp
    end

    rhos = wms*pa/(rr*ts)

    return SLAB_Release_Gas_Props(wms,cps,ts,rhos,tbp,cmed0,cpsl,dhe,rhosl,spa,spb,spc)
end

function _slab_init_met(z0::F,za::F,ua::F,ta::F,rh::F,stab::F,ala::F) where F <: AbstractFloat

    #c ================================
    #c  calculate stability parameters
    #c ================================
    al1 = 0.0081/(z0^0.3044)
    if z0 > 0.0111
       al2 = 0.0385/(z0^0.1715)
       al3 = 0.0875/(z0^0.1028)
    else
       al2 = al1 + 0.0137/(z0^0.1715) + 0.0218
       al3 = al2 + 0.0557
    end

    en2 = log(al2/al1)/log(2.)
    en3 = log(al3/al1)/log(3.)
    eni = en3 + en3 - en2
    dln = en2 - eni
    alm = al1*(3.5^( eni + dln/3.25 ))

    if stab == zero(F)
        aal = abs(ala)       
        if aal < al2
            astb = (aal/al1)^(1/en2)
        elseif aal < alm
            ral = (aal-al2)/(al3-al2)
            en = eni + dln/(1+ral*ral)
            astb = (aal/al1)^(1/en)
        else
            astb = 3.5
        end

        if ala >= zero(F)
            stb = astb
        else
            stb = -astb
        end

        stab = 4.0 + stb

    else
        stb = stab - 4.0
        astb = abs(stb)
        if astb < 2.0
            aal = al1*(astb^en2)
        elseif astb < 3.5
            en = eni + dln/(1+(astb-2)*(astb-2))
            aal = al1*(astb^en)
        else
            aal = alm
        end

        if stb >= zero(F)
            ala = aal
        else
            ala = -aal
        end
    end

    #c  atmospheric and physical constants
    rpwa = 0.01*rh*exp(spaw-spbw/ta)
    cmwa = wmw*rpwa/(wma+(wmw-wma)*rpwa)
    cmdaa = 1-cmwa
    cpaa = cmdaa*cpa + cmwa*cpwv
    wmae = wma*wmw/(wmw+(wma-wmw)*cmwa)

    #c  calculate ambient meteorological values
    rhoa = wmae*pa/(rr*ta)
    
    if stb < zero(F)
        zl = exp(-0.8*stb)
    else
        zl = 1.0 + 0.8*stb
    end

    if za > 3.0
        z = za
    else
        z = 3.0
    end

    ala0 = ala*(1 + z/zl)
    hmx = 130.0*(2.0^(3.0 - stb))

    phimi, phgam = 0.0, 0.0
    if stb < zero(F)
        phimi = 1.0/sqrt(sqrt(1.0 - 16.0*zl*ala0))
        phgam = -8.0*ala0/(1.0 - phimi)
    end

    #c  initialize velocity function
    zt = 2.71828183*z0
    cu1 = 1.0/zt
    cu2 = 0.0
    uf = _slab_uafn(zt,z0,ala0,zl,hmx,zt,cu1,cu2)

    if ala0 < zero(F)
        phmi = 1/sqrt(sqrt(1 - 16*zl*ala0))
        gu = -8*ala0/(1 - phmi)
        phm = phmi + (1 - phmi)/sqrt(1 + gu*zt)
    else
        phm = 1 + 5*ala0*zt/(1 + zt/zl)
    end
    
    ufp = phm*(1 - zt/hmx)/zt
    cu1 = (2*uf - zt*ufp)/zt
    cu2 = (zt*ufp - uf)/(zt*zt)

    _wp = SLAB_Wind_Profile(z0,ala0,zl,hmx,zt,cu1,cu2)
  
    #c  calculatre friction velocity
    uastr = xk*ua/_slab_uafn(za,_wp)

    _met = SLAB_Ambient_Met_Props(wmae,cpaa,rhoa,za,pa,ua,ta,rh,uastr,stab,ala,z0,stb,
                                  phimi,phgam,cmwa,cmdaa)

    return _met, _wp
end

function _slab_init_hjet(idspl::I,ncalc::I,msfm::I,mnfm::I,mffm::I,wms::F,cps::F,tbp::F,
                         cmed0::F,dhe::F,cpsl::F,rhosl::F,spb::F,spc::F,ts::F,qs::F,as::F,
                         tsd::F,qtis::F,hs::F,tav::F,xffm::F,zp::AbstractVector{F},z0::F,
                         za::F,ua::F,ta::F,rh::F,stab::F,ala::F) where { I <: Integer, 
                         F <: AbstractFloat}
    
    # intialize the release gas properties
    _rgps = _slab_init_rgps(wms,cps,tbp,cmed0,dhe,cpsl,rhosl,spb,spc,ts)

    # initialize the ambient meteorological parameters
    if z0 < zero(F) 
        @error "z0 must be positive"
    end
    _met, _wp = _slab_init_met(z0,za,ua,ta,rh,stab,ala)
    
    # initialize the spill characteristics
    qtcs = qs*tsd
    qtis = 0.0
    ws = 0.0
    bs = 0.5*sqrt(as)

    cmev0 = 1.0 - cmed0
    rhoa = _met.rhoa
    alfm = (_met.wmae/wms)*cmev0
    betm = (rhoa/rhosl)*cmed0
    rho = rhoa*ta/(alfm*ts+betm*ta)
    us = qs/(rho*as)

    _spl = SLAB_Spill_Chars(idspl,qs,tsd,qtcs,qtis,as,ws,bs,hs,us)

    # initialize the field parameters

    hmx = _wp.hmx
    if hs > (hmx-bs)
        @error "input source height hs is greater than the calculated mixing layer height hmx minus the stack half width bs.  program terminated."
    end

    uastr = _met.uastr
    uaz2 = (uastr/xk)*_slab_uafn(2.0,_wp)
    tffm = xffm/uaz2
    _fld = SLAB_Field_Params(tav,hmx,xffm,tffm,zp)

    # initialize additional parameters
    nssm = 3*ncalc
    _aps = SLAB_Additional_Params(ncalc,nssm,grav,rr,xk)

    # initialize other constants
    tgon = 1.0
    bse = 1.0
    hrf = 4.0
    urf = _slab_uafn(hrf,_wp)
    cf0 = uastr/uaz2
    rcf = sqrt(cf0/cf00)
    tau0 = 10.0
    at0 = 0.0
    afa = 0.08*(((at0+tau0*exp(-at0/tau0))/tav0)^0.2)
       
    # intialize main storage arrays
    vecs = SLAB_Vecs(F,mffm)
        
    alfg = 0.0
    bb = bs
    h = bb+bb
    hh = 0.5*h
    htp = hs+hh
    if hs <= hh
        alfg = 0.25
        bb = 0.5*(sqrt(hs*hs + as+as) - hs)
        h = bb + hs
        htp = h
    end

    if htp > hmx
        htp = hmx
    end

    if h > hmx
        bb = h*bb/hmx
        h = hmx
    end

    vecs.t[1] = ts
    vecs.cm[1] = 1.0
    vecs.cv[1] = 1.0
    vecs.cmw[1] = 0.0
    vecs.cmwv[1] = 0.0

    cmev0 = 1 - cmed0
    vecs.cmev[1] = cmev0
    vecs.cmda[1] = 0.0
    cp0 = cmev0*cps + cmed0*cpsl

    zc = hs
    zcmx = hmx - 0.5*h

    if zc > zcmx
        zc = zcmx
    end
    vecs.zc[1] = zc
    xcc0 = 1.0
    bxs0 = 0.0
    tim = tsd

    vecs.x[1] = 1.0
    vecs.xccp[1] = 1.0
    vecs.tccp[1] = 0.0
    vecs.tim[1] = tsd

    vecs.vg[1] = 0.0

    vecs.wc[1] = 0.0
    vecs.ug[1] = 0.0

    fu = 0.0
    fv = 0.0
    fw = 0.0
    ft = 0.0
    vecs.qint[1] = 0.0

    vecs.bb[1] = bb
    bbv0 = bb
    b = 0.9*bb
    vecs.b[1] = b
    bv0 = b

    beta = sqrt(bb*bb-b*b)/sr3
    vecs.beta[1] = beta
    vecs.rho[1] = rho

    vecs.h[1] = h
    htp0 = htp

    uab = 0.0
    dh = htp/5.0
    h1 = -0.5*dh
    h2 = 0.0

    u1, u2 = 0, 0
    for i in 1:5
        h1 = h1 + dh
        h2 = h2 + dh
        u1 = _slab_uafn(h1,_wp)
        u2 = _slab_uafn(h2,_wp)
        uab = uab + 4*u1 + 2*u2
    end

    uastr = _met.uastr
    uab = (uastr/xk)*(uab-u2)/30
    vecs.uab[1] = uab

    u = us
    vecs.u[1] = u

    r0 = rho*u*bb*h
    srug = .5*alfg*grav*(rho-rhoa)*bb*h*h
    sru0 = r0*u + srug

    delu = (rhoa/rho)*abs(uab-us)
    if htp > h
        ubs20 = uastr^2
        ustr2 = uastr^2 + cf1*delu*delu
    else
        umgs = (uastr/uab)*us
        ubs20 = umgs*umgs
        ustr2 = ubs20 + cf1*delu*delu
    end
    uah = _slab_uafn(htp,_wp)
    w = (urf/uah)* aa*xk*sqrt(ustr2)
    vecs.w[1] = w

    if ala < 0
        stby = 1 - rcf*10*ala
    else
        stby = 1/(1 + rcf*10*ala)
    end

    atot = afa*stby
    rab = 0.5*sigb/atot
    va = atot*uab/(1.0 + rab*bb/sr3)
    vjp = aa*xk*delu
    v = sqrt(va*va + cf1*vjp*vjp)
    vecs.v[1] = v
    vecs.vx[1] = 0.0

    # additional constants to initialize the integration
    msfm = 1
    mnfm = mffm - 1
    nxi = msfm + 1
  
    vars = SLAB_Loop_Init(nxi,msfm,mnfm,mffm,tgon,bse,urf,cf0,rcf,afa,ft,fu,fv,fw,bbv0,
                          bv0,r0,cp0,alfg,sru0,htp0,ubs20,xcc0,bxs0)
    
    params = SLAB_Params(_rgps,_spl,_fld,_met,_aps,_wp)

    return vecs,vars,params


end


