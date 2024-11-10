function _slab_init_vjet(idspl::I,ncalc::I,msfm::I,mnfm::I,mffm::I,wms::F,cps::F,tbp::F,
                        cmed0::F,dhe::F,cpsl::F,rhosl::F,spb::F,spc::F,ts::F,qs::F,as::F,
                        tsd::F,qtis::F,hs::F,tav::F,xffm::F,zp::AbstractVector{F},z0::F,
                        za::F,ua::F,ta::F,rh::F,stab::F,ala::F) where { I <: Integer, 
                        F <: AbstractFloat}

# initialize additional parameters
nssm = 3*ncalc
_aps = SLAB_Additional_Params(ncalc,nssm,grav,rr,xk)

# intialize the release gas properties
_rgps = _slab_init_rgps(wms,cps,tbp,cmed0,dhe,cpsl,rhosl,spb,spc,ts)

# initialize the ambient meteorological parameters
if z0 < zero(F) 
    @error "z0 must be positive"
end
_met, _wp = _slab_init_met(z0,za,ua,ta,rh,stab,ala)

# initialize the field parameters
hmx = _wp.hmx
bs = 0.5*sqrt(as)
if hs > ( hmx - bs )
    @error "input source height hs is greater than the calculated mixing layer height hmx minus the stack half width bs.  program terminated."
end

# calculate friction velocity
uastr = _met.uastr
uaz2 = (uastr/xk)*_slab_uafn(2.0,_wp)
tffm = xffm/uaz2
cf0 = uastr/uaz2
rcf = sqrt(cf0/cf00)
hrf = 4.0
urf = _slab_uafn(hrf,_wp)

if xffm < 2.0
    xffm = 2.0
end

_fld = SLAB_Field_Params(tav,hmx,xffm,tffm,zp)

# initialize other constants
tau0 = 10.0
at0 = 0.0
afa = 0.08*(((at0+tau0*exp(-at0/tau0))/tav0)^0.2)

# line 150
nxtr = mffm + 1
idpf = 0
qtcs = qs*tsd
tgon = 1.0


# c ======================
# c  vertical jet release
# c ======================

t0 = ts
cmev0 = 1.0 - cmed0
rhoa = _met.rhoa
alfm = (_met.wmae/wms)*cmev0
betm = (rhoa/rhosl)*cmed0
rho = rhoa*ta/(alfm*ts+betm*ta)

ws = qs/(rho*as)
us = 0.0
bse = 1.0
qtis = 0.0
xcc0 = 1.0
bxs0 = 0.0
#tim = tsd
_othr = SLAB_Other_Params(tgon,bse,hrf,urf,cf0,rcf,tau0,at0,afa)

# c plume rise calculation
hsu = max(hs,1.0)
uahs = (uastr/xk)*_slab_uafn(hsu,_wp)
betu = 0.4 + 1.2*uahs/ws
hprm = 1.0363*ws*bs/(betu*sqrt(uahs*uastr))
h67 = hprm^0.85714

for i in 1:3
    hprt = h67*((hs+hprm)^0.14286)
    hprm = (hs + 0.85714*hprm)*hprt/(hs + hprm - 0.14286*hprt)
end

# c dense gas, momentum, or buoyant jet
if rho > rhoa
    #c dense gas jet plume rise
    ds = 4.0*bs/√(π)
    fr2 = ws^2/(grav*ds*(rho-rhoa)/rhoa)
    sgs = rho/rhoa
    rwus = ws/uahs
    hprd = 1.32*ds*((rwus*sgs*fr2)^0.3333)
    hpr = hprm*hprd/sqrt(hprd*hprd + hprm*hprm)

    if (hs + hpr) > hmx
        hpr = hmx - hs
    end

    xpr = 0.435*hpr^3/(sgs*(ds*rwus)^2)
    cvpk = 1.69*rwus/((hpr/ds)^1.85)
    cvpk = max(0.99,cvpk)
    cv = (0.5236*(1.0 - cvpk) + cvpk)*cvpk
    cm = wms*cv/(_met.wmae+(wms-_met.wmae)*cv)
    r0 = 0.5*qs
    r = r0/cm
    cm0 = 1.0
    cmev0 = 1.0 - cmed0
    cmw0 = 0.0
    cmwv0 = 0.0
    cp0 = cmev0*cps + cmed0*cpsl
    xn = 1.0 + xpr
    sft = 0.0

    _spl = SLAB_Spill_Chars(idspl,qs,tsd,qtcs,qtis,as,ws,bs,hs,us)
    params = SLAB_Params(_rgps,_spl,_fld,_met,_aps,_wp,_othr)

    # call thermo
    timn, rmi = 0.0, 0.0 # these need to be initialized for thermo, but not used
    cm,cv,cmw,cmwv,cmev,t,rho,_cp = _slab_sub_thermo(params,idpf,xn,timn,rmi,
                                     t0,cmev0,cm0,cmw0,cmwv0,cp0,r,r0,sft,bse)
    

    # intialize main storage arrays
    vecs = SLAB_Vecs(F,mffm)
    
    # go to plume rise interpolation - line 500
    
    # c ========================================
    # c  vertical jet plume rise interpolations
    # c ========================================
    msfm = 2 + xpr/(2bs)
    msfm = min(round(Int,msfm),11)
    dxpr = xpr/(msfm-1.0)

    zcr = hs + hpr
    uab = (uastr/xk)*_slab_uafn(zcr,_wp)
    cfpr = 3.0*(2.0 + (xpr/bs))*cf1
    cmpr = sqrt((1.0 - cm)^2 + 4.0*cm*cfpr)
    cupr = 0.5*uab/(1.0 - cfpr)
    u = cupr * (1.0 - cm - 2cfpr + cmpr)

    ar = qs/(rho*u*cm)
    am = qs/(rho*u)
    hft = 2.4 + 1.6/(1.0 + xpr^2/(100bs^2))

    bb = sqrt(ar/hft)
    bb = max(bb, bs) # if (bb .lt. bs)  bb = bs

    b = .9*sqrt(am/hft)
    bsm = .9*bs
    b = max(b, bsm) # if (b .lt. bsm)  b = bsm

    h = ar/(bb+bb)
    hh = 0.5*h
    htp = zcr + hh

    alfg = 0.0
    if hh > zcr
       bb = h*bb/htp
       b = h*b/htp
       h = htp
       alfg = .25
    end

    htp = min(htp, hmx) # if (htp .gt. hmx) htp=hmx

    if h > hmx
       bb = h*bb/hmx
       b = h*b/hmx
       h = hmx
    end

    bb0 = bb
    bbv = bb
    bbv0 = bb
    h0 = h
    htp0 = htp
    zc = zcr
    zcmx = hmx - 0.5*h
    zcmx = min(zc, zcmx) # if (zc .gt. zcmx) zc = zcmx
    zc0 = zc

    b0 = b
    bv = b
    bv0 = b

    rho0 = rho
    t0 = t
    hs0 = 0.0
    if ws > u
       hs0 = bs+bs
    else
       hs0 = (bs+bs)*ws/u
    end

    # initializing variables that are called before they are set
    qint = 0.5*qtcs
    qint0 = qint

    for i in 1:msfm
        vecs.uab[i] = uab
        vecs.x[i] = 1.0 + (i - 1.0)*dxpr
        xir = (vecs.x[i] - 1.0)/xpr
        xiz = (vecs.x[i] - 1.0 - xpr)/xpr
        
        vecs.u[i] = u*sqrt(xir)
        vecs.wc[i] = ws*(u - vecs.u[i])/u

        vecs.h[i] = hs0 + (h-hs0)*xir
        vecs.zc[i] = hs + hpr*sqrt(1.0 - xiz*xiz)
        zcmx = hmx - 0.5*vecs.h[i]
        vecs.zc[i] = min(vecs.zc[i], zcmx) # if (zcp(i) .gt. zcmx) zcp(i) = zcmx
        
        vecs.bb[i] = bs + (bb-bs)*xir
        vecs.b[i] = bsm + (b-bsm)*xir
        vecs.beta[i] = sqrt(vecs.bb[i]^2 - vecs.b[i]^2)/√(3)

        vecs.tccp[i] = 2.0*(vecs.x[i] - 1.0)/u
        vecs.xccp[i] = vecs.x[i]

        if vecs.tccp[i] < tsd
            qint = .5*qs*vecs.tccp[i]
            vecs.qint[i] = qint
            qint0 = qint
        else
            nxtr = min(nxtr,i)
            qint = 0.5*qtcs
            qint0 = qint
            vecs.tim[i] = vecs.tccp[i]
            vecs.bbx[i] = 0.5*u*tsd
            vecs.bx[i] = .9999*vecs.bbx[i]
            vecs.betax[i] = sqrt(vecs.bbx[i]^2 - vecs.bx[i]^2)/√(3)
        end

        rtcm = (1.0 - cm)/cm
        vecs.cm[i] = 1.0/(1.0 + rtcm*xir^2)
        vecs.cv[i] = _met.wmae*vecs.cm[i]/(wms + (_met.wmae - wms)*vecs.cm[i])
        vecs.cmda[i] = (1.0 - vecs.cm[i])*_met.cmdaa
        vecs.cmw[i] = (1.0 - vecs.cm[i])*_met.cmwa
        vecs.rho[i] = -1.0
        vecs.t[i] = -1.0
        vecs.ug[i] = 0.0
        vecs.vg[i] = 0.0
        vecs.cmev[i] = -1.0
        vecs.cmwv[i] = -1.0
        vecs.w[i] = -1.0
        vecs.v[i] = -1.0
        vecs.vx[i] = 0.0
    end

    vecs.rho[msfm] = rho
    vecs.t[msfm] = t
    vecs.cmwv[msfm] = cmwv
    vecs.cmev[msfm] = cmev

    cm0 = cm
    cmev0 = cmev
    cmw0 = cmw
    cmwv0 = cmwv
    cp0 = _cp

    mnfm = mffm-msfm
    x = vecs.x[msfm]
    # wse = 0.0 this doesn't appear to be used anywhere

    vg = 0.0
    vg0 = 0.0
    wc = 0.0
    wc0 = 0.0
    ug = 0.0
    ug0 = 0.0

    fu = 0.0
    fv = 0.0
    fw = 0.0
    ft = 0.0
    fug = 0.0

    r0 = r
    srug = 0.5*alfg*grav*(rho - rhoa)*bb*h^2
    sru0 = r*(u - (1.0 - cm)*uab) + srug
    delu = (rhoa/rho)*abs(uab - u)
    if htp > h
       ubs20 = uastr^2
       ustr2 = uastr^2 + cf1*delu^2
    else
       umgs = (uastr/uab)*u
       ubs20 = umgs^2
       ustr2 = ubs20 + cf1*delu^2
    end

    uah = _slab_uafn(htp,_wp)
    w = (urf/uah)*aa*xk*sqrt(ustr2)
    vecs.w[msfm] = w

    if ala < 0.0
       stby = 1.0 - 10rcf*ala
    else
       stby = 1.0/(1.0 + 10rcf*ala)
    end
    atot = afa*stby
    rab = 0.5*sigb/atot
    va = atot*uab/(1.0 + rab*bb/√(3))
    vjp = aa*xk*delu
    v = sqrt(va^2 + cf1*vjp^2)
    vecs.v[msfm] = v

    # c  check for entry to plume or puff mode
    if qint < 0.5*qtcs
    # c  go to nearfield region calculation
    #        go to 600
        nxi = msfm + 1
        dt = 1.0 # this does nothing

        angam = log10((xffm/vecs.x[msfm]) - 1) + 1
        nstp = nssm*mnfm
        anstp = nstp
        gam = 10^(angam/anstp)

        # initialize some variables that are (at this point) undefined
        bx, bbx, bbvx0, bvx0 = zeros(F,4)
        vars = SLAB_Loop_Init(nxi,msfm,mnfm,mffm,gam,ft,fu,fv,fw,fug,bbv0,bv0,r0,cp0,alfg,sru0,
                            htp0,ubs20,rmi,bx,bbx,bbvx0,bvx0,xcc0,bxs0)

        return idpf,nxtr,vecs,vars,params,dt
    else
        idpf = 2
        nxi = msfm + 1
        rho0 = rho
        r = 0.25*qs*tsd/cm
        r0 = r
        rmi = 0.0
        bbx = r/(rho*bb*h)
        bbx0 = bbx
        bbvx = bbx
        bbvx0 = bbx
        bx = .9999*bbx
        bx0 = bx
        bvx = bx
        bvx0 = bx
        sru0 = r*(u - (1.0 - cm)*uab)
        fug = 0.0
        ug = 0.0
        ug0 = 0.0
        cm0 = cm
        cmw0 = cmw
        cmwv0 = cmwv
        cmev0 = cmev
        cp0 = _cp
        t0 = t

        if xffm < (2.0*vecs.x[msfm])
            xffm = 2.0*vecs.v[msfm]
        end
        _fld = SLAB_Field_Params(tav,hmx,xffm,tffm,zp)

        ngam = log10((xffm/vecs.v[msfm]) - 1.0) + 1.0
        nstp = nssm*mnfm
        anstp = nstp
        gam = 10.0^(ngam/anstp)
        dx = (gam - 1.0)*(xffm - vecs.v[msfm])/((gam^nstp) - 1.0)
        dt = dx/u

    # c  go to transient puff dispersion mode
    #        go to 730
        vars = SLAB_Loop_Init(nxi,msfm,mnfm,mffm,gam,ft,fu,fv,fw,fug,bbv0,bv0,r0,cp0,alfg,sru0,
                              htp0,ubs20,rmi,bx,bbx,bbvx0,bvx0,xcc0,bxs0)

        params = SLAB_Params(_rgps,_spl,_fld,_met,_aps,_wp,_othr)

        return idpf,nxtr,vecs,vars,params,dt
    end
else
    #c momentum or buoyant jet plume rise
    hprb = 1.7254*grav*(rhoa-rho)*ws*bs*bs/(rhoa*uahs*uastr*uastr)
    h35 = hprb^0.6
    
    for i in 1:3
        hprt = h35*((hs + hprb)^0.4)
        hprb = (hs + 0.6*hprb)*hprt/(hs + hprb - 0.4*hprt)
    end
    
    hpr = sqrt(hprb*hprb + hprm*hprm)
    hs = hs + hpr
    hs = min(hs, hmx - bs) # if (hs .gt. (hmx-bs))  hs = hmx-bs

    us = (uastr/xk)*_slab_uafn(hs,_wp)
    ws = 0.0
    rho = rhoa
    rhos = rhoa # need to update rhos in _rgps
    as = qs/(rho*us)
    alfg = 0.0
    bb = 0.5*sqrt(as)
    bb = max(bb, bs) # if (bb .lt. bs) bb = bs
    h = as/(bb+bb)
    hh = .5*h
    htp = hs+hh
    
    if (hs < hh)
        alfg = .25
        bb = h*bb/htp
        h = htp
    end

    htp = min(htp, hmx) # if (htp .gt. hmx) htp=hmx
         
    if (h > hmx)
        bb = h*bb/hmx
        h = hmx
    end

    # intialize main storage arrays
    vecs = SLAB_Vecs(F,mffm)

    ts = (wms/_met.wmae)*ta
    t0 = ts
    vecs.t[1] = t0

    # this is deranged, it redefines basic material properties
    # I don't understand what this is supposed to be doing.
    # tbp = 0.99*min(ts,ta) # wtf? redefining the boiling point!?
    # spa = spb/(tbp+spc) # !?
    # spaw = spbw/tbp # wtf? this is a material constant!

    vecs.cm[1] = 1.0
    vecs.cv[1] = 1.0
    vecs.cmw[1] = 0.0
    vecs.cmwv[1] = 0.0
    vecs.cmev[1] = 1.0
    vecs.cmda[1] = 0.0
    cp0 = (_met.wmae/wms)*cpa

    zc = hs
    zcmx = hmx - 0.5*h
    zc = min(zc, zcmx) # if (zc .gt. zcmx) zc = zcmx
    vecs.zc[1] = zc

    tgon = 0.

    _spl = SLAB_Spill_Chars(idspl,qs,tsd,qtcs,qtis,as,ws,bs,hs,us)
       
    # go to horizontal jet - line 400
    msfm,mnfm,ft,fu,fv,fw,fug,bbv0,bv0,r0,sru0,htp0,ubs20,bx,bbx,bbvx0,bvx0 = _slab_init_hjet_400!(vecs,_met,_wp,_othr,_spl,mffm,bb,rho,rhoa,h,htp,alfg)

    # go to near field region calculation - line 600
    nxi = msfm + 1
    dt = 1.0 # this does nothing

    angam = log10((xffm/vecs.x[msfm]) - 1) + 1
    nstp = nssm*mnfm
    anstp = nstp
    gam = 10^(angam/anstp)

    vars = SLAB_Loop_Init(nxi,msfm,mnfm,mffm,gam,ft,fu,fv,fw,fug,bbv0,bv0,r0,cp0,alfg,sru0,
                          htp0,ubs20,rmi,bx,bbx,bbvx0,bvx0,xcc0,bxs0)

    params = SLAB_Params(_rgps,_spl,_fld,_met,_aps,_wp,_othr)

    return idpf,nxtr,vecs,vars,params,dt

end


end
