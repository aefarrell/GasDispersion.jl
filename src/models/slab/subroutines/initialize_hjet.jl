function _slab_init_hjet_400!(vecs::SLAB_Vecs, _met::SLAB_Ambient_Met_Props, _wp::SLAB_Wind_Profile, 
                              _othr::SLAB_Other_Params, _spl::SLAB_Spill_Chars,
                              mffm, bb, rho, rhoa, h, htp, alfg)

# c ======================================
# c  horizontal jet source initialization
# c ======================================
    msfm = 1
    mnfm = mffm - 1
    
    vecs.x[1] = 1.0
    vecs.xccp[1] = 1.0
    vecs.tccp[1] = 0.0
    vecs.tim[1] = _spl.tsd
    # wse = 0.0 this doesn't appear to be used anywhere

    vecs.vg[1] = 0.0

    vecs.wc[1] = 0.0
    vecs.ug[1] = 0.0

    fu = 0.0
    fv = 0.0
    fw = 0.0
    ft = 0.0
    fug = 0.0
    vecs.qint[1] = 0.0

    vecs.bb[1] = bb
    bbv0 = bb
    b = 0.9*bb
    vecs.b[1] = b
    bv0 = b
    bx = 0.0
    bbx = 0.0
    bbvx0 = 0.0
    bvx0 = 0.0

    beta = sqrt(bb*bb-b*b)/√(3)
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

    u = _spl.us
    vecs.u[1] = u

    r0 = rho*u*bb*h
    srug = .5*alfg*grav*(rho-rhoa)*bb*h*h
    sru0 = r0*u + srug

    delu = (rhoa/rho)*abs(uab-_spl.us)
    if htp > h
        ubs20 = uastr^2
        ustr2 = uastr^2 + cf1*delu*delu
    else
        umgs = (uastr/uab)*_spl.us
        ubs20 = umgs*umgs
        ustr2 = ubs20 + cf1*delu*delu
    end
    uah = _slab_uafn(htp,_wp)
    w = (_othr.urf/uah)* aa*xk*sqrt(ustr2)
    vecs.w[1] = w

    if _met.ala < 0
        stby = 1 - _othr.rcf*10*_met.ala
    else
        stby = 1/(1 + _othr.rcf*10*_met.ala)
    end

    atot = _othr.afa*stby
    rab = 0.5*sigb/atot
    va = atot*uab/(1.0 + rab*bb/√(3))
    vjp = aa*xk*delu
    v = sqrt(va*va + cf1*vjp*vjp)
    vecs.v[1] = v
    vecs.vx[1] = 0.0
    
    return (msfm,mnfm,ft,fu,fv,fw,fug,bbv0,bv0,r0,sru0,
            htp0,ubs20,bx,bbx,bbvx0,bvx0)
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

if xffm < 2.0
    xffm = 2.0
end

_fld = SLAB_Field_Params(tav,hmx,xffm,tffm,zp)

# initialize additional parameters
nssm = 3*ncalc
_aps = SLAB_Additional_Params(ncalc,nssm,grav,rr,xk)

# initialize other constants
rmi = 0.0
tgon = 1.0
bse = 1.0
hrf = 4.0
urf = _slab_uafn(hrf,_wp)
cf0 = uastr/uaz2
rcf = sqrt(cf0/cf00)
tau0 = 10.0
at0 = 0.0
afa = 0.08*(((at0+tau0*exp(-at0/tau0))/tav0)^0.2)
_othr = SLAB_Other_Params(tgon,bse,hrf,urf,cf0,rcf,tau0,at0,afa)

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

htp = min(htp, hmx) # if (htp .gt. hmx) htp=hmx

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
zc = min(zc, zcmx) # if (zc .gt. zcmx) zc = zcmx
vecs.zc[1] = zc

xcc0 = 1.0
bxs0 = 0.0
# tim = tsd

# go to 400
msfm,mnfm,ft,fu,fv,fw,fug,bbv0,bv0,r0,sru0,htp0,ubs20,bx,bbx,bbvx0,bvx0 = _slab_init_hjet_400!(vecs,_met,_wp,_othr,_spl,mffm,bb,rho,rhoa,h,htp,alfg)

# go to 600
# additional constants to initialize the integration
nxi = msfm + 1
nxtr = mffm + 1
idpf = 0
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