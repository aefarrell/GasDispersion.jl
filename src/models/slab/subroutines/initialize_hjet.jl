function _slab_init_hjet_400!(vecs::SLAB_Vecs, _met::SLAB_Ambient_Met_Props, _wp::SLAB_Wind_Profile, 
                              _othr::SLAB_Other_Params, _spl::SLAB_Spill_Chars, mffm::I, bb::F, 
                              rho::F, rhoa::F, h::F, htp::F, alfg::F, t0::F, cm0::F, cv0::F, cmw0::F, 
                              cmwv0::F, cmev0::F, zc::F) where { I <: Integer, F <: AbstractFloat}

# c ======================================
# c  horizontal jet source initialization
# c ======================================
    msfm = 1
    mnfm = mffm - 1
    
    x = one(F)
    # vecs.x[1] = 1.0 # handled by the store! call
    # vecs.xccp[1] = 1.0 # handled by the store! call
    vecs.tccp[1] = zero(F)
    # vecs.tim[1] = _spl.tsd # handled by the store! call
    tim = _spl.tsd
    # wse = 0.0 this doesn't appear to be used anywhere

    vg = zero(F)
    # vecs.vg[1] = 0.0 # handled by the store! call

    wc = zero(F)
    # vecs.wc[1] = 0.0 # handled by the store! call

    ug = zero(F)
    # vecs.ug[1] = 0.0 # handled by the store! call

    fu = zero(F)
    fv = zero(F)
    fw = zero(F)
    ft = zero(F)
    fug = zero(F)
    qint = zero(F)
    # vecs.qint[1] = 0.0 # handled by the store! call

    # vecs.bb[1] = bb # handled by the store! call
    bbv0 = bb
    b = 0.9*bb
    # vecs.b[1] = b # handled by the store! call
    bv0 = b
    bx = zero(F)
    bbx = zero(F)
    bbvx0 = zero(F)
    bvx0 = zero(F)

    beta = sqrt(bb^2 - b^2)/√(3)
    # vecs.beta[1] = beta # handled by the store! call
    # vecs.rho[1] = rho # handled by the store! call

    # vecs.h[1] = h # handled by the store! call
    htp0 = htp

    uab = zero(F)
    dh = htp/5.0
    h1 = -0.5*dh
    h2 = zero(F)

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
    # vecs.uab[1] = uab # handled by the store! call

    u = _spl.us
    # vecs.u[1] = u # handled by the store! call

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
    # vecs.w[1] = w # handled by the store! call

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
    # vecs.v[1] = v # handled by the store! call

    vx = zero(F)
    # vecs.vx[1] = 0.0 # handled by the store! call

    cmdaa, betax = zeros(F,2)
    _slab_sub_store!(vecs,1,x,bb,b,vg,cm0,t0,rho,u,h,cv0,
        beta,w,v,cmdaa,cmw0,cmwv0,cmev0,uab,wc,zc,qint,tim,bbx,bx,betax,
        ug,vx)
    
    # this stores the initial state for the integration step to retrieve
    # the fortran code has all the program state variables as globals
    # to avoid that, they are stored at the starting position of the integration
    _slab_sub_store!(vecs,2,x,bb,b,vg,cm0,t0,rho,u,h,cv0,
        beta,w,v,cmdaa,cmw0,cmwv0,cmev0,uab,wc,zc,qint,tim,bbx,bx,betax,
        ug,vx)
    
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
qtis = zero(F)
ws = zero(F)
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
rmi = zero(F)
tgon = one(F)
bse = one(F)
hrf = 4.0
urf = _slab_uafn(hrf,_wp)
cf0 = uastr/uaz2
rcf = sqrt(cf0/cf00)
tau0 = 10.0
at0 = zero(F)
afa = 0.08*(((at0+tau0*exp(-at0/tau0))/tav0)^0.2)
_othr = SLAB_Other_Params(tgon,bse,hrf,urf,cf0,rcf,tau0,at0,afa)

# intialize main storage arrays
vecs = SLAB_Vecs(F,mffm)

alfg = zero(F)
bb = bs
h = 2bb
htp = hs+0.5h
if hs <= 0.5h
    alfg = 0.25
    bb = 0.5*(sqrt(hs^2 + 2as) - hs)
    h = bb + hs
    htp = h
end

htp = min(htp, hmx) # if (htp .gt. hmx) htp=hmx

if h > hmx
    bb = h*bb/hmx
    h = hmx
end

t0 = ts
# vecs.t[1] = t0 # passed to line 400, where initial state is stored
cm0 = one(F)
# vecs.cm[1] = cm0 # passed to line 400, where initial state is stored
cv0 = one(F)
# vecs.cv[1] = cv0 # passed to line 400, where initial state is stored
cmw0 = zero(F)
# vecs.cmw[1] = cmw0 # passed to line 400, where initial state is stored
cmwv0 = zero(F)
# vecs.cmwv[1] = cmwv0 # passed to line 400, where initial state is stored

cmev0 = 1 - cmed0
# vecs.cmev[1] = cmev0 # passed to line 400, where initial state is stored
# vecs.cmda[1] = 0.0
cp0 = cmev0*cps + cmed0*cpsl

zc = hs
zcmx = hmx - 0.5*h
zc = min(zc, zcmx) # if (zc .gt. zcmx) zc = zcmx
# vecs.zc[1] = zc # passed to line 400, where initial state is stored

xcc0 = one(F)
bxs0 = zero(F)
# tim = tsd

# go to 400
msfm,mnfm,ft,fu,fv,fw,fug,bbv0,bv0,r0,sru0,htp0,ubs20,bx,bbx,bbvx0,bvx0 = _slab_init_hjet_400!(vecs,_met,_wp,_othr,_spl,mffm,bb,rho,rhoa,h,htp,alfg,t0,cm0,cv0,cmw0,cmwv0,cmev0,zc)
vecs.cmda[1] = zero(F)

# go to 600
# additional constants to initialize the integration
nxi = msfm + 1
nxtr = mffm + 1
idpf = zero(I)
dt = one(F) # this does nothing

angam = log10((xffm/vecs.x[msfm]) - 1) + 1
nstp = nssm*mnfm
anstp = nstp
gam = 10^(angam/anstp)


vars = SLAB_Loop_Init(nxi,msfm,mnfm,mffm,gam,ft,fu,fv,fw,fug,bbv0,bv0,r0,cp0,alfg,sru0,
                      htp0,ubs20,rmi,bx,bbx,bbvx0,bvx0,xcc0,bxs0)

params = SLAB_Params(_rgps,_spl,_fld,_met,_aps,_wp,_othr)

return idpf,nxtr,vecs,vars,params,dt


end