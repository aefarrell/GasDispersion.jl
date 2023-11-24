__precompile__()

module slab

export SLAB_Input, SLAB_Output
export slab_main

# defining structs, how the data is passed into and out of SLAB
include("structs.jl")

# declare global variables and constants
include("globals.jl")


# code is organized differently than in SLAB, the functions and subroutines
# are defined first, the main program is the function _slab_main
include("functions.jl")
include("subroutines.jl")




"""
    slab_main(inp::SLAB_Input)

to do doc string
"""
slab_main(inp::SLAB_Input) = slab_main(inp.idspl,inp.ncalc,inp.wms,inp.cps,inp.tbp,inp.cmed0,
                                       inp.dhe,inp.cpsl,inp.rhosl,inp.spb,inp.spc,inp.ts,inp.qs,
                                       inp.as,inp.tsd,inp.qtis,inp.hs,inp.tav,inp.xffm,inp.zp,
                                       inp.z0,inp.za,inp.ua,inp.ta,inp.rh,inp.stab,inp.ala)

function slab_main(idspl::I,ncalc::I,wms::F,cps::F,tbp::F,cmed0::F,dhe::F,cpsl::F,rhosl::F,
                   spb::F,spc::F,ts::F,qs::F,as::F,tsd::F,qtis::F,hs::F,tav::F,xffm::F,
                   zp::AbstractVector{F},z0::F,za::F,ua::F,ta::F,rh::F,stab::F,
                   ala::F) where {I <: Integer, F <: AbstractFloat}

    # intialize main storage arrays
    _ic = SLAB_Inst_Cloud(F,61)

    # some program constants
    msfm = 11
    mnfm = 50
    mffm = 61

    @assert z0 ≥ zero(F) "z0 must be positive"

    #stabin = stab assigned but never used
    
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


    #c  ===============
    #c   set constants
    #c  ===============
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

    #c  initial droplet concentration
    cmed = cmed0

    #c  atmospheric and physical constants
    rpwa = 0.01*rh*exp(spaw-spbw/ta)
    cmwa = wmw*rpwa/(wma+(wmw-wma)*rpwa)
    cmdaa = 1-cmwa
    cpaa = cmdaa*cpa + cmwa*cpwv
    wmae = wma*wmw/(wmw+(wma-wmw)*cmwa)

    #c  number of integration steps
    nssm = ncalc + ncalc + ncalc

    #c  number of zp values
    nzpm = 1
    for i in 2:4
        if zp[i] == zero(F)
            break
        else
            nzpm = i
        end
    end
    
    #c   ====================
    #c   calculate parameters
    #c   ====================
    #c  calculate ambient meteorological values
    rhoa = wmae*pa/(rr*ta)

    tau0 = 10.0
    at0 = 0.0
    afa = 0.08*(((at0+tau0*exp(-at0/tau0))/tav0)^0.2)

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
    bs = 0.5*sqrt(as)

    #c  check if source below mixing layer height
    if hs > (hmx-bs)
        @error "input source height hs is greater than the calculated mixing layer height hmx minus the stack half width bs.  program terminated."
    end
    
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
    uaz2 = (uastr/xk)*_slab_uafn(2.,_wp)
    tffm = xffm/uaz2
    cf0 = uastr/uaz2
    rcf = sqrt(cf0/cf00)
    hrf = 4.0
    urf = _slab_uafn(hrf,_wp)

    #c
    #c  *     *    ***    ***  *     *           *****     ***    *         *****
    #c  **   **   *   *    *   **    *          *     *   *   *   *        *     *
    #c  * * * *  *     *   *   * *   *          *        *     *  *        *
    #c  *  *  *  *******   *   *  *  *          *        *******  *        *
    #c  *     *  *     *   *   *   * *          *        *     *  *        *
    #c  *     *  *     *   *   *    **          *     *  *     *  *        *     *
    #c  *     *  *     *  ***  *     *           *****   *     *  *******   *****
    #c
    #c
    #c   main calculation
    
    
    #c ====================================
    #c ====================================
    #c  select source type and write input
    #c ====================================
    #c ====================================

    rhos = wms*pa/(rr*ts)
    qtcs = qs*tsd
    hqtcs = .5*qtcs
    nxtr = mffm + 1
    idpf = 0
    tgon = 1.0

    #c ========================
    #c  horizontal jet release
    #c ========================
    t0 = ts
    cmev0 = 1 - cmed
    alfm = (wmae/wms)*cmev0
    betm = (rhoa/rhosl)*cmed
    rho = rhoa*ta/(alfm*t0+betm*ta)

    us = qs/(rho*as)
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

    ws = 0.0
    qtis = 0.0
    
    t = t0
    _ic.t[1] = t0

    cm0 = 1.0
    cm = cm0
    _ic.cm[1] = cm0

    cv = 1.0
    _ic.cv[1] = cv
    
    cmw0 = 0.0
    cmw = cmw0
    _ic.cmw[1] = cmw0
    
    cmwv0 = 0.0
    cmwv = cmwv0
    _ic.cmwv[1] = cmwv0

    cmev = cmev0
    _ic.cmev[1] = cmev0
    _ic.cmda[1] = 0.0
    cp0 = cmev0*cps + cmed*cpsl
    _cp = cp0

    bse = 1.0
    zc = hs
    zcmx = hmx - 0.5*h

    if zc > zcmx
        zc = zcmx
    end
    zc0 = zc
    _ic.zc[1] = zc
    xcc0 = 1.0
    bxs0 = 0.0
    tim = tsd

    #call editin
    _prms = SLAB_Params(
        SLAB_Release_Gas_Props(wms,cps,ts,rhos,tbp,cmed0,
                               cpsl,dhe,rhosl,spa,spb,spc),
        SLAB_Spill_Chars(idspl,qs,tsd,qtcs,qtis,as,ws,bs,
                         hs,us),
        SLAB_Field_Params(tav,hmx,xffm,zp),
        SLAB_Ambient_Met_Props(wmae,cpaa,rhoa,za,pa,ua,ta,
                               rh,uastr,stab,ala ,z0),
        SLAB_Additional_Params(ncalc,nssm,grav,rr,xk),
        _wp,
        SLAB_Other_Consts(tgon,bse,urf,rcf,afa,stb,phimi,phgam,cmwa)
    )
    
    #go to 400

    #c ======================================
    #c  horizontal jet source initialization
    #c ======================================

    msfm = 1
    mnfm = mffm - 1
    x = 1.0
    _ic.x[1] = 1.0
    xccp = zeros(F,61)
    tccp = zeros(F,61)
    xccp[1] = 1.0
    tccp[1] = 0.0
    _ic.tim[1] = tsd
    wse = 0.0

    vg = 0.0
    vg0 = 0.0
    _ic.vg[1] = 0.0

    wc = 0.0
    wc0 = 0.0
    _ic.wc[1] = 0.0
    ug = 0.0
    ug0 = 0.0
    _ic.ug[1] = 0.0

    fu = 0.0
    fv = 0.0
    fw = 0.0
    ft = 0.0
    qint = 0.0
    qint0 = 0.0

    bb0 = bb
    _ic.bb[1] = bb
    bbv = bb
    bbv0 = bb
    b = 0.9*bb
    b0 = b
    _ic.b[1] = b
    bv = b
    bv0 = b

    beta = sqrt(bb*bb-b*b)/sr3
    _ic.beta[1] = beta
    _ic.rho[1] = rho
    rho0 = rho

    _ic.h[1] = h
    h0 = h
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

    uab = (uastr/xk)*(uab-u2)/30
    uab0 = uab
    _ic.uab[1] = uab

    u = us
    u0 = u
    _ic.u[1] = u

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
    ubs2 = ubs20
    uah = _slab_uafn(htp,_wp)
    w = (urf/uah)* aa*xk*sqrt(ustr2)
    _ic.w[1] = w

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
    _ic.v[1] = v

    vx = 0.0
    _ic.vx[1] = 0.0
    
    #c  go to near field region calculation
    #      go to 600
    #
    #c =====================================
    #c =====================================
    #c  steady state plume dispersion mode
    #c  near field region calculation
    #c =====================================
    #c =====================================
    #
    #c =============================
    #c  pool, jet, or stack release
    #c =============================
    #
    #c ================
    #c  initialization
    #c ================

    nxi = msfm + 1
    wss = 0.0
    
    if xffm < 2*_ic.x[msfm]
        xffm = 2*_ic.x[msfm]
    end

    angam = log10((xffm/_ic.x[msfm]) - 1) + 1
    nstp = nssm*mnfm
    anstp = nstp
    gam = 10^(angam/anstp)
    dx = (gam - 1) * (xffm - _ic.x[msfm])/((gam^nstp) - 1)

    #c ===============================================
    #c  integrate steady state conservation equations
    #c ===============================================
    
    # initializing here due to scope issues later on
    dt,xn,rmi,timn,xstr = zeros(F,5)
    r,g,gw,sft,sfu,sfy,sfz,fug = zeros(F,8)
    betax,bx,bx0,bbx,bbx0,bvx,bvx0,bbvx,bbvx0 = zeros(F,9)

    # initializing arrays for integration
    dxxi = zeros(F,3)
    dxrk = zeros(F,4)
    _sum = zeros(F,11)
    dy = zeros(F,11)
    f = zeros(F,11)
    
    #do 675 nx=nxi,mffm
    for nx in nxi:mffm
        #do 660 ns = 1,nssm
        for ns in 1:nssm
            #do 636 i=1,11
            #636 sum(i) = 0.
            for i in 1:11
                _sum[i] = 0.0
            end

            dxxi[1] = 0.5*dx
            dxxi[2] = 0.5*dx
            dxxi[3] = dx

            dxrk[1] = dx/6
            dxrk[2] = dx/3
            dxrk[3] = dx/3
            dxrk[4] = dx/6

            #do 655 k=1,4
            for k in 1:4
                #call slope
                _slab_sub_slope!(f,_prms,rho,x,h,v,w,b,bb,vg,u,wc,cm,ft,fu,fv,fw)

                #do 650 j=1,11
                for j in 1:11
                    _sum[j] = _sum[j] + dxrk[k] * f[j]
                    if k == 4
                        dy[j] = _sum[j]
                    else
                        dy[j] = dxxi[k] * f[j]
                        xn = x + dxxi[k]
                    end
                #650 continue
                end

                #call solve
                bbv,bv,qint,zc,r,g,gw,sft,sfu,sfy,sfz = _slab_sub_solve(_prms,dy,bbv0,bv0,
                                                         zc0,r0,h,qint0)

                #call thermo
                cm,cv,cmw,cmwv,cmev,t,rho,_cp = _slab_sub_thermo(_prms,idpf,xn,timn,rmi,
                                                 t0,cmev0,cm0,cmw0,cmwv0,cp0,r,r0,sft)

                #call eval
                u,uab,b,bb,beta,h,zc,vg,vg0,wc,htp = _slab_sub_eval(_prms,xn,alfg,sru0,
                                                      zc,h0,u0,uab0,b0,bb0,r,r0,bv,bv0,
                                                      bbv,bbv0,rho,rho0,vg0,wc0,cm,htp0,
                                                      b,htp,uab,beta,vg,wc,h,u,bb,sfu,
                                                      sfz,sfy,g,gw)

                #call entran
                w,v,vx,ubs2,fug,ft,fu,fv,fw = _slab_sub_entran(_prms,idpf,xn,timn,wss,
                                               xstr,ubs20,u,ug,vg,uab,rho,zc,t,h,htp,bb,
                                               bbx,wc,_cp)

            #655 continue
            end

            x = xn
            htp0 = htp
            r0 = r
            bb0 = bb
            b0 = b
            bbv0 = bbv
            bv0 = bv
            rho0 = rho
            h0 = h
            vg0 = vg
            wc0 = wc
            zc0 = zc
            t0 = t
            u0 = u
            uab0 = uab
            ubs20 = ubs2

            if htp > h
                # lfg = 0.0 assigned but never used
                srug = 0.0
            else
                if rho > rhoa
                    alfg = 0.25
                    srug = 0.5*alfg*grav*(rho-rhoa)*bb*h*h
                else
                    alfg = 0.0
                    srug = 0.0
                end
            end

            sru0 = r*u - r*(1 - cm)*uab + srug
            qint0 = qint

            cm0 = cm
            cmw0 = cmw
            cmwv0 = cmwv
            cmev0 = cmev
            cp0 = _cp

            dx = gam*dx

        #660 continue
        end

        xccp[nx] = x
        tccp[nx] = (qint+qint)/qs

        #call store
        _slab_sub_store!(_ic,nx,x,bb,b,vg,cm,t,rho,u,h,cv,beta,w,v,cmdaa,cmw,cmwv,
                         cmev,uab,wc,zc,qint,tim,bbx,bx,betax,ug,vx)

        if qint < hqtcs
            continue
        else
            idpf = 2
            nxtr = nx
            nxi = nx+1
            dt = dx/u
            rho0 = rho
            r = 0.25*qs*tsd/cm
            r0 = r
            rmi = 0.0
            bbx = r/(rho*bb*h)
            bbx0 = bbx
            _ic.bbx[nx] = bbx
            bbvx = bbx
            bbvx0 = bbx
            bbx = 0.9999*bbx
            bx0 = bbx
            _ic.bx[nx] = bbx
            bvx = bbx
            bvx0 = bbx
            _ic.betax[nx] = sqrt(bbx*bbx-bbx*bbx)/sr3
            sru0 = r*(u - (1 - cm)*uab)
            fv = bbx*fv
            fu = bbx*fu
            fw = bbx*fw
            fug = bbx*fug
            ft = bbx*ft

            ug = (bb/bbx)*vg
            ug0 = ug
            _ic.ug[nx] = ug

            cm0 = cm
            cmw0 = cmw
            cmwv0 = cmwv
            cmev0 = cmev
            cp0 = _cp
            t0 = t

            #c  go to transient puff dispersion mode
            #go to 730
            break
        end

    # 675
    end

    #c ================================
    #c ================================
    #c  transient puff dispersion mode
    #c ================================
    #c ================================
    #
    #c ============================================
    #c  pool, jet, stack, or instantaneous release
    #c ============================================
    #
    #
    #c ========================================
    #c   integrate puff conservation equations
    #c ========================================

    # initializing arrays for integration
    dxxi = zeros(F,3)
    dxrk = zeros(F,4)
    _sum = zeros(F,15)
    dy = zeros(F,15)
    f = zeros(F,15)

    #730 do 775 nx=nxi,mffm
    for nx in nxi:mffm

        #do 765 ns=1,nssm
        for ns in 1:nssm
            if tim < tsd
                wss = ws
                xstr = h/cf0
            else
                wss = 0.0
            end

            #do 735 i=1,15
            for i in 1:15
                _sum[i] = 0.0
            
            #735 continue
            end

            dxxi[1] = 0.5*dt
            dxxi[2] = 0.5*dt
            dxxi[3] = dt

            dxrk[1] = dt/6.0
            dxrk[2] = dt/3.0
            dxrk[3] = dt/3.0
            dxrk[4] = dt/6.0

            #do 745 k=1,4
            for k in 1:4
                #call slopepf
                if tim >= tsd
                    rqs = 0.0
                    wss = 0.0
                    atau = 0.0
                else
                    rqs = 0.25*qs
                    wss = ws
                    if r >= 0.0
                        atau = rqs/r
                    else
                        atau = 0.0
                    end
                end
                _slab_sub_slopepf!(f,_prms,rqs,atau,x,u,ug,v,vg,vx,w,wc,b,bx,bb,
                                   bbx,h,rho,fug,ft,fu,fv,fw)

                #do 740 j=1,15
                for j in 1:15

                    _sum[j] = _sum[j] + dxrk[k]*f[j]

                    if k == 4
                        dy[j] = _sum[j]
                    else
                        timn = tim + dxxi[k]
                        dy[j] = dxxi[k]*f[j]
                    end

                #740 continue
                end

                #call solvepf
                xn,bbv,bv,bbvx,bvx,zc,r,g,gw,gx,sft,sfu,sfx,sfy,sfz = _slab_sub_solvepf(_prms,dy,x,bbv0,bv0,bbvx0,bvx0,zc0,r0,h)
                
                #call thermo
                cm,cv,cmw,cmwv,cmev,t,rho,_cp = _slab_sub_thermo(_prms,idpf,xn,timn,rmi,
                                                 t0,cmev0,cm0,cmw0,cmwv0,cp0,r,r0,sft)
                
                #call evalpf
                bb,b,bbx,bx,beta,betax,h,zc,htp,u,uab,vg,vg0,ug,ug0,wc,wc0 = _slab_sub_evalpf(_prms,sru0,zc,h0,rho,rho0,bb0,b0,bbx0,bx0,r,r0,bbv,bbv0,bv,bv0,bbvx,bbvx0,bvx,bvx0,u0,uab0,cm,vg0,ug0,wc0,htp,htp0,sfu,sfx,sfy,sfz,g,gw,gx)

                #call entran
                w,v,vx,ubs2,fug,ft,fu,fv,fw = _slab_sub_entran(_prms,idpf,xn,timn,wss,
                                               xstr,ubs20,u,ug,vg,uab,rho,zc,t,h,htp,bb,
                                               bbx,wc,_cp)


            #745 continue
            end

            x = xn
            tim = timn
            htp0 = htp
            r0 = r
      
            bb0 = bb
            b0 = b
            bbv0 = bbv
            bv0 = bv
            bbx0 = bbx
            bbvx0 = bbvx
            bx0 = bx
            bvx0 = bvx
            rho0 = rho
      
            h0 = h
            vg0 = vg
            ug0 = ug
            wc0 = wc
            zc0 = zc
            t0 = t
            u0 = u
            uab0 = uab
            ubs20 = ubs2
            sru0 = r*(u - (1.0 - cm)*uab)

            cm0 = cm
            cmw0 = cmw
            cmwv0 = cmwv
            cmev0 = cmev
            cp0 = _cp

            dt = gam*dt

        #765 continue
        end

        xccp[nx] = x
        tccp[nx] = tim

        #call store
        _slab_sub_store!(_ic,nx,x,bb,b,vg,cm,t,rho,u,h,cv,beta,w,v,cmdaa,cmw,cmwv,
                         cmev,uab,wc,zc,qint,tim,bbx,bx,betax,ug,vx)

    #775 continue
    end


    #900 continue

    #c =============
    #c =============
    #c  time averaged volume concentration
    #c =============
    #c =============

    #c   steady state calc of timp

    if nxtr > 61
        xptr = (hqtcs/qint)*(_ic.x[61]-_ic.x[1]) + _ic.x[1]
        bxtr = bxs0 + xptr - xcc0
        itr = 61
    else
        xptr = _ic.x[nxtr]
        bxtr = _ic.bbx[nxtr]
        itr = nxtr - 1
    end

    if itr > 0
        txt = xptr - xcc0 - xcc0
        txb = txt + xptr
        bxr = (bxtr-bxs0)/(xptr-xcc0)
        for i in 1:itr
            _ic.tim[i] = tsd*(_ic.x[i]+txt)/txb
            _ic.bbx[i] = bxs0 + bxr*(xccp[i]-xcc0)
            _ic.bx[i] = .9999*_ic.bbx[i]
            _ic.betax[i] = sqrt(_ic.bbx[i]*_ic.bbx[i]-_ic.bx[i]*_ic.bx[i])/sr3
        end
    end

    if idspl == 1
        for i in 1:5
            _ic.tim[i] = _ic.tim[12-i]
        end
    end

    return SLAB_Output(_prms, _ic)
end

end
