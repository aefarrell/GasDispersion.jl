
function _slab_int_steady_state!(vecs::SLAB_Vecs{F,A},vars::SLAB_Loop_Init{I,F},
                                  params::SLAB_Params{I,F,A},idpf::I,nxtr::I) where {
                                  I <: Integer, F <: AbstractFloat, A <: AbstractVector{F}}
                             
    # unpack parameters
    cmdaa = params.met.cmdaa
    qs = params.spl.qs
    rhoa = params.met.rhoa
    tsd = params.spl.tsd
    qtcs = params.spl.qtcs
    tgon = params.othr.tgon
    bse = params.othr.bse
    urf = params.othr.urf
    cf0 = params.othr.cf0
    rcf = params.othr.rcf
    afa = params.othr.afa

    # unpack intial loop variables
    wss = zero(F)
    ft = vars.ft
    fu = vars.fu
    fv = vars.fv
    fw = vars.fw
    fug = vars.fug
    bbv = bbv0 = vars.bbv0
    bv = bv0 = vars.bv0
    r0 = vars.r0
    _cp = cp0 = vars.cp0
    alfg = vars.alfg
    sru0 = vars.sru0
    htp = htp0 = vars.htp0
    ubs2 = ubs20 = vars.ubs20
    rmi = vars.rmi
    bx = bx0 = vars.bx
    bbx = bbx0 = vars.bbx
    bvx = bvx0 = vars.bvx0
    bbvx = bbvx0 = vars.bbvx0
    bxs0 = vars.bxs0
    xcc0 = vars.xcc0

    # initialize other loop variables to zero
    xn = timn = xstr = zero(F)
    r = g = gw = sft = sfu = sfy = sfz = zero(F)
    betax = zero(F)

    # intialize boundaries and step size
    msfm = vars.msfm
    mnfm = vars.mnfm
    mffm = vars.mffm
    nxi = vars.nxi
    gam = vars.gam
    nssm = params.xtra.nssm
    xffm = params.fld.xffm
    nstp = nssm*mnfm
    dx = (gam - 1) * (xffm - vecs.x[msfm])/((gam^nstp) - 1)
    
    # initial state variable location
    n = max(1, nxi)

    # initialize state variables
    x = vecs.x[n]
    zc = zc0 = vecs.zc[n]
    h = h0 = vecs.h[n]
    bb = bb0 = vecs.bb[n]
    b = b0 = vecs.b[n]
    cv = vecs.cv[n]
    rho = rho0 = vecs.rho[n]
    t = t0 = vecs.t[n]
    u = u0 = vecs.u[n]
    uab = uab0 = vecs.uab[n]
    cm = cm0 = vecs.cm[n]
    cmev = cmev0 = vecs.cmev[n]
    cmw = cmw0 = vecs.cmw[n]
    cmwv = cmwv0 = vecs.cmwv[n]
    wc = wc0 = vecs.wc[n]
    vg = vg0 = vecs.vg[n]
    ug = ug0 = vecs.ug[n]
    w = vecs.w[n]
    v = vecs.v[n]
    vx = vecs.vx[n]
    tim = vecs.tim[n]
    beta = vecs.beta[n]
    qint = qint0 = vecs.qint[n]

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
                _slab_sub_slope!(f,params,rho,x,h,v,w,b,bb,vg,u,wc,cm,ft,fu,fv,fw,bse)

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
                bbv,bv,qint,zc,r,g,gw,sft,sfu,sfy,sfz = _slab_sub_solve(params,dy,bbv0,bv0,
                                                         zc0,r0,h,qint0)

                #call thermo
                cm,cv,cmw,cmwv,cmev,t,rho,_cp = _slab_sub_thermo(params,idpf,xn,timn,rmi,
                                                 t0,cmev0,cm0,cmw0,cmwv0,cp0,r,r0,sft,bse)

                #call eval
                u,uab,b,bb,beta,h,zc,vg,vg0,wc,htp = _slab_sub_eval(params,xn,alfg,sru0,
                                                      zc,h0,u0,uab0,b0,bb0,r,r0,bv,bv0,
                                                      bbv,bbv0,rho,rho0,vg0,wc0,cm,htp0,
                                                      b,htp,uab,beta,vg,wc,h,u,bb,sfu,
                                                      sfz,sfy,g,gw,bse)

                #call entran
                w,v,vx,ubs2,fug,ft,fu,fv,fw = _slab_sub_entran(params,idpf,xn,timn,wss,
                                               xstr,ubs20,u,ug,vg,uab,rho,zc,t,h,htp,bb,
                                               bbx,wc,_cp,tgon,bse,urf,rcf,afa)

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

        _slab_sub_store!(vecs,nx,x,bb,b,vg,cm,t,rho,u,h,cv,beta,w,v,cmdaa,cmw,cmwv,
                         cmev,uab,wc,zc,qint,tim,bbx,bx,betax,ug,vx)
        vecs.tccp[nx] = (qint+qint)/qs

        if qint < 0.5*qtcs
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
            vecs.bbx[nx] = bbx
            bbvx = bbx
            bbvx0 = bbx
            bbx = 0.9999*bbx
            bx0 = bbx
            vecs.bx[nx] = bbx
            bvx = bbx
            bvx0 = bbx
            vecs.betax[nx] = sqrt(bbx*bbx-bbx*bbx)/√(3)
            sru0 = r*(u - (1 - cm)*uab)
            fv = bbx*fv
            fu = bbx*fu
            fw = bbx*fw
            fug = bbx*fug
            ft = bbx*ft

            ug = (bb/bbx)*vg
            ug0 = ug
            vecs.ug[nx] = ug

            cm0 = cm
            cmw0 = cmw
            cmwv0 = cmwv
            cmev0 = cmev
            cp0 = _cp
            t0 = t

            tvars = SLAB_Loop_Init(nxi,msfm,mnfm,mffm,gam,ft,fu,fv,fw,fug,bbv0,bv0,r0,
                                   cp0,alfg,sru0,htp0,ubs20,rmi,bx,bbx,bbvx0,bvx0,xcc0,bxs0)

            #c  go to transient puff dispersion mode
            #go to 730
            _slab_int_transient!(vecs,tvars,params,idpf,nxtr,dt)
            break
        end
    end

    #c   steady state calc of timp
    if nxtr > length(vecs.x)
        xptr = (0.5*params.spl.qtcs/qint)*(vecs.x[end]-vecs.x[1]) + vecs.x[1]
        bxtr = bxs0 + xptr - xcc0
        itr = length(vecs.x)
    else
        xptr = vecs.x[nxtr]
        bxtr = vecs.bbx[nxtr]
        itr = nxtr - 1
    end
    
    if itr > 0
        txt = xptr - xcc0 - xcc0
        txb = txt + xptr
        bxr = (bxtr-bxs0)/(xptr-xcc0)
        for i in 1:itr
            vecs.tim[i] = tsd*(vecs.x[i]+txt)/txb
            vecs.bbx[i] = bxs0 + bxr*(vecs.xccp[i]-xcc0)
            vecs.bx[i] = .9999*vecs.bbx[i]
            vecs.betax[i] = sqrt(vecs.bbx[i]*vecs.bbx[i]-vecs.bx[i]*vecs.bx[i])/√(3)
        end
    end

    # this only runs for an evaporating pool release
    # if idspl == 1
    #     for i in 1:5
    #         vecs.tim[i] = vecs.tim[12-i]
    #     end
    # end
end