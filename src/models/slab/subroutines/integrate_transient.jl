function _slab_int_transient!(vecs::SLAB_Vecs{F,A},vars::SLAB_Transient_Loop_Init{I,F},
                              params::SLAB_Params{I,F,A},idpf::I,nxtr::I) where {
                              I <: Integer, F <: AbstractFloat, A <: AbstractVector{F}}
# fucking SLAB_Output
    # unpack parameters
    cmdaa = params.met.cmdaa
    qs = params.spl.qs
    tsd = params.spl.tsd
    ws = params.spl.ws

    # unpack intial loop variables
    wss = zero(F)
    tgon = vars.tgon
    bse = vars.bse
    urf = vars.urf
    cf0 = vars.cf0
    rcf = vars.rcf
    afa = vars.afa
    ft = vars.ft
    fu = vars.fu
    fv = vars.fv
    fw = vars.fw
    fug = vars.fug
    bbv = bbv0 = vars.bbv0
    bv = bv0 = vars.bv0
    bbvx = bbvx0 = vars.bbvx0
    bvx = bvx0 = vars.bvx0
    r = r0 = vars.r0
    rmi = vars.rmi
    bx = vars.bx
    bbx = vars.bbx
    _cp = cp0 = vars.cp0
    alfg = vars.alfg
    sru0 = vars.sru0
    htp = htp0 = vars.htp0
    ubs2 = ubs20 = vars.ubs20

    # initialize other loop variables to zero
    xn = timn = xstr = zero(F)
    g = gw = sft = sfu = sfy = sfz = zero(F)
    betax = zero(F)

    # initialize boundaries and step size
    nxi = vars.nxi
    msfm = vars.msfm
    mnfm = vars.mnfm
    mffm = vars.mffm
    nssm = params.xtra.nssm
    gam = vars.gam
    dt = vars.dt

    # initialize state variables
    x = vecs.x[nxtr]
    zc = zc0 = vecs.zc[nxtr]
    h = h0 = vecs.h[nxtr]
    bb = bb0 = vecs.bb[nxtr]
    b = b0 = vecs.b[nxtr]
    bbx0 = vecs.bbx[nxtr]
    bx0 = vecs.bx[nxtr]
    cv = vecs.cv[nxtr]
    rho = rho0 = vecs.rho[nxtr]
    t = t0 = vecs.t[nxtr]
    u = u0 = vecs.u[nxtr]
    uab = uab0 = vecs.uab[nxtr]
    cm = cm0 = vecs.cm[nxtr]
    cmev = cmev0 = vecs.cmev[nxtr]
    cmw = cmw0 = vecs.cmw[nxtr]
    cmwv = cmwv0 = vecs.cmwv[nxtr]
    wc = wc0 = vecs.wc[nxtr]
    vg = vg0 = vecs.vg[nxtr]
    ug = ug0 = vecs.ug[nxtr]
    w = vecs.w[nxtr]
    v = vecs.v[nxtr]
    vx = vecs.vx[nxtr]
    tim = vecs.tim[nxtr]
    beta = vecs.beta[nxtr]
    qint = qint0 = vecs.qint[nxtr] 
    
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
                _slab_sub_slopepf!(f,params,rqs,atau,x,u,ug,v,vg,vx,w,wc,b,bx,bb,
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
                xn,bbv,bv,bbvx,bvx,zc,r,g,gw,gx,sft,sfu,sfx,sfy,sfz = _slab_sub_solvepf(params,dy,x,bbv0,bv0,bbvx0,bvx0,zc0,r0,h)
                
                #call thermo
                cm,cv,cmw,cmwv,cmev,t,rho,_cp = _slab_sub_thermo(params,idpf,xn,timn,rmi,
                                                 t0,cmev0,cm0,cmw0,cmwv0,cp0,r,r0,sft,bse)
                
                #call evalpf
                bb,b,bbx,bx,beta,betax,h,zc,htp,u,uab,vg,vg0,ug,ug0,wc,wc0 = _slab_sub_evalpf(params,sru0,zc,h0,rho,rho0,bb0,b0,bbx0,bx0,r,r0,bbv,bbv0,bv,bv0,bbvx,bbvx0,bvx,bvx0,u0,uab0,cm,vg0,ug0,wc0,htp,htp0,sfu,sfx,sfy,sfz,g,gw,gx)

                #call entran
                w,v,vx,ubs2,fug,ft,fu,fv,fw = _slab_sub_entran(params,idpf,xn,timn,wss,
                                               xstr,ubs20,u,ug,vg,uab,rho,zc,t,h,htp,bb,
                                               bbx,wc,_cp,tgon,bse,urf,rcf,afa)


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

        _slab_sub_store!(vecs,nx,x,bb,b,vg,cm,t,rho,u,h,cv,beta,w,v,cmdaa,cmw,cmwv,
                         cmev,uab,wc,zc,qint,tim,bbx,bx,betax,ug,vx)
        vecs.tccp[nx] = tim

    #775 continue
    end

end