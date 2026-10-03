# subroutine editcc
# calculates the time averaged volume concentration

function _slab_editcc_point(params::SLAB_Params, x, x0, zc, h, b, beta, uab, cm,
                            cv, bx, bbx, tcld)
    rhos = params.rgp.rhos
    idspl = params.spl.idspl
    us = params.spl.us
    tav = params.fld.tav
    ala = params.met.ala
    rhoa = params.met.rhoa
    tau0 = params.othr.tau0
    rcf = params.othr.rcf
    afa = params.othr.afa

    if ala < 0
        stby = 1 - 10*rcf*ala
    else
        stby = 1/(1 + 10*rcf*ala)
    end
    asig = 2*afa*stby/sigb
    tav = tav == 0 ? 1.0 : tav
    tmdr = min(tcld, tav)
    afata = 0.08*((tmdr + tau0*exp(-tmdr/tau0))/tav0)^0.2
    rav2 = (afata/afa)^2
    sig0 = asig*(sqrt(1 + sigb*(x-x0)) - 1)

    rjt = 1.0
    if idspl == 2
        rm3 = ((rhoa*uab*(1-cm))/(rhos*us*cm))^(1/3)
        rlg = log((1-rm3+rm3^2)/((1+rm3)*(1+rm3)))
        rit = atan((2*rm3-1)/sqrt(3)) + atan(1/sqrt(3))
        if rm3 < 0.1
            rjt = 0.0
        else
            rjt = 1 - (2/(3*rm3^3))*(0.5*rlg + sqrt(3)*rit)
        end
    end

    sigm2 = (rav2-1)*(sig0*rjt)^2
    betac2 = beta^2 + sigm2
    betac = sqrt(betac2)
    bbc = sqrt(b^2 + 3*betac2)
    hhf = 0.5*h
    htpp = zc > hhf ? zc + hhf : h
    sig = sqrt((htpp-zc)^2/3)
    voln = bbx*bbc*h*cv
    vold = 4*sqrt(2*pi)*bx*b*sig
    cc = vold != 0 ? voln/vold : 0.0
    return (; cc,betac,sig,bbc)
end

function editcc(vecs::SLAB_Vecs{F,A}, params::SLAB_Params{I,F,A}, mffm::I) where {I <: Integer, F <: AbstractFloat, A <: AbstractVector}
    
    # unpack parameters
    rhos = params.rgp.rhos
    idspl = params.spl.idspl
    tsd = params.spl.tsd
    us = params.spl.us
    tav = params.fld.tav
    ala = params.met.ala
    rhoa = params.met.rhoa
    tau0 = params.othr.tau0
    rcf = params.othr.rcf
    afa = params.othr.afa

    # tcld - cloud duration
    tcld = zeros(F,mffm)
    
    tcmx = 2.0*vecs.bbx[end]/vecs.u[end]
    tcld[end] = max(tcmx,tsd)
    for i in 2:(mffm-1)
        in = mffm+1-i
        tcmx = 2.0*vecs.bbx[in]/vecs.u[in]
        tcld[in] = max(tsd,min(tcld[in+1],tcmx))
    end
    tcld[1] = tcld[2]

    if ala < zero(F)
        stby = 1.0 - 10.0*rcf*ala
    else
        stby = 1.0/(1.0 + 10.0*rcf*ala)
    end
    asig = 2.0*afa*stby/sigb

    if vecs.u[1] == zero(F)
        vecs.u[1] = 0.001
    end

    if tav == zero(F)
        tav = 1.0
    end

    bbcp = zeros(F,mffm) # effective half width with meander
    betacp = zeros(F,mffm) # effective beta with meander
    sigp = zeros(F,mffm) # dispersion
    ccp = zeros(F,mffm) # centerline concentration
    for i in 1:mffm
        point = _slab_editcc_point(params,vecs.x[i],vecs.x[1],vecs.zc[i],vecs.h[i],
            vecs.b[i],vecs.beta[i],vecs.uab[i],vecs.cm[i],vecs.cv[i],vecs.bx[i],
            vecs.bbx[i],tcld[i])
        ccp[i],betacp[i],sigp[i],bbcp[i] = point.cc,point.betac,point.sig,point.bbc

    end

    return SLAB_CC_Vecs(vecs.x,ccp,vecs.b,betacp,vecs.zc,sigp,vecs.tccp,vecs.xccp,
                        vecs.bx,vecs.bbx,vecs.betax,vecs.tim,tcld,bbcp)
end