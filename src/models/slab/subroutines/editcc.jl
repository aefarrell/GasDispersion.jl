# subroutine editcc
# calculates the time averaged volume concentration

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
        tmdr = min(tcld[i],tav)
        afata = 0.08*((tmdr + tau0*exp(-tmdr/tau0))/tav0)^0.2
        rav2 = (afata/afa)^2
        sig0 = asig*(√(1.0 + sigb*(vecs.x[i]-vecs.x[1])) - 1.0)

        rjt = 1.0
        if idspl == 2
            rm3 = ((rhoa*vecs.uab[i]*(1.0 - vecs.cm[i]))/(rhos*us*vecs.cm[i]))^(1.0/3.0)
            rlg = log((1.0 - rm3 + rm3^2)/((1.0 + rm3)*(1.0 + rm3)))
            rit = atan(((2.0*rm3)-1.)/√(3)) + atan(1.0/√(3))
            if rm3 < 0.1
                rjt = 0.0
            else
                rjt = 1.0 - (2.0/(3.0*rm3^3))*(0.5*rlg + √(3)*rit)
            end
        end

        sigm2 = (rav2  - 1.0)*(sig0*rjt)^2
        betac2 = vecs.beta[i]^2 + sigm2
        betacp[i] = √(betac2)
        bbcp[i] = √(vecs.b[i]^2 + 3.0*betac2)

        hhf = 0.5*vecs.h[i]
        if vecs.zc[i] > hhf
            htpp = vecs.zc[i] + hhf
        else
            htpp = vecs.h[i]
        end
        sig2 = ((htpp - vecs.zc[i])^2)/3
        sigp[i] = √(sig2)

        if sigp[i] == zero(F)
            hdsig = √(3)
        else
            hdsig = vecs.h[i]/sigp[i]
        end

        voln = vecs.bbx[i]*bbcp[i]*vecs.h[i]*vecs.cv[i]
        vold = 4.0*√(2π)*vecs.bx[i]*vecs.b[i]*sigp[i]
  
        if vold != zero(F)
            ccp[i] = voln/vold
        else
            ccp[i] = 0.0
        end

    end

    return SLAB_CC_Vecs(vecs.x,ccp,vecs.b,betacp,vecs.zc,sigp,vecs.tccp,vecs.xccp,
                        vecs.bx,vecs.betax,vecs.tim,tcld,bbcp)
end