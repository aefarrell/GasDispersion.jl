function _slab_sub_thermo(p::SLAB_Params,idpf,xn,timn,rmi,t0,cmev0,cm0,cmw0,cmwv0,cp0,
                          r,r0,sft)
    # unpacking parameters
    qs = p.spl.qs
    tsd = p.spl.tsd
    bse = p.ocs.bse
    wmae = p.met.wmae
    wms = p.rgp.wms
    dhe = p.rgp.dhe
    tbp = p.rgp.tbp
    cpsl = p.rgp.cpsl
    cps = p.rgp.cps
    cmwa = p.ocs.cmwa
    ta = p.met.ta
    ts = p.rgp.ts
    cpaa = p.met.cpaa
    spa = p.rgp.spa
    spb = p.rgp.spb
    spc = p.rgp.spc
    rhoa = p.met.rhoa
    rhosl = p.rgp.rhosl
    
    #subroutine thermo
    #c  calculates the thermodynamic properties
    if idpf == 2
        if timn > tsd
            sqs = 0.25*qs*tsd
        else
            sqs = 0.25*qs*timn
        end
        cm = (rmi+sqs)/r
    else
        if xn > bse
            sqs = 0.5*qs
        else
            sqs = 0.25*qs*(bse+xn)/bse
        end
        cm = sqs/r
    end

    cv = (wmae*cm)/(wms+(wmae-wms)*cm)

    cmw = (1 - cm)*cmwa
    cmda = (1 - cm)*(1 - cmwa)
    cmwvt = cmw + (r0/r)*(cmwv0-cmw0)
    cmevt = cm + (r0/r)*(cmev0-cm0)
    cmwv = cmwvt
    cmev = cmevt
    cmwd = cmw-cmwv
    cmed = cm-cmev

    dhe0 = dhe + tbp*(cpsl-cps)
    dhw0 = dhw + 298.2*(cpwl-cpwv)

    _cp = cmda*cpa + cmwv*cpwv + cmwd*cpwl + cmev*cps + cmed*cpsl

    if sft == 0
        sfte = 0.0
    else
        sfte = -cp0*(t0-ta)*abs(sft)/(abs(sft) + r0*cp0*abs(t0-ta))
    end

    etrn = (1 - cm)*cpaa*ta + cm*cps*ts + (r0/r)*(cp0*t0-(1 - cm0)*cpaa*ta-cm0*cps*ts + sfte)
    t = etrn/_cp

    alfm = wmae*(cmda/wma + cmwv/wmw + cmev/wms)
    fft = exp(spa-spb/(t+spc))
    ggt = exp(spaw-spbw/t)
    pw = (wmae*cmw)/(alfm*wmw)
    pe = (wmae*cm)/(alfm*wms)

    ie = 0
    iw = 0
    if cmev0 == cm0
        if fft > pe
            ie = 1
        end
    end

    if cmwv0 == cmw0
        if ggt > pw
            iw = 1
        end
    end

    if (ie+iw) == 2
        #go to 70
    else
        if ie == 1
            #go to 5
        else
            ffpr = (spb*fft)/((t+spc)*(t+spc))

            if fft > 1
                fft = 1 + spa*(t-tbp)/(tbp+spc)
                ffpr = spa/(tbp+spc)
            end

            alfs = wms*(cmda/wma + cmwv/wmw)
            ab2 = 0.5*((1 - fft)/ffpr + (dhe0/_cp)*(alfs+cmevt))
            ac = dhe0*(alfs - (1 - fft)*cmevt)/(_cp*ffpr)
            dlti = ab2 - sqrt(ab2*ab2 + ac)
            dltm = dhe0*(cmevt-cm)/_cp
            dlti = max(dltm,dlti)
            t = t + dlti
            fft = exp(spa - spb/(t+spc))
            ggt = exp(spaw - spbw/t)

        end

        #5

        #do 65 irep = 1,2
        for irep in 1:2

            tstr = t
            cmwvs = cmwv
            cmevs = cmev
            icyc = 0

            #10 do 30 i=1,15
            for i in 1:15
                if ie == 0
                    cmwvm = min(cmwvs,cmw)
                    alfs = wms*(cmda/wma + cmwvm/wmw)
                    fftt = cm/(cm+alfs)

                    if fft < fftt
                        cmevs = alfs*fft/(1 - fft)
                        cmevspr = spb*cmevs/((1 - fft)*(tstr+spc)*(tstr+spc))
                    else
                        tesrt = spb/(spa-log(fftt)) - spc
                        cmevspr = spb*cm/((1 - fftt)*(tesrt+spc)*(tesrt+spc))
                        cmevs = cmevspr*(tstr-tesrt) + cm
                    end
                else
                    cmevs = cm
                    cmevspr = 0.0
                end
            
                #20
                if iw == 0
                    cmevm = min(cmevs,cm)
                    alfw = wmw*(cmda/wma + cmevm/wms)
                    ggtt = cmw/(cmw+alfw)
                    
                    if ggt < ggtt
                        cmwvs = alfw*ggt/(1 - ggt)
                        cmwvspr = spbw*cmwvs/((1 - ggt)*tstr*tstr)
                    else
                        twsrt = spbw/(spaw-log(ggtt))
                        cmwvspr = spbw*cmw/((1 - ggtt)*twsrt*twsrt)
                        cmwvs = cmwvspr*(tstr-twsrt) + cmw
                    end
                else
                    cmwvs = cmw
                    cmwvspr = 0.0
                end
            
                cmwvm = min(cmwvs,cmw)
                cmevm = min(cmevs,cm)
            
                cmwds = cmw-cmwvm
                cmeds = cm-cmevm
                cpstr = cmda*cpa + cmwvm*cpwv + cmwds*cpwl + cmevm*cps + cmeds*cpsl
                tstra = (etrn + dhw0*(cmwvt-cmwvs) + dhe0*(cmevt-cmevs))/cpstr
            
                gfnc = tstr - tstra
                dhwt = dhw0 + tstra*(cpwv-cpwl)
                dhet = dhe0 + tstra*(cps-cpsl)
                gfncp = 1 + (dhwt*cmwvspr + dhet*cmevspr)/cpstr
                tstr1 = tstr - gfnc/gfncp
            
                dtstr = abs(tstr1-tstra)
                tstr = tstr1
            
                fft = exp(spa-spb/(tstr+spc))
                ggt = exp(spaw-spbw/tstr)
            
                if dtstr < 0.001
                    #go to 40
                    break
                end
            #30 continue
            end
            
            #40
            cmwvm = min(cmwvs,cmw)
            cmevm = min(cmevs,cm)
            alfm = wmae*(cmda/wma+cmwvm/wmw+cmevm/wms)

            if ie == 1
                cmev = cm
            else
                cmevs = alfm*wms*fft/wmae
                if cmevs > cm
                    cmev = cm
                    ie = 1
                    
                    if iw == 0
                        icyc = 1
                    end
                else
                    cmev = cmevs
                end
            end

            if iw == 1
                cmwv = cmw
            else
                cmwvs = alfm*wmw*ggt/wmae
                if cmwvs > cmw
                    cmwv = cmw
                    iw = 1
                    if ie == 0
                        icyc = 1
                    end
                else
                    cmwv = cmwvs
                end
            end

            cmwd = cmw-cmwv
            cmed = cm-cmev
            alfm = (cmda/wma + cmwv/wmw + cmev/wms)*wmae
            _cp = cmda*cpa + cmwv*cpwv + cmwd*cpwl + cmev*cps + cmed*cpsl
            t = (etrn+dhw0*(cmwvt-cmwv) + dhe0*(cmevt-cmev))/_cp

            if ie == 1
                fft = exp(spa-spb/(t+spc))
                cmevsc = (alfm*wms*fft)/wmae
                if cmevsc < cm
                    ie = 0
                    icyc = 1
                end
            end

            if iw == 1
                ggt = exp(spaw-spbw/t)
                cmwvsc = (alfm*wmw*ggt)/wmae
                if cmwvsc < cmw
                    iw = 0
                    icyc = 1
                end
            end

            if icyc == 0
                #go to 70
                break
            end

        #65 continue
        end
    end
    
    # 70
    betm = (rhoa/rhowl)*cmwd + (rhoa/rhosl)*cmed
    rho = rhoa*ta/(alfm*t + betm*ta)

    #end thermo
    return cm,cv,cmw,cmwv,cmev,t,rho,_cp
end