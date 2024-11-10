# routines common to initializing all release types

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