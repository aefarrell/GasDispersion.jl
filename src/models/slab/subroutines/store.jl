function _slab_sub_store!(vecs::SLAB_Inst_Cloud,nx::Integer,x,bb,b,vg,cm,t,rho,u,h,cv,
                          beta,w,v,cmdaa,cmw,cmwv,cmev,uab,wc,zc,qint,tim,bbx,bx,betax,
                          ug,vx)
    #subroutine store

    #c  stores the output variables
    vecs.x[nx] = x
    vecs.bb[nx] = bb
    vecs.b[nx] = b
    vecs.vg[nx] = vg
    vecs.cm[nx] = cm
    vecs.t[nx] = t
    vecs.rho[nx] = rho
    vecs.u[nx] = u
    vecs.h[nx] = h
    vecs.cv[nx] = cv
    vecs.beta[nx] = beta
    vecs.w[nx] = w
    vecs.v[nx] = v
    vecs.cmda[nx] = (1 - cm)*cmdaa
    vecs.cmw[nx] = cmw
    vecs.cmwv[nx] = cmwv
    vecs.cmev[nx] = cmev
    vecs.uab[nx] = uab
    vecs.wc[nx] = wc
    vecs.zc[nx] = zc
    vecs.qint[nx] = qint
    vecs.tim[nx] = tim
    vecs.bbx[nx] = bbx
    vecs.bx[nx] = bx
    vecs.betax[nx] = betax
    vecs.ug[nx] = ug
    vecs.vx[nx] = vx
    
    #end store
end