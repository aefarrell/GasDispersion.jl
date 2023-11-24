function _slab_sub_solve(p::SLAB_Params,dy,bbv0,bv0,zc0,r0,h,qint0)
    # unpacking parameters
    hmx = p.wps.hmx
   
    #subroutine solve
    #c  solves the steady state plume conservation equations
    r = r0 + dy[1]
    bbv = bbv0 + dy[2]
    bv = bv0 + dy[3]
    g = dy[4]
    sft = dy[5]
    sfu = dy[6]
    gw = dy[7]
    zc = zc0 + dy[8]
    zcmx = hmx - 0.5*h

    if zc > zcmx
        zc = zcmx
    end

    if zc < 0
        zc = 0.0
    end

    qint = qint0 + dy[9]
    sfy = dy[10]
    sfz = dy[11]
    #end solve
    
    return bbv,bv,qint,zc,r,g,gw,sft,sfu,sfy,sfz
end