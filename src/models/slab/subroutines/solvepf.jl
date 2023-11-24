function _slab_sub_solvepf(p::SLAB_Params,dy,x,bbv0,bv0,bbvx0,bvx0,zc0,r0,h)
    # unpacking parameters
    hmx = p.wps.hmx

    #subroutine solvepf

    #c  solves the transient puff conservation equations
    r = r0 + dy[1]
    xn = x + dy[2]
    bbv = bbv0 + dy[3]
    bv = bv0 + dy[4]
    bbvx = bbvx0 + dy[5]
    bvx = bvx0 + dy[6]
    g = dy[7]
    gx = dy[8]
    sft = dy[9]
    sfu = dy[10]
    gw = dy[11]
    zc = zc0 + dy[12]
    zcmx = hmx - .5*h
    if zc > zcmx
        zc = zcmx
    end
    if zc < 0
        zc = 0.0
    end
    sfy = dy[13]
    sfx = dy[14]
    sfz = dy[15]

    return xn,bbv,bv,bbvx,bvx,zc,r,g,gw,gx,sft,sfu,sfx,sfy,sfz
end