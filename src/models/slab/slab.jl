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
                   ala::F;msfm::I=11,mnfm::I=50,mffm::I=61) where {I <: Integer, F <: AbstractFloat}

    #c  number of zp values
    # nzpm = 1
    # for i in 2:4
    #     if zp[i] == zero(F)
    #         break
    #     else
    #         nzpm = i
    #     end
    # end

   # select appropriate release type
    if idspl == 3
        # vertical jet
        idpf,nxtr,vecs,vars,params,dt = _slab_init_vjet(3,ncalc,msfm,mnfm,mffm,wms,cps,tbp,cmed0,
                                                    dhe,cpsl,rhosl,spb,spc,ts,qs,as,tsd,qtis,hs,tav,
                                                    xffm,zp,z0,za,ua,ta,rh,stab,ala)
        if idpf < 2
            _slab_int_steady_state!(vecs,vars,params,idpf,nxtr)
        else
            _slab_int_transient!(vecs,vars,params,idpf,nxtr,dt)
        end
    else
        # default is a horizontal jet
        idpf,nxtr,vecs,vars,params,dt = _slab_init_hjet(2,ncalc,msfm,mnfm,mffm,wms,cps,tbp,cmed0,
                                                    dhe,cpsl,rhosl,spb,spc,ts,qs,as,tsd,qtis,hs,tav,
                                                    xffm,zp,z0,za,ua,ta,rh,stab,ala)
        if idpf < 2
            _slab_int_steady_state!(vecs,vars,params,idpf,nxtr)
        else
            _slab_int_transient!(vecs,vars,params,idpf,nxtr,dt)
        end
    end

    cc_vecs = editcc(vecs,params,mffm)

    return SLAB_Output(params,vecs,cc_vecs)

end

end