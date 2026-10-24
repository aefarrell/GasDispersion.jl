# this is a little script to log 
# the internal state of the slab thermo subroutine
# because it is large and hard to unit test as a big blob of code
# this only exists for testing purposes and should be removed
# along with trace.jl and the `checkpoint!` function calls

using GasDispersion

import GasDispersion.slab: TraceLog, with_trace

trace_log = TraceLog("test/test_data/slab_inpr2_thermo_data.jls")

inp = GasDispersion.SLAB_Input(idspl =  2, ncalc =  1, wms   =  0.017031, cps   =  2045.90, tbp   =  239.57, cmed0 =  0.81, dhe   =  1170000.0, cpsl  =  4611.80, rhosl =  603.00, spb   =  2976.01, spc   =  0.00, ts    =  239.57, qs    =  107.87, as    =  0.93, tsd   =  381.0, qtis  =  0.00, hs    =  1.00, tav   =  10.00, xffm  =  2800.00, zp    =  [0.00, 1.00, 0.00, 0.00], z0    =  0.003000, za    =  2.00, ua    =  4.50, ta    =  306.20, rh    =  21.30, stab  =  0.00, ala   =  0.0221)

with_trace(trace_log) do 
    GasDispersion.slab.slab_main(inp)    
end

close(trace_log)