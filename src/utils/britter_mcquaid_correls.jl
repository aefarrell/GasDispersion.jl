# plume correlations
function _bm_pl_c_1(α)
    # corresponds to (Cm/C0) = 0.10
    if α ≤ -0.55
        return 1.75
    elseif -0.55 < α ≤ -0.14
        return 0.24*α + 1.88
    elseif -0.14 < α ≤ 1.0
        return -0.5*α + 1.78
    else
        @warn "α= $α is out of range"
        return -0.5*α + 1.78
    end
end

function _bm_pl_c_05(α)
    # corresponds to (Cm/C0) = 0.05
    if α ≤ -0.68
        return 1.92
    elseif -0.68 < α ≤ -0.29
        return 0.36*α + 2.16
    elseif -0.29 < α ≤ -0.18
        return 2.06
    elseif -0.18 < α ≤ 1.0
        return -0.56*α + 1.96
    else
        @warn "α= $α is out of range"
        return -0.56*α + 1.96
    end
end

function _bm_pl_c_02(α)
    # corresponds to (Cm/C0) = 0.02
    if α ≤ -0.69
        return 2.08
    elseif -0.69 < α ≤ -0.31
        return 0.45*α + 2.39
    elseif -0.31 < α ≤ -0.16
        return 2.25
    elseif -0.16 < α ≤ 1.0
        return -0.54*α + 2.16
    else
        @warn "α= $α is out of range"
        return -0.54*α + 2.16
    end
end

function _bm_pl_c_01(α)
    # corresponds to (Cm/C0) = 0.01
    if α ≤ -0.70
        return 2.25
    elseif -0.70 < α ≤ -0.29
        return 0.49*α + 2.59
    elseif -0.29 < α ≤ -0.20
        return 2.45
    elseif -0.20 < α ≤ 1.0
        return -0.52*α + 2.35
    else
        @warn "α= $α is out of range"
        return -0.52*α + 2.35
    end
end

function _bm_pl_c_005(α)
    # corresponds to (Cm/C0) = 0.005
    if α ≤ -0.67
        return 2.40
    elseif -0.67 < α ≤ -0.28
        return 0.59*α + 2.80
    elseif -0.28 < α ≤ -0.15
        return 2.63
    elseif -0.15 < α ≤ 1.0
        return -0.48*α + 2.56
    else
        @warn "α= $α is out of range"
        return -0.48*α + 2.56
    end
end

function _bm_pl_c_002(α)
    # corresponds to (Cm/C0) = 0.002
    if α ≤ -0.69
        return 2.60
    elseif -0.69 < α ≤ -0.25
        return 0.39*α + 2.87
    elseif -0.25 < α ≤ -0.13
        return 2.77
    elseif -0.13 < α ≤ 1.0
        return -0.50*α + 2.71
    else
        @warn "α= $α is out of range"
        return -0.50*α + 2.71
    end
end

"""
    _bm_pl_c(α)

Britter-McQuaid plume correlations

# References
+ CCPS, *Guidelines for Consequence Analysis of Chemical Releases*, American Institute of Chemical Engineers, New York (1999)

"""
function _bm_pl_c(α)
    concs = [0.10, 0.05, 0.02, 0.01, 0.005, 0.002]
    βs = [
            _bm_pl_c_1(α),   # (C_m/C_0) = 0.1
            _bm_pl_c_05(α),  # (C_m/C_0) = 0.05
            _bm_pl_c_02(α),  # (C_m/C_0) = 0.02
            _bm_pl_c_01(α),  # (C_m/C_0) = 0.01
            _bm_pl_c_005(α), # (C_m/C_0) = 0.005
            _bm_pl_c_002(α)  # (C_m/C_0) = 0.002
    ]
    return concs, βs
end


###############################################################################
# puff correlations

function _bm_pf_c_1(α)
    # corresponds to (Cm/C0) = 0.10
    if α ≤ -0.44
        return 0.70
    elseif -0.44 < α ≤ 0.43
        return 0.26*α + 0.81
    elseif 0.43 < α ≤ 1.0
        return 0.93
    else
        @warn "α= $α is out of range"
        return 0.93
    end
end

function _bm_pf_c_05(α)
    # corresponds to (Cm/C0) = 0.05
    if α ≤ -0.56
        return 0.85
    elseif -0.56 < α ≤ 0.31
        return 0.26*α + 1.00
    elseif 0.31 < α ≤ 1.0
        return -0.12*α + 1.12
    else
        @warn "α= $α is out of range"
        α = min(α,10.0)
        return -0.12*α + 1.12
    end
end

function _bm_pf_c_02(α)
    # corresponds to (Cm/C0) = 0.02
    if α ≤ -0.66
        return 0.95
    elseif -0.66 < α ≤ 0.32
        return 0.36*α + 1.19
    elseif 0.32 < α ≤ 1.0
        return -0.26*α + 1.38
    else
        @warn "α= $α is out of range"
        α = min(α,10.0)
        return -0.26*α + 1.38
    end
end

function _bm_pf_c_01(α)
    # corresponds to (Cm/C0) = 0.01
    if α ≤ -0.71
        return 1.15
    elseif -0.71 < α ≤ 0.37
        return 0.34*α + 1.39
    elseif 0.37 < α ≤ 1.0
        return -0.38*α + 1.66
    else
        @warn "α= $α is out of range"
        α = min(α,10.0)
        return -0.38*α + 1.66
    end
end

function _bm_pf_c_005(α)
    # corresponds to (Cm/C0) = 0.005
    if α ≤ -0.52
        return 1.48
    elseif -0.52 < α ≤ 0.24
        return 0.26*α + 1.62
    elseif 0.24 < α ≤ 1.0
        return -0.30*α + 1.75
    else
        @warn "α= $α is out of range"
        α = min(α,10.0)
        return -0.30*α + 1.75
    end
end

function _bm_pf_c_002(α)
    # corresponds to (Cm/C0) = 0.002
    if α ≤ 0.27
        return 1.83
    elseif 0.27 < α ≤ 1.0
        return -0.32*α + 1.92
    else
        @warn "α= $α is out of range"
        α = min(α,10.0)
        return -0.32*α + 1.92
    end
end

function _bm_pf_c_001(α)
    # corresponds to (Cm/C0) = 0.001
    if α ≤ -0.10
        return 2.075
    elseif -0.10 < α ≤ 1.0
        return -0.27*α + 2.05
    else
        @warn "α= $α is out of range"
        α = min(α,10.0)
        return -0.27*α + 2.05
    end
end

"""
    _bm_pf_c(α)

Britter-McQuaid puff correlations

# References
+ CCPS, *Guidelines for Consequence Analysis of Chemical Releases*, American Institute of Chemical Engineers, New York (1999)

"""
function _bm_pf_c(α)
    concs = [0.10, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001]
    βs = [
            _bm_pf_c_1(α),   # (C_m/C_0) = 0.1
            _bm_pf_c_05(α),  # (C_m/C_0) = 0.05
            _bm_pf_c_02(α),  # (C_m/C_0) = 0.02
            _bm_pf_c_01(α),  # (C_m/C_0) = 0.01
            _bm_pf_c_005(α), # (C_m/C_0) = 0.005
            _bm_pf_c_002(α), # (C_m/C_0) = 0.002
            _bm_pf_c_001(α)  # (C_m/C_0) = 0.001
    ]
    return concs, βs
end
