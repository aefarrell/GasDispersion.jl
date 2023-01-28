
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

Britter-McQuaid plume correlations as digtized in TSCREEN

# References
+ EPA, *User's Guide to TSCREEN* U.S. Environmental Protection Agency EPA-454/B-94-023 (1994)

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
