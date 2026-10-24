### A Pluto.jl notebook ###
# v1.0.4

using Markdown
using InteractiveUtils

# ╔═╡ 3c2d24f1-3688-4da2-b59f-94735eb66842
# ╠═╡ skip_as_script = true
#=╠═╡
include("../functions.jl")
  ╠═╡ =#

# ╔═╡ 302d2f4d-5c3e-4393-9a08-775abb86c03a
# ╠═╡ skip_as_script = true
#=╠═╡
include("../structs.jl")
  ╠═╡ =#

# ╔═╡ 7f38c6d2-c4ce-11f1-9d8d-ebbb8a18134b
md"""
This is a re-implementation of the original SLAB `thermo` subroutine. The intention is for it to be more modular and align more closely with the SLAB User's Guide model description in section 2.5
"""

# ╔═╡ 0f3358c6-99ee-4fd7-8962-5aa211c09960
md"""
# Concentrations
"""

# ╔═╡ b403b73b-7515-4222-bb09-74b228546eb8
md"""
## Cross-wind averaged mass fractions

The cross-wind averaged mass concentration, ``m``, is the total mass fraction of the released substance within the cloud, averaged over the cloud cross sectional area.

This is given in the SLAB User's Guide as (equation 1a, page 28)
```math
m = \begin{cases} 
          q_s \frac{B_{s} + x}{4B_{s} R} & x \leq B_{s} \\
          \frac{q_{s}}{2R} & x \gt B_{s}
       \end{cases}
```

Where ``q_s`` is the source strength and ``B_{s}`` is the effective source radius.
"""

# ╔═╡ 5bee5899-54ef-4a31-a535-7eafe3408cd7
@inline function total_mass_fraction(downwind_distance,
									 curr::SLAB_Plume_Integrated, 
									 p::SLAB_Release_Parameters)
	Bₛ = p.source.bse
	qₛ  = p.source.qs
	R   = curr.r
	x   = downwind_distance
	if x > Bₛ
		return qₛ/(2R)
	else
		return qₛ*(Bₛ + x)/(4Bₛ*R)
	end
end

# ╔═╡ 793d5460-5d70-422f-857b-1d0483bc9cd4
md"""
The cross-wind averaged vapour mass fraction, ``m_{ev}`` is given in equation 40d (page 29) in the SLAB User's Guide

```math
m_{ev} = m + \left( \frac{R_0}{R} \right) \left( m_{ev0} - m_{0} \right) + \delta m_{epc}
```

where subscript 0 indicates the value for the previous accepted integration step and ``\delta m_{epc}`` is the change in vapour mass fraction due to phase change during the current step.
"""

# ╔═╡ 23306e15-17d6-4538-b218-336528c95a92
@inline function vapour_mass_fraction(total_mass_fraction,
									  curr::SLAB_Plume_Integrated,
									  last::SLAB_Plume_State; δm=0)
	m = total_mass_fraction
    Rₒ = last.integrated.r
    R = curr.r
    mₑᵥₒ = last.thermo.cmev
    mₒ = last.thermo.cm
    return m + (Rₒ/R)*(mₑᵥₒ - mₒ) + δm
end

# ╔═╡ 09586e1e-28c1-452b-8c99-bbb970d82ece
md"""
## Cross-wind averaged water mass fractions

The cross-wind averaged mass fraction of water in the plume, $m_w$ is given by equation 40a (page 28) in the SLAB User's Guide

```math
m_{w} = (1-m) m_{wa}
```

Where $m$ is the total mass fraction of the released substance and $m_{wa}$ is the total mass fraction of water in ambient air.

The plume is assumed to be made of three components: the released substance, with total mass fraction $m$, entrained dry air, and entrained water that originally came from ambient air. The dry air and water has total mass fraction $1-m$ and has the same ratio of dry air/water as ambient air because that is where it came from.

"""

# ╔═╡ 3fe16300-7c9d-42e2-846b-7aee12baf839
@inline function water_mass_fraction(total_mass_fraction,
									 p::SLAB_Release_Parameters)
	m_wa = p.ambient_air.cmwa # mass fraction water in ambient air
	m = total_mass_fraction
	return (1-m)*m_wa
end

# ╔═╡ 31676b7c-3822-4a27-abab-881da22df212
md"""
The cross-wind averaged mass fraction of water which is in the vapour phase, $m_{wv}$ is given by equation 40c (page 29) in the SLAB User's Guide

```math
m_{wv} = m_w + \frac{R_0}{R} \left(m_{wv0} - m_{w0} \right) + \delta m_{wpc}
```

where subscript 0 indicates the value for the previous accepted integration step and $\delta m_{wpc}$ is the change in water vapour mass fraction due to phase change during the current step.
"""

# ╔═╡ f30ea944-57ec-4d13-9c7b-18aecc6c1aca
@inline function water_vapour_mass_fraction(total_water_mass_fraction,
											curr::SLAB_Plume_Integrated,
											last::SLAB_Plume_State; δm=0)
	m_w = total_water_mass_fraction
	R₀ = last.integrated.r
    R = curr.r
    m_wv₀ = last.thermo.cmwv # mass fraction water vapour, previous step
    m_w₀ = last.thermo.cmw # total mass fraction previous step
    return m_w + (R₀/R)*(m_wv₀ - m_w₀) + δm
end

# ╔═╡ c5635f10-04b9-4f14-bde5-f3be2d05ed03
md"""
## Cross-wind averaged dry air mass fraction

The cross-wind averaged mass fraction of dry air in the plume, $m_{da}$ is given by equation 40b (page 29) in the SLAB User's Guide

```math
m_{da} = \left( 1 -m \right) \left( 1 - m_{wa} \right)
```

The plume is assumed to be made of three components: the released substance, with total mass fraction $m$, entrained dry air, $m_{da}$, and entrained water that originally came from ambient air. The dry air and water has total mass fraction $1-m$ and has the same ratio of dry air/water as ambient air because that is where it came from.
"""

# ╔═╡ 994492cf-6d83-46f7-b712-d6547d453414
@inline function dry_air_mass_fraction(total_mass_fraction,
									   p::SLAB_Release_Parameters)
	m_wa = p.ambient_air.cmwa # mass fraction water in ambient air
	m = total_mass_fraction
	return (1-m)*(1-m_wa)
end

# ╔═╡ 96bd7796-4465-45cb-b7b2-102e30c9a6a4
md"""
## Volume Fractions

The volume fraction is calculated from equation 12 (page 15) of the SLAB User's Guide

```math
C = { {M_{ae} m } \over {M_s + \left( M_{ae} - M_s \right) m }}
```

Where $M_{ae}$ is the molar mass of ambient air, $M_s$ the molar mass of the released substance, and $m$ the total mass fraction of released substance.

Note: the User's Guide uses $M_a$ for the molar mass of ambient air in equation 12, the SLAB fortran code uses $M_{ae}$ which is replicated here.
"""

# ╔═╡ 9c751efa-7ea4-4d1e-8ead-6ce17f029f1e
@inline function volume_fraction(total_mass_fraction,
								 p::SLAB_Release_Parameters)
	Mₐₑ = p.ambient_air.wmae # molar mass ambient air
	Mₛ  = p.source.wms  # molar mass released substance
	m   = total_mass_fraction
	return (Mₐₑ*m)/(Mₛ + (Mₐₑ - Mₛ)*m)
end

# ╔═╡ ee454b66-c076-4024-ba6d-dee2adccf9ec
md"""
# Enthalpies and Specific Heats
"""

# ╔═╡ a3361184-40ab-4374-9d1d-94a5635ae620
md"""
## Total Cloud Specific Heat
"""

# ╔═╡ 2794e095-3884-479b-b4e6-d2d7252aaee4
@inline function total_specific_heat(dry_air_mass_fraction,
									 water_vapour_mass_fraction,
									 water_liquid_mass_fraction,
									 vapour_mass_fraction,
									 liquid_mass_fraction,
									 p::SLAB_Release_Parameters)
	m_da = dry_air_mass_fraction
	cpₐ  = p.ambient_air.cpa # dry air specific heat
	
	m_wv  = water_vapour_mass_fraction
	m_wl  = water_liquid_mass_fraction
	cp_wv = p.ambient_air.cpwv # water vapour specific heat
	cp_wl = p.ambient_air.cpwl # water liquid specific heat
	
	mₑᵥ  = vapour_mass_fraction
	mₑₗ  = liquid_mass_fraction
	cpₛ  = p.source.cps # released substance vapour specific heat
	cpₛₗ = p.source.cpsl # released substance liquid specific heat

	return m_da*cpₐ + m_wv*cp_wv + m_wl*cp_wl + mₑᵥ*cpₛ + mₑₗ*cpₛₗ
end
	

# ╔═╡ ce7cdd1b-a5c9-41b2-9170-34e17b784e41
md"""
# Transported Energy

The transported energy is calculated from equation 3c (page 21) of the SLAB User's Guide, and comes from the slab thermo routine line 2276-9 where $T = \frac{E_{trn}}{C_p}$

```math
E_{trn} = { {R_0 C_{p0} T_0 + E - E_0 + F_{pc} + F_t} \over {R}}
```

where $C_p$ is the total specific heat of the plume, $F_{pc} = \int_{x_0}^x f_{pc} dx$ is the net phase change energy, and $F_t = \int_{x_0}^x f_{t} dx$ is the ned ground heat flux and 0 indicates the previous accepted integration step. The introduced variable $E$ is defined as

```math
E = R \left( 1 - m \right) C_{pa} T_a + R m C_{ps} T_s
```

Substituting and collecting like terms

```math
E_{trn} = \left( 1 - m \right) C_{pa} T_a + m C_{ps} T_s + \frac{R_0}{R} 
\left( C_{p0} T_0 - \left( 1 - m_0 \right) C_{pa} T_a - m_0 C_{ps} T_s + F_{pc} F_{t} \right)
```

This is used to generate the initial guess of the plume temperature and phase behaviour, which is refined using subsequent Newton iteration.

"""

# ╔═╡ 1f715ead-8bbe-4735-9b21-0e9e0a0666a8
@inline function transported_energy(total_mass_fraction, 
                                    curr::SLAB_Plume_Integrated, 
                                    last::SLAB_Plume_State, 
                                    p::SLAB_Release_Parameters)
    m = total_mass_fraction
    R = curr.r

    m₀  = last.thermo.cm
    R₀  = last.integrated.r
    Cₚ₀ = last.thermo.cp
    T₀  = last.thermo.t
    
    Cₚₐ = p.ambient_air.cpaa
    Tₐ  = p.ambient_air.ta

    Cₚₛ = p.source.cps
    Tₛ  = p.source.ts
    
	# this adjusts the step $F_t$ to ensure conservation of energy
    Fₜ = curr.sft
    Fₜₑ = iszero(Fₜ) ? zero(Fₜ) :
           -Cₚ₀*(T₀ - Tₐ)*abs(Fₜ)/
               (abs(Fₜ) + R₀*Cₚ₀*abs(T₀ - Tₐ))
    
    base = (1-m)*Cₚₐ*Tₐ + m*Cₚₛ*Tₛ
    excess = Cₚ₀*T₀ - (1-m₀)*Cₚₐ*Tₐ - m₀*Cₚₛ*Tₛ + Fₜₑ

    return base + (R₀/R)*excess
end

# ╔═╡ a8ade6fc-3ec6-4da4-a67b-8cd6a4bcf05a
md"""
# Initial Phase Check

SLAB checks the phase of the released substance and water in the plume, and uses this to decide whether a more detailed solve of the phase partition is required. It does this by calculating two properties:

1. The partial pressure for a pure substance at the initial temperature, i.e. the ratio of the saturation pressure to the actual pressure. This initial temperature is the temperature calculated from the transported energy, from the previous accepted integration step.
2. The partial pressure from the current integration step, assuming $\delta m_{pc} =0$

If the material was entirely in the vapour phase in the previous accepted integration step *and* the partial pressure calculated in 1. is greater than in 2. (i.e. the thermodynamics favours an all vapour phase), then the material is assumed to be entirely a vapour.

This check is done independently for the released substance and water. If both are entirely vapour, all the phase partition calculations are skipped. If one of them potentially has a liquid phase present, more detailed routines are used (described below) to determine the actual phase partition.
"""

# ╔═╡ bc0b92ec-4112-4cff-badf-2ffd1909244a
@inline function vapour_partial_pressure(mass_fraction,
										molar_mass,
										vapour_mass_fraction,
										water_vapour_mass_fraction,
										dry_air_mass_fraction,
										p::SLAB_Release_Parameters)

    Mₐ  = p.ambient_air.wma  # Molar mass dry air
    M_w = p.ambient_air.wmw  # Molar mass water
    Mₛ  = p.source.wms    # Molar mass released substance

    m  = mass_fraction
    M  = molar_mass
    mₑᵥ = vapour_mass_fraction
    m_wv= water_vapour_mass_fraction
    m_da = dry_air_mass_fraction

    return m/M/(m_da/Mₐ + m_wv/M_w + mₑᵥ/Mₛ)
end

# ╔═╡ 4f00285b-9609-4a26-8af5-11abf30998ed
@inline function vapour_check_substance(total_mass_fraction,
                                        vapour_mass_fraction,
                                        water_vapour_mass_fraction,
                                        dry_air_mass_fraction,
                                        temperature,
                                        last::SLAB_Plume_State, 
                                        p::SLAB_Release_Parameters)

    if !isapprox(last.thermo.cm, last.thermo.cmev)
        return false
    end
    
    m  = total_mass_fraction
    Mₛ  = p.source.wms    # Molar mass released substance

    y₁ = _slab_psat(temperature, p.source.antoine)
    y₂ = vapour_partial_pressure(m, Mₛ, vapour_mass_fraction,
										water_vapour_mass_fraction,
										dry_air_mass_fraction, p)

    return y₁ > y₂
end

# ╔═╡ d23524d9-5e6b-4157-9898-3760853b2b8f
@inline function vapour_check_water(water_mass_fraction,
                                        water_vapour_mass_fraction,
                                        vapour_mass_fraction,
                                        dry_air_mass_fraction,
                                        temperature,
                                        last::SLAB_Plume_State, 
                                        p::SLAB_Release_Parameters)

    if !isapprox(last.thermo.cmw, last.thermo.cmwv)
        return false
    end
    
    m_w  = water_mass_fraction
    M_w = p.ambient_air.wmw  # Molar mass water

    y₁ = _slab_psat(temperature, p.ambient_air.antoine)
    y₂ = vapour_partial_pressure(m_w, M_w, vapour_mass_fraction,
										water_vapour_mass_fraction,
										dry_air_mass_fraction, p)

    return y₁ > y₂
end

# ╔═╡ 02f18a27-7032-4d09-8820-3641f31c70ce
md"""
# Final Subroutines
"""

# ╔═╡ 1f0c08f7-2651-4832-8e8a-6af6fcb70554
function plume_thermodynamics(x, 
                              curr::SLAB_Plume_Integrated, 
                              last::SLAB_Plume_State, 
                              p::SLAB_Release_Parameters)
	# calculate the mass fractions and volume fraction of the plume
    cm = total_mass_fraction(x, curr, p)
    cv = volume_fraction(cm,p)
    cmw = water_mass_fraction(cm, p)
    cmda = dry_air_mass_fraction(cm, p)

    # initial guess for vapour phase mass concentrations
    cmwvt = water_vapour_mass_fraction(cmw, curr, last)
    cmevt = vapour_mass_fraction(cm, curr, last)

	# initial guess for heat capacity
	cp = total_specific_heat(cmda,cmwvt,cmw-cmwvt,cmevt,cm-cmevt,p)
    
    # temperature from transported energy
    etrn = transported_energy(cm, curr, last, p)
    t = etrn/last.thermo.cp

	cmwv = cmwvt
	cmwd = cmw - cmwv
	cmev = cmevt
	cmed = cm - cmev
	rho = zero(eltype(cm))
	phase = nothing

    return SLAB_Thermo_State(cm,cv,cmw,cmda,cmwv,cmwd,cmev,cmed,cp,t,rho,phase,etrn)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.13.2"
manifest_format = "2.1"
project_hash = "71853c6197a6a7f222db0f1978c7cb232b87c5ee"

[deps]
"""

# ╔═╡ Cell order:
# ╟─7f38c6d2-c4ce-11f1-9d8d-ebbb8a18134b
# ╠═3c2d24f1-3688-4da2-b59f-94735eb66842
# ╠═302d2f4d-5c3e-4393-9a08-775abb86c03a
# ╟─0f3358c6-99ee-4fd7-8962-5aa211c09960
# ╟─b403b73b-7515-4222-bb09-74b228546eb8
# ╠═5bee5899-54ef-4a31-a535-7eafe3408cd7
# ╟─793d5460-5d70-422f-857b-1d0483bc9cd4
# ╠═23306e15-17d6-4538-b218-336528c95a92
# ╟─09586e1e-28c1-452b-8c99-bbb970d82ece
# ╠═3fe16300-7c9d-42e2-846b-7aee12baf839
# ╟─31676b7c-3822-4a27-abab-881da22df212
# ╠═f30ea944-57ec-4d13-9c7b-18aecc6c1aca
# ╟─c5635f10-04b9-4f14-bde5-f3be2d05ed03
# ╠═994492cf-6d83-46f7-b712-d6547d453414
# ╟─96bd7796-4465-45cb-b7b2-102e30c9a6a4
# ╠═9c751efa-7ea4-4d1e-8ead-6ce17f029f1e
# ╟─ee454b66-c076-4024-ba6d-dee2adccf9ec
# ╟─a3361184-40ab-4374-9d1d-94a5635ae620
# ╠═2794e095-3884-479b-b4e6-d2d7252aaee4
# ╟─ce7cdd1b-a5c9-41b2-9170-34e17b784e41
# ╠═1f715ead-8bbe-4735-9b21-0e9e0a0666a8
# ╟─a8ade6fc-3ec6-4da4-a67b-8cd6a4bcf05a
# ╠═bc0b92ec-4112-4cff-badf-2ffd1909244a
# ╠═4f00285b-9609-4a26-8af5-11abf30998ed
# ╠═d23524d9-5e6b-4157-9898-3760853b2b8f
# ╟─02f18a27-7032-4d09-8820-3641f31c70ce
# ╠═1f0c08f7-2651-4832-8e8a-6af6fcb70554
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
