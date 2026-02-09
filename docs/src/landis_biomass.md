# LANDIS Biomass Module

## Overview

The LANDIS biomass module models forest growth, mortality, and dead woody biomass
decomposition for individual species-age cohorts. It integrates aboveground net
primary productivity (ANPP), biomass-related mortality, age-related mortality, and
dead wood decomposition into a cohort-level biomass dynamics model suitable for
landscape-scale forest simulations.

The model tracks two state variables per cohort: living aboveground biomass and
dead woody biomass. Growth follows a peaked function that increases at low biomass
and decreases as the cohort approaches its maximum potential biomass. Mortality
combines a logistic biomass-dependent term (increasing as the site fills up) and
an exponential age-dependent term (increasing sharply as the cohort approaches its
maximum lifespan). Dead biomass decays exponentially.

**Reference**: Scheller, R.M. and Mladenoff, D.J. (2004). A forest growth and
biomass module for a landscape simulation model, LANDIS: design, validation,
and application. *Ecological Modelling*, 180, 211-229.
doi:10.1016/j.ecolmodel.2004.01.022

```@docs
LANDISBiomass
```

## Implementation

### State Variables

```@example landis
using DataFrames, ModelingToolkit, Symbolics, DynamicQuantities
using Vegetation

sys = LANDISBiomass()
compiled = mtkcompile(sys)
vars = unknowns(compiled)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars],
)
```

### Parameters

```@example landis
params = parameters(compiled)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [dimension(ModelingToolkit.get_unit(p)) for p in params],
    :Default => [ModelingToolkit.hasdefault(p) ? ModelingToolkit.getdefault(p) : missing for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params],
)
```

### Equations

```@example landis
eqs = equations(sys)
```

### Table 1: Acronyms, Definitions, and Units

The following table reproduces Table 1 from Scheller and Mladenoff (2004),
showing acronyms, definitions, and units for the biomass model variables.

```@example landis
DataFrame(
    :Acronym => [
        "ANPP_MAX", "ANPP_ACT", "B", "B_AP",
        "B_MAX", "B_POT", "B_PM", "D_wood",
        "M_BIO", "M_AGE",
    ],
    :Definition => [
        "Maximum ANPP for species",
        "Actual ANPP for species, age cohort",
        "Actual biomass for species, age cohort",
        "Ratio of actual to potential biomass for species, age cohort",
        "Maximum possible biomass for species",
        "Potential biomass — the limit to biomass based on growing space available",
        "Ratio of potential to maximum biomass for species, age cohort",
        "Dead biomass, including coarse woody debris (CWD) and snags",
        "Mortality rate as a function of total stand biomass",
        "Mortality rate as a function of species cohort age and maximum species longevity",
    ],
    Symbol("Units (paper)") => [
        "Mg ha⁻¹ year⁻¹", "Mg ha⁻¹ year⁻¹", "Mg ha⁻¹", "dimensionless",
        "Mg ha⁻¹", "Mg ha⁻¹", "dimensionless", "Mg ha⁻¹",
        "Mg ha⁻¹ year⁻¹", "Mg ha⁻¹ year⁻¹",
    ],
)
```

### Table 2: Species Life History Attributes

The following table reproduces Table 2 from Scheller and Mladenoff (2004),
showing life history attributes and maximum ANPP for 16 species characteristic
of northern Wisconsin. These values are used to parameterize the biomass module
for individual species.

```@example landis
species_data = DataFrame(
    :Species => [
        "Abies balsamea", "Acer rubrum", "Acer saccharum",
        "Betula alleghaniensis", "Betula papyrifera", "Fraxinus americana",
        "Picea glauca", "Pinus banksiana", "Pinus resinosa",
        "Pinus strobus", "Populus tremuloides", "Quercus ellipsoidalis",
        "Quercus rubra", "Thuja occidentalis", "Tilia americana",
        "Tsuga canadensis",
    ],
    Symbol("Longevity (yr)") => [
        200, 150, 400, 350, 120, 300,
        300, 70, 250, 450, 120, 300,
        250, 400, 250, 640,
    ],
    Symbol("Shade tol.") => [5, 4, 5, 4, 2, 4, 3, 1, 2, 3, 1, 2, 3, 4, 4, 5],
    Symbol("Fire tol.") => [1, 1, 1, 2, 2, 1, 2, 3, 4, 3, 2, 5, 3, 1, 2, 3],
    Symbol("Max ANPP (Mg/ha/yr)") => [
        8.09, 7.46, 7.45, 7.64, 6.72, 7.00,
        6.50, 5.77, 5.01, 10.19, 7.46, 7.73,
        7.58, 4.93, 8.48, 6.27,
    ],
)
```

## Analysis

### Fig. 3a: ANPP and Mortality vs. Actual/Potential Biomass Ratio

This figure reproduces Fig. 3a from the paper, showing the correlation between
ANPP (as a fraction of ANPP\_MAX) and biomass-related mortality (M\_BIO as a
fraction of ANPP\_MAX) as a function of B\_AP (the ratio of actual to potential
biomass). The ANPP curve peaks at B\_AP = 1 while mortality increases
logistically.

```@example landis
using Plots
using OrdinaryDiffEqDefault

compiled = mtkcompile(LANDISBiomass())

yr_to_s = 3.15576e7
Mg_ha_to_kg_m2 = 0.1

# Sweep B values to vary B_AP from ~0 to ~1
# B_MAX = ANPP_MAX * 30 yr = 7.45 * 0.1 / yr_to_s * 30 * yr_to_s = 22.35 kg/m²
B_MAX_val = 7.45 * Mg_ha_to_kg_m2 / yr_to_s * 30.0 * yr_to_s
ANPP_MAX_val = 7.45 * Mg_ha_to_kg_m2 / yr_to_s

B_AP_vals = Float64[]
anpp_frac = Float64[]
mort_frac = Float64[]
for frac in 0.01:0.01:1.0
    B_val = frac * B_MAX_val
    prob = ODEProblem(compiled,
        Dict(compiled.B => B_val, compiled.D_wood => 0.0),
        (0.0, 1.0))
    sol = solve(prob)
    push!(B_AP_vals, sol[compiled.B_AP][1])
    push!(anpp_frac, sol[compiled.ANPP_ACT][1] / ANPP_MAX_val)
    push!(mort_frac, sol[compiled.M_BIO][1] / ANPP_MAX_val)
end

p = plot(B_AP_vals, anpp_frac, label = "ANPP", linewidth = 2,
    xlabel = "B_AP (Actual biomass / Potential biomass)",
    ylabel = "Fraction ANPP_MAX",
    title = "Fig. 3a: Growth and Mortality vs. B_AP",
    legend = :right, ylim = (0, 1.05))
plot!(p, B_AP_vals, mort_frac, label = "Mortality", linewidth = 2, linestyle = :dash)
p
```

### Fig. 3b: Age-Related Mortality

This figure reproduces Fig. 3b from the paper, showing the fraction of biomass
removed by age-related mortality as a function of the fraction of species
lifespan. With the default shape parameter d = 10, mortality begins to increase
noticeably around 50% of lifespan and reaches 100% at the maximum lifespan.

```@example landis
# Sweep age_init to vary age fraction from 0 to 1
B_ref = 10.0 * Mg_ha_to_kg_m2  # Reference biomass for computing M_AGE fraction
max_age_val = 400.0 * yr_to_s

age_frac_vals = Float64[]
age_mort_frac = Float64[]
for frac in 0.0:0.01:1.0
    age_val = frac * max_age_val
    prob = ODEProblem(compiled,
        Dict(compiled.B => B_ref, compiled.D_wood => 0.0,
            compiled.age_init => age_val),
        (0.0, 1.0))
    sol = solve(prob)
    # M_AGE fraction = M_AGE / (B / one_yr) where one_yr = yr_to_s
    m_age = sol[compiled.M_AGE][1]
    push!(age_frac_vals, frac)
    push!(age_mort_frac, m_age / (B_ref / yr_to_s))
end

p = plot(age_frac_vals, age_mort_frac, linewidth = 2, label = nothing,
    xlabel = "Fraction of species' lifespan",
    ylabel = "Fraction of biomass removed",
    title = "Fig. 3b: Age-Related Mortality",
    ylim = (0, 1.05))
p
```

### Fig. 6a: A. saccharum Single-Species Growth

This figure reproduces Fig. 6a from the paper, showing the growth trajectory of
a single *Acer saccharum* (sugar maple) cohort with default parameters
(ANPP\_MAX = 7.45 Mg/ha/yr, max\_age = 400 years). Living biomass reaches an
asymptote near 240 Mg/ha, while dead woody biomass accumulates over time. The
paper reports that the mortality rate does not exceed the decomposition rate
until around year 50, after which dead biomass remains roughly constant near
25 Mg/ha.

```@example landis
tspan_s = 200.0 * yr_to_s
prob = ODEProblem(compiled, [], (0.0, tspan_s))
sol = solve(prob)

years = 1:200
B_vals = [sol(yr * yr_to_s; idxs = compiled.B) / Mg_ha_to_kg_m2 for yr in years]
D_vals = [sol(yr * yr_to_s; idxs = compiled.D_wood) / Mg_ha_to_kg_m2 for yr in years]

p = plot(years, B_vals, label = "A. saccharum", linewidth = 2,
    xlabel = "Simulation year",
    ylabel = "Biomass (Mg ha⁻¹)",
    title = "Fig. 6a: Single Species Growth (A. saccharum)",
    legend = :right, ylim = (0, 300))
plot!(p, years, D_vals, label = "Dead biomass (D)", linewidth = 2, linestyle = :dash)
p
```

### Fig. 6b: P. tremuloides Single-Species Growth

This figure reproduces Fig. 6b from the paper, showing the growth trajectory of
a single *Populus tremuloides* (aspen) cohort. Because aspen is a shade-intolerant
species with a short lifespan (120 years), its biomass exhibits cyclical
behavior: it rises rapidly, then declines due to age-related mortality. The
paper reports that aboveground living biomass of *P. tremuloides* ranges from
approximately 94 to 124 Mg/ha, with empirical measurements indicating a range
of 95--125 Mg/ha for mature aspen stands.

```@example landis
# P. tremuloides parameters from Table 2
prob_poptre = ODEProblem(
    compiled,
    Dict(
        compiled.B => 5.0 * Mg_ha_to_kg_m2, compiled.D_wood => 0.0,
        compiled.max_age => 120.0 * yr_to_s,
        compiled.ANPP_MAX => 7.46 * Mg_ha_to_kg_m2 / yr_to_s
    ),
    (0.0, tspan_s)
)
sol_poptre = solve(prob_poptre)

B_poptre = [sol_poptre(yr * yr_to_s; idxs = compiled.B) / Mg_ha_to_kg_m2 for yr in years]
D_poptre = [sol_poptre(yr * yr_to_s; idxs = compiled.D_wood) / Mg_ha_to_kg_m2 for yr in years]

p = plot(years, B_poptre, label = "P. tremuloides", linewidth = 2,
    xlabel = "Simulation year",
    ylabel = "Biomass (Mg ha⁻¹)",
    title = "Fig. 6b: Single Species Growth (P. tremuloides)",
    legend = :right, ylim = (0, 200))
plot!(p, years, D_poptre, label = "Dead biomass (D)", linewidth = 2, linestyle = :dash)
p
```

### Fig. 6c: Multi-Species Competition

This figure reproduces Fig. 6c from the paper, showing a site initialized with
two mature mid-tolerant species (*Pinus strobus* and *Betula alleghaniensis*)
and two young shade-tolerant species (*Acer saccharum* and *Tsuga canadensis*).
In the paper, the shade-tolerant species eventually overtake as the mature
cohorts senesce. Here, each species is simulated independently with its
respective parameters from Table 2 to illustrate their individual growth
trajectories. Note: true multi-cohort competition with dynamic `B_other` would
require coupling multiple `LANDISBiomass` systems.

```@example landis
# Species parameters from Table 2
# P. strobus: Longevity=450, ANPP_MAX=10.19, mature (age_init=200 yr)
prob_pinstro = ODEProblem(
    compiled,
    Dict(
        compiled.B => 50.0 * Mg_ha_to_kg_m2, compiled.D_wood => 0.0,
        compiled.max_age => 450.0 * yr_to_s,
        compiled.ANPP_MAX => 10.19 * Mg_ha_to_kg_m2 / yr_to_s,
        compiled.age_init => 200.0 * yr_to_s
    ),
    (0.0, tspan_s)
)
sol_pinstro = solve(prob_pinstro)

# B. alleghaniensis: Longevity=350, ANPP_MAX=7.64, mature (age_init=200 yr)
prob_betall = ODEProblem(
    compiled,
    Dict(
        compiled.B => 50.0 * Mg_ha_to_kg_m2, compiled.D_wood => 0.0,
        compiled.max_age => 350.0 * yr_to_s,
        compiled.ANPP_MAX => 7.64 * Mg_ha_to_kg_m2 / yr_to_s,
        compiled.age_init => 200.0 * yr_to_s
    ),
    (0.0, tspan_s)
)
sol_betall = solve(prob_betall)

# A. saccharum: Longevity=400, ANPP_MAX=7.45, young (age_init=10 yr, default)
prob_acesac = ODEProblem(compiled, [], (0.0, tspan_s))
sol_acesac = solve(prob_acesac)

# T. canadensis: Longevity=640, ANPP_MAX=6.27, young (age_init=10 yr)
prob_tsuga = ODEProblem(
    compiled,
    Dict(
        compiled.B => 5.0 * Mg_ha_to_kg_m2, compiled.D_wood => 0.0,
        compiled.max_age => 640.0 * yr_to_s,
        compiled.ANPP_MAX => 6.27 * Mg_ha_to_kg_m2 / yr_to_s
    ),
    (0.0, tspan_s)
)
sol_tsuga = solve(prob_tsuga)

B_pinstro = [sol_pinstro(yr * yr_to_s; idxs = compiled.B) / Mg_ha_to_kg_m2 for yr in years]
B_betall = [sol_betall(yr * yr_to_s; idxs = compiled.B) / Mg_ha_to_kg_m2 for yr in years]
B_acesac = [sol_acesac(yr * yr_to_s; idxs = compiled.B) / Mg_ha_to_kg_m2 for yr in years]
B_tsuga = [sol_tsuga(yr * yr_to_s; idxs = compiled.B) / Mg_ha_to_kg_m2 for yr in years]
D_total = [sol_acesac(yr * yr_to_s; idxs = compiled.D_wood) / Mg_ha_to_kg_m2 for yr in years]

p = plot(years, B_pinstro, label = "P. strobus", linewidth = 2,
    xlabel = "Simulation year",
    ylabel = "Biomass (Mg ha⁻¹)",
    title = "Fig. 6c: Multi-Species Growth",
    legend = :right, ylim = (0, 350))
plot!(p, years, B_betall, label = "B. alleghaniensis", linewidth = 2, linestyle = :dash)
plot!(p, years, B_acesac, label = "A. saccharum", linewidth = 2, linestyle = :dot)
plot!(p, years, B_tsuga, label = "T. canadensis", linewidth = 2, linestyle = :dashdot)
plot!(p, years, D_total, label = "Dead biomass (D)", linewidth = 2, color = :black, linestyle = :dash)
p
```

### Species Comparison: Short-Lived vs. Long-Lived Species

This analysis compares the growth trajectories of two contrasting species: a
short-lived, shade-intolerant species (*P. banksiana*, max\_age = 70 years,
ANPP\_MAX = 5.77 Mg/ha/yr) and a long-lived, shade-tolerant species
(*A. saccharum*, max\_age = 400 years, ANPP\_MAX = 7.45 Mg/ha/yr). The
short-lived species shows a rapid rise and decline in biomass due to age-related
mortality, while the long-lived species achieves higher peak biomass over a
longer time horizon.

```@example landis
# Short-lived species: P. banksiana
prob_short = ODEProblem(
    compiled,
    Dict(
        compiled.B => 5.0 * Mg_ha_to_kg_m2, compiled.D_wood => 0.0,
        compiled.max_age => 70.0 * yr_to_s,
        compiled.ANPP_MAX => 5.77 * Mg_ha_to_kg_m2 / yr_to_s
    ),
    (0.0, tspan_s)
)
sol_short = solve(prob_short)

# Long-lived species: A. saccharum (defaults)
prob_long = ODEProblem(compiled, [], (0.0, tspan_s))
sol_long = solve(prob_long)

B_short = [sol_short(yr * yr_to_s; idxs = compiled.B) / Mg_ha_to_kg_m2 for yr in years]
B_long = [sol_long(yr * yr_to_s; idxs = compiled.B) / Mg_ha_to_kg_m2 for yr in years]

p = plot(years, B_long, label = "A. saccharum (400 yr)", linewidth = 2,
    xlabel = "Simulation year",
    ylabel = "Living Biomass (Mg ha⁻¹)",
    title = "Species Longevity Comparison",
    legend = :topright, ylim = (0, 250))
plot!(p, years, B_short, label = "P. banksiana (70 yr)", linewidth = 2, linestyle = :dash)
p
```

### Competition Effect

This analysis demonstrates the effect of inter-cohort competition on growth. When
other cohorts occupy a significant fraction of the site's maximum biomass capacity
(B\_MAX\_site = 500 Mg/ha), the focal cohort's potential biomass (B\_POT) is
reduced, leading to lower ANPP and lower peak biomass. Competition only has an
effect when B\_other > B\_MAX\_site - B\_MAX (i.e., when the site is sufficiently
crowded that the growing space available to the cohort is less than its species
maximum).

```@example landis
tspan_comp = 150.0 * yr_to_s

# No competition
prob_no_comp = ODEProblem(
    compiled,
    Dict(
        compiled.B => 5.0 * Mg_ha_to_kg_m2, compiled.D_wood => 0.0,
        compiled.B_other => 0.0
    ),
    (0.0, tspan_comp)
)
sol_no_comp = solve(prob_no_comp)

# High competition: 400 Mg/ha of other species
prob_comp = ODEProblem(
    compiled,
    Dict(
        compiled.B => 5.0 * Mg_ha_to_kg_m2, compiled.D_wood => 0.0,
        compiled.B_other => 400.0 * Mg_ha_to_kg_m2
    ),
    (0.0, tspan_comp)
)
sol_comp = solve(prob_comp)

years_comp = 1:150
B_no_comp = [sol_no_comp(yr * yr_to_s; idxs = compiled.B) / Mg_ha_to_kg_m2 for yr in years_comp]
B_comp = [sol_comp(yr * yr_to_s; idxs = compiled.B) / Mg_ha_to_kg_m2 for yr in years_comp]

p = plot(years_comp, B_no_comp, label = "No competition", linewidth = 2,
    xlabel = "Simulation year",
    ylabel = "Living Biomass (Mg ha⁻¹)",
    title = "Effect of Competition on Growth",
    legend = :topright)
plot!(p, years_comp, B_comp, label = "B_other = 400 Mg/ha", linewidth = 2, linestyle = :dash)
p
```

## Limitations

The following aspects of the paper are not implemented in this module:

- **Multi-cohort and multi-species interactions**: This implementation models a single cohort. The paper describes landscape-scale simulations with multiple species-age cohorts competing within sites.
- **Disturbance modules**: Fire, wind, and harvesting disturbance processes (Sections 4.2.1, 4.2.2) are not included.
- **Shade calculation**: The percent full sunlight calculation (Eq. 8) and shade tolerance classes are not implemented.
- **Leaf litter partitioning**: The fine/coarse dead biomass partitioning described in Section 3.2.4 is not modeled.
- **PnET-II coupling**: The ANPP\_MAX parameterization via PnET-II (Section 4.3) is external to this module.
