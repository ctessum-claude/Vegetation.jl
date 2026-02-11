# Stage (1973) Prognosis Model

## Overview

The Prognosis Model for Stand Development is an individual-tree growth model
originally developed for lodgepole pine (*Pinus contorta*) stands. It predicts
periodic changes in tree diameter, height, crown ratio, and survival as functions
of tree size, site quality, stand density, and competitive position. The model is
notable for its stochastic treatment of diameter growth, where random normal errors
are applied to the logarithm of basal area increment to capture the
deviation-amplifying nature of the growth process.

This model is the foundation of the USDA Forest Service's Forest Vegetation
Simulator (FVS), one of the most widely used forest growth models in North America.

**Reference**: Stage, A.R. (1973). Prognosis model for stand development.
USDA Forest Service Research Paper INT-137.
Intermountain Forest and Range Experiment Station, Ogden, Utah. 32 p.

```@docs
StagePrognosis
```

## Implementation

The model tracks four state variables for each individual tree record:
- `Dsq`: Squared diameter at breast height (DBH²) — the primary growth variable
- `HT`: Total tree height
- `HCB`: Height to crown base
- `N_trees`: Tree density (trees per m²)

The key innovation is that diameter growth is modeled in diameter-squared space
(proportional to basal area), with stochastic noise applied multiplicatively.
This is implemented as a stochastic differential equation (SDE) using
ModelingToolkit's `@brownians` macro.

### State Variables

```@example stage_prognosis
using DataFrames, ModelingToolkit, Symbolics, DynamicQuantities
using Vegetation

sys = StagePrognosis()
compiled_sys = mtkcompile(sys)

vars = unknowns(compiled_sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape=false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

### Parameters

```@example stage_prognosis
params = parameters(compiled_sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape=false)) for p in params],
    :Units => [dimension(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

### Equations

```@example stage_prognosis
eqs = equations(compiled_sys)
```

## Analysis

### Deterministic Growth Trajectory

The following figure shows the deterministic growth trajectory (no stochastic noise)
for a lodgepole pine tree starting from conditions similar to Figure 1 in Stage (1973):
DBH = 6.0 inches, height = 60 feet, crown ratio = 0.33, 500 trees/acre,
site index = 80 feet, elevation = 4300 feet.

```@example stage_prognosis
using StochasticDiffEq
using Plots

one_year = 3.15576e7
one_inch = 0.0254
one_foot = 0.3048

prob = SDEProblem(compiled_sys, [], (0.0, 50.0 * one_year))
prob_det = remake(prob; p=[compiled_sys.σ_growth => 0.0])
sol = solve(prob_det, EM(), dt=one_year/10)

years = sol.t ./ one_year
dbh_in = sol[compiled_sys.DBH] ./ one_inch
ht_ft = sol[compiled_sys.HT] ./ one_foot
cr = sol[compiled_sys.CR]
n_acre = sol[compiled_sys.N_trees] .* 4046.86

p1 = plot(years, dbh_in, xlabel="Years", ylabel="DBH (inches)",
    title="Diameter Growth", legend=false, linewidth=2)
p2 = plot(years, ht_ft, xlabel="Years", ylabel="Height (feet)",
    title="Height Growth", legend=false, linewidth=2)
p3 = plot(years, cr, xlabel="Years", ylabel="Crown Ratio",
    title="Crown Development", legend=false, linewidth=2)
p4 = plot(years, n_acre, xlabel="Years", ylabel="Trees/Acre",
    title="Stand Density", legend=false, linewidth=2)

p = plot(p1, p2, p3, p4, layout=(2, 2), size=(800, 600),
    plot_title="Stage (1973) Prognosis Model - Deterministic")
p
```

### Stochastic Ensemble

The model includes multiplicative stochastic noise on the basal area increment,
as described in Stage (1973, p. 12). The following figure shows 20 stochastic
realizations alongside the deterministic trajectory, demonstrating the
deviation-amplifying behavior of the growth process.

```@example stage_prognosis
n_runs = 20
p1 = plot(xlabel="Years", ylabel="DBH (inches)", title="Stochastic Diameter Growth")
p2 = plot(xlabel="Years", ylabel="Height (feet)", title="Stochastic Height Growth")

for i in 1:n_runs
    sol_s = solve(prob, EM(), dt=one_year/10, seed=i)
    plot!(p1, sol_s.t ./ one_year, sol_s[compiled_sys.DBH] ./ one_inch,
        alpha=0.3, color=:blue, label="")
    plot!(p2, sol_s.t ./ one_year, sol_s[compiled_sys.HT] ./ one_foot,
        alpha=0.3, color=:blue, label="")
end

# Add deterministic trajectory
plot!(p1, years, dbh_in, color=:red, linewidth=2, label="Deterministic")
plot!(p2, years, ht_ft, color=:red, linewidth=2, label="Deterministic")

p = plot(p1, p2, layout=(1, 2), size=(900, 400))
p
```

### Response to Site Index

The diameter growth equation includes site index (SI) as a predictor variable,
with coefficient 0.4143 on ln(SI). Higher site index indicates better growing
conditions and should produce faster diameter growth.

```@example stage_prognosis
si_values = [50.0, 60.0, 70.0, 80.0, 90.0, 100.0]

p = plot(xlabel="Years", ylabel="DBH (inches)",
    title="Diameter Growth by Site Index")

for si in si_values
    prob_si = remake(prob_det; p=[compiled_sys.σ_growth => 0.0,
        compiled_sys.SI => si * one_foot])
    sol_si = solve(prob_si, EM(), dt=one_year/10)
    plot!(p, sol_si.t ./ one_year, sol_si[compiled_sys.DBH] ./ one_inch,
        label="SI = $(Int(si)) ft", linewidth=2)
end
p
```

### Response to Stand Density

The crown competition factor (CCF) reduces diameter growth through the coefficient
-0.3781 on ln(CCF). Denser stands experience more competition, leading to slower
individual tree growth but potentially higher total stand yield.

```@example stage_prognosis
ccf_values = [50.0, 75.0, 100.0, 150.0, 200.0]

p = plot(xlabel="Years", ylabel="DBH (inches)",
    title="Diameter Growth by Crown Competition Factor")

for ccf in ccf_values
    prob_ccf = remake(prob_det; p=[compiled_sys.σ_growth => 0.0,
        compiled_sys.CCF => ccf])
    sol_ccf = solve(prob_ccf, EM(), dt=one_year/10)
    plot!(p, sol_ccf.t ./ one_year, sol_ccf[compiled_sys.DBH] ./ one_inch,
        label="CCF = $(Int(ccf))", linewidth=2)
end
p
```

### Limitations

1. **Height increment coefficients**: The specific regression coefficients for the
   height increment equation were referenced to an unpublished manuscript
   (Cole and Stage, in preparation) and are not provided in Stage (1973). The
   default values used here are calibrated to approximate the growth trajectories
   shown in Figure 1B of the paper.

2. **Mortality rates**: The base endemic mortality rates from Lee (1971) are
   described qualitatively but not tabulated in Stage (1973). The implementation
   uses a quadratic approximation with minimum at 10.6 inches mean DBH.

3. **Beetle-induced mortality**: The mountain pine beetle mortality model
   (Figure 3 in the paper) is not implemented, as the specific sub-functions
   are not provided in sufficient detail.

4. **Self-calibration**: The self-calibration feature described in Stage (1973)
   for adjusting growth predictions based on observed past growth is not
   implemented in this SDE formulation.
