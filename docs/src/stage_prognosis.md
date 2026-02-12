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

```@docs
StagePrognosisHCB
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

### Predictor Variables (Stage 1973, p. 15)

The following table summarizes the tree growth components, their predictor variables,
and data sources as described in Stage (1973):

| Tree growth component | Predictor variables | Data source |
|---|---|---|
| Annual basal area increment (b.a.i.) | D.b.h., relative stand density, site, elevation, habitat type, percentile in basal area distribution, crown ratio | Increment cores, remeasured plots |
| Height increment | Radial increment, habitat type, d.b.h., height | Stem analyses |
| Crown dimensions | Relative density, percentile in basal area distribution, d.b.h. | Temporary plots |
| Bark ratio | Same as b.a.i. | Temporary plots or tree samples |
| Mortality rates | Same as b.a.i. and radial increment plus pest population models where applicable | Remeasured plots, "last *n* years mortality", "years since death" |

### Component Equations

The model implements the following equations from Stage (1973):

1. **Diameter Change Model** (p. 15, Eq. 1): Basal area increment (DDS) is predicted from
   a log-linear regression on site index, elevation, crown competition factor, crown ratio,
   DBH, and basal area percentile.

2. **Diameter Growth** (p. 15): The actual diameter increment is computed as
   DG = √(DBH² + DDS) - DBH, matching the paper's exact formulation.

3. **Height Increment Model** (p. 15-16, Eq. 2): Height growth is predicted from DG,
   DBH, and total height using species-specific coefficients from Cole and Stage
   (in preparation, referenced in Stage 1973).

4. **Crown Ratio Development** (p. 16, Eq. 3): Crown base rises at 1/5 of the height
   increment rate when CCF < 125, and at 0.61 times the height increment rate when
   CCF ≥ 125.

5. **Endemic Mortality** (p. 16-17, Eq. 4): Based on Lee (1971), mortality rates depend
   on mean stand DBH with minimum at 10.6 inches. Mortality is distributed among trees
   by the factor [0.25 + 1.5·(1 - PCT/100)].

6. **HCB Prediction** (p. 16): A static equation for estimating crown base height when
   field crown measurements are unavailable.

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

### Comparison with Figure 1B (Stage, 1973)

Figure 1B from Stage (1973) shows the development of five sample trees at the
boundaries of the initial diameter distribution for a lodgepole pine stand.
The following table compares the paper's values for the 50th percentile (median) tree
with the model predictions using default parameters.

The paper reports for the 50th percentile tree:
- 1969: DBH = 6.8 in, Height = 67.0 ft, Crown Ratio = 35%
- 1980: DBH = 7.2 in, Height = 74.2 ft, Crown Ratio = 34%
- 1990: DBH = 7.5 in, Height = 80.6 ft, Crown Ratio = 34%
- 2000: DBH = 8.2 in, Height = 87.3 ft, Crown Ratio = 37%

```@example stage_prognosis
using StochasticDiffEq
using Plots

one_year = 3.15576e7
one_inch = 0.0254
one_foot = 0.3048

# Set initial conditions to match the 0.5 percentile tree in Figure 1B
# DBH = 6.8 in, HT = 67.0 ft, crown ratio = 0.35, 509 trees/acre, SI = 80 ft
cr0 = 0.35
ht0 = 67.0 * one_foot
hcb0 = ht0 * (1 - cr0)

prob = SDEProblem(compiled_sys, [
    compiled_sys.Dsq => (6.8 * one_inch)^2,
    compiled_sys.HT => ht0,
    compiled_sys.HCB => hcb0,
    compiled_sys.N_trees => 509.0 / 4046.86,
], (0.0, 51.0 * one_year))
prob_det = remake(prob; p=[compiled_sys.σ_growth => 0.0, compiled_sys.PCT => 50.0])
sol = solve(prob_det, EM(), dt=one_year/10)

# Extract values at 10-year intervals
years = sol.t ./ one_year
idx_0 = 1
idx_11 = findfirst(y -> y >= 11.0, years)  # 1969→1980
idx_21 = findfirst(y -> y >= 21.0, years)  # 1969→1990
idx_31 = findfirst(y -> y >= 31.0, years)  # 1969→2000
idx_51 = findfirst(y -> y >= 51.0, years)  # 1969→2020

paper_years = [1969, 1980, 1990, 2000, 2020]
paper_dbh = [6.8, 7.2, 7.5, 8.2, 8.8]  # inches
paper_ht = [67.0, 74.2, 80.6, 87.3, 99.5]  # feet

model_idx = [idx_0, idx_11, idx_21, idx_31, idx_51]
model_dbh = [sol[compiled_sys.DBH][i] / one_inch for i in model_idx]
model_ht = [sol[compiled_sys.HT][i] / one_foot for i in model_idx]

DataFrame(
    :Year => paper_years,
    Symbol("Paper DBH (in)") => paper_dbh,
    Symbol("Model DBH (in)") => round.(model_dbh, digits=1),
    Symbol("Paper HT (ft)") => paper_ht,
    Symbol("Model HT (ft)") => round.(model_ht, digits=1),
)
```

### Deterministic Growth Trajectory

The following figure shows the deterministic growth trajectory (no stochastic noise)
for a lodgepole pine tree starting from conditions similar to Figure 1 in Stage (1973):
DBH = 6.0 inches, height = 60 feet, crown ratio = 0.33, 500 trees/acre,
site index = 80 feet, elevation = 4300 feet.

```@example stage_prognosis
prob2 = SDEProblem(compiled_sys, [], (0.0, 50.0 * one_year))
prob2_det = remake(prob2; p=[compiled_sys.σ_growth => 0.0])
sol2 = solve(prob2_det, EM(), dt=one_year/10)

years2 = sol2.t ./ one_year
dbh_in = sol2[compiled_sys.DBH] ./ one_inch
ht_ft = sol2[compiled_sys.HT] ./ one_foot
cr = sol2[compiled_sys.CR]
n_acre = sol2[compiled_sys.N_trees] .* 4046.86

p1 = plot(years2, dbh_in, xlabel="Years", ylabel="DBH (inches)",
    title="Diameter Growth", legend=false, linewidth=2)
p2 = plot(years2, ht_ft, xlabel="Years", ylabel="Height (feet)",
    title="Height Growth", legend=false, linewidth=2)
p3 = plot(years2, cr, xlabel="Years", ylabel="Crown Ratio",
    title="Crown Development", legend=false, linewidth=2)
p4 = plot(years2, n_acre, xlabel="Years", ylabel="Trees/Acre",
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
    sol_s = solve(prob2, EM(), dt=one_year/10, seed=i)
    plot!(p1, sol_s.t ./ one_year, sol_s[compiled_sys.DBH] ./ one_inch,
        alpha=0.3, color=:blue, label="")
    plot!(p2, sol_s.t ./ one_year, sol_s[compiled_sys.HT] ./ one_foot,
        alpha=0.3, color=:blue, label="")
end

# Add deterministic trajectory
plot!(p1, years2, dbh_in, color=:red, linewidth=2, label="Deterministic")
plot!(p2, years2, ht_ft, color=:red, linewidth=2, label="Deterministic")

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
    prob_si = remake(prob2_det; p=[compiled_sys.σ_growth => 0.0,
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
    prob_ccf = remake(prob2_det; p=[compiled_sys.σ_growth => 0.0,
        compiled_sys.CCF => ccf])
    sol_ccf = solve(prob_ccf, EM(), dt=one_year/10)
    plot!(p, sol_ccf.t ./ one_year, sol_ccf[compiled_sys.DBH] ./ one_inch,
        label="CCF = $(Int(ccf))", linewidth=2)
end
p
```

### Endemic Mortality Rate (Stage 1973, p. 16-17)

The endemic mortality model is based on Lee (1971) yield tables for lodgepole pine
in Alberta. Annual mortality rates depend on mean stand DBH, declining with increasing
mean DBH to a minimum at 10.6 inches and then increasing. Within the stand, mortality
is distributed by the percentile factor [0.25 + 1.5·(1 - PCT/100)], giving the lowest
rate to the largest tree (PCT = 100) and the highest to the smallest (PCT = 0).

```@example stage_prognosis
# Mortality rate as a function of mean stand DBH
dbh_range = 4.0:0.1:16.0  # inches

mort_a0 = 0.0536
mort_a1 = -0.0088
mort_a2 = 0.000415

base_rate = [mort_a0 + mort_a1 * d + mort_a2 * d^2 for d in dbh_range]

p1 = plot(collect(dbh_range), base_rate .* 100,
    xlabel="Mean Stand DBH (inches)", ylabel="Annual Mortality Rate (%)",
    title="Base Mortality Rate (Lee 1971)",
    legend=false, linewidth=2)
vline!(p1, [10.6], linestyle=:dash, color=:gray, label="")

# Percentile adjustment factor
pct_range = 0.0:1.0:100.0
pct_factor = [0.25 + 1.5 * (1.0 - p / 100.0) for p in pct_range]

p2 = plot(collect(pct_range), pct_factor,
    xlabel="Percentile in Basal Area Distribution",
    ylabel="Mortality Factor",
    title="Percentile Adjustment (Stage 1973, p. 17)",
    legend=false, linewidth=2)

p = plot(p1, p2, layout=(1, 2), size=(900, 400))
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
