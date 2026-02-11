"""
    StagePrognosis(; name=:StagePrognosis)

A stochastic individual-tree growth model for lodgepole pine (*Pinus contorta*)
based on the Prognosis Model for Stand Development by Stage (1973).

The model predicts periodic changes in tree diameter (DBH), height, crown ratio,
and survival for individual trees in a forest stand. It is the foundation of the
USDA Forest Service's Forest Vegetation Simulator (FVS).

The core growth equation predicts basal area increment as a function of tree size,
site quality, stand density, competitive position, and crown development. Stochastic
variation is represented by multiplicative noise on the basal area increment rate,
implemented as a Brownian motion term in the SDE formulation.

The model tracks the squared diameter (Dsq = DBH²) as a state variable because
basal area increment is most naturally linear in diameter-squared space.

Units: The original paper uses imperial units (inches for DBH, feet for height and
elevation, trees per acre). This implementation uses SI units (meters for lengths,
m⁻² for tree density) with explicit conversion constants.

**Reference**: Stage, A.R. (1973). Prognosis model for stand development.
USDA Forest Service Research Paper INT-137.
Intermountain Forest and Range Experiment Station, Ogden, Utah. 32 p.

$(SIGNATURES)
"""
@component function StagePrognosis(; name = :StagePrognosis)

    # ===== Unit conversion constants =====
    # The regression coefficients were fit using imperial units.
    # We define reference quantities so that dividing a SI quantity by the
    # reference yields a dimensionless number equal to the imperial value.
    @constants begin
        one_inch = 0.0254, [description = "Reference length: 1 inch", unit = u"m"]
        one_foot = 0.3048, [description = "Reference length: 1 foot", unit = u"m"]
        hundred_feet = 30.48, [description = "Reference length: 100 feet", unit = u"m"]
        one_inch_sq = 0.0254^2, [description = "Reference area: 1 square inch", unit = u"m^2"]
        one_year = 3.15576e7, [description = "Reference time: 1 year", unit = u"s"]
        # 1 acre = 4046.86 m², so trees/acre * (1 acre / 4046.86 m²) = trees/m²
        one_acre = 4046.86, [description = "Reference area: 1 acre", unit = u"m^2"]
    end

    # ===== Parameters =====
    @parameters begin
        # Site parameters
        SI = 80.0 * 0.3048, [description = "Site index (dominant height at base age, converted from 80 ft)", unit = u"m"]
        EL = 43.0 * 30.48, [description = "Elevation above sea level (converted from 43 hundreds of feet)", unit = u"m"]

        # Stand density parameters
        CCF = 100.0, [description = "Crown competition factor (dimensionless)"]
        PCT = 50.0, [description = "Percentile in basal area distribution (0-100, dimensionless)"]
        RMSQD = 7.0 * 0.0254, [description = "Diameter of tree of mean basal area", unit = u"m"]

        # Tree parameters
        BKR = 1.0, [description = "Bark ratio (dob/dib)² (dimensionless)"]

        # Habitat type parameter for HCB equation
        HAB = 0.0, [description = "Habitat type adjustment for HCB (dimensionless, in feet-space)"]
        # HAB values: 0.0 for Abies/Xerophyllum, -4.24 for Abies/Vaccinium,
        #            -3.86 for Abies/Pachistima, -5.47 for Pseudotsuga/Calamagrostis

        # Height increment coefficients (species-specific, lodgepole pine estimates)
        # From Cole and Stage (in preparation), referenced in Stage (1973).
        # These defaults are calibrated to reproduce ~0.65 ft/yr height growth
        # for a 6-inch DBH, 60-foot tall lodgepole pine (consistent with Figure 1B).
        c1_ht = 0.0, [description = "Height increment intercept coefficient (dimensionless)"]
        c2_ht = 0.5, [description = "Height increment coefficient on ln(DG+0.05) (dimensionless)"]
        c3_ht = -0.1, [description = "Height increment coefficient on ln(DBH) (dimensionless)"]
        c4_ht = 0.2, [description = "Height increment coefficient on ln(HT) (dimensionless)"]

        # Endemic mortality parameters
        # Approximation of Lee (1971) mortality curve for lodgepole pine
        # Base annual mortality rate as function of mean DBH (in inches):
        #   rate = mort_a0 + mort_a1*DBH_mean_in + mort_a2*DBH_mean_in²
        # Minimum at DBH = 10.6 inches (0.269 m)
        mort_a0 = 0.0536, [description = "Mortality intercept coefficient (dimensionless)"]
        mort_a1 = -0.0088, [description = "Mortality linear coefficient (dimensionless)"]
        mort_a2 = 0.000415, [description = "Mortality quadratic coefficient (dimensionless)"]
        mean_DBH = 7.0 * 0.0254, [description = "Stand mean DBH for mortality calculation", unit = u"m"]

        # Stochastic parameters
        # The paper specifies σ ≈ 0.3 in log-space per year. Since the time unit is
        # seconds, we scale: σ_growth = σ_paper * sqrt(seconds_per_year)
        # = 0.3 * sqrt(3.15576e7) ≈ 1685.4. This ensures the noise over one year
        # gives the correct log-space standard deviation.
        σ_growth = 0.3 * sqrt(3.15576e7), [description = "Scaled std dev of stochastic growth noise (dimensionless)"]
    end

    # ===== State Variables =====
    @variables begin
        Dsq(t) = (6.0 * 0.0254)^2, [description = "Squared diameter at breast height (DBH²)", unit = u"m^2"]
        HT(t) = 60.0 * 0.3048, [description = "Total tree height", unit = u"m"]
        HCB(t) = 40.0 * 0.3048, [description = "Height to crown base", unit = u"m"]
        N_trees(t) = 500.0 / 4046.86, [description = "Tree density (trees per m²)", unit = u"m^-2"]
    end

    # ===== Derived variables (algebraic) =====
    @variables begin
        DBH(t), [description = "Diameter at breast height", unit = u"m"]
        CR(t), [description = "Crown ratio (dimensionless)"]
        DDS_rate(t), [description = "Rate of change of diameter squared", unit = u"m^2/s"]
        HTGF_rate(t), [description = "Rate of height growth", unit = u"m/s"]
        crown_recession_rate(t), [description = "Rate of crown base rise", unit = u"m/s"]
        mortality_rate(t), [description = "Mortality rate", unit = u"s^-1"]
        mean_DBH_in(t), [description = "Stand mean DBH in inches (dimensionless ratio)"]
    end

    # ===== Brownian motion for stochastic diameter growth =====
    @brownians B_growth

    eqs = [
        # --- Derived quantities ---

        # DBH from squared diameter
        DBH ~ sqrt(Dsq),  # Auxiliary

        # Crown ratio from height and crown base height
        CR ~ max(0.01, (HT - HCB) / HT),  # Auxiliary, bounded away from zero

        # Mean DBH in inches (dimensionless) for mortality calculation
        mean_DBH_in ~ mean_DBH / one_inch,

        # --- Diameter Change Model (Stage 1973, p. 15) ---
        # DDS = BKR * FINT * exp(-1.66955 + 0.4143*ln(SI) - 0.004388*EL
        #        - 0.3781*ln(CCF) + 0.4879*ln(CR) + 0.9948*ln(DBH) + 0.006141*PCT)
        # Original equation uses imperial units. We convert: DBH/one_inch gives
        # dimensionless value numerically equal to DBH in inches, etc.
        # DDS has units of in²/year, converted back to m²/s.
        DDS_rate ~ (BKR * one_inch_sq / one_year) * exp(
            -1.66955
            + 0.4143 * log(SI / one_foot)       # SI in feet
            - 0.004388 * (EL / hundred_feet)     # EL in hundreds of feet
            - 0.3781 * log(CCF)                  # CCF dimensionless
            + 0.4879 * log(CR)                   # CR dimensionless
            + 0.9948 * log(DBH / one_inch)       # DBH in inches
            + 0.006141 * PCT                     # PCT dimensionless
        ),  # Eq. 1 - Basal area increment rate

        # --- Height Increment Model (Stage 1973, p. 15-16) ---
        # ln(HTGF) = c1 + c2*ln(DG+0.05) + c3*ln(DBH) + c4*ln(HT)
        # DG is periodic diameter increment in inches.
        # We approximate DG from the instantaneous DDS_rate:
        #   DG ≈ DDS_rate * one_year / (2*DBH) (small increment approximation)
        #   converted to inches: DG_in = DG / one_inch
        # HTGF is periodic height growth in feet/year, converted to m/s
        HTGF_rate ~ (one_foot / one_year) * exp(
            c1_ht
            + c2_ht * log(max(0.01,
                (DDS_rate * one_year / one_inch_sq) / (2.0 * DBH / one_inch)
            ) + 0.05)
            + c3_ht * log(DBH / one_inch)        # DBH in inches
            + c4_ht * log(HT / one_foot)          # HT in feet
        ),  # Eq. 2 - Height growth rate

        # --- Crown Ratio Development (Stage 1973, p. 16) ---
        # Crown base recedes (rises) at a fraction of height increment rate
        # If CCF < 125: rate = 1/5 * height growth rate
        # If CCF >= 125: rate = 0.61 * height growth rate
        crown_recession_rate ~ ifelse(CCF < 125.0,
            0.2 * HTGF_rate,    # 1/5 of height growth
            0.61 * HTGF_rate    # 0.61 of height growth
        ),  # Eq. 3 - Crown recession

        # --- Endemic Mortality Model (Stage 1973, p. 16-17) ---
        # Based on Lee (1971): rates depend on mean DBH with minimum at 10.6 inches
        # Factor distributes mortality by percentile: [0.25 + 1.5*(1 - PCT/100)]
        mortality_rate ~ (1.0 / one_year) * max(0.0,
            (mort_a0 + mort_a1 * mean_DBH_in + mort_a2 * mean_DBH_in * mean_DBH_in)
            * (0.25 + 1.5 * (1.0 - PCT / 100.0))
        ),  # Eq. 4 - Mortality rate

        # ===== Differential Equations =====

        # Stochastic diameter squared growth (Stage 1973, p. 12, 15)
        # Random error added to ln(basal area increment) → multiplicative noise.
        # The Brownian term B_growth is dimensionless in MTK's unit system;
        # σ_growth is pre-scaled by sqrt(seconds_per_year) so that over one year
        # the integrated noise has the correct standard deviation.
        D(Dsq) ~ DDS_rate + σ_growth * abs(DDS_rate) * B_growth,

        # Height growth
        D(HT) ~ HTGF_rate,  # Deterministic height growth

        # Crown base height growth
        D(HCB) ~ crown_recession_rate,  # Crown recession

        # Survival (tree density decreases with mortality)
        D(N_trees) ~ -mortality_rate * N_trees,  # Mortality
    ]

    return System(eqs, t; name)
end
