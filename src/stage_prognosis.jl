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
        one_acre = 4046.86, [description = "Reference area: 1 acre", unit = u"m^2"]
    end

    # ===== Diameter growth regression constants (Stage 1973, p. 15) =====
    # DDS = BKR * FINT * exp(b0 + b1*ln(SI) + b2*EL + b3*ln(CCF) + b4*ln(CR) + b5*ln(DBH) + b6*PCT)
    # Regression fitted on ln(DDS) in imperial units (in²/yr);
    # all predictor variables are non-dimensionalized by reference values, so
    # these coefficients are dimensionless.
    @constants begin
        b_DDS_0 = -1.66955, [description = "DDS regression intercept (dimensionless)"]
        b_DDS_SI = 0.4143, [description = "DDS coefficient on ln(SI in feet) (dimensionless)"]
        b_DDS_EL = -0.004388, [description = "DDS coefficient on EL in hundreds of feet (dimensionless)"]
        b_DDS_CCF = -0.3781, [description = "DDS coefficient on ln(CCF) (dimensionless)"]
        b_DDS_CR = 0.4879, [description = "DDS coefficient on ln(CR) (dimensionless)"]
        b_DDS_DBH = 0.9948, [description = "DDS coefficient on ln(DBH in inches) (dimensionless)"]
        b_DDS_PCT = 0.006141, [description = "DDS coefficient on PCT (dimensionless)"]
    end

    # ===== Crown recession constants (Stage 1973, p. 16) =====
    @constants begin
        crown_recess_low = 0.2, [description = "Crown recession fraction when CCF < 125 (dimensionless)"]
        crown_recess_high = 0.61, [description = "Crown recession fraction when CCF >= 125 (dimensionless)"]
        CCF_threshold = 125.0, [description = "CCF threshold for crown recession rate (dimensionless)"]
    end

    # ===== Endemic mortality constants (Lee 1971, via Stage 1973, p. 16-17) =====
    # Base annual mortality rate as function of mean stand DBH (in inches):
    #   rate_per_year = mort_a0 + mort_a1 * DBH_mean_in + mort_a2 * DBH_mean_in²
    # Minimum at DBH = 10.6 inches (0.269 m). All inputs are non-dimensionalized,
    # so these coefficients are dimensionless annual rates.
    @constants begin
        mort_a0 = 0.0536, [description = "Mortality intercept, annual rate (dimensionless)"]
        mort_a1 = -0.0088, [description = "Mortality linear coeff on DBH in inches (dimensionless)"]
        mort_a2 = 0.000415, [description = "Mortality quadratic coeff on DBH² in inches² (dimensionless)"]
    end

    # ===== Literal constants used in equations =====
    @constants begin
        CR_min = 0.01, [description = "Minimum crown ratio floor (dimensionless)"]
        DG_min = 0.001, [description = "Minimum diameter growth increment floor (dimensionless)"]
        DG_offset = 0.05, [description = "Offset for DG in height growth log transform (dimensionless)"]
        zero_rate = 0.0, [unit = u"s^-1", description = "Zero rate for mortality floor"]
        one_rate_per_year = 1.0 / 3.15576e7, [unit = u"s^-1", description = "Unit annual rate converted to per-second"]
        mort_base_factor = 0.25, [description = "Base mortality distribution factor (dimensionless)"]
        mort_scale_factor = 1.5, [description = "Mortality distribution scale factor (dimensionless)"]
        mort_pct_one = 1.0, [description = "Unity constant for PCT normalization (dimensionless)"]
        mort_pct_max = 100.0, [description = "Maximum percentile value (dimensionless)"]
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

        # Height increment coefficients (species-specific, lodgepole pine estimates)
        # From Cole and Stage (in preparation), referenced in Stage (1973).
        # These defaults are calibrated to reproduce ~0.65 ft/yr height growth
        # for a 6-inch DBH, 60-foot tall lodgepole pine (consistent with Figure 1B).
        c1_ht = 0.0, [description = "Height increment intercept coefficient (dimensionless)"]
        c2_ht = 0.5, [description = "Height increment coefficient on ln(DG+0.05) (dimensionless)"]
        c3_ht = -0.1, [description = "Height increment coefficient on ln(DBH) (dimensionless)"]
        c4_ht = 0.2, [description = "Height increment coefficient on ln(HT) (dimensionless)"]

        # Mortality parameters
        mean_DBH = 7.0 * 0.0254, [description = "Stand mean DBH for mortality calculation", unit = u"m"]

        # Stochastic parameters
        # The paper specifies σ ≈ 0.3 in log-space per year (Stage 1973, p. 12).
        # In the SDE formulation D(Dsq) ~ drift + σ_growth * |drift| * B, the
        # relative noise over one year has std dev = σ_growth / √(one_year).
        # Setting σ_growth = 0.3 * √(one_year) gives the correct 0.3 relative
        # std dev over one year.
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
        DG(t), [description = "Periodic diameter increment (dimensionless, in inches-space)"]
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
        CR ~ max(CR_min, (HT - HCB) / HT),  # Auxiliary, bounded away from zero

        # Mean DBH in inches (dimensionless) for mortality calculation
        mean_DBH_in ~ mean_DBH / one_inch,

        # --- Diameter Change Model (Stage 1973, p. 15) ---
        # DDS = BKR * FINT * exp(b0 + b1*ln(SI) + b2*EL + b3*ln(CCF)
        #        + b4*ln(CR) + b5*ln(DBH) + b6*PCT)
        # Original equation uses imperial units. We convert: DBH/one_inch gives
        # dimensionless value numerically equal to DBH in inches, etc.
        # DDS has units of in²/year, converted back to m²/s.
        DDS_rate ~ (BKR * one_inch_sq / one_year) * exp(
            b_DDS_0
                + b_DDS_SI * log(SI / one_foot)       # SI in feet
                + b_DDS_EL * (EL / hundred_feet)       # EL in hundreds of feet
                + b_DDS_CCF * log(CCF)                 # CCF dimensionless
                + b_DDS_CR * log(CR)                   # CR dimensionless
                + b_DDS_DBH * log(DBH / one_inch)      # DBH in inches
                + b_DDS_PCT * PCT                      # PCT dimensionless
        ),  # Eq. 1 - Basal area increment rate

        # --- Diameter Growth Increment (Stage 1973, p. 15) ---
        # DG = sqrt(DBH² + DDS) - DBH, all in inches
        # DDS_rate has units m²/s, multiply by one_year to get m²/yr,
        # then divide by one_inch_sq to get in²/yr (dimensionless).
        # DBH/one_inch gives DBH in inches (dimensionless).
        DG ~ max(DG_min, sqrt((DBH / one_inch)^2 + DDS_rate * one_year / one_inch_sq) - DBH / one_inch),

        # --- Height Increment Model (Stage 1973, p. 15-16) ---
        # ln(HTGF) = c1 + c2*ln(DG+0.05) + c3*ln(DBH) + c4*ln(HT)
        # where DG = periodic diameter increment in inches (dimensionless)
        # HTGF is periodic height growth in feet/year, converted to m/s
        HTGF_rate ~ (one_foot / one_year) * exp(
            c1_ht
                + c2_ht * log(DG + DG_offset)          # DG in inches (dimensionless)
                + c3_ht * log(DBH / one_inch)          # DBH in inches
                + c4_ht * log(HT / one_foot)           # HT in feet
        ),  # Eq. 2 - Height growth rate

        # --- Crown Ratio Development (Stage 1973, p. 16) ---
        # Crown base recedes (rises) at a fraction of height increment rate
        # If CCF < 125: rate = 1/5 * height growth rate
        # If CCF >= 125: rate = 0.61 * height growth rate
        crown_recession_rate ~ ifelse(
            CCF < CCF_threshold,
            crown_recess_low * HTGF_rate,    # 1/5 of height growth
            crown_recess_high * HTGF_rate    # 0.61 of height growth
        ),  # Eq. 3 - Crown recession

        # --- Endemic Mortality Model (Stage 1973, p. 16-17) ---
        # Based on Lee (1971): rates depend on mean DBH with minimum at 10.6 inches
        # Factor distributes mortality by percentile: [0.25 + 1.5*(1 - PCT/100)]
        mortality_rate ~ max(
            zero_rate,
            one_rate_per_year
                * (mort_a0 + mort_a1 * mean_DBH_in + mort_a2 * mean_DBH_in * mean_DBH_in)
                * (mort_base_factor + mort_scale_factor * (mort_pct_one - PCT / mort_pct_max))
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

"""
    StagePrognosisHCB(; name=:StagePrognosisHCB)

Height to crown base (HCB) prediction model for lodgepole pine from Stage (1973, p. 16).

This static model predicts the height to the base of the live crown as a function of
tree height, stand density, elevation, relative tree size, and habitat type. It is used
to estimate initial crown dimensions when field crown measurements are not available.

The equation is:

    HCB = -29.26 + 0.61*HT + 9.178*ln(CCF) - 0.222*EL - 5.80*DBH/RMSQD + HAB

where all variables are in imperial units (feet for heights/elevation, inches for diameters).
This implementation converts to SI units.

**Reference**: Stage, A.R. (1973). Prognosis model for stand development.
USDA Forest Service Research Paper INT-137, p. 16.

$(SIGNATURES)
"""
@component function StagePrognosisHCB(; name = :StagePrognosisHCB)
    # ===== Unit conversion constants =====
    @constants begin
        one_inch = 0.0254, [description = "Reference length: 1 inch", unit = u"m"]
        one_foot = 0.3048, [description = "Reference length: 1 foot", unit = u"m"]
        hundred_feet = 30.48, [description = "Reference length: 100 feet", unit = u"m"]
    end

    # ===== HCB regression constants (Stage 1973, p. 16) =====
    # HCB = b0 + b1*HT + b2*ln(CCF) + b3*EL + b4*DBH/RMSQD + HAB
    # All in feet-space; result converted to meters.
    @constants begin
        b_HCB_0 = -29.26, [description = "HCB regression intercept in feet (dimensionless)"]
        b_HCB_HT = 0.61, [description = "HCB coefficient on height in feet (dimensionless)"]
        b_HCB_CCF = 9.178, [description = "HCB coefficient on ln(CCF) (dimensionless)"]
        b_HCB_EL = -0.222, [description = "HCB coefficient on elevation in hundreds of feet (dimensionless)"]
        b_HCB_DBH_RMSQD = -5.8, [description = "HCB coefficient on DBH/RMSQD ratio (dimensionless)"]
    end

    @parameters begin
        HT_input = 60.0 * 0.3048, [description = "Total tree height", unit = u"m"]
        DBH_input = 6.0 * 0.0254, [description = "Diameter at breast height", unit = u"m"]
        CCF = 100.0, [description = "Crown competition factor (dimensionless)"]
        EL = 43.0 * 30.48, [description = "Elevation above sea level", unit = u"m"]
        RMSQD = 7.0 * 0.0254, [description = "Diameter of tree of mean basal area", unit = u"m"]
        HAB = 0.0, [description = "Habitat type adjustment (dimensionless, in feet-space)"]
        # HAB values: 0.0 for Abies/Xerophyllum, -4.24 for Abies/Vaccinium,
        #            -3.86 for Abies/Pachistima, -5.47 for Pseudotsuga/Calamagrostis
    end

    @variables begin
        HCB_pred(t), [description = "Predicted height to crown base", unit = u"m"]
    end

    # Stage (1973, p. 16): HCB = -29.26 + 0.61*HT + 9.178*ln(CCF) - 0.222*EL - 5.80*DBH/RMSQD + HAB
    # All terms in feet-space; convert result to meters
    eqs = [
        HCB_pred ~ one_foot * (
            b_HCB_0
                + b_HCB_HT * (HT_input / one_foot)
                + b_HCB_CCF * log(CCF)
                + b_HCB_EL * (EL / hundred_feet)
                + b_HCB_DBH_RMSQD * (DBH_input / one_inch) / (RMSQD / one_inch)
                + HAB
        ),  # HCB prediction (Stage 1973, p. 16)
    ]

    return System(eqs, t; name)
end
