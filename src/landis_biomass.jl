"""
    LANDISBiomass(; name=:LANDISBiomass)

A forest growth and biomass module based on the LANDIS landscape simulation model.

This component implements the biomass module from Scheller and Mladenoff (2004),
which calculates aboveground net primary productivity (ANPP), mortality
(biomass-related and age-related), and dead woody biomass decomposition for a
single species-age cohort.

The model tracks living biomass and dead woody biomass for a species-age cohort
using three core processes: growth (ANPP), mortality, and decomposition.
Growth follows a peaked function of the ratio of actual to potential biomass.
Mortality combines a logistic biomass-dependent term and an exponential
age-dependent term. Dead biomass decays exponentially.

Units: The paper uses Mg/ha and years. This implementation uses SI base units
(kg/m² for biomass density, seconds for time). Conversion: 1 Mg/ha = 0.1 kg/m²,
1 year = 3.15576e7 seconds. All parameter values are specified in SI base units.

**Reference**: Scheller, R.M. and Mladenoff, D.J. (2004). A forest growth and
biomass module for a landscape simulation model, LANDIS: design, validation,
and application. *Ecological Modelling*, 180, 211-229.
doi:10.1016/j.ecolmodel.2004.01.022

$(SIGNATURES)
"""
@component function LANDISBiomass(; name = :LANDISBiomass)

    # Conversion factors:
    # 1 Mg/ha = 0.1 kg/m²
    # 1 year = 3.15576e7 seconds
    # 1 Mg/ha/yr = 0.1 / 3.15576e7 kg/m²/s ≈ 3.169e-9 kg/m²/s

    # --- Constants ---
    @constants begin
        e1 = exp(1), [description = "Euler's number (dimensionless)"]
        # 30 years in seconds: ratio of B_MAX to ANPP_MAX (Eq. 2)
        B_MAX_ratio = 30.0 * 3.15576e7, [description = "Ratio of maximum biomass to maximum ANPP (Eq. 2)", unit = u"s"]
        # 1 year in seconds: for converting per-year rates to per-second
        one_yr = 3.15576e7, [description = "One year in seconds for rate conversion", unit = u"s"]
    end

    # --- Parameters ---
    @parameters begin
        # ANPP_MAX = 7.45 Mg/ha/yr = 7.45 * 0.1 / 3.15576e7 kg/m²/s
        ANPP_MAX = 7.45 * 0.1 / 3.15576e7, [description = "Maximum ANPP for species (7.45 Mg/ha/yr for A. saccharum)", unit = u"kg/m^2/s"]
        # max_age = 400 years = 400 * 3.15576e7 seconds
        max_age = 400.0 * 3.15576e7, [description = "Maximum species longevity (400 years for A. saccharum)", unit = u"s"]
        r = 0.08, [description = "Mortality rate growth parameter (dimensionless)"]
        y0 = 0.01, [description = "Initial mortality rate, rescaled 0-1 (dimensionless)"]
        d = 10.0, [description = "Age-related mortality shape parameter (dimensionless)"]
        # k = 0.03/yr = 0.03 / 3.15576e7 s^-1
        k = 0.03 / 3.15576e7, [description = "Dead biomass decomposition rate (0.03/yr)", unit = u"s^-1"]
        # B_MAX_site = 500 Mg/ha = 500 * 0.1 = 50 kg/m²
        B_MAX_site = 500.0 * 0.1, [description = "Maximum site biomass over all cohorts (500 Mg/ha)", unit = u"kg/m^2"]
        B_other = 0.0, [description = "Biomass of all other cohorts at the site", unit = u"kg/m^2"]
        # Initial cohort age = 10 years = 10 * 3.15576e7 seconds
        age_init = 10.0 * 3.15576e7, [description = "Initial cohort age (10 years)", unit = u"s"]
    end

    # --- Variables ---
    @variables begin
        # B = 5.0 Mg/ha = 5.0 * 0.1 = 0.5 kg/m² (viable initial biomass for positive net growth)
        B(t) = 5.0 * 0.1, [description = "Aboveground living biomass of cohort", unit = u"kg/m^2"]
        D_wood(t) = 0.0, [description = "Dead woody biomass", unit = u"kg/m^2"]
        cohort_age(t), [description = "Cohort age (age_init + t)", unit = u"s"]
        B_MAX(t), [description = "Maximum possible biomass for species", unit = u"kg/m^2"]
        B_POT(t), [description = "Potential biomass (available growing space)", unit = u"kg/m^2"]
        B_AP(t), [description = "Ratio of actual to potential biomass (dimensionless)"]
        B_PM(t), [description = "Ratio of potential to maximum biomass (dimensionless)"]
        ANPP_ACT(t), [description = "Actual ANPP", unit = u"kg/m^2/s"]
        M_BIO(t), [description = "Biomass-related mortality rate", unit = u"kg/m^2/s"]
        M_AGE(t), [description = "Age-related mortality rate", unit = u"kg/m^2/s"]
        M_total(t), [description = "Total mortality rate", unit = u"kg/m^2/s"]
    end

    eqs = [
        # Cohort age = initial age + elapsed time
        cohort_age ~ age_init + t,

        # Eq. 2 - Maximum biomass from maximum ANPP
        B_MAX ~ ANPP_MAX * B_MAX_ratio,

        # Eq. 3 - Potential biomass: available growing space minus space occupied by others
        B_POT ~ max(B, min(B_MAX, B_MAX_site - B_other)),

        # Ratio of actual to potential biomass (Eq. 4 auxiliary)
        B_AP ~ B / B_POT,

        # Ratio of potential to maximum biomass (Eq. 4 auxiliary)
        B_PM ~ B_POT / B_MAX,

        # Eq. 4 - Actual ANPP: peaked growth function
        ANPP_ACT ~ ANPP_MAX * (e1 * B_AP * exp(-B_AP)) * B_PM,

        # Eq. 5 - Biomass-related mortality (logistic)
        # Note: The exponent in the paper is r*y0^(-1)*B_AP (r divided by y0)
        M_BIO ~ ANPP_MAX * (y0 / (y0 + (1 - y0) * exp(-r / y0 * B_AP))) * B_PM,

        # Eq. 6 - Age-related mortality (exponential increase with age)
        M_AGE ~ (B / one_yr) * exp((cohort_age / max_age) * d) / exp(d),

        # Total mortality
        M_total ~ M_BIO + M_AGE,

        # Eq. 1 - Biomass dynamics: dB/dt = ANPP_ACT - M_total
        D(B) ~ ANPP_ACT - M_total,

        # Eq. 7 - Dead biomass: gains from mortality, losses from decomposition
        D(D_wood) ~ M_total - k * D_wood,
    ]

    return System(eqs, t; name)
end
