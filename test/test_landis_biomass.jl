@testsnippet LANDISSetup begin
    using Test
    using ModelingToolkit
    using ModelingToolkit: t, D
    using OrdinaryDiffEqDefault
    using OrdinaryDiffEqDefault: SciMLBase
    using Vegetation

    # Unit conversion constants
    const Mg_ha_to_kg_m2 = 0.1        # 1 Mg/ha = 0.1 kg/m²
    const yr_to_s = 3.15576e7          # 1 year = 3.15576e7 seconds
end

@testitem "LANDISBiomass: Structural Verification" setup = [LANDISSetup] tags = [:landis] begin
    sys = LANDISBiomass()
    @test sys isa ModelingToolkit.System
    @test nameof(sys) == :LANDISBiomass

    vars = unknowns(sys)
    eqs = equations(sys)

    # 2 state variables (B, D_wood) + 9 algebraic variables = 11 unknowns
    @test length(vars) == 11
    # 11 equations (2 ODEs + 9 algebraic)
    @test length(eqs) == 11

    # Verify key variable names exist
    var_names = [string(v) for v in vars]
    for expected in [
            "B(t)", "D_wood(t)", "cohort_age(t)", "B_MAX(t)",
            "B_POT(t)", "B_AP(t)", "B_PM(t)", "ANPP_ACT(t)",
            "M_BIO(t)", "M_AGE(t)", "M_total(t)",
        ]
        @test any(n -> contains(n, expected), var_names)
    end

    # Verify compilation succeeds
    compiled = mtkcompile(sys)
    @test compiled !== nothing
    @test length(unknowns(compiled)) == 2  # B and D_wood are the only ODE states
end

@testitem "LANDISBiomass: Parameter Defaults" setup = [LANDISSetup] tags = [:landis] begin
    sys = LANDISBiomass()
    compiled = mtkcompile(sys)

    params = parameters(compiled)
    pdict = Dict(
        Symbol(p) => ModelingToolkit.getdefault(p) for p in params
            if ModelingToolkit.hasdefault(p)
    )

    # Verify ANPP_MAX = 7.45 Mg/ha/yr in SI
    @test pdict[:ANPP_MAX] ≈ 7.45 * Mg_ha_to_kg_m2 / yr_to_s rtol = 1.0e-6

    # Verify max_age = 400 years in seconds
    @test pdict[:max_age] ≈ 400.0 * yr_to_s rtol = 1.0e-6

    # Verify dimensionless parameters
    @test pdict[:r] ≈ 0.08
    @test pdict[:y0] ≈ 0.01
    @test pdict[:d] ≈ 10.0

    # Verify k = 0.03/yr in s^-1
    @test pdict[:k] ≈ 0.03 / yr_to_s rtol = 1.0e-6

    # Verify B_MAX_site = 500 Mg/ha in kg/m²
    @test pdict[:B_MAX_site] ≈ 500.0 * Mg_ha_to_kg_m2 rtol = 1.0e-6

    # Verify age_init = 10 years in seconds
    @test pdict[:age_init] ≈ 10.0 * yr_to_s rtol = 1.0e-6
end

@testitem "LANDISBiomass: Eq. 2 - Maximum Biomass" setup = [LANDISSetup] tags = [:landis] begin
    # B_MAX = ANPP_MAX * 30 years
    sys = LANDISBiomass()
    compiled = mtkcompile(sys)

    prob = ODEProblem(compiled, [], (0.0, 1.0))
    sol = solve(prob)

    # B_MAX should be 7.45 * 30 = 223.5 Mg/ha = 22.35 kg/m²
    B_MAX = sol[compiled.B_MAX][1]
    @test B_MAX ≈ 7.45 * 30.0 * Mg_ha_to_kg_m2 rtol = 1.0e-4
end

@testitem "LANDISBiomass: Eq. 4 - ANPP Growth Function Shape" setup = [LANDISSetup] tags = [:landis] begin
    # Verify the growth function peaks at B_AP = 1 by running the actual model
    sys = LANDISBiomass()
    compiled = mtkcompile(sys)

    # Set B = B_POT (so B_AP = 1, B_PM = 1) to test peak ANPP
    # B_MAX = ANPP_MAX * 30 yr = 7.45 * 0.1 / yr_to_s * 30 * yr_to_s = 22.35 kg/m²
    B_MAX_val = 7.45 * Mg_ha_to_kg_m2 / yr_to_s * 30.0 * yr_to_s
    prob_peak = ODEProblem(
        compiled,
        Dict(compiled.B => B_MAX_val, compiled.D_wood => 0.0),
        (0.0, 1.0)
    )
    sol_peak = solve(prob_peak)
    # At B_AP = 1, B_PM = 1: ANPP_ACT should equal ANPP_MAX * e * 1 * exp(-1) * 1 = ANPP_MAX
    ANPP_at_peak = sol_peak[compiled.ANPP_ACT][1]
    ANPP_MAX_val = 7.45 * Mg_ha_to_kg_m2 / yr_to_s
    @test ANPP_at_peak ≈ ANPP_MAX_val rtol = 1.0e-4

    # At lower B (B_AP < 1), ANPP should be less than at peak
    prob_low = ODEProblem(
        compiled,
        Dict(compiled.B => 0.5 * B_MAX_val, compiled.D_wood => 0.0),
        (0.0, 1.0)
    )
    sol_low = solve(prob_low)
    ANPP_at_low = sol_low[compiled.ANPP_ACT][1]
    @test ANPP_at_low < ANPP_at_peak
    # e * 0.5 * exp(-0.5) ≈ 0.8244
    @test ANPP_at_low / ANPP_MAX_val ≈ 0.8244 rtol = 0.01
end

@testitem "LANDISBiomass: Eq. 5 - Mortality Logistic Shape" setup = [LANDISSetup] tags = [:landis] begin
    # Verify the mortality logistic shape by running the actual model
    sys = LANDISBiomass()
    compiled = mtkcompile(sys)

    B_MAX_val = 7.45 * Mg_ha_to_kg_m2 / yr_to_s * 30.0 * yr_to_s
    ANPP_MAX_val = 7.45 * Mg_ha_to_kg_m2 / yr_to_s

    # At B_AP ≈ 1.0 (B = B_MAX, B_PM = 1): M_BIO/ANPP_MAX ≈ 0.968
    prob_full = ODEProblem(
        compiled,
        Dict(compiled.B => B_MAX_val, compiled.D_wood => 0.0),
        (0.0, 1.0)
    )
    sol_full = solve(prob_full)
    M_BIO_full = sol_full[compiled.M_BIO][1]
    @test M_BIO_full / ANPP_MAX_val ≈ 0.968 rtol = 0.01

    # At B_AP ≈ 0.4 (B = 0.4 * B_MAX, B_PM = 1): M_BIO/ANPP_MAX ≈ 0.199
    prob_mid = ODEProblem(
        compiled,
        Dict(compiled.B => 0.4 * B_MAX_val, compiled.D_wood => 0.0),
        (0.0, 1.0)
    )
    sol_mid = solve(prob_mid)
    M_BIO_mid = sol_mid[compiled.M_BIO][1]
    @test M_BIO_mid / ANPP_MAX_val ≈ 0.199 rtol = 0.01

    # At very low B_AP, M_BIO should approach y0 * ANPP_MAX * B_PM
    prob_low = ODEProblem(
        compiled,
        Dict(compiled.B => 0.001 * B_MAX_val, compiled.D_wood => 0.0),
        (0.0, 1.0)
    )
    sol_low = solve(prob_low)
    M_BIO_low = sol_low[compiled.M_BIO][1]
    @test M_BIO_low / ANPP_MAX_val ≈ 0.01 rtol = 0.1
end

@testitem "LANDISBiomass: Eq. 6 - Age-related Mortality Shape" setup = [LANDISSetup] tags = [:landis] begin
    # Verify M_AGE shape by running the actual model at different ages
    sys = LANDISBiomass()
    compiled = mtkcompile(sys)

    B_init = 10.0 * Mg_ha_to_kg_m2  # 10 Mg/ha

    # At age = 0.5 * max_age (200 years): M_AGE fraction ≈ exp(-d/2) ≈ 0.00674
    prob_half = ODEProblem(
        compiled,
        Dict(
            compiled.B => B_init, compiled.D_wood => 0.0,
            compiled.age_init => 200.0 * yr_to_s
        ),
        (0.0, 1.0)
    )
    sol_half = solve(prob_half)
    M_AGE_half = sol_half[compiled.M_AGE][1]
    # M_AGE = B/yr * exp(0.5*d)/exp(d) = B/yr * exp(-d/2)
    expected_half = B_init / yr_to_s * exp(-10.0 / 2)
    @test M_AGE_half ≈ expected_half rtol = 1.0e-4

    # At age = 0.8 * max_age (320 years): fraction ≈ 0.135 (Fig. 3b)
    prob_80 = ODEProblem(
        compiled,
        Dict(
            compiled.B => B_init, compiled.D_wood => 0.0,
            compiled.age_init => 320.0 * yr_to_s
        ),
        (0.0, 1.0)
    )
    sol_80 = solve(prob_80)
    M_AGE_80 = sol_80[compiled.M_AGE][1]
    expected_80 = B_init / yr_to_s * exp(0.8 * 10.0) / exp(10.0)
    @test M_AGE_80 ≈ expected_80 rtol = 1.0e-4
    # Fraction should be ≈ 0.135
    @test exp(0.8 * 10.0) / exp(10.0) ≈ 0.135 rtol = 0.01

    # M_AGE should increase with age
    @test M_AGE_80 > M_AGE_half
end

@testitem "LANDISBiomass: Single Species Growth (A. saccharum)" setup = [LANDISSetup] tags = [:landis] begin
    sys = LANDISBiomass()
    compiled = mtkcompile(sys)

    # Default parameters: A. saccharum (ANPP_MAX=7.45, max_age=400, d=10, r=0.08, y0=0.01)
    tspan_s = 200.0 * yr_to_s
    prob = ODEProblem(compiled, [], (0.0, tspan_s))
    sol = solve(prob)

    @test sol.retcode == SciMLBase.ReturnCode.Success

    # At ~100 years, biomass should be between 150-220 Mg/ha
    B_100 = sol(100.0 * yr_to_s; idxs = compiled.B) / Mg_ha_to_kg_m2
    @test 150.0 < B_100 < 220.0

    # Biomass should peak and then decline due to age-related mortality
    B_150 = sol(150.0 * yr_to_s; idxs = compiled.B) / Mg_ha_to_kg_m2
    B_200 = sol(200.0 * yr_to_s; idxs = compiled.B) / Mg_ha_to_kg_m2
    @test B_200 < B_150

    # Peak biomass should be around 190-210 Mg/ha
    B_vals = [sol(yr * yr_to_s; idxs = compiled.B) / Mg_ha_to_kg_m2 for yr in 10:10:200]
    @test 180.0 < maximum(B_vals) < 240.0
end

@testitem "LANDISBiomass: Eq. 3 - Potential Biomass" setup = [LANDISSetup] tags = [:landis] begin
    sys = LANDISBiomass()
    compiled = mtkcompile(sys)

    # Without competition, B_POT should equal B_MAX
    prob = ODEProblem(compiled, [], (0.0, 1.0))
    sol = solve(prob)
    B_POT_val = sol[compiled.B_POT][1]
    B_MAX_val = sol[compiled.B_MAX][1]
    @test B_POT_val ≈ B_MAX_val rtol = 1.0e-6

    # With heavy competition (B_other = 490 Mg/ha), B_POT should be reduced
    # B_MAX_site - B_other = 500 - 490 = 10 Mg/ha = 1.0 kg/m²
    # B_POT = max(B, min(B_MAX, B_MAX_site - B_other))
    #        = max(0.5, min(22.35, 1.0)) = max(0.5, 1.0) = 1.0 kg/m² = 10 Mg/ha
    prob_comp = ODEProblem(
        compiled,
        Dict(
            compiled.B => 5.0 * Mg_ha_to_kg_m2, compiled.D_wood => 0.0,
            compiled.B_other => 490.0 * Mg_ha_to_kg_m2
        ),
        (0.0, 1.0)
    )
    sol_comp = solve(prob_comp)
    B_POT_comp = sol_comp[compiled.B_POT][1]
    @test B_POT_comp ≈ 1.0 rtol = 1.0e-4  # 10 Mg/ha = 1.0 kg/m²

    # B_PM should reflect the reduced potential
    B_PM_comp = sol_comp[compiled.B_PM][1]
    @test B_PM_comp < 1.0  # B_POT < B_MAX implies B_PM < 1
end

@testitem "LANDISBiomass: Eq. 7 - Dead Biomass Decomposition" setup = [LANDISSetup] tags = [:landis] begin
    sys = LANDISBiomass()
    compiled = mtkcompile(sys)

    tspan_s = 200.0 * yr_to_s
    prob = ODEProblem(compiled, [], (0.0, tspan_s))
    sol = solve(prob)

    @test sol.retcode == SciMLBase.ReturnCode.Success

    # Paper (p. 221): Dead woody biomass near 25 Mg/ha after year 50
    D_50 = sol(50.0 * yr_to_s; idxs = compiled.D_wood) / Mg_ha_to_kg_m2
    @test 15.0 < D_50 < 40.0  # ~25 Mg/ha per paper

    # Increasing k should reduce steady-state dead biomass
    prob_fast_decay = ODEProblem(
        compiled,
        Dict(compiled.k => 0.06 / yr_to_s),
        (0.0, tspan_s)
    )
    sol_fast = solve(prob_fast_decay)
    D_fast = sol_fast(100.0 * yr_to_s; idxs = compiled.D_wood) / Mg_ha_to_kg_m2
    D_default = sol(100.0 * yr_to_s; idxs = compiled.D_wood) / Mg_ha_to_kg_m2
    @test D_fast < D_default
end

@testitem "LANDISBiomass: Dead Biomass Accumulation" setup = [LANDISSetup] tags = [:landis] begin
    sys = LANDISBiomass()
    compiled = mtkcompile(sys)

    tspan_s = 200.0 * yr_to_s
    prob = ODEProblem(compiled, [], (0.0, tspan_s))
    sol = solve(prob)

    @test sol.retcode == SciMLBase.ReturnCode.Success

    # Dead biomass should accumulate over time
    D_50 = sol(50.0 * yr_to_s; idxs = compiled.D_wood) / Mg_ha_to_kg_m2
    D_100 = sol(100.0 * yr_to_s; idxs = compiled.D_wood) / Mg_ha_to_kg_m2
    D_200 = sol(200.0 * yr_to_s; idxs = compiled.D_wood) / Mg_ha_to_kg_m2

    @test D_50 > 0
    @test D_100 > D_50
    @test D_200 > D_100
    @test D_200 < 500.0  # Sanity bound
end

@testitem "LANDISBiomass: Biomass Positivity" setup = [LANDISSetup] tags = [:landis] begin
    sys = LANDISBiomass()
    compiled = mtkcompile(sys)

    tspan_s = 200.0 * yr_to_s
    prob = ODEProblem(compiled, [], (0.0, tspan_s))
    sol = solve(prob)

    @test sol.retcode == SciMLBase.ReturnCode.Success

    # B should remain positive at all solution points
    B_vals = sol[compiled.B]
    @test all(B_vals .> 0)

    # D_wood should remain non-negative
    D_vals = sol[compiled.D_wood]
    @test all(D_vals .>= 0)
end

@testitem "LANDISBiomass: Competition via B_other" setup = [LANDISSetup] tags = [:landis] begin
    sys = LANDISBiomass()
    compiled = mtkcompile(sys)

    tspan_s = 100.0 * yr_to_s

    # No competition case
    prob_no_comp = ODEProblem(
        compiled,
        Dict(
            compiled.B => 5.0 * Mg_ha_to_kg_m2, compiled.D_wood => 0.0,
            compiled.B_other => 0.0
        ),
        (0.0, tspan_s)
    )
    sol_no_comp = solve(prob_no_comp)

    # High competition case: 400 Mg/ha of other species
    # B_other must exceed B_MAX_site - B_MAX = 500 - 223.5 = 276.5 Mg/ha
    # for competition to reduce B_POT below B_MAX
    prob_comp = ODEProblem(
        compiled,
        Dict(
            compiled.B => 5.0 * Mg_ha_to_kg_m2, compiled.D_wood => 0.0,
            compiled.B_other => 400.0 * Mg_ha_to_kg_m2
        ),
        (0.0, tspan_s)
    )
    sol_comp = solve(prob_comp)

    @test sol_no_comp.retcode == SciMLBase.ReturnCode.Success
    @test sol_comp.retcode == SciMLBase.ReturnCode.Success

    # With competition, biomass should be lower
    B_no_comp = sol_no_comp(100.0 * yr_to_s; idxs = compiled.B) / Mg_ha_to_kg_m2
    B_comp = sol_comp(100.0 * yr_to_s; idxs = compiled.B) / Mg_ha_to_kg_m2

    @test B_comp < B_no_comp
end

@testitem "LANDISBiomass: Species Longevity Effect" setup = [LANDISSetup] tags = [:landis] begin
    sys = LANDISBiomass()
    compiled = mtkcompile(sys)

    tspan_s = 150.0 * yr_to_s

    # Short-lived species: P. banksiana (70 years)
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

    # Long-lived species: A. saccharum (400 years, default)
    prob_long = ODEProblem(compiled, [], (0.0, tspan_s))
    sol_long = solve(prob_long)

    @test sol_short.retcode == SciMLBase.ReturnCode.Success
    @test sol_long.retcode == SciMLBase.ReturnCode.Success

    # Short-lived species should have lower biomass at 100 years
    B_short_100 = sol_short(100.0 * yr_to_s; idxs = compiled.B) / Mg_ha_to_kg_m2
    B_long_100 = sol_long(100.0 * yr_to_s; idxs = compiled.B) / Mg_ha_to_kg_m2

    @test B_short_100 < B_long_100
end

@testitem "LANDISBiomass: Fig. 6b - P. tremuloides Cyclical Biomass" setup = [LANDISSetup] tags = [:landis] begin
    sys = LANDISBiomass()
    compiled = mtkcompile(sys)

    tspan_s = 200.0 * yr_to_s

    # P. tremuloides parameters from Table 2
    prob = ODEProblem(
        compiled,
        Dict(
            compiled.B => 5.0 * Mg_ha_to_kg_m2, compiled.D_wood => 0.0,
            compiled.max_age => 120.0 * yr_to_s,
            compiled.ANPP_MAX => 7.46 * Mg_ha_to_kg_m2 / yr_to_s
        ),
        (0.0, tspan_s)
    )
    sol = solve(prob)
    @test sol.retcode == SciMLBase.ReturnCode.Success

    # Paper reports P. tremuloides biomass ranges ~94-124 Mg/ha
    B_vals = [sol(yr * yr_to_s; idxs = compiled.B) / Mg_ha_to_kg_m2 for yr in 10:10:200]
    peak_B = maximum(B_vals)
    @test 80.0 < peak_B < 160.0  # Peak should be in reasonable range

    # Biomass should decline after the species reaches its max_age
    # At 120+ years, age-related mortality should cause decline
    B_100 = sol(100.0 * yr_to_s; idxs = compiled.B) / Mg_ha_to_kg_m2
    B_150 = sol(150.0 * yr_to_s; idxs = compiled.B) / Mg_ha_to_kg_m2
    @test B_150 < B_100  # Biomass declines past max_age
end

@testitem "LANDISBiomass: Fig. 6c - Multi-Species Relative Growth" setup = [LANDISSetup] tags = [:landis] begin
    sys = LANDISBiomass()
    compiled = mtkcompile(sys)

    tspan_s = 200.0 * yr_to_s

    # Mature P. strobus (age_init=200 yr, Longevity=450, ANPP_MAX=10.19)
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

    # Young A. saccharum (age_init=10 yr, default)
    prob_acesac = ODEProblem(compiled, [], (0.0, tspan_s))
    sol_acesac = solve(prob_acesac)

    # Young T. canadensis (age_init=10 yr, Longevity=640, ANPP_MAX=6.27)
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

    @test sol_pinstro.retcode == SciMLBase.ReturnCode.Success
    @test sol_acesac.retcode == SciMLBase.ReturnCode.Success
    @test sol_tsuga.retcode == SciMLBase.ReturnCode.Success

    # A. saccharum has higher ANPP_MAX so it should grow faster initially
    B_acesac_50 = sol_acesac(50.0 * yr_to_s; idxs = compiled.B) / Mg_ha_to_kg_m2
    B_tsuga_50 = sol_tsuga(50.0 * yr_to_s; idxs = compiled.B) / Mg_ha_to_kg_m2
    @test B_acesac_50 > B_tsuga_50  # Higher ANPP gives faster early growth

    # P. strobus starts mature but should eventually decline due to age
    # (age at end = 200+200 = 400 yr, longevity = 450)
    B_pinstro_50 = sol_pinstro(50.0 * yr_to_s; idxs = compiled.B) / Mg_ha_to_kg_m2
    B_pinstro_200 = sol_pinstro(200.0 * yr_to_s; idxs = compiled.B) / Mg_ha_to_kg_m2
    @test B_pinstro_200 < B_pinstro_50  # Mature species declines over time
end
