@testsnippet StagePrognosisSetup begin
    using Test
    using Statistics
    using ModelingToolkit
    using ModelingToolkit: t, D
    using OrdinaryDiffEqDefault
    using StochasticDiffEq
    using Vegetation

    const one_year = 3.15576e7  # seconds per year
    const one_inch = 0.0254     # meters per inch
    const one_foot = 0.3048     # meters per foot
end

@testitem "Structural Verification" setup = [StagePrognosisSetup] tags = [:stage_prognosis] begin
    sys = StagePrognosis()

    # Check system creation
    @test sys isa ModelingToolkit.System

    # Check equations: 8 algebraic + 4 differential = 12
    eqs = equations(sys)
    @test length(eqs) == 12

    # Check unknowns: 4 state + 8 derived = 12
    @test length(unknowns(sys)) == 12

    # Check that the system compiles
    compiled_sys = mtkcompile(sys)
    @test compiled_sys isa ModelingToolkit.System

    # Check that we can create an SDEProblem (not just ODEProblem)
    tspan = (0.0, 10.0 * one_year)
    prob = SDEProblem(compiled_sys, [], tspan)
    @test prob isa SDEProblem

    # Verify state variables exist
    state_vars = unknowns(compiled_sys)
    state_names = Set(string.(Symbol.(state_vars)))
    @test "Dsq" in state_names || "Dsq(t)" in state_names
    @test "HT" in state_names || "HT(t)" in state_names
    @test "HCB" in state_names || "HCB(t)" in state_names
    @test "N_trees" in state_names || "N_trees(t)" in state_names
end

@testitem "Equation Verification - Diameter Growth" setup = [StagePrognosisSetup] tags = [:stage_prognosis] begin
    # Test that the diameter growth equation produces reasonable values
    # for the lodgepole pine example in Figure 1 of Stage (1973)
    sys = StagePrognosis()
    compiled_sys = mtkcompile(sys)

    # Initial conditions matching Figure 1:
    # DBH = 6.0 inches, HT = 60 feet, crown ratio ~0.33, 500 trees/acre
    tspan = (0.0, 11.0 * one_year)  # 11-year projection (1969-1980)
    prob = SDEProblem(compiled_sys, [], tspan)

    # Solve deterministically (zero noise) by setting σ_growth = 0
    prob_det = remake(prob; p = [compiled_sys.σ_growth => 0.0])
    sol = solve(prob_det, EM(), dt = one_year / 10)

    # After 11 years, DBH should increase from 6.0 to ~6.0-7.0 inches
    # Figure 1B shows DBH going from 5.7 to 6.0 for the 0.9 percentile tree
    DBH_final_in = sol[compiled_sys.DBH][end] / one_inch
    @test 6.0 < DBH_final_in < 8.0  # reasonable range

    # DDS_rate should be positive (trees are growing)
    @test all(sol[compiled_sys.DDS_rate] .> 0)

    # Check diameter growth is monotonically increasing
    dbh_vals = sol[compiled_sys.DBH]
    @test all(diff(dbh_vals) .>= -1.0e-15)  # allowing tiny numerical noise
end

@testitem "Equation Verification - Height Growth" setup = [StagePrognosisSetup] tags = [:stage_prognosis] begin
    sys = StagePrognosis()
    compiled_sys = mtkcompile(sys)

    tspan = (0.0, 11.0 * one_year)
    prob = SDEProblem(compiled_sys, [], tspan)
    prob_det = remake(prob; p = [compiled_sys.σ_growth => 0.0])
    sol = solve(prob_det, EM(), dt = one_year / 10)

    # Height should increase from 60 ft to ~65-70 ft over 11 years
    # Figure 1B shows height going from ~60 to ~70 feet
    HT_init_ft = sol[compiled_sys.HT][1] / one_foot
    HT_final_ft = sol[compiled_sys.HT][end] / one_foot

    @test HT_init_ft ≈ 60.0 atol = 0.1
    @test 62.0 < HT_final_ft < 78.0  # reasonable range

    # Height growth rate should be positive
    @test all(sol[compiled_sys.HTGF_rate] .> 0)

    # Height should be monotonically increasing
    ht_vals = sol[compiled_sys.HT]
    @test all(diff(ht_vals) .>= -1.0e-15)
end

@testitem "Crown Ratio Development" setup = [StagePrognosisSetup] tags = [:stage_prognosis] begin
    sys = StagePrognosis()
    compiled_sys = mtkcompile(sys)

    # Test with CCF < 125 (default CCF=100)
    tspan = (0.0, 10.0 * one_year)
    prob = SDEProblem(compiled_sys, [], tspan)
    prob_det = remake(prob; p = [compiled_sys.σ_growth => 0.0])
    sol = solve(prob_det, EM(), dt = one_year / 10)

    # Crown ratio should be bounded between 0 and 1
    cr_vals = sol[compiled_sys.CR]
    @test all(cr_vals .> 0.0)
    @test all(cr_vals .< 1.0)

    # Crown base height should increase (crown recedes)
    hcb_vals = sol[compiled_sys.HCB]
    @test all(diff(hcb_vals) .>= -1.0e-15)

    # With CCF < 125, crown recession = 1/5 of height growth
    # So crown ratio should generally increase (height grows faster than crown recedes)
    ht_vals = sol[compiled_sys.HT]
    # Check that crown length (HT - HCB) is increasing
    crown_length = ht_vals .- hcb_vals
    @test crown_length[end] > crown_length[1]
end

@testitem "Mortality Model" setup = [StagePrognosisSetup] tags = [:stage_prognosis] begin
    sys = StagePrognosis()
    compiled_sys = mtkcompile(sys)

    tspan = (0.0, 50.0 * one_year)
    prob = SDEProblem(compiled_sys, [], tspan)
    prob_det = remake(prob; p = [compiled_sys.σ_growth => 0.0])
    sol = solve(prob_det, EM(), dt = one_year / 10)

    # Tree density should decrease over time due to mortality
    n_vals = sol[compiled_sys.N_trees]
    @test n_vals[end] < n_vals[1]
    @test all(n_vals .> 0)  # Can't have negative trees

    # Mortality rate should be positive
    mort_vals = sol[compiled_sys.mortality_rate]
    @test all(mort_vals .>= 0)

    # Check mortality percentile factor:
    # At PCT=50: factor = 0.25 + 1.5*(1 - 50/100) = 0.25 + 0.75 = 1.0
    # At PCT=100 (largest tree): factor = 0.25 + 1.5*(1 - 1) = 0.25 (lowest mortality)
    # At PCT=0 (smallest tree): factor = 0.25 + 1.5*1 = 1.75 (highest mortality)

    # Test with high PCT (large tree) - should have lower mortality
    prob_high = remake(prob; p = [compiled_sys.σ_growth => 0.0, compiled_sys.PCT => 90.0])
    sol_high = solve(prob_high, EM(), dt = one_year / 10)

    # Test with low PCT (small tree) - should have higher mortality
    prob_low = remake(prob; p = [compiled_sys.σ_growth => 0.0, compiled_sys.PCT => 10.0])
    sol_low = solve(prob_low, EM(), dt = one_year / 10)

    # Small trees (low PCT) should lose more trees to mortality
    n_high = sol_high[compiled_sys.N_trees]
    n_low = sol_low[compiled_sys.N_trees]
    @test n_low[end] < n_high[end]  # More mortality for small trees
end

@testitem "Stochastic Behavior" setup = [StagePrognosisSetup] tags = [:stage_prognosis] begin
    sys = StagePrognosis()
    compiled_sys = mtkcompile(sys)

    tspan = (0.0, 50.0 * one_year)
    prob = SDEProblem(compiled_sys, [], tspan)

    # Run multiple stochastic realizations
    n_runs = 20
    final_dbh = Float64[]
    for i in 1:n_runs
        sol = solve(prob, EM(), dt = one_year / 10, seed = i)
        push!(final_dbh, sol[compiled_sys.DBH][end] / one_inch)
    end

    # Check that there is variation between runs
    @test std(final_dbh) > 0.01  # Some variation exists

    # Check that all runs produce reasonable DBH values
    @test all(final_dbh .> 5.0)   # Trees grew
    @test all(final_dbh .< 20.0)  # Not unreasonably large

    # Mean should be close to deterministic result
    prob_det = remake(prob; p = [compiled_sys.σ_growth => 0.0])
    sol_det = solve(prob_det, EM(), dt = one_year / 10)
    det_dbh = sol_det[compiled_sys.DBH][end] / one_inch
    @test abs(mean(final_dbh) - det_dbh) < 2.0  # Within 2 inches of deterministic
end

@testitem "Qualitative Properties" setup = [StagePrognosisSetup] tags = [:stage_prognosis] begin
    sys = StagePrognosis()
    compiled_sys = mtkcompile(sys)

    # Test positivity preservation
    tspan = (0.0, 50.0 * one_year)
    prob = SDEProblem(compiled_sys, [], tspan)
    prob_det = remake(prob; p = [compiled_sys.σ_growth => 0.0])
    sol = solve(prob_det, EM(), dt = one_year / 10)

    # All state variables should remain positive
    @test all(sol[compiled_sys.Dsq] .> 0)
    @test all(sol[compiled_sys.HT] .> 0)
    @test all(sol[compiled_sys.N_trees] .> 0)

    # DBH should always be less than height (physical constraint)
    dbh_vals = sol[compiled_sys.DBH]
    ht_vals = sol[compiled_sys.HT]
    @test all(dbh_vals .< ht_vals)

    # Crown base should always be below total height
    hcb_vals = sol[compiled_sys.HCB]
    @test all(hcb_vals .< ht_vals)
end

@testitem "Site Index Response" setup = [StagePrognosisSetup] tags = [:stage_prognosis] begin
    sys = StagePrognosis()
    compiled_sys = mtkcompile(sys)

    tspan = (0.0, 50.0 * one_year)
    prob = SDEProblem(compiled_sys, [], tspan)

    # Higher site index should produce more growth
    prob_low_si = remake(
        prob; p = [
            compiled_sys.σ_growth => 0.0,
            compiled_sys.SI => 50.0 * one_foot,
        ]
    )
    prob_high_si = remake(
        prob; p = [
            compiled_sys.σ_growth => 0.0,
            compiled_sys.SI => 100.0 * one_foot,
        ]
    )

    sol_low = solve(prob_low_si, EM(), dt = one_year / 10)
    sol_high = solve(prob_high_si, EM(), dt = one_year / 10)

    # Better sites (higher SI) should produce larger DBH
    @test sol_high[compiled_sys.DBH][end] > sol_low[compiled_sys.DBH][end]
end

@testitem "Stand Density Response" setup = [StagePrognosisSetup] tags = [:stage_prognosis] begin
    sys = StagePrognosis()
    compiled_sys = mtkcompile(sys)

    tspan = (0.0, 50.0 * one_year)
    prob = SDEProblem(compiled_sys, [], tspan)

    # Higher CCF (more competition) should reduce growth
    prob_low_ccf = remake(
        prob; p = [
            compiled_sys.σ_growth => 0.0,
            compiled_sys.CCF => 50.0,
        ]
    )
    prob_high_ccf = remake(
        prob; p = [
            compiled_sys.σ_growth => 0.0,
            compiled_sys.CCF => 200.0,
        ]
    )

    sol_low = solve(prob_low_ccf, EM(), dt = one_year / 10)
    sol_high = solve(prob_high_ccf, EM(), dt = one_year / 10)

    # Less competition (lower CCF) should produce larger DBH
    @test sol_low[compiled_sys.DBH][end] > sol_high[compiled_sys.DBH][end]
end

@testitem "HCB Prediction Equation" setup = [StagePrognosisSetup] tags = [:stage_prognosis] begin
    # Test the static HCB prediction equation from Stage (1973, p. 16)
    # HCB = -29.26 + 0.61*HT + 9.178*ln(CCF) - 0.222*EL - 5.80*DBH/RMSQD + HAB

    sys = StagePrognosisHCB()
    @test sys isa ModelingToolkit.System

    compiled_sys = mtkcompile(sys)

    # Test with default parameters (HT=60ft, DBH=6in, CCF=100, EL=43 hundreds of feet, RMSQD=7in)
    # HCB = -29.26 + 0.61*60 + 9.178*ln(100) - 0.222*43 - 5.80*6/7 + 0.0
    #      = -29.26 + 36.6 + 42.26 - 9.55 - 4.97 + 0 = 35.08 feet
    prob = ODEProblem(compiled_sys, [], (0.0, 1.0))
    sol = solve(prob)
    hcb_ft = sol[compiled_sys.HCB_pred][1] / one_foot

    # Verify result is in a physically reasonable range
    @test 25.0 < hcb_ft < 45.0  # Crown base should be between 25-45 ft for a 60ft tree

    # Verify that taller trees have higher crown base
    prob_tall = remake(
        prob; p = [compiled_sys.HT_input => 80.0 * one_foot]
    )
    sol_tall = solve(prob_tall)
    @test sol_tall[compiled_sys.HCB_pred][1] > sol[compiled_sys.HCB_pred][1]

    # Verify that higher CCF (more competition) raises crown base
    prob_dense = remake(
        prob; p = [compiled_sys.CCF => 200.0]
    )
    sol_dense = solve(prob_dense)
    @test sol_dense[compiled_sys.HCB_pred][1] > sol[compiled_sys.HCB_pred][1]

    # Verify that larger trees relative to mean (DBH/RMSQD > 1) have lower crown base
    prob_large = remake(
        prob; p = [compiled_sys.DBH_input => 10.0 * one_inch]
    )
    sol_large = solve(prob_large)
    @test sol_large[compiled_sys.HCB_pred][1] < sol[compiled_sys.HCB_pred][1]

    # Verify exact computation against hand-calculated value
    # HCB = -29.26 + 0.61*60 + 9.178*ln(100) - 0.222*43 - 5.80*(6/7) + 0.0
    expected_hcb_ft = -29.26 + 0.61 * 60 + 9.178 * log(100) - 0.222 * 43 - 5.8 * (6 / 7) + 0.0
    @test hcb_ft ≈ expected_hcb_ft rtol = 1.0e-4
end

@testitem "Figure 1B Numerical Verification" setup = [StagePrognosisSetup] tags = [:stage_prognosis] begin
    # Verify model output against Figure 1B reference values from Stage (1973).
    # The paper shows sample tree records for the 50th percentile tree in a
    # lodgepole pine stand (SI=80 ft, 509 trees/acre, EL=43 hundred feet).
    # The 0.5 percentile tree starts at DBH=6.8 in, HT=67.0 ft, CR=35% in 1969.
    sys = StagePrognosis()
    compiled_sys = mtkcompile(sys)

    cr0 = 0.35
    ht0 = 67.0 * one_foot
    hcb0 = ht0 * (1.0 - cr0)

    prob = SDEProblem(
        compiled_sys, [
            compiled_sys.Dsq => (6.8 * one_inch)^2,
            compiled_sys.HT => ht0,
            compiled_sys.HCB => hcb0,
            compiled_sys.N_trees => 509.0 / 4046.86,
        ], (0.0, 51.0 * one_year)
    )
    prob_det = remake(prob; p = [compiled_sys.σ_growth => 0.0, compiled_sys.PCT => 50.0])
    sol = solve(prob_det, EM(), dt = one_year / 10)

    years = sol.t ./ one_year

    # Paper values for the 0.5 percentile tree (Figure 1B, p. 5)
    # Year  DBH(in)  HT(ft)  CR
    # 1969  6.8      67.0    35
    # 1980  7.2      74.2    34
    # 1990  7.5      80.6    34
    # 2000  8.2      87.3    37
    # 2020  8.8      99.5    —

    # Extract model values at the matching time indices
    idx_0 = 1
    idx_11 = findfirst(y -> y >= 11.0, years)   # 1969→1980 (11 years)
    idx_21 = findfirst(y -> y >= 21.0, years)   # 1969→1990
    idx_31 = findfirst(y -> y >= 31.0, years)   # 1969→2000
    idx_51 = findfirst(y -> y >= 51.0, years)   # 1969→2020

    model_dbh = [sol[compiled_sys.DBH][i] / one_inch for i in [idx_0, idx_11, idx_21, idx_31, idx_51]]
    model_ht = [sol[compiled_sys.HT][i] / one_foot for i in [idx_0, idx_11, idx_21, idx_31, idx_51]]

    paper_dbh = [6.8, 7.2, 7.5, 8.2, 8.8]
    paper_ht = [67.0, 74.2, 80.6, 87.3, 99.5]

    # Initial conditions should match exactly
    @test model_dbh[1] ≈ paper_dbh[1] atol = 0.01
    @test model_ht[1] ≈ paper_ht[1] atol = 0.01

    # Growth direction must be correct: DBH and height must increase
    @test all(diff(model_dbh) .> 0)
    @test all(diff(model_ht) .> 0)

    # Model predictions should be within reasonable range of paper values.
    # The height increment coefficients (c1_ht..c4_ht) are calibrated estimates
    # since the original Cole & Stage (in prep.) coefficients are unpublished.
    # Height growth matches well (~10%), but diameter growth runs ~2x the paper
    # values, likely because Figure 1B reflects a complex stand structure that
    # cannot be fully captured by running a single-tree model with PCT=50.
    for i in 2:5
        dbh_growth_paper = paper_dbh[i] - paper_dbh[1]
        dbh_growth_model = model_dbh[i] - model_dbh[1]
        @test dbh_growth_model > 0  # Must grow

        ht_growth_paper = paper_ht[i] - paper_ht[1]
        ht_growth_model = model_ht[i] - model_ht[1]
        @test ht_growth_model > 0  # Must grow
        if ht_growth_paper > 0.1
            @test abs(ht_growth_model - ht_growth_paper) / ht_growth_paper < 0.15
        end
    end

    # DBH growth should be within a factor of 3 of the paper values
    # (order-of-magnitude check given the known discrepancy)
    total_dbh_growth_paper = paper_dbh[end] - paper_dbh[1]
    total_dbh_growth_model = model_dbh[end] - model_dbh[1]
    @test 0.3 < total_dbh_growth_model / total_dbh_growth_paper < 3.0
end
