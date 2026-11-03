* =========================================================
* 01_baseline_ppml.do  Baseline PPML: adaptation finance on lagged MA3 centrality
* Origin : Code/Reg_model.do, section 1.0 (lines 93-129)
* Output : output/tables/
*            regression_results_5_centrality_lagged_MA3_PPML.csv
* =========================================================

do "$ANALYSIS/00_prepare.do"

* ---------------------------------------------------------
* 1.0 Baseline Regression: PPML
* ---------------------------------------------------------

local centralities in_degree out_degree closeness eigenvector
local controls gdppc fdi gov_effectiveness disaster_total_impact_pc agri_va rural_pop co2 total_erp_pct trade

eststo clear

local i = 1
foreach var in `centralities' {
    eststo m`i': ppmlhdfe total_adaptation_finance `var'_ma3 `controls', absorb(region_id year) vce(cluster iso3_id) d

	local pseudo_r2 = e(r2_p)
    local nobs = e(N)

	margins, dydx(`var'_ma3) post

	eststo m`i'
	estadd scalar r2_p = `pseudo_r2'
    estadd scalar N_obs = `nobs'
    estadd local Region_FE "Yes"
    estadd local Year_FE "Yes"
    estadd local Controls "Yes"

    local i = `i' + 1
}

esttab m1 m2 m3 m4 ///
    using "$TABLES/regression_results_5_centrality_lagged_MA3_PPML.csv", ///
    replace se star(* 0.1 ** 0.05 *** 0.01) ///
    b(3) se(3) ///
    compress nogaps ///
    nodepvars nomtitles ///
    keep(*_ma3) ///
    scalars("Controls Controls" "Region_FE Region FE" "Year_FE Year FE" "r2_p Pseudo R-squared" "N Observations") ///
    label
