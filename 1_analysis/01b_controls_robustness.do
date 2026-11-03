* =========================================================
* 01b_controls_robustness.do  Baseline PPML with stepwise control sets
* Response to reviewer: are results driven by mitigation-related controls?
*   (1) FE only
*   (2) + general economic / institutional controls
*   (3) + adaptation-need controls incl. ND-GAIN   <- proposed main spec
*   (4) + mitigation controls (recap, co2)
* All columns are estimated on a common sample (non-missing for (4)).
* Output : output/tables/
*            regression_results_5b_controls_stepwise_PPML.csv
* =========================================================

do "$ANALYSIS/00_prepare.do"

* gain_index is not standardised in 00_prepare.do
egen temp_gain = std(gain_index)
replace gain_index = temp_gain
drop temp_gain
label var gain_index "GAIN Index (Standardized)"

local centralities in_degree out_degree closeness eigenvector
local ctrl_general    gdppc fdi trade gov_effectiveness total_erp_pct
local ctrl_adaptation disaster_total_impact_pc agri_va rural_pop gain_index
local ctrl_mitigation recap co2

local set1
local set2 `ctrl_general'
local set3 `ctrl_general' `ctrl_adaptation'
local set4 `ctrl_general' `ctrl_adaptation' `ctrl_mitigation'

local outfile "$TABLES/regression_results_5b_controls_stepwise_PPML.csv"
capture erase "`outfile'"

foreach var in `centralities' {

    eststo clear

    * Common sample: non-missing for the fullest specification
    capture drop _insample
    quietly ppmlhdfe total_adaptation_finance `var'_ma3 `set4', ///
        absorb(region_id year) vce(cluster iso3_id) d
    gen byte _insample = e(sample)

    forvalues s = 1/4 {
        eststo m`s': ppmlhdfe total_adaptation_finance `var'_ma3 `set`s'' ///
            if _insample, absorb(region_id year) vce(cluster iso3_id) d

        local pseudo_r2 = e(r2_p)
        local nobs = e(N)

        margins, dydx(`var'_ma3) post

        eststo m`s'
        estadd scalar r2_p = `pseudo_r2'
        estadd scalar N_obs = `nobs'
        estadd local Region_FE "Yes"
        estadd local Year_FE "Yes"
        estadd local General   = cond(`s' >= 2, "Yes", "No")
        estadd local Adaptation = cond(`s' >= 3, "Yes", "No")
        estadd local Mitigation = cond(`s' >= 4, "Yes", "No")
    }

    local mode = cond("`var'" == "in_degree", "replace", "append")

    esttab m1 m2 m3 m4, ///
        se star(* 0.1 ** 0.05 *** 0.01) b(3) se(3) ///
        compress nogaps nodepvars nomtitles keep(*_ma3) ///
        scalars("General General controls" "Adaptation Adaptation controls" ///
                "Mitigation Mitigation controls" "r2_p Pseudo R-squared" "N Observations") ///
        title("`var'_ma3")

    esttab m1 m2 m3 m4 using "`outfile'", `mode' ///
        se star(* 0.1 ** 0.05 *** 0.01) b(3) se(3) ///
        compress nogaps nodepvars nomtitles keep(*_ma3) ///
        scalars("General General controls" "Adaptation Adaptation controls" ///
                "Mitigation Mitigation controls" "Region_FE Region FE" "Year_FE Year FE" ///
                "r2_p Pseudo R-squared" "N Observations") ///
        title("`var'_ma3") label

    drop _insample
}
