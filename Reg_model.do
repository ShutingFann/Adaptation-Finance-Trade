clear all
set more off

cd "/Users/fanshuting/Desktop/UCL-CSC/1_SupplyChain"

import delimited "regression_panel_treated_non_high_income.csv", clear

keep if inrange(year, 2000, 2024)

label variable in_degree 			"In-Degree Centrality"
label variable gdppc 				"GDP per capita"
label variable fdi 					"Foreign Direct Investment"
label variable gain_index 			"GAIN Index"
label variable gov_effectiveness 	"Government Effectiveness"
label variable disaster_affected 	"Disaster Affected"
label variable total_erp_pct 		"Total Equity Risk Premium (%)"
label variable out_degree    		"Out-Degree Centrality"
label variable closeness     		"Closeness Centrality"
label variable betweenness   		"Betweenness Centrality"
label variable katz   				"katz Centrality"
label variable eigenvector   		"Eigenvector Centrality"
label variable recap   				"Renewable Capacity Installed"
label variable agri_va        		"Agricultural Value Added (%)"
label variable rural_pop     		"Rural Population (%)"
label variable co2           		"CO2 Emissions per capita"
label variable trade           		"Trade (% of GDP)"

local control_vars gdppc fdi gov_effectiveness recap agri_va rural_pop co2 total_erp_pct trade ///
		disaster_total_impact_pc idealpoint_align iv_seaborne_trade_total iv_seaborne_trade_total_ch iv_seaborne_trade_total_gr trade_loaded trade_discharged trade_total ln_iv_seaborne_trade_total

foreach var of local control_vars {
    egen temp_`var' = std(`var')
    replace `var' = temp_`var'
    drop temp_`var'
    local var_label : variable label `var'
    label var `var' "`var_label' (Standardized)"
}

ds, has(type numeric)
local numvars `r(varlist)'

foreach var of local numvars {
    
    local lower_var = strlower("`var'")
    
    if inlist("`lower_var'", "in_degree", "out_degree", "in_degree_iv", "out_degree_iv") {
        replace `var' = ln(`var')
       
        egen tem_`var' = std(`var')
        replace `var' = tem_`var'
        drop tem_`var'
    }
    
	else if inlist("`lower_var'", "betweenness", "eigenvector", "betweenness_iv", "eigenvector_iv") {
        replace `var' = ln(`var' * 1000000 + 1)
        
        egen tem_`var' = std(`var')
        replace `var' = tem_`var'
        drop tem_`var'
    }
	
    else if inlist("`lower_var'", "closeness", "closeness_iv") {
        egen tem_`var' = std(`var')
        replace `var' = tem_`var'
        drop tem_`var'
    }
}

encode iso3, gen(iso3_id)
encode region, gen(region_id)
xtset iso3_id year

tab year, gen(yd_)
drop yd_1


local centralities in_degree out_degree closeness betweenness eigenvector

foreach var of local centralities {
    rangestat (mean) `var', interval(year -3 -1) by(iso3_id)
    rename `var'_mean `var'_ma3
    label var `var'_ma3 "`: var label `var''"
}

local centralities in_degree out_degree closeness betweenness eigenvector

foreach var of local centralities {
    rangestat (mean) `var', interval(year -5 -1) by(iso3_id)
    rename `var'_mean `var'_ma5
    label var `var'_ma5 "`: var label `var''"
}

* ---------------------------------------------------------
* 1.0 Baseline Regression: PPML
* ---------------------------------------------------------

local centralities in_degree out_degree closeness eigenvector
local controls gdppc fdi gov_effectiveness disaster_total_impact_pc recap agri_va rural_pop co2 total_erp_pct trade

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
    using "regression_results_5_centrality_lagged_MA3_PPML.csv", ///
    replace se star(* 0.1 ** 0.05 *** 0.01) ///
    b(3) se(3) ///
    compress nogaps ///
    nodepvars nomtitles ///
    keep(*_ma3) ///
    scalars("Controls Controls" "Region_FE Region FE" "Year_FE Year FE" "r2_p Pseudo R-squared" "N Observations") ///
    label	

* =========================================================
* 1.1 IV Regression: Gravity IV -- Control Function Approach (CF) -- Adaptation Financing
* =========================================================

global controls gdppc fdi gov_effectiveness disaster_total_impact_pc recap agri_va rural_pop co2 total_erp_pct trade
capture xtset, clear
capture tsset, clear
eststo clear

capture program drop my_cf_boot
program my_cf_boot, eclass
    version 14.0 
    capture drop v_resid_boot
    
    quietly reghdfe $current_endog $current_iv $controls, absorb(region_id year) resid
    quietly predict v_resid_boot, resid
    quietly ppmlhdfe total_adaptation_finance $current_endog $controls v_resid_boot, absorb(region_id year) d
	quietly margins, dydx($current_endog v_resid_boot) post
end

local centralities in_degree out_degree closeness eigenvector

foreach var in `centralities' {
    
    global current_endog "`var'_ma3"
    global current_iv "`var'_iv"
    
    quietly reghdfe $current_endog $current_iv $controls, absorb(region_id year) vce(cluster iso3_id)
    test $current_iv
    local f_stat = r(F)
    
    capture drop new_iso3 
    set seed 123456
    quietly bootstrap, reps(1000) cluster(iso3_id) idcluster(new_iso3): my_cf_boot
    
    estadd scalar F_first = `f_stat'
    estadd local Controls "Yes"
    estadd local FE "Yes"
    eststo iv_`var'
}

esttab iv_in_degree iv_out_degree iv_closeness iv_eigenvector ///
    using "regression_results_5_centrality_GravityIV_CF_Bootstrap_Adaptation.csv", replace ///
    keep(*_ma3 v_resid_boot) /// 
	order(in_degree_ma3 out_degree_ma3 closeness_ma3 eigenvector_ma3 v_resid_boot) ///
    b(3) se(3) star(* 0.1 ** 0.05 *** 0.01) ///
    mtitles("In-Degree" "Out-Degree" "Closeness" "Eigenvector") ///
    scalars( ///
        "F_first First-Stage F-stat" ///
        "Controls Controls" ///
        "FE Region & Year FE" ///
        "N Observations" ///
    ) ///
    compress nogaps label	
	
	
* =========================================================
* 1.2 IV Regression: Seaborne Trade Shift-Share IV -- Control Function Approach (CF) -- Adatation Financing
* =========================================================

global controls gdppc fdi gov_effectiveness disaster_total_impact_pc recap agri_va rural_pop co2 total_erp_pct trade
capture xtset, clear
capture tsset, clear
eststo clear

capture program drop my_cf_boot
program my_cf_boot, eclass
    version 14.0 
    capture drop v_resid_boot
    
    quietly reghdfe $current_endog $current_iv $controls, absorb(region_id year) resid
    quietly predict v_resid_boot, resid
    quietly ppmlhdfe total_adaptation_finance $current_endog $controls v_resid_boot, absorb(region_id year) d
	quietly margins, dydx($current_endog v_resid_boot) post

end

local centralities in_degree out_degree closeness eigenvector

foreach var in `centralities' {
    
    global current_endog "`var'_ma3"
    global current_iv ln_iv_seaborne_trade_total
    
    reghdfe $current_endog $current_iv $controls, absorb(region_id year) vce(cluster iso3_id)
    quietly test $current_iv
    local f_stat = r(F)
    
    capture drop new_iso3 
    set seed 123
    bootstrap, reps(1000) cluster(iso3_id) idcluster(new_iso3): my_cf_boot
    
    estadd scalar F_first = `f_stat'
    estadd local Controls "Yes"
    estadd local FE "Yes"
    eststo iv_`var'
}

esttab iv_in_degree iv_out_degree iv_closeness iv_eigenvector ///
    using "regression_results_5_centrality_BartikIV_CF_Bootstrap_Adaptation.csv", replace ///
    keep(*_ma3 v_resid_boot) ///
	order(in_degree_ma3 out_degree_ma3 closeness_ma3 eigenvector_ma3 v_resid_boot) ///
    b(3) se(3) star(* 0.1 ** 0.05 *** 0.01) ///
    mtitles("In-Degree" "Out-Degree" "Closeness" "Eigenvector") ///
    scalars( ///
        "F_first First-Stage F-stat" ///
        "Controls Controls" ///
        "FE Region & Year FE" ///
        "N Observations" ///
    ) ///
    compress nogaps label	

* =========================================================
* 2. Continuous CF-PPML:
* 4x4 Quartile-Entry-Threshold Counterfactual Matrix
*
* First stage:  continuous centrality ~ IV + controls + FE
* Second stage: PPML(finance ~ continuous centrality
*                    + CF residual + controls + FE)
*
* Counterfactual:
* - Baseline: each country-year retains its actual centrality
* - Target: move centrality to the nearest quartile-entry
*           threshold in the direction of movement
*
* Examples:
*   Q1 -> Q4: set centrality to annual P75
*   Q4 -> Q1: set centrality to annual P25
*   Q3 -> Q2: set centrality to annual P50
*	Near-threshold random draws
* =========================================================

clear matrix
set more off

* ---------------------------------------------------------
* 0. Global settings
* ---------------------------------------------------------
global controls gdppc fdi gov_effectiveness disaster_total_impact_pc ///
    recap agri_va rural_pop co2 total_erp_pct trade

capture xtset, clear
capture tsset, clear
eststo clear

local centralities closeness eigenvector in_degree out_degree
// local centralities closeness eigenvector
local reps 1000
global mc_reps 50

set seed 12345


* ---------------------------------------------------------
* 1. Bootstrap subprogram
*
* In each bootstrap replication:
*   1. Re-estimate continuous first stage;
*   2. Construct CF residual;
*   3. Estimate continuous CF-PPML;
*   4. Predict baseline finance at actual centrality;
*   5. Move centrality to the relevant quartile boundary;
*   6. Return predicted levels and changes.
* ---------------------------------------------------------
capture program drop my_cont_entry_cf_boot

program define my_cont_entry_cf_boot, rclass
    version 16.0

    tempvar fsample ssample vhat ppml_d
    tempvar c_actual mu_actual

    * -----------------------------------------------------
    * First stage: continuous centrality
    * -----------------------------------------------------
    capture quietly reghdfe ///
        $current_endog ///
        $current_iv ///
        $controls ///
        if $current_sample == 1, ///
        absorb(region_id year) resid

    if _rc {
        forvalues k = 1/4 {
            return scalar m`k' = .
            return scalar d`k' = .
        }
        exit
    }

    gen byte `fsample' = e(sample)

    capture quietly predict double `vhat' if `fsample' == 1, resid

    if _rc {
        forvalues k = 1/4 {
            return scalar m`k' = .
            return scalar d`k' = .
        }
        exit
    }


    * -----------------------------------------------------
    * Second stage: continuous CF-PPML
    * -----------------------------------------------------
    capture quietly ppmlhdfe ///
        total_adaptation_finance ///
        c.$current_endog ///
        $controls ///
        `vhat' ///
        if `fsample' == 1, ///
        absorb(region_id year) ///
        d(`ppml_d')

    if _rc {
        forvalues k = 1/4 {
            return scalar m`k' = .
            return scalar d`k' = .
        }
        exit
    }

    gen byte `ssample' = e(sample)

    * Check that this bootstrap draw contains observations
    * from the relevant original source quartile.
    count if `ssample' == 1 & is_base == 1

    if r(N) == 0 {
        forvalues k = 1/4 {
            return scalar m`k' = .
            return scalar d`k' = .
        }
        exit
    }


    * -----------------------------------------------------
    * Baseline prediction:
    * retain every country-year's observed centrality
    * -----------------------------------------------------
    gen double `c_actual' = $current_endog

    capture quietly predict double `mu_actual' ///
        if `ssample' == 1, mu

    if _rc {
        forvalues k = 1/4 {
            return scalar m`k' = .
            return scalar d`k' = .
        }
        exit
    }

    quietly summarize `mu_actual' ///
        if `ssample' == 1 & is_base == 1, meanonly

    if r(N) == 0 {
        forvalues k = 1/4 {
            return scalar m`k' = .
            return scalar d`k' = .
        }
        exit
    }

    local m_actual = r(mean)


	* -----------------------------------------------------
	* Counterfactual predictions:
	* Near-threshold random draws
	*
	* Instead of assigning observations to a single
	* quartile-entry threshold, draw random values from
	* a near-threshold band of the destination quartile.
	*
	* Upward movements:
	*   target Q2: P25   to P31.25
	*   target Q3: P50   to P56.25
	*   target Q4: P75   to P81.25
	*
	* Downward movements:
	*   target Q1: P18.75 to P25
	*   target Q2: P43.75 to P50
	*   target Q3: P68.75 to P75
	* -----------------------------------------------------
	forvalues t = 1/4 {

		* Restore actual continuous centrality first
		replace $current_endog = `c_actual'

		* Diagonal cell:
		* Maintain actual centrality.
		if `t' == $base_q {

			local m`t' = `m_actual'
			local d`t' = 0
		}

		* Off-diagonal cell:
		* Randomly draw from a near-threshold band
		* of the destination quartile.
		else {

			local lower
			local upper

			* -------------------------------------------------
			* Upward movement:
			* Draw from the lower segment of the destination
			* quartile.
			* -------------------------------------------------
			if `t' > $base_q {

				if `t' == 2 {
					local lower "$current_p25"
					local upper "$current_p3125"
				}

				if `t' == 3 {
					local lower "$current_p50"
					local upper "$current_p5625"
				}

				if `t' == 4 {
					local lower "$current_p75"
					local upper "$current_p8125"
				}
			}

			* -------------------------------------------------
			* Downward movement:
			* Draw from the upper segment of the destination
			* quartile.
			* -------------------------------------------------
			if `t' < $base_q {

				if `t' == 1 {
					local lower "$current_p1875"
					local upper "$current_p25"
				}

				if `t' == 2 {
					local lower "$current_p4375"
					local upper "$current_p50"
				}

				if `t' == 3 {
					local lower "$current_p6875"
					local upper "$current_p75"
				}
			}

			* Safety check
			if "`lower'" == "" | "`upper'" == "" {
				replace $current_endog = `c_actual'

				forvalues k = 1/4 {
					return scalar m`k' = .
					return scalar d`k' = .
				}
				exit
			}

			* -------------------------------------------------
			* Monte Carlo draws within each bootstrap replication
			* -------------------------------------------------
			local m_sum = 0
			local d_sum = 0
			local ok_draws = 0

			forvalues r = 1/$mc_reps {

				* Restore actual centrality before each draw
				replace $current_endog = `c_actual'

				* Assign random values from the near-threshold band.
				* Only source-group observations need to be changed,
				* because predictions are averaged over is_base == 1.
				replace $current_endog = ///
					`lower' + runiform() * (`upper' - `lower') ///
					if `ssample' == 1 & is_base == 1

				tempvar mu_cf

				capture quietly predict double `mu_cf' ///
					if `ssample' == 1, mu

				if !_rc {

					quietly summarize `mu_cf' ///
						if `ssample' == 1 & is_base == 1, meanonly

					if r(N) > 0 & r(mean) < . {

						local m_r = r(mean)
						local d_r = `m_r' - `m_actual'

						local m_sum = `m_sum' + `m_r'
						local d_sum = `d_sum' + `d_r'
						local ok_draws = `ok_draws' + 1
					}
				}

				capture drop `mu_cf'
			}

			* If no valid MC draw, return missing
			if `ok_draws' == 0 {

				replace $current_endog = `c_actual'

				forvalues k = 1/4 {
					return scalar m`k' = .
					return scalar d`k' = .
				}
				exit
			}

			* Average across Monte Carlo draws
			local m`t' = `m_sum' / `ok_draws'
			local d`t' = `d_sum' / `ok_draws'
		}
	}

    * Restore observed continuous centrality before exit.
    replace $current_endog = `c_actual'


    * -----------------------------------------------------
    * Return results to bootstrap
    *
    * m1-m4: predicted financing under target Q1-Q4
    * d1-d4: target prediction minus actual prediction
    * -----------------------------------------------------
    forvalues t = 1/4 {
        return scalar m`t' = `m`t''
        return scalar d`t' = `d`t''
    }

end


* ---------------------------------------------------------
* 2. Main loop: construct one 4x4 matrix per centrality
* ---------------------------------------------------------
foreach var in `centralities' {

    display "==================================================="
    display "Continuous CF-PPML entry-threshold matrix: `var'"
    display "==================================================="


    * -----------------------------------------------------
    * Clean variables created in a previous run
    * -----------------------------------------------------
    capture drop sample_`var'
    capture drop q_`var'

    capture drop p25_`var'
    capture drop p50_`var'
    capture drop p75_`var'

    capture drop is_base


    * -----------------------------------------------------
    * Construct common analytical sample
    * -----------------------------------------------------
    gen byte sample_`var' = !missing( ///
        total_adaptation_finance, ///
        `var'_ma3, ///
        `var'_iv, ///
        gdppc, ///
        fdi, ///
        gov_effectiveness, ///
        disaster_total_impact_pc, ///
        recap, ///
        agri_va, ///
        rural_pop, ///
        co2, ///
        total_erp_pct, ///
        trade, ///
        region_id, ///
        year, ///
        iso3_id ///
    )


    * -----------------------------------------------------
    * Identify the final PPML estimation sample.
    *
    * Quartiles and threshold anchors are defined on this
    * common final sample rather than a larger raw sample.
    * -----------------------------------------------------
    tempvar fs0 v0 ss0 d0

    quietly reghdfe ///
        `var'_ma3 ///
        `var'_iv ///
        $controls ///
        if sample_`var' == 1, ///
        absorb(region_id year) resid

    gen byte `fs0' = e(sample)

    quietly predict double `v0' if `fs0' == 1, resid

    quietly ppmlhdfe ///
        total_adaptation_finance ///
        c.`var'_ma3 ///
        $controls ///
        `v0' ///
        if `fs0' == 1, ///
        absorb(region_id year) ///
        d(`d0')

    gen byte `ss0' = e(sample)

    replace sample_`var' = `ss0'

    drop `fs0' `v0' `ss0' `d0'


    * -----------------------------------------------------
    * Define within-year original quartiles.
    * These define each country-year's source group.
    * -----------------------------------------------------
    astile q_`var' = `var'_ma3 ///
        if sample_`var' == 1, ///
        nq(4) by(year)


    * -----------------------------------------------------
    * Define annual quartile-entry thresholds:
    *
    * P25 = boundary between Q1 and Q2
    * P50 = boundary between Q2 and Q3
    * P75 = boundary between Q3 and Q4
    * -----------------------------------------------------
    bysort year: egen double p25_`var' = ///
        pctile(`var'_ma3) ///
        if sample_`var' == 1, p(25)

    bysort year: egen double p50_`var' = ///
        pctile(`var'_ma3) ///
        if sample_`var' == 1, p(50)

    bysort year: egen double p75_`var' = ///
        pctile(`var'_ma3) ///
        if sample_`var' == 1, p(75)

	* Near-threshold percentile anchors
	bysort year: egen double p1875_`var' = ///
		pctile(`var'_ma3) ///
		if sample_`var' == 1, p(18.75)

	bysort year: egen double p3125_`var' = ///
		pctile(`var'_ma3) ///
		if sample_`var' == 1, p(31.25)

	bysort year: egen double p4375_`var' = ///
		pctile(`var'_ma3) ///
		if sample_`var' == 1, p(43.75)

	bysort year: egen double p5625_`var' = ///
		pctile(`var'_ma3) ///
		if sample_`var' == 1, p(56.25)

	bysort year: egen double p6875_`var' = ///
		pctile(`var'_ma3) ///
		if sample_`var' == 1, p(68.75)

	bysort year: egen double p8125_`var' = ///
		pctile(`var'_ma3) ///
		if sample_`var' == 1, p(81.25)
	
    * -----------------------------------------------------
    * Diagnostic: every estimation observation must have
    * all three annual threshold values.
    * -----------------------------------------------------
	count if sample_`var' == 1 & ///
		(missing(p1875_`var') | missing(p25_`var') | ///
		 missing(p3125_`var') | missing(p4375_`var') | ///
		 missing(p50_`var') | missing(p5625_`var') | ///
		 missing(p6875_`var') | missing(p75_`var') | ///
		 missing(p8125_`var'))
	 
    if r(N) > 0 {
        display as error ///
            "Some estimation observations have missing annual quartile thresholds for `var'."
        display as error ///
            "Check whether all years contain enough usable observations."
        exit 459
    }


    * -----------------------------------------------------
    * Globals used inside bootstrap program
    * -----------------------------------------------------
    global current_endog "`var'_ma3"
    global current_iv "`var'_iv"
    global current_sample "sample_`var'"

	global current_p1875 "p1875_`var'"
	global current_p25   "p25_`var'"
	global current_p3125 "p3125_`var'"

	global current_p4375 "p4375_`var'"
	global current_p50   "p50_`var'"
	global current_p5625 "p5625_`var'"

	global current_p6875 "p6875_`var'"
	global current_p75   "p75_`var'"
	global current_p8125 "p8125_`var'"

    eststo clear


    * -----------------------------------------------------
    * Source-group loop:
    * original Q1 / Q2 / Q3 / Q4
    * -----------------------------------------------------
    forvalues b = 1/4 {

        display ">>> Source group: actual Q`b' <<<"

        global base_q = `b'

        capture drop is_base

        gen byte is_base = ///
            (q_`var' == `b' & sample_`var' == 1)


        * -------------------------------------------------
        * Country-clustered bootstrap
        *
        * m1-m4 = predicted financing under target Q1-Q4
        * d1-d4 = target financing minus actual financing
        * -------------------------------------------------
         bootstrap ///
            m1 = r(m1) ///
            m2 = r(m2) ///
            m3 = r(m3) ///
            m4 = r(m4) ///
            d1 = r(d1) ///
            d2 = r(d2) ///
            d3 = r(d3) ///
            d4 = r(d4), ///
            reps(`reps') ///
            cluster(iso3_id): ///
            my_cont_entry_cf_boot


        * -------------------------------------------------
        * Extract point estimates and bootstrap SEs
        * -------------------------------------------------
        matrix B = e(b)
        matrix V = e(V)

        forvalues t = 1/4 {

            local pos_m = colnumb(B, "m`t'")
            local pos_d = colnumb(B, "d`t'")

            local margin`t' = B[1, `pos_m']
            local diff`t'   = B[1, `pos_d']
            local se`t'     = sqrt(V[`pos_d', `pos_d'])

            * Diagonal: retained actual centrality by definition
            if `t' == `b' {
                local diff`t' = 0
                local se`t'   = .
                local pval`t' = .
            }
            else {
                if (`diff`t'' < . & `se`t'' < . & `se`t'' > 0) {
                    local z = `diff`t'' / `se`t''
                    local pval`t' = 2 * normal(-abs(`z'))
                }
                else {
                    local pval`t' = .
                }
            }
        }


        * -------------------------------------------------
        * Store results
        * -------------------------------------------------
        local margin_actual = `margin`b''

        estadd scalar margin_actual = `margin_actual'

        estadd scalar margin_target_q1 = `margin1'
        estadd scalar margin_target_q2 = `margin2'
        estadd scalar margin_target_q3 = `margin3'
        estadd scalar margin_target_q4 = `margin4'

        forvalues t = 1/4 {
            estadd scalar diff_to_q`t' = `diff`t''
            estadd scalar se_to_q`t'   = `se`t''
            estadd scalar pval_to_q`t' = `pval`t''
        }

        eststo `var'_sourceQ`b'
    }


    * -----------------------------------------------------
    * Export complete 4x4 entry-threshold matrix
    * -----------------------------------------------------
    esttab ///
        `var'_sourceQ1 ///
        `var'_sourceQ2 ///
        `var'_sourceQ3 ///
        `var'_sourceQ4 ///
        using "`var'_continuous_entry_threshold_4x4_matrix.csv", ///
        replace ///
        cells(none) ///
        mtitles("Original = Q1" "Original = Q2" ///
                "Original = Q3" "Original = Q4") ///
        scalars( ///
            "margin_actual Actual expected finance" ///
            "margin_target_q1 Predicted finance at target Q1 boundary" ///
            "margin_target_q2 Predicted finance at target Q2 boundary" ///
            "margin_target_q3 Predicted finance at target Q3 boundary" ///
            "margin_target_q4 Predicted finance at target Q4 boundary" ///
            "diff_to_q1 Change: target Q1 boundary minus actual" ///
            "se_to_q1 SE: target Q1 boundary minus actual" ///
            "pval_to_q1 P-value: target Q1 boundary minus actual" ///
            "diff_to_q2 Change: target Q2 boundary minus actual" ///
            "se_to_q2 SE: target Q2 boundary minus actual" ///
            "pval_to_q2 P-value: target Q2 boundary minus actual" ///
            "diff_to_q3 Change: target Q3 boundary minus actual" ///
            "se_to_q3 SE: target Q3 boundary minus actual" ///
            "pval_to_q3 P-value: target Q3 boundary minus actual" ///
            "diff_to_q4 Change: target Q4 boundary minus actual" ///
            "se_to_q4 SE: target Q4 boundary minus actual" ///
            "pval_to_q4 P-value: target Q4 boundary minus actual" ///
        ) ///
        sfmt(%9.4f) ///
        compress nogaps label
}


display "==================================================="
display "All continuous CF-PPML entry-threshold matrices completed."
display "==================================================="


* =========================================================
* 3.1 IV Regression: Gravity IV -- CF Approach -- Bilateral & Multilateral Adaptation
* =========================================================

global controls gdppc fdi gov_effectiveness disaster_total_impact_pc recap agri_va rural_pop co2 total_erp_pct trade
capture xtset, clear
capture tsset, clear
eststo clear

capture program drop my_cf_boot
program my_cf_boot, eclass
    version 14.0 
    capture drop v_resid_boot
    
    quietly reghdfe $current_endog $current_iv $controls, absorb(region_id year) resid
    quietly predict v_resid_boot, resid
    
    quietly ppmlhdfe $current_dv $current_endog $controls v_resid_boot, absorb(region_id year) d
	
	quietly margins, dydx($current_endog v_resid_boot) post

end

local dvs "bilateral_adaptation multilateral_adaptation"
local centralities "in_degree out_degree closeness eigenvector"

foreach dv in `dvs' {
    
    if "`dv'" == "bilateral_adaptation" {
        local prefix "bi"
        local dv_label "Bilateral"
    }
    else {
        local prefix "multi"
        local dv_label "Multilateral"
    }
    
    global current_dv "`dv'"
        
    foreach var in `centralities' {
                
        global current_endog "`var'_ma3"
        global current_iv "`var'_iv"
        
        quietly reghdfe $current_endog $current_iv $controls, absorb(region_id year) vce(cluster iso3_id)
        quietly test $current_iv
		local f_stat = r(F)
        
        capture drop new_iso3 
        set seed 123456
        quietly bootstrap, reps(1000) cluster(iso3_id) idcluster(new_iso3): my_cf_boot
        
        estadd scalar F_first = `f_stat'
        estadd local Controls "Yes"
        estadd local FE "Yes"
        estadd local DV_Type "`dv_label'"
        
        eststo `prefix'_`var' 
    }
}

esttab bi_in_degree bi_out_degree bi_closeness bi_eigenvector ///
       multi_in_degree multi_out_degree multi_closeness multi_eigenvector ///
    using "regression_results_Bi_Multi_GravityIV_CF.csv", replace ///
    keep(*_ma3 v_resid_boot) ///
    order(in_degree_ma3 out_degree_ma3 closeness_ma3 eigenvector_ma3 v_resid_boot) ///
    b(3) se(3) star(* 0.1 ** 0.05 *** 0.01) ///
    mtitles("Bi: In" "Bi: Out" "Bi: Close" "Bi: Eigen" "Multi: In" "Multi: Out" "Multi: Close" "Multi: Eigen") ///
    scalars( ///
        "DV_Type Dependent Var" /// 
        "F_first First-Stage F-stat" ///
        "Controls Controls" ///
        "FE Region & Year FE" ///
        "N Observations" ///
    ) ///
    compress nogaps label
	
	
	
* =========================================================
* 3.2 IV Regression: Bartik IV -- CF Approach -- Bilateral & Multilateral Adaptation
* =========================================================
global controls gdppc fdi gov_effectiveness disaster_total_impact_pc recap agri_va rural_pop co2 total_erp_pct trade
global current_endog closeness_ma3
global current_iv iv_seaborne_trade_total

capture xtset, clear
capture tsset, clear
eststo clear

capture program drop my_cf_boot
program my_cf_boot, eclass
    version 14.0 
    capture drop v_resid_boot
    
    quietly reghdfe $current_endog $current_iv $controls, absorb(region_id year) resid
    quietly predict v_resid_boot, resid
    
    quietly ppmlhdfe $current_dv $current_endog $controls v_resid_boot, absorb(region_id year) d
	
	quietly margins, dydx($current_endog v_resid_boot) post
end

local dvs "bilateral_adaptation multilateral_adaptation"
local centralities "in_degree out_degree closeness eigenvector"

foreach dv in `dvs' {
    
    if "`dv'" == "bilateral_adaptation" {
        local prefix "bi"
        local dv_label "Bilateral"
    }
    else {
        local prefix "multi"
        local dv_label "Multilateral"
    }
    
    global current_dv "`dv'"
        
    foreach var in `centralities' {
        
        global current_endog "`var'_ma3"
        global current_iv ln_iv_seaborne_trade_total
        
        quietly reghdfe $current_endog $current_iv $controls, absorb(region_id year) vce(cluster iso3_id)
        quietly test $current_iv
		local f_stat = r(F)
        
        capture drop new_iso3 
        set seed 1234
        quietly bootstrap, reps(1000) cluster(iso3_id) idcluster(new_iso3): my_cf_boot
        		
        estadd scalar F_first = `f_stat'
        estadd local Controls "Yes"
        estadd local FE "Yes"
        estadd local DV_Type "`dv_label'"
        		
        eststo `prefix'_`var' 
    }
}

esttab bi_in_degree bi_out_degree bi_closeness bi_eigenvector ///
       multi_in_degree multi_out_degree multi_closeness multi_eigenvector ///
    using "regression_results_Bi_Multi_BartikIV_CF.csv", replace ///
    keep(*_ma3 v_resid_boot) ///
    order(in_degree_ma3 out_degree_ma3 closeness_ma3 eigenvector_ma3 v_resid_boot) ///
    b(3) se(3) star(* 0.1 ** 0.05 *** 0.01) ///
    mtitles("Bi: In" "Bi: Out" "Bi: Close" "Bi: Eigen" "Multi: In" "Multi: Out" "Multi: Close" "Multi: Eigen") ///
	title("Average Marginal Effects (AME) from PPML with CF Approach") /// 
    scalars( ///
        "DV_Type Dependent Var" /// 
        "F_first First-Stage F-stat" ///
        "Controls Controls" ///
        "FE Region & Year FE" ///
        "N Observations" ///
    ) ///
    compress nogaps label

	
	
* =========================================================
* 4.1 IV Regression: Gravity IV -- CF Approach -- by Income Group
* =========================================================

global controls gdppc fdi gov_effectiveness disaster_total_impact_pc recap agri_va rural_pop co2 total_erp_pct trade
global depvar total_adaptation_finance

capture xtset, clear
capture tsset, clear
eststo clear

capture program drop my_cf_boot_group
program my_cf_boot_group, eclass
    version 14.0 
    capture drop v_resid_boot
    
    quietly reghdfe $current_endog $current_iv $controls, absorb(region_id year) resid
    quietly predict v_resid_boot, resid
    
    quietly ppmlhdfe $depvar $current_endog $controls v_resid_boot, absorb(region_id year) d
	
	quietly margins, dydx($current_endog v_resid_boot) post

end

levelsof incomegroup, local(groups)
local centralities in_degree out_degree closeness  eigenvector

foreach g of local groups {
    local clean_g = strtoname("`g'")
    
    preserve
    keep if incomegroup == "`g'"
    
    local i = 1
    
    foreach var in `centralities' {

	global current_endog "`var'_ma3"
        global current_iv "`var'_iv"
        
        quietly reghdfe $current_endog $current_iv $controls, absorb(region_id year) vce(cluster iso3_id)
        quietly test $current_iv
        local f_stat = r(F) 
        
        capture drop new_iso3 
        set seed 123
        quietly bootstrap, reps(1000) cluster(iso3_id) idcluster(new_iso3): my_cf_boot_group
        
        estadd scalar F_first = `f_stat'
        estadd local Controls "Yes"
        estadd local FE "Yes"
        eststo m`i'_`clean_g'
        
        local i = `i' + 1
    }
    
    restore 
}

esttab ///
    m1_Low_income m2_Low_income m3_Low_income m4_Low_income ///
    m1_Lower_middle_income m2_Lower_middle_income m3_Lower_middle_income m4_Lower_middle_income ///
    m1_Upper_middle_income m2_Upper_middle_income m3_Upper_middle_income m4_Upper_middle_income ///
    using "regression_results_GravityIV_CF_Bootstrap_by_income.csv", /// 
    replace ///
    b(3) se(3) ///
    star(* 0.1 ** 0.05 *** 0.01) ///
    compress nogaps ///
    nodepvars nomtitles ///
    keep(*_ma3 v_resid_boot) ///
    order(in_degree_ma3 out_degree_ma3 closeness_ma3 eigenvector_ma3 v_resid_boot) ///
    scalars( ///
        "F_first First-Stage F-stat" ///
        "Controls Controls" ///
        "FE Region & Year FE" ///
        "N Observations" ///
    ) ///
    label
	
	
* =========================================================
* 4.2 IV Regression: Seaborne Trade Shift-Share IV -- CF Approach -- by Income Group
* =========================================================

global controls gdppc fdi gov_effectiveness disaster_affected recap agri_va rural_pop co2 total_erp_pct trade
global depvar total_adaptation_finance

capture xtset, clear
capture tsset, clear
eststo clear

capture program drop my_cf_boot_bartik
program my_cf_boot_bartik, eclass
    version 14.0 
    capture drop v_resid_boot
    
    quietly reghdfe $current_endog $current_iv $controls, absorb(region_id year) resid
    quietly predict v_resid_boot, resid
    
    quietly ppmlhdfe $depvar $current_endog $controls v_resid_boot, absorb(region_id year) d
	quietly margins, dydx($current_endog v_resid_boot) post

end

levelsof incomegroup, local(groups)
local centralities in_degree out_degree closeness eigenvector

foreach g of local groups {
    local clean_g = strtoname("`g'")
    
    preserve
    keep if incomegroup == "`g'"
    
    local i = 1 
    
    foreach var in `centralities' {
        
        global current_endog "`var'_ma3"
        global current_iv "ln_iv_seaborne_trade_total" 
        
        quietly reghdfe $current_endog $current_iv $controls, absorb(region_id year) vce(cluster iso3_id)
        quietly test $current_iv
        local f_stat = r(F) 
        
        capture drop new_iso3 
        set seed 123
        quietly bootstrap, reps(1000) cluster(iso3_id) idcluster(new_iso3): my_cf_boot_bartik
        
        estadd scalar F_first = `f_stat'
        estadd local Controls "Yes"
        estadd local FE "Yes"
        eststo m`i'_`clean_g'
        
        local i = `i' + 1
    }
    
    restore 
}

esttab ///
    m1_Low_income m2_Low_income m3_Low_income m4_Low_income ///
    m1_Lower_middle_income m2_Lower_middle_income m3_Lower_middle_income m4_Lower_middle_income ///
    m1_Upper_middle_income m2_Upper_middle_income m3_Upper_middle_income m4_Upper_middle_income ///
    using "regression_results_BartikIV_CF_Bootstrap_by_income.csv", /// 
    replace ///
    b(3) se(3) ///
    star(* 0.1 ** 0.05 *** 0.01) ///
    compress nogaps ///
    nodepvars nomtitles ///
    keep(*_ma3 v_resid_boot) ///
    order(in_degree_ma3 out_degree_ma3 closeness_ma3 eigenvector_ma3 v_resid_boot) ///
    scalars( ///
        "F_first First-Stage F-stat" ///
        "Controls Controls" ///
        "FE Region & Year FE" ///
        "N Observations" ///
    ) ///
    label	

	
* =========================================================
* 5.1 IV Regression Moderation Effects: Gravity IV -- CF Approach -- disaster
* =========================================================

global controls gdppc fdi gov_effectiveness recap agri_va rural_pop co2 total_erp_pct trade
global moderator disaster_total_impact_pc

capture xtset, clear
capture tsset, clear
eststo clear

capture program drop my_cf_reg_boot
program my_cf_reg_boot, eclass
    version 14.0
    syntax [if] [in]
    
    marksample touse
    capture drop v_resid1_boot v_resid2_boot
    
    quietly reghdfe $current_endog ///
        $current_iv $current_iv_X_mod $moderator $controls ///
        if `touse', absorb(region_id year) resid
    
    quietly predict v_resid1_boot if e(sample), resid
    
    quietly reghdfe $current_endog_X_mod ///
        $current_iv $current_iv_X_mod $moderator $controls ///
        if `touse', absorb(region_id year) resid
    
    quietly predict v_resid2_boot if e(sample), resid
    
    quietly ppmlhdfe total_adaptation_finance ///
        c.${current_endog}##c.${moderator} ///
        $controls v_resid1_boot v_resid2_boot ///
        if `touse', absorb(region_id year) d
end


capture program drop my_cf_margins_boot
program my_cf_margins_boot, eclass
    version 14.0
    syntax [if] [in]
    
    marksample touse
    capture drop v_resid1_boot v_resid2_boot

    quietly reghdfe $current_endog ///
        $current_iv $current_iv_X_mod $moderator $controls ///
        if `touse', absorb(region_id year) resid
    
    quietly predict v_resid1_boot if e(sample), resid
    
    quietly reghdfe $current_endog_X_mod ///
        $current_iv $current_iv_X_mod $moderator $controls ///
        if `touse', absorb(region_id year) resid
    
    quietly predict v_resid2_boot if e(sample), resid
    
    quietly ppmlhdfe total_adaptation_finance ///
        c.${current_endog}##c.${moderator} ///
        $controls v_resid1_boot v_resid2_boot ///
        if `touse', absorb(region_id year) d
    
    quietly margins, dydx(${current_endog}) ///
        at(${moderator}=($p50 $p75 $p90 $p95 $p99)) post
end

capture program drop my_cf_me_diff_boot
program my_cf_me_diff_boot, rclass
    version 14.0
    syntax [if] [in]
    
    marksample touse
    capture drop v_resid1_boot v_resid2_boot
    
    quietly reghdfe $current_endog ///
        $current_iv $current_iv_X_mod $moderator $controls ///
        if `touse', absorb(region_id year) resid
    
    quietly predict v_resid1_boot if e(sample), resid
    
    quietly reghdfe $current_endog_X_mod ///
        $current_iv $current_iv_X_mod $moderator $controls ///
        if `touse', absorb(region_id year) resid
    
    quietly predict v_resid2_boot if e(sample), resid
    
    quietly ppmlhdfe total_adaptation_finance ///
        c.${current_endog}##c.${moderator} ///
        $controls v_resid1_boot v_resid2_boot ///
        if `touse', absorb(region_id year) d
    
    quietly margins, dydx(${current_endog}) ///
        at(${moderator}=($p50 $p75 $p90 $p95 $p99)) post
    
    return scalar me_p50 = _b[1._at]
    return scalar me_p75 = _b[2._at]
    return scalar me_p90 = _b[3._at]
    return scalar me_p95 = _b[4._at]
    return scalar me_p99 = _b[5._at]
    
    return scalar diff_p75_p50 = _b[2._at] - _b[1._at]
    return scalar diff_p90_p50 = _b[3._at] - _b[1._at]
    return scalar diff_p95_p50 = _b[4._at] - _b[1._at]
    return scalar diff_p99_p50 = _b[5._at] - _b[1._at]
    return scalar diff_p90_p75 = _b[3._at] - _b[2._at]
    return scalar diff_p95_p90 = _b[4._at] - _b[3._at]
	return scalar diff_p95_p75 = _b[4._at] - _b[2._at]
    return scalar diff_p99_p95 = _b[5._at] - _b[4._at]
end

local centralities in_degree out_degree closeness eigenvector

foreach var in `centralities' {
    
    global current_endog "`var'_ma3"
    global current_iv "`var'_iv"
    
    capture drop `var'_X_mod `var'_iv_X_mod
    quietly gen `var'_X_mod = ${current_endog} * ${moderator}
    quietly gen `var'_iv_X_mod = ${current_iv} * ${moderator}
    
    global current_endog_X_mod "`var'_X_mod"
    global current_iv_X_mod "`var'_iv_X_mod"
    
    capture drop estimation_sample
    
    quietly reghdfe total_adaptation_finance ///
        ${current_endog} ${current_endog_X_mod} ///
        ${current_iv} ${current_iv_X_mod} ///
        ${moderator} $controls, ///
        absorb(region_id year)
    
    gen byte estimation_sample = e(sample)
    
    quietly summarize $moderator if estimation_sample == 1, detail
    
    global p50 = r(p50)
    global p75 = r(p75)
    global p90 = r(p90)
    global p95 = r(p95)
    global p99 = r(p99)
    
    display "P50 = $p50"
    display "P75 = $p75"
    display "P90 = $p90"
    display "P95 = $p95"
    display "P99 = $p99"
    
    capture drop new_iso3_reg
    set seed 123456
    
    quietly bootstrap, reps(1000) cluster(iso3_id) idcluster(new_iso3_reg): ///
        my_cf_reg_boot if estimation_sample == 1
    
    estadd local Controls "Yes"
    estadd local FE "Yes"
    eststo reg_`var'
    
    capture drop new_iso3_mar
    set seed 123456
    
    quietly bootstrap, reps(1000) cluster(iso3_id) idcluster(new_iso3_mar): ///
        my_cf_margins_boot if estimation_sample == 1
    
    eststo margins_`var'
	
    capture drop new_iso3_diff
    set seed 123456
    
    quietly bootstrap ///
        me_p50=r(me_p50) ///
        me_p75=r(me_p75) ///
        me_p90=r(me_p90) ///
        me_p95=r(me_p95) ///
        me_p99=r(me_p99) ///
        diff_p75_p50=r(diff_p75_p50) ///
        diff_p90_p50=r(diff_p90_p50) ///
        diff_p95_p50=r(diff_p95_p50) ///
        diff_p99_p50=r(diff_p99_p50) ///
        diff_p90_p75=r(diff_p90_p75) ///
        diff_p95_p75=r(diff_p95_p75) ///		
        diff_p95_p90=r(diff_p95_p90) ///
        diff_p99_p95=r(diff_p99_p95), ///
        reps(1000) cluster(iso3_id) idcluster(new_iso3_diff): ///
        my_cf_me_diff_boot if estimation_sample == 1
    
    eststo mediff_`var'
}


esttab reg_in_degree reg_out_degree reg_closeness reg_eigenvector ///
    using "Moderation_disaster_Table1_Main_Regression_CF.csv", replace ///
    keep(*_ma3 ${moderator} *#* v_resid1_boot v_resid2_boot) ///
    b(3) se(3) star(* 0.1 ** 0.05 *** 0.01) ///
    mtitles("In-Degree" "Out-Degree" "Closeness" "Eigenvector") ///
    scalars("Controls Controls" "FE Region & Year FE") ///
    compress nogaps label ///
    title("Second-Stage PPML Regression Results")


esttab margins_in_degree margins_out_degree margins_closeness margins_eigenvector ///
    using "Moderation_disaster_Table2_Marginal_Effects_CF.csv", replace ///
    b(3) se(3) star(* 0.1 ** 0.05 *** 0.01) ///
    mtitles("In-Degree" "Out-Degree" "Closeness" "Eigenvector") ///
    coeflabels( ///
        1._at "ME at P50" ///
        2._at "ME at P75" ///
        3._at "ME at P90" ///
        4._at "ME at P95" ///
        5._at "ME at P99" ///
    ) ///
    compress nogaps label ///
    title("Marginal Effects at Different Percentiles")

esttab mediff_in_degree mediff_out_degree mediff_closeness mediff_eigenvector ///
    using "Moderation_disaster_Table3_Marginal_Effects_Differences_CF.csv", replace ///
    keep(diff_p75_p50 diff_p90_p50 diff_p95_p50 diff_p99_p50 ///
         diff_p90_p75 diff_p95_p75 diff_p95_p90 diff_p99_p95) ///
    b(3) se(3) star(* 0.1 ** 0.05 *** 0.01) ///
    mtitles("In-Degree" "Out-Degree" "Closeness" "Eigenvector") ///
    coeflabels( ///
        diff_p75_p50 "ME(P75) - ME(P50)" ///
        diff_p90_p50 "ME(P90) - ME(P50)" ///
        diff_p95_p50 "ME(P95) - ME(P50)" ///
        diff_p99_p50 "ME(P99) - ME(P50)" ///
        diff_p90_p75 "ME(P90) - ME(P75)" ///
		diff_p95_p75 "ME(P95) - ME(P75)" ///
        diff_p95_p90 "ME(P95) - ME(P90)" ///
        diff_p99_p95 "ME(P99) - ME(P95)" ///
    ) ///
    compress nogaps label ///
    title("Differences in Marginal Effects across Disaster Percentiles")	


	
* =========================================================
* 5.2 IV Regression Moderation Effects: Gravity IV -- CF Approach -- political alignment
* =========================================================

global controls gdppc fdi gov_effectiveness recap agri_va rural_pop co2 total_erp_pct trade
global moderator idealpoint_align

capture xtset, clear
capture tsset, clear
eststo clear

capture program drop my_cf_reg_boot
program my_cf_reg_boot, eclass
    version 14.0
    syntax [if] [in]
    
    marksample touse
    capture drop v_resid1_boot v_resid2_boot
    
    quietly reghdfe $current_endog ///
        $current_iv $current_iv_X_mod $moderator $controls ///
        if `touse', absorb(region_id year) resid
    
    quietly predict v_resid1_boot if e(sample), resid
    
    quietly reghdfe $current_endog_X_mod ///
        $current_iv $current_iv_X_mod $moderator $controls ///
        if `touse', absorb(region_id year) resid
    
    quietly predict v_resid2_boot if e(sample), resid
    
    quietly ppmlhdfe total_adaptation_finance ///
        c.${current_endog}##c.${moderator} ///
        $controls v_resid1_boot v_resid2_boot ///
        if `touse', absorb(region_id year) d
end


capture program drop my_cf_margins_boot
program my_cf_margins_boot, eclass
    version 14.0
    syntax [if] [in]
    
    marksample touse
    capture drop v_resid1_boot v_resid2_boot
    
    quietly reghdfe $current_endog ///
        $current_iv $current_iv_X_mod $moderator $controls ///
        if `touse', absorb(region_id year) resid
    
    quietly predict v_resid1_boot if e(sample), resid
    
    quietly reghdfe $current_endog_X_mod ///
        $current_iv $current_iv_X_mod $moderator $controls ///
        if `touse', absorb(region_id year) resid
    
    quietly predict v_resid2_boot if e(sample), resid
    
    quietly ppmlhdfe total_adaptation_finance ///
        c.${current_endog}##c.${moderator} ///
        $controls v_resid1_boot v_resid2_boot ///
        if `touse', absorb(region_id year) d
    
    quietly margins, dydx(${current_endog}) ///
        at(${moderator}=($p10 $p25 $p50 $p75 $p90)) post
end

capture program drop my_cf_me_diff_boot
program my_cf_me_diff_boot, rclass
    version 14.0
    syntax [if] [in]
    
    marksample touse
    capture drop v_resid1_boot v_resid2_boot
    
    quietly reghdfe $current_endog ///
        $current_iv $current_iv_X_mod $moderator $controls ///
        if `touse', absorb(region_id year) resid
    
    quietly predict v_resid1_boot if e(sample), resid
    
    quietly reghdfe $current_endog_X_mod ///
        $current_iv $current_iv_X_mod $moderator $controls ///
        if `touse', absorb(region_id year) resid
    
    quietly predict v_resid2_boot if e(sample), resid
    
    quietly ppmlhdfe total_adaptation_finance ///
        c.${current_endog}##c.${moderator} ///
        $controls v_resid1_boot v_resid2_boot ///
        if `touse', absorb(region_id year) d
    
    quietly margins, dydx(${current_endog}) ///
        at(${moderator}=($p10 $p25 $p50 $p75 $p90)) post
    
    return scalar me_p10 = _b[1._at]
    return scalar me_p25 = _b[2._at]
    return scalar me_p50 = _b[3._at]
    return scalar me_p75 = _b[4._at]
    return scalar me_p90 = _b[5._at]
    
    return scalar diff_p10_p50 = _b[1._at] - _b[3._at]
    return scalar diff_p25_p50 = _b[2._at] - _b[3._at]
    return scalar diff_p75_p50 = _b[4._at] - _b[3._at]
    return scalar diff_p90_p50 = _b[5._at] - _b[3._at]
    return scalar diff_p75_p25 = _b[4._at] - _b[2._at]
	return scalar diff_p75_p10 = _b[4._at] - _b[1._at]
    return scalar diff_p90_p10 = _b[5._at] - _b[1._at]
    return scalar diff_p25_p10 = _b[2._at] - _b[1._at]
    return scalar diff_p50_p25 = _b[3._at] - _b[2._at]
    return scalar diff_p75_p50_adj = _b[4._at] - _b[3._at]
    return scalar diff_p90_p75 = _b[5._at] - _b[4._at]
end

local centralities in_degree out_degree closeness eigenvector

foreach var in `centralities' {
    
    global current_endog "`var'_ma3"
    global current_iv "`var'_iv"
    
    capture drop `var'_X_mod `var'_iv_X_mod
    quietly gen `var'_X_mod = ${current_endog} * ${moderator}
    quietly gen `var'_iv_X_mod = ${current_iv} * ${moderator}
    
    global current_endog_X_mod "`var'_X_mod"
    global current_iv_X_mod "`var'_iv_X_mod"
    
    capture drop estimation_sample
    
    quietly reghdfe total_adaptation_finance ///
        ${current_endog} ${current_endog_X_mod} ///
        ${current_iv} ${current_iv_X_mod} ///
        ${moderator} $controls, ///
        absorb(region_id year)
    
    gen byte estimation_sample = e(sample)
    
    quietly summarize $moderator if estimation_sample == 1, detail
    
	global p10 = r(p10)   
	global p25 = r(p25)   
    global p50 = r(p50)
    global p75 = r(p75)
    global p90 = r(p90)
	
    display "P10 = $p10"
    display "P25 = $p25"
    display "P50 = $p50"
    display "P75 = $p75"
    display "P90 = $p90"
    
    capture drop new_iso3_reg
    set seed 123456
    
    quietly bootstrap, reps(1000) cluster(iso3_id) idcluster(new_iso3_reg): ///
        my_cf_reg_boot if estimation_sample == 1
    
    estadd local Controls "Yes"
    estadd local FE "Yes"
    eststo reg_`var'
    
    capture drop new_iso3_mar
    set seed 123456
    
    quietly bootstrap, reps(1000) cluster(iso3_id) idcluster(new_iso3_mar): ///
        my_cf_margins_boot if estimation_sample == 1
    
    eststo margins_`var'
	
    capture drop new_iso3_diff
    set seed 123456
    
    quietly bootstrap ///
        me_p10=r(me_p10) ///
        me_p25=r(me_p25) ///
        me_p50=r(me_p50) ///
        me_p75=r(me_p75) ///
        me_p90=r(me_p90) ///
        diff_p10_p50=r(diff_p10_p50) ///
        diff_p25_p50=r(diff_p25_p50) ///
        diff_p75_p50=r(diff_p75_p50) ///
        diff_p90_p50=r(diff_p90_p50) ///
        diff_p75_p25=r(diff_p75_p25) ///
        diff_p75_p10=r(diff_p75_p10) ///
        diff_p90_p10=r(diff_p90_p10) ///
        diff_p25_p10=r(diff_p25_p10) ///
        diff_p50_p25=r(diff_p50_p25) ///
        diff_p75_p50_adj=r(diff_p75_p50_adj) ///
        diff_p90_p75=r(diff_p90_p75), ///
        reps(1000) cluster(iso3_id) idcluster(new_iso3_diff): ///
        my_cf_me_diff_boot if estimation_sample == 1
    
    eststo mediff_`var'
}


esttab reg_in_degree reg_out_degree reg_closeness reg_eigenvector ///
    using "Moderation_idealpoint_Table1_Main_Regression_CF.csv", replace ///
    keep(*_ma3 ${moderator} *#* v_resid1_boot v_resid2_boot) ///
    b(3) se(3) star(* 0.1 ** 0.05 *** 0.01) ///
    mtitles("In-Degree" "Out-Degree" "Closeness" "Eigenvector") ///
    scalars("Controls Controls" "FE Region & Year FE") ///
    compress nogaps label ///
    title("Second-Stage PPML Regression Results")


esttab margins_in_degree margins_out_degree margins_closeness margins_eigenvector ///
    using "Moderation_idealpoint_Table2_Main_Regression_CF.csv", replace ///
    b(3) se(3) star(* 0.1 ** 0.05 *** 0.01) ///
    mtitles("In-Degree" "Out-Degree" "Closeness" "Eigenvector") ///
    coeflabels( ///
        1._at "ME at P10" ///
        2._at "ME at P25" ///
        3._at "ME at P50" ///
        4._at "ME at P75" ///
        5._at "ME at P90" ///
    ) ///
    compress nogaps label ///
    title("Marginal Effects at Different Percentiles")

esttab mediff_in_degree mediff_out_degree mediff_closeness mediff_eigenvector ///
    using "Moderation_idealpoint_Table3_Marginal_Effects_Differences_CF.csv", replace ///
    keep(diff_p10_p50 diff_p25_p50 diff_p75_p50 diff_p90_p50 ///
         diff_p75_p25 diff_p75_p10 diff_p90_p10 ///
         diff_p25_p10 diff_p50_p25 diff_p75_p50_adj diff_p90_p75) ///
    b(3) se(3) star(* 0.1 ** 0.05 *** 0.01) ///
    mtitles("In-Degree" "Out-Degree" "Closeness" "Eigenvector") ///
    coeflabels( ///
        diff_p10_p50 "ME(P10) - ME(P50)" ///
        diff_p25_p50 "ME(P25) - ME(P50)" ///
        diff_p75_p50 "ME(P75) - ME(P50)" ///
        diff_p90_p50 "ME(P90) - ME(P50)" ///
        diff_p75_p25 "ME(P75) - ME(P25)" ///
        diff_p75_p10 "ME(P75) - ME(P10)" ///
        diff_p90_p10 "ME(P90) - ME(P10)" ///
        diff_p25_p10 "ME(P25) - ME(P10)" ///
        diff_p50_p25 "ME(P50) - ME(P25)" ///
        diff_p75_p50_adj "ME(P75) - ME(P50)" ///
        diff_p90_p75 "ME(P90) - ME(P75)" ///
    ) ///
    compress nogaps label ///
    title("Differences in Marginal Effects across Ideal-Point Alignment Percentiles")
	
	
* =========================================================
* 5.3 IV Regression Moderation Effects: Bartik IV -- CF Approach -- disaster
* =========================================================

global controls gdppc fdi gov_effectiveness recap agri_va rural_pop co2 total_erp_pct trade
global moderator disaster_total_impact_pc

capture xtset, clear
capture tsset, clear
eststo clear

capture program drop my_cf_reg_boot
program my_cf_reg_boot, eclass
    version 14.0
    syntax [if] [in]
    
    marksample touse
    capture drop v_resid1_boot v_resid2_boot
    
    quietly reghdfe $current_endog ///
        $current_iv $current_iv_X_mod $moderator $controls ///
        if `touse', absorb(region_id year) resid
    
    quietly predict v_resid1_boot if e(sample), resid
    
    quietly reghdfe $current_endog_X_mod ///
        $current_iv $current_iv_X_mod $moderator $controls ///
        if `touse', absorb(region_id year) resid
    
    quietly predict v_resid2_boot if e(sample), resid
    
    quietly ppmlhdfe total_adaptation_finance ///
        c.${current_endog}##c.${moderator} ///
        $controls v_resid1_boot v_resid2_boot ///
        if `touse', absorb(region_id year) d
end


capture program drop my_cf_margins_boot
program my_cf_margins_boot, eclass
    version 14.0
    syntax [if] [in]
    
    marksample touse
    capture drop v_resid1_boot v_resid2_boot
    
    quietly reghdfe $current_endog ///
        $current_iv $current_iv_X_mod $moderator $controls ///
        if `touse', absorb(region_id year) resid
    
    quietly predict v_resid1_boot if e(sample), resid
    
    quietly reghdfe $current_endog_X_mod ///
        $current_iv $current_iv_X_mod $moderator $controls ///
        if `touse', absorb(region_id year) resid
    
    quietly predict v_resid2_boot if e(sample), resid
    
    quietly ppmlhdfe total_adaptation_finance ///
        c.${current_endog}##c.${moderator} ///
        $controls v_resid1_boot v_resid2_boot ///
        if `touse', absorb(region_id year) d
    
    quietly margins, dydx(${current_endog}) ///
        at(${moderator}=($p50 $p75 $p90 $p95 $p99)) post
end

capture program drop my_cf_me_diff_boot
program my_cf_me_diff_boot, rclass
    version 14.0
    syntax [if] [in]
    
    marksample touse
    capture drop v_resid1_boot v_resid2_boot
    
    quietly reghdfe $current_endog ///
        $current_iv $current_iv_X_mod $moderator $controls ///
        if `touse', absorb(region_id year) resid
    
    quietly predict v_resid1_boot if e(sample), resid
    
    quietly reghdfe $current_endog_X_mod ///
        $current_iv $current_iv_X_mod $moderator $controls ///
        if `touse', absorb(region_id year) resid
    
    quietly predict v_resid2_boot if e(sample), resid
    
    quietly ppmlhdfe total_adaptation_finance ///
        c.${current_endog}##c.${moderator} ///
        $controls v_resid1_boot v_resid2_boot ///
        if `touse', absorb(region_id year) d
    
    quietly margins, dydx(${current_endog}) ///
        at(${moderator}=($p50 $p75 $p90 $p95 $p99)) post
    
    return scalar me_p50 = _b[1._at]
    return scalar me_p75 = _b[2._at]
    return scalar me_p90 = _b[3._at]
    return scalar me_p95 = _b[4._at]
    return scalar me_p99 = _b[5._at]
    
    return scalar diff_p75_p50 = _b[2._at] - _b[1._at]
    return scalar diff_p90_p50 = _b[3._at] - _b[1._at]
    return scalar diff_p95_p50 = _b[4._at] - _b[1._at]
    return scalar diff_p99_p50 = _b[5._at] - _b[1._at]
    return scalar diff_p90_p75 = _b[3._at] - _b[2._at]
	return scalar diff_p95_p75 = _b[4._at] - _b[2._at]
    return scalar diff_p95_p90 = _b[4._at] - _b[3._at]
    return scalar diff_p99_p95 = _b[5._at] - _b[4._at]
end

local centralities in_degree out_degree closeness eigenvector

foreach var in `centralities' {
    
    global current_endog "`var'_ma3"
    global current_iv "ln_iv_seaborne_trade_total"
    
    capture drop `var'_X_mod `var'_iv_X_mod
    quietly gen `var'_X_mod = ${current_endog} * ${moderator}
    quietly gen `var'_iv_X_mod = ${current_iv} * ${moderator}
    
    global current_endog_X_mod "`var'_X_mod"
    global current_iv_X_mod "`var'_iv_X_mod"
    
    capture drop estimation_sample
    
    quietly reghdfe total_adaptation_finance ///
        ${current_endog} ${current_endog_X_mod} ///
        ${current_iv} ${current_iv_X_mod} ///
        ${moderator} $controls, ///
        absorb(region_id year)
    
    gen byte estimation_sample = e(sample)
    
    quietly summarize $moderator if estimation_sample == 1, detail
    
    global p50 = r(p50)
    global p75 = r(p75)
    global p90 = r(p90)
    global p95 = r(p95)
    global p99 = r(p99)
    
    display "P50 = $p50"
    display "P75 = $p75"
    display "P90 = $p90"
    display "P95 = $p95"
    display "P99 = $p99"
    
    capture drop new_iso3_reg
    set seed 123456
    
    quietly bootstrap, reps(1000) cluster(iso3_id) idcluster(new_iso3_reg): ///
        my_cf_reg_boot if estimation_sample == 1
    
    estadd local Controls "Yes"
    estadd local FE "Yes"
    eststo reg_`var'
    
    capture drop new_iso3_mar
    set seed 123456
    
    quietly bootstrap, reps(1000) cluster(iso3_id) idcluster(new_iso3_mar): ///
        my_cf_margins_boot if estimation_sample == 1
    
    eststo margins_`var'
	
    capture drop new_iso3_diff
    set seed 123456
    
    quietly bootstrap ///
        me_p50=r(me_p50) ///
        me_p75=r(me_p75) ///
        me_p90=r(me_p90) ///
        me_p95=r(me_p95) ///
        me_p99=r(me_p99) ///
        diff_p75_p50=r(diff_p75_p50) ///
        diff_p90_p50=r(diff_p90_p50) ///
        diff_p95_p50=r(diff_p95_p50) ///
        diff_p99_p50=r(diff_p99_p50) ///
        diff_p90_p75=r(diff_p90_p75) ///
		diff_p95_p75=r(diff_p95_p75) ///
        diff_p95_p90=r(diff_p95_p90) ///
        diff_p99_p95=r(diff_p99_p95), ///
        reps(1000) cluster(iso3_id) idcluster(new_iso3_diff): ///
        my_cf_me_diff_boot if estimation_sample == 1
    
    eststo mediff_`var'

}


esttab reg_in_degree reg_out_degree reg_closeness reg_eigenvector ///
    using "Moderation_disaster_Table1_Main_Regression_CF_Bartik.csv", replace ///
    keep(*_ma3 ${moderator} *#* v_resid1_boot v_resid2_boot) ///
    b(3) se(3) star(* 0.1 ** 0.05 *** 0.01) ///
    mtitles("In-Degree" "Out-Degree" "Closeness" "Eigenvector") ///
    scalars("Controls Controls" "FE Region & Year FE") ///
    compress nogaps label ///
    title("Second-Stage PPML Regression Results")


esttab margins_in_degree margins_out_degree margins_closeness margins_eigenvector ///
    using "Moderation_disaster_Table2_Marginal_Effects_CF_Bartik.csv", replace ///
    b(3) se(3) star(* 0.1 ** 0.05 *** 0.01) ///
    mtitles("In-Degree" "Out-Degree" "Closeness" "Eigenvector") ///
    coeflabels( ///
        1._at "ME at P50" ///
        2._at "ME at P75" ///
        3._at "ME at P90" ///
        4._at "ME at P95" ///
        5._at "ME at P99" ///
    ) ///
    compress nogaps label ///
    title("Marginal Effects at Different Percentiles")

esttab mediff_in_degree mediff_out_degree mediff_closeness mediff_eigenvector ///
    using "Moderation_disaster_Table3_Marginal_Effects_Differences_CF_Bartik.csv", replace ///
    keep(diff_p75_p50 diff_p90_p50 diff_p95_p50 diff_p99_p50 ///
         diff_p90_p75 diff_p95_p75 diff_p95_p90 diff_p99_p95) ///
    b(3) se(3) star(* 0.1 ** 0.05 *** 0.01) ///
    mtitles("In-Degree" "Out-Degree" "Closeness" "Eigenvector") ///
    coeflabels( ///
        diff_p75_p50 "ME(P75) - ME(P50)" ///
        diff_p90_p50 "ME(P90) - ME(P50)" ///
        diff_p95_p50 "ME(P95) - ME(P50)" ///
        diff_p99_p50 "ME(P99) - ME(P50)" ///
        diff_p90_p75 "ME(P90) - ME(P75)" ///
		diff_p95_p75 "ME(P95) - ME(P75)" ///
        diff_p95_p90 "ME(P95) - ME(P90)" ///
        diff_p99_p95 "ME(P99) - ME(P95)" ///
    ) ///
    compress nogaps label ///
    title("Differences in Marginal Effects across Disaster Percentiles")	

	
* =========================================================
* 5.4 IV Regression Moderation Effects: Bartik IV -- CF Approach -- political alignment
* =========================================================

global controls gdppc fdi gov_effectiveness recap agri_va rural_pop co2 total_erp_pct trade
global moderator idealpoint_align

capture xtset, clear
capture tsset, clear
eststo clear

capture program drop my_cf_reg_boot
program my_cf_reg_boot, eclass
    version 14.0
    syntax [if] [in]
    
    marksample touse
    capture drop v_resid1_boot v_resid2_boot
    
    quietly reghdfe $current_endog ///
        $current_iv $current_iv_X_mod $moderator $controls ///
        if `touse', absorb(region_id year) resid
    
    quietly predict v_resid1_boot if e(sample), resid
    
    quietly reghdfe $current_endog_X_mod ///
        $current_iv $current_iv_X_mod $moderator $controls ///
        if `touse', absorb(region_id year) resid
    
    quietly predict v_resid2_boot if e(sample), resid
    
    quietly ppmlhdfe total_adaptation_finance ///
        c.${current_endog}##c.${moderator} ///
        $controls v_resid1_boot v_resid2_boot ///
        if `touse', absorb(region_id year) d
end


capture program drop my_cf_margins_boot
program my_cf_margins_boot, eclass
    version 14.0
    syntax [if] [in]
    
    marksample touse
    capture drop v_resid1_boot v_resid2_boot
    
    quietly reghdfe $current_endog ///
        $current_iv $current_iv_X_mod $moderator $controls ///
        if `touse', absorb(region_id year) resid
    
    quietly predict v_resid1_boot if e(sample), resid
    
    quietly reghdfe $current_endog_X_mod ///
        $current_iv $current_iv_X_mod $moderator $controls ///
        if `touse', absorb(region_id year) resid
    
    quietly predict v_resid2_boot if e(sample), resid
    
    quietly ppmlhdfe total_adaptation_finance ///
        c.${current_endog}##c.${moderator} ///
        $controls v_resid1_boot v_resid2_boot ///
        if `touse', absorb(region_id year) d
    
    quietly margins, dydx(${current_endog}) ///
        at(${moderator}=($p10 $p25 $p50 $p75 $p90)) post
end

capture program drop my_cf_me_diff_boot
program my_cf_me_diff_boot, rclass
    version 14.0
    syntax [if] [in]
    
    marksample touse
    capture drop v_resid1_boot v_resid2_boot
    
    quietly reghdfe $current_endog ///
        $current_iv $current_iv_X_mod $moderator $controls ///
        if `touse', absorb(region_id year) resid
    
    quietly predict v_resid1_boot if e(sample), resid
    
    quietly reghdfe $current_endog_X_mod ///
        $current_iv $current_iv_X_mod $moderator $controls ///
        if `touse', absorb(region_id year) resid
    
    quietly predict v_resid2_boot if e(sample), resid
    
    quietly ppmlhdfe total_adaptation_finance ///
        c.${current_endog}##c.${moderator} ///
        $controls v_resid1_boot v_resid2_boot ///
        if `touse', absorb(region_id year) d
    
    quietly margins, dydx(${current_endog}) ///
        at(${moderator}=($p10 $p25 $p50 $p75 $p90)) post
    
    return scalar me_p10 = _b[1._at]
    return scalar me_p25 = _b[2._at]
    return scalar me_p50 = _b[3._at]
    return scalar me_p75 = _b[4._at]
    return scalar me_p90 = _b[5._at]
    
    return scalar diff_p10_p50 = _b[1._at] - _b[3._at]
    return scalar diff_p25_p50 = _b[2._at] - _b[3._at]
    return scalar diff_p75_p50 = _b[4._at] - _b[3._at]
    return scalar diff_p90_p50 = _b[5._at] - _b[3._at]
    
    return scalar diff_p75_p25 = _b[4._at] - _b[2._at]
    return scalar diff_p75_p10 = _b[4._at] - _b[1._at]
    return scalar diff_p90_p10 = _b[5._at] - _b[1._at]
    
    return scalar diff_p25_p10 = _b[2._at] - _b[1._at]
    return scalar diff_p50_p25 = _b[3._at] - _b[2._at]
    return scalar diff_p90_p75 = _b[5._at] - _b[4._at]
end

local centralities in_degree out_degree closeness eigenvector

foreach var in `centralities' {
    
    global current_endog "`var'_ma3"
    global current_iv "ln_iv_seaborne_trade_total"
    
    capture drop `var'_X_mod `var'_iv_X_mod
    quietly gen `var'_X_mod = ${current_endog} * ${moderator}
    quietly gen `var'_iv_X_mod = ${current_iv} * ${moderator}
    
    global current_endog_X_mod "`var'_X_mod"
    global current_iv_X_mod "`var'_iv_X_mod"
    
    capture drop estimation_sample
    
    quietly reghdfe total_adaptation_finance ///
        ${current_endog} ${current_endog_X_mod} ///
        ${current_iv} ${current_iv_X_mod} ///
        ${moderator} $controls, ///
        absorb(region_id year)
    
    gen byte estimation_sample = e(sample)
    
    quietly summarize $moderator if estimation_sample == 1, detail
    
	global p10 = r(p10)   
	global p25 = r(p25)   
    global p50 = r(p50)
    global p75 = r(p75)
    global p90 = r(p90)
	
    display "P10 = $p10"
    display "P25 = $p25"
    display "P50 = $p50"
    display "P75 = $p75"
    display "P90 = $p90"
    
    capture drop new_iso3_reg
    set seed 123456
    
    quietly bootstrap, reps(1000) cluster(iso3_id) idcluster(new_iso3_reg): ///
        my_cf_reg_boot if estimation_sample == 1
    
    estadd local Controls "Yes"
    estadd local FE "Yes"
    eststo reg_`var'
    
    capture drop new_iso3_mar
    set seed 123456
    
    quietly bootstrap, reps(1000) cluster(iso3_id) idcluster(new_iso3_mar): ///
        my_cf_margins_boot if estimation_sample == 1
    
    eststo margins_`var'

    capture drop new_iso3_diff
    set seed 123456
    
    quietly bootstrap ///
        me_p10=r(me_p10) ///
        me_p25=r(me_p25) ///
        me_p50=r(me_p50) ///
        me_p75=r(me_p75) ///
        me_p90=r(me_p90) ///
        diff_p10_p50=r(diff_p10_p50) ///
        diff_p25_p50=r(diff_p25_p50) ///
        diff_p75_p50=r(diff_p75_p50) ///
        diff_p90_p50=r(diff_p90_p50) ///
        diff_p75_p25=r(diff_p75_p25) ///
        diff_p75_p10=r(diff_p75_p10) ///
        diff_p90_p10=r(diff_p90_p10) ///
        diff_p25_p10=r(diff_p25_p10) ///
        diff_p50_p25=r(diff_p50_p25) ///
        diff_p90_p75=r(diff_p90_p75), ///
        reps(1000) cluster(iso3_id) idcluster(new_iso3_diff): ///
        my_cf_me_diff_boot if estimation_sample == 1
    
    eststo mediff_`var'

}


esttab reg_in_degree reg_out_degree reg_closeness reg_eigenvector ///
    using "Moderation_idealpoint_Table1_Main_Regression_CF_Bartik.csv", replace ///
    keep(*_ma3 ${moderator} *#* v_resid1_boot v_resid2_boot) ///
    b(3) se(3) star(* 0.1 ** 0.05 *** 0.01) ///
    mtitles("In-Degree" "Out-Degree" "Closeness" "Eigenvector") ///
    scalars("Controls Controls" "FE Region & Year FE") ///
    compress nogaps label ///
    title("Second-Stage PPML Regression Results")


esttab margins_in_degree margins_out_degree margins_closeness margins_eigenvector ///
    using "Moderation_idealpoint_Table2_Main_Regression_CF_Bartik.csv", replace ///
    b(3) se(3) star(* 0.1 ** 0.05 *** 0.01) ///
    mtitles("In-Degree" "Out-Degree" "Closeness" "Eigenvector") ///
    coeflabels( ///
        1._at "ME at P10" ///
        2._at "ME at P25" ///
        3._at "ME at P50" ///
        4._at "ME at P75" ///
        5._at "ME at P90" ///
    ) ///
    compress nogaps label ///
    title("Marginal Effects at Different Percentiles")

esttab mediff_in_degree mediff_out_degree mediff_closeness mediff_eigenvector ///
    using "Moderation_idealpoint_Table3_Marginal_Effects_Differences_CF_Bartik.csv", replace ///
    keep(diff_p10_p50 diff_p25_p50 diff_p75_p50 diff_p90_p50 ///
         diff_p75_p25 diff_p75_p10 diff_p90_p10 ///
         diff_p25_p10 diff_p50_p25 diff_p75_p50_adj diff_p90_p75) ///
    b(3) se(3) star(* 0.1 ** 0.05 *** 0.01) ///
    mtitles("In-Degree" "Out-Degree" "Closeness" "Eigenvector") ///
    coeflabels( ///
        diff_p10_p50 "ME(P10) - ME(P50)" ///
        diff_p25_p50 "ME(P25) - ME(P50)" ///
        diff_p75_p50 "ME(P75) - ME(P50)" ///
        diff_p90_p50 "ME(P90) - ME(P50)" ///
        diff_p75_p25 "ME(P75) - ME(P25)" ///
        diff_p75_p10 "ME(P75) - ME(P10)" ///
        diff_p90_p10 "ME(P90) - ME(P10)" ///
        diff_p25_p10 "ME(P25) - ME(P10)" ///
        diff_p50_p25 "ME(P50) - ME(P25)" ///
        diff_p75_p50_adj "ME(P75) - ME(P50)" ///
        diff_p90_p75 "ME(P90) - ME(P75)" ///
    ) ///
    compress nogaps label ///
    title("Differences in Marginal Effects across Ideal-Point Alignment Percentiles")

	
* =========================================================
* 6.1 IV Regression Robustness Tests -- Control Function Approach (CF)
* =========================================================

global controls gdppc fdi gov_effectiveness disaster_affected recap agri_va rural_pop co2 total_erp_pct trade
global centralities in_degree out_degree closeness eigenvector

foreach var in $centralities {
    capture drop `var'_ma5
    rangestat (mean) `var', interval(year -5 -1) by(iso3_id)
    rename `var'_mean `var'_ma5
    label var `var'_ma5 "MA5 of `: var label `var''"
}

capture xtset iso3_id year 
global lag_controls ""
foreach var in $controls {
    capture drop l_`var'
    gen l_`var' = L.`var'
    global lag_controls "$lag_controls l_`var'"
}

local vars_to_winsor "total_adaptation_finance"
foreach var in $centralities {
    local vars_to_winsor "`vars_to_winsor' `var'_ma3 `var'_iv"
}
capture drop *_w 
winsor2 `vars_to_winsor', cuts(1 99) suffix(_w)

* ---------------------------------------------------------
capture xtset, clear
capture tsset, clear
eststo clear

capture program drop my_cf_boot_robust
program my_cf_boot_robust, eclass
    version 14.0 
    capture drop v_resid_boot
    
    quietly reghdfe $current_endog $current_iv $current_controls $if_cond, absorb(region_id year) resid
    quietly predict v_resid_boot, resid
    
    quietly ppmlhdfe $depvar $current_endog $current_controls v_resid_boot $if_cond, absorb(region_id year) d(`ppml_d')
	quietly margins, dydx($current_endog v_resid_boot) post
end

* ---------------------------------------------------------

* --- Group 1: Alternative Measure (MA5) ---
foreach var in $centralities {
    display "Group 1 (MA5 IV): `var'"
    global depvar "total_adaptation_finance"
    global current_endog "`var'_ma5"

    global current_iv "`var'_iv"  
    global current_controls "$controls"
    global if_cond ""
    
    quietly reghdfe $current_endog $current_iv $current_controls, absorb(region_id year) vce(cluster iso3_id)
    quietly test $current_iv
	local f_stat = r(F)
    
    capture drop new_iso3 
    set seed 123456
    quietly bootstrap, reps(1000) cluster(iso3_id) idcluster(new_iso3): my_cf_boot_robust
    
    estadd scalar F_first = `f_stat'
    estadd local Region_FE "Yes"
    estadd local Year_FE "Yes"
    eststo G1_`var'
}

* --- Group 2: Winsorized ---
foreach var in $centralities {
    display "Group 2 (Winsorized IV): `var'"
    global depvar "total_adaptation_finance_w"
    global current_endog "`var'_ma3_w"
    global current_iv "`var'_iv_w" 
    global current_controls "$controls"
    global if_cond ""
    
    quietly reghdfe $current_endog $current_iv $current_controls, absorb(region_id year) vce(cluster iso3_id)
	quietly test $current_iv
	local f_stat = r(F)  
	
    capture drop new_iso3 
    set seed 123456
    quietly bootstrap, reps(1000) cluster(iso3_id) idcluster(new_iso3): my_cf_boot_robust
    
    estadd scalar F_first = `f_stat'
    estadd local Region_FE "Yes"
    estadd local Year_FE "Yes"
    eststo G2_`var'
}

* --- Group 3: Lagged Controls ---
foreach var in $centralities {
    display "Group 3 (Lagged Controls IV): `var'"
    global depvar "total_adaptation_finance"
    global current_endog "`var'_ma3"
    global current_iv "`var'_iv" 
    global current_controls "$lag_controls"  //
    global if_cond ""
    
    quietly reghdfe $current_endog $current_iv $current_controls, absorb(region_id year) vce(cluster iso3_id)
    quietly test $current_iv
	local f_stat = r(F)

    capture drop new_iso3 
    set seed 123456
    quietly bootstrap, reps(1000) cluster(iso3_id) idcluster(new_iso3): my_cf_boot_robust
    
    estadd scalar F_first = `f_stat'
    estadd local Region_FE "Yes"
    estadd local Year_FE "Yes"
    eststo G3_`var'
}

* --- Group 4: Excluding COVID-19 ---
foreach var in $centralities {
    display "Group 4 (Ex-Covid IV): `var'"
    global depvar "total_adaptation_finance"
    global current_endog "`var'_ma3"
    global current_iv "`var'_iv" 
    global current_controls "$controls"
    global if_cond "if year != 2020 & year != 2021" //
    
    quietly reghdfe $current_endog $current_iv $current_controls $if_cond, absorb(region_id year) vce(cluster iso3_id)
    quietly test $current_iv
	local f_stat = r(F)
    
    capture drop new_iso3 
    set seed 123456
    quietly bootstrap, reps(1000) cluster(iso3_id) idcluster(new_iso3): my_cf_boot_robust
    
    estadd scalar F_first = `f_stat'
    estadd local Region_FE "Yes"
    estadd local Year_FE "Yes"
    eststo G4_`var'
}

* ---------------------------------------------------------
esttab G1_* G2_* G3_* G4_* ///
    using "regression_results_IV_CF_robustness_all.csv", replace ///
    b(3) se(3) star(* 0.1 ** 0.05 *** 0.01) compress nogaps nodepvars ///
    mtitles("MA5" "MA5" "MA5" "MA5" ///
            "Winsor" "Winsor" "Winsor" "Winsor" ///
            "Lag-Ctrl" "Lag-Ctrl" "Lag-Ctrl" "Lag-Ctrl" ///
            "Ex-Covid" "Ex-Covid" "Ex-Covid" "Ex-Covid") ///
    rename(in_degree_ma5 in_degree_ma3 out_degree_ma5 out_degree_ma3 closeness_ma5 closeness_ma3 eigenvector_ma5 eigenvector_ma3 ///
           in_degree_ma3_w in_degree_ma3 out_degree_ma3_w out_degree_ma3 closeness_ma3_w closeness_ma3 eigenvector_ma3_w eigenvector_ma3) /// 
    keep(*_ma3 v_resid_boot) ///
    order(in_degree_ma3 out_degree_ma3 closeness_ma3 eigenvector_ma3 v_resid_boot) ///
    scalars( ///
        "F_first First-Stage F-stat" ///
        "N Observations" ///
        "Region_FE Region FE" ///
        "Year_FE Year FE" ///
    ) label



* =========================================================
* 6.2 IV Regression Robustness Tests: Bartik IV -- CF Approach
* =========================================================

global controls gdppc fdi gov_effectiveness disaster_affected recap agri_va rural_pop co2 total_erp_pct trade
global centralities in_degree out_degree closeness eigenvector

foreach var in $centralities {
    capture drop `var'_ma5
    quietly rangestat (mean) `var', interval(year -5 -1) by(iso3_id)
    rename `var'_mean `var'_ma5
    label var `var'_ma5 "MA5 of `: var label `var''"
}

capture xtset iso3_id year 
global lag_controls ""
foreach var in $controls {
    capture drop l_`var'
    gen l_`var' = L.`var'
    global lag_controls "$lag_controls l_`var'"
}

local vars_to_winsor "total_adaptation_finance ln_iv_seaborne_trade_total"
foreach var in $centralities {
    local vars_to_winsor "`vars_to_winsor' `var'_ma3"
}
capture drop *_w 
winsor2 `vars_to_winsor', cuts(1 99) suffix(_w)

* ---------------------------------------------------------
capture xtset, clear
capture tsset, clear
eststo clear

capture program drop my_cf_boot_bartik_rob
program my_cf_boot_bartik_rob, eclass
    version 14.0 
    capture drop v_resid_boot
    
    quietly reghdfe $current_endog $current_iv $current_controls $if_cond, absorb(region_id year) resid
    quietly predict v_resid_boot, resid
    
	quietly ppmlhdfe $depvar $current_endog $current_controls v_resid_boot $if_cond, absorb(region_id year) d(`ppml_d')
	quietly margins, dydx($current_endog v_resid_boot) post

end

* ---------------------------------------------------------

* --- Group 1: Alternative Measure (MA5) ---
foreach var in $centralities {
    display "Group 1 (MA5 | Bartik IV): `var'"
    global depvar "total_adaptation_finance"
    global current_endog "`var'_ma5"
    global current_iv "ln_iv_seaborne_trade_total" 
    global current_controls "$controls"
    global if_cond ""
    
    quietly reghdfe $current_endog $current_iv $current_controls, absorb(region_id year) vce(cluster iso3_id)
    quietly test $current_iv
	local f_stat = r(F)
    
    capture drop new_iso3 
    set seed 123456
    quietly bootstrap, reps(1000) cluster(iso3_id) idcluster(new_iso3): my_cf_boot_bartik_rob
    
    estadd scalar F_first = `f_stat'
    estadd local Region_FE "Yes"
    estadd local Year_FE "Yes"
    eststo G1_`var'
}

* --- Group 2: Winsorized ---
foreach var in $centralities {
    display "Group 2 (Winsorized | Bartik IV): `var'"
    global depvar "total_adaptation_finance_w"
    global current_endog "`var'_ma3_w"
    global current_iv "ln_iv_seaborne_trade_total_w"
    global current_controls "$controls"
    global if_cond ""
    
    quietly reghdfe $current_endog $current_iv $current_controls, absorb(region_id year) vce(cluster iso3_id)
    quietly test $current_iv
	local f_stat = r(F)
    
    capture drop new_iso3 
    set seed 123456
    quietly bootstrap, reps(1000) cluster(iso3_id) idcluster(new_iso3): my_cf_boot_bartik_rob
    
    estadd scalar F_first = `f_stat'
    estadd local Region_FE "Yes"
    estadd local Year_FE "Yes"
    eststo G2_`var'
}

* --- Group 3: Lagged Controls ---
foreach var in $centralities {
    display "Group 3 (Lagged Controls | Bartik IV): `var'"
    global depvar "total_adaptation_finance"
    global current_endog "`var'_ma3"
    global current_iv "iv_seaborne_trade_total"
    global current_controls "$lag_controls"
    global if_cond ""
    
    quietly reghdfe $current_endog $current_iv $current_controls, absorb(region_id year) vce(cluster iso3_id)
    quietly test $current_iv
	local f_stat = r(F)
    
    capture drop new_iso3 
    set seed 123456
    quietly bootstrap, reps(1000) cluster(iso3_id) idcluster(new_iso3): my_cf_boot_bartik_rob
    
    estadd scalar F_first = `f_stat'
    estadd local Region_FE "Yes"
    estadd local Year_FE "Yes"
    eststo G3_`var'
}

* --- Group 4: Excluding COVID-19 ---
foreach var in $centralities {
    display "Group 4 (Ex-Covid | Bartik IV): `var'"
    global depvar "total_adaptation_finance"
    global current_endog "`var'_ma3"
    global current_iv "iv_seaborne_trade_total"
    global current_controls "$controls"
    global if_cond "if year != 2020 & year != 2021" 
    
    quietly reghdfe $current_endog $current_iv $current_controls $if_cond, absorb(region_id year) vce(cluster iso3_id)
    quietly test $current_iv
	local f_stat = r(F)
    
    capture drop new_iso3 
    set seed 123456
    quietly bootstrap, reps(1000) cluster(iso3_id) idcluster(new_iso3): my_cf_boot_bartik_rob
    
    estadd scalar F_first = `f_stat'
    estadd local Region_FE "Yes"
    estadd local Year_FE "Yes"
    eststo G4_`var'
}

* ---------------------------------------------------------
esttab G1_* G2_* G3_* G4_* ///
    using "regression_results_BartikIV_CF_robustness_all.csv", replace ///
    b(3) se(3) star(* 0.1 ** 0.05 *** 0.01) compress nogaps nodepvars ///
    mtitles("MA5" "MA5" "MA5" "MA5" ///
            "Winsor" "Winsor" "Winsor" "Winsor" ///
            "Lag-Ctrl" "Lag-Ctrl" "Lag-Ctrl" "Lag-Ctrl" ///
            "Ex-Covid" "Ex-Covid" "Ex-Covid" "Ex-Covid" ) ///
    rename(in_degree_ma5 in_degree_ma3 out_degree_ma5 out_degree_ma3 closeness_ma5 closeness_ma3 eigenvector_ma5 eigenvector_ma3 ///
           in_degree_ma3_w in_degree_ma3 out_degree_ma3_w out_degree_ma3 closeness_ma3_w closeness_ma3 eigenvector_ma3_w eigenvector_ma3) /// 
    keep(*_ma3 v_resid_boot) ///
    order(in_degree_ma3 out_degree_ma3 closeness_ma3 eigenvector_ma3 v_resid_boot) ///
    scalars( ///
        "F_first First-Stage F-stat" ///
        "N Observations" ///
        "Region_FE Region FE" ///
        "Year_FE Year FE" ///
    ) label
