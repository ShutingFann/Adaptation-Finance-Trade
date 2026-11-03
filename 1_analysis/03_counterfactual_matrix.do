* =========================================================
* 03_counterfactual_matrix.do  Counterfactual 4x4 quartile-entry-threshold matrices (continuous CF-PPML)
* Origin : Code/Reg_model.do, section 2 (lines 243-881)
* Output : output/tables/
*            <centrality>_continuous_entry_threshold_4x4_matrix.csv  (x4)
* =========================================================

do "$ANALYSIS/00_prepare.do"

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
    agri_va rural_pop co2 total_erp_pct trade

capture xtset, clear
capture tsset, clear
eststo clear

local centralities closeness eigenvector in_degree out_degree
// local centralities closeness eigenvector
local reps $BOOT_REPS
global mc_reps $MC_REPS

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
        `var'_iv_ma3, ///
        gdppc, ///
        fdi, ///
        gov_effectiveness, ///
        disaster_total_impact_pc, ///
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
        `var'_iv_ma3 ///
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
    global current_iv "`var'_iv_ma3"
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
        using "$TABLES/`var'_continuous_entry_threshold_4x4_matrix.csv", ///
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
