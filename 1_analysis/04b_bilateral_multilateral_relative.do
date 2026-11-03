* =========================================================
* 04b_bilateral_multilateral_relative.do  Relative sensitivity: bilateral vs multilateral
* Output : output/tables/
*            regression_results_Bi_Multi_relative_GravityIV_CF.csv
*            regression_results_Bi_Multi_relative_BartikIV_CF.csv
*
* Module 04 reports average marginal effects in USD million, which mostly reflect
* the size of each channel (multilateral ~86% of adaptation finance). Here we
* compare the PPML coefficients instead. With fixed effects, PPML fitted values
* average to the sample mean of y, so AME / mean(y) = beta: beta is the change in
* finance per 1 SD of centrality relative to the channel's mean, and
* 100*(exp(beta)-1) is the percentage change in expected finance.
*
* Both channels are estimated inside the same bootstrap replicate (same first
* stage, same resampled countries), so beta_diff = b_bi - b_multi gets a valid
* bootstrap SE. Inference: country-clustered bootstrap SEs, as in module 04.
* =========================================================

do "$ANALYSIS/00_prepare.do"

global controls gdppc fdi gov_effectiveness disaster_total_impact_pc agri_va rural_pop co2 total_erp_pct trade
capture xtset, clear
capture tsset, clear

capture program drop my_cf_rel_boot
program my_cf_rel_boot, rclass
    version 14.0
    capture drop v_resid_boot

    quietly reghdfe $current_endog $current_iv $controls, absorb(region_id year) resid
    quietly predict v_resid_boot, resid

    quietly ppmlhdfe bilateral_adaptation $current_endog $controls v_resid_boot, absorb(region_id year)
    local b_bi = _b[$current_endog]
    quietly ppmlhdfe multilateral_adaptation $current_endog $controls v_resid_boot, absorb(region_id year)
    local b_multi = _b[$current_endog]

    return scalar b_bi      = `b_bi'
    return scalar b_multi   = `b_multi'
    return scalar beta_diff = `b_bi' - `b_multi'
end

local centralities "in_degree out_degree closeness eigenvector"

foreach iv in gravity bartik {
    eststo clear

    foreach var in `centralities' {
        global current_endog "`var'_ma3"
        if "`iv'" == "gravity" global current_iv "`var'_iv_ma3"
        else                   global current_iv ln_iv_seaborne_trade_total_ma3

        quietly reghdfe $current_endog $current_iv $controls, absorb(region_id year) vce(cluster iso3_id)
        quietly test $current_iv
        local f_stat = r(F)

        capture drop new_iso3
        set seed 20261022
        quietly bootstrap b_bi=r(b_bi) b_multi=r(b_multi) beta_diff=r(beta_diff), ///
            reps($BOOT_REPS) cluster(iso3_id) idcluster(new_iso3): my_cf_rel_boot

        display as text "`iv' IV, `var': b_bi = " %6.3f _b[b_bi] ", b_multi = " %6.3f _b[b_multi] ///
            ", diff = " %6.3f _b[beta_diff] " (SE " %6.3f _se[beta_diff] ")"

        estadd scalar F_first = `f_stat'
        estadd local Controls "Yes"
        estadd local FE "Yes"
        eststo r_`var'
    }

    local ivname = cond("`iv'" == "gravity", "GravityIV", "BartikIV")
    esttab r_in_degree r_out_degree r_closeness r_eigenvector ///
        using "$TABLES/regression_results_Bi_Multi_relative_`ivname'_CF.csv", replace ///
        coeflabels(b_bi "Bilateral (PPML beta)" b_multi "Multilateral (PPML beta)" ///
                   beta_diff "Difference (bilateral - multilateral)") ///
        b(3) se(3) star(* 0.1 ** 0.05 *** 0.01) ///
        mtitles("In-Degree" "Out-Degree" "Closeness" "Eigenvector") ///
        title("PPML coefficients with CF approach, `ivname' (beta = AME / mean of dependent variable)") ///
        scalars( ///
            "F_first First-Stage F-stat" ///
            "Controls Controls" ///
            "FE Region & Year FE" ///
            "N Observations" ///
        ) ///
        compress nogaps
}
