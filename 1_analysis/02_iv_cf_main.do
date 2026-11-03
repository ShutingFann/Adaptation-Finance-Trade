* =========================================================
* 02_iv_cf_main.do  Main IV results: control-function PPML with gravity IV and Bartik IV
* Origin : Code/Reg_model.do, sections 1.1-1.2 (lines 131-241)
* Output : output/tables/
*            regression_results_5_centrality_GravityIV_CF_Bootstrap_Adaptation.csv
*            regression_results_5_centrality_BartikIV_CF_Bootstrap_Adaptation.csv
* =========================================================

do "$ANALYSIS/00_prepare.do"

* =========================================================
* 1.1 IV Regression: Gravity IV -- Control Function Approach (CF) -- Adaptation Financing
* =========================================================

global controls gdppc fdi gov_effectiveness disaster_total_impact_pc agri_va rural_pop co2 total_erp_pct trade
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
    global current_iv "`var'_iv_ma3"

    quietly reghdfe $current_endog $current_iv $controls, absorb(region_id year) vce(cluster iso3_id)
    test $current_iv
    local f_stat = r(F)

    capture drop new_iso3
    set seed 123456
    quietly bootstrap, reps($BOOT_REPS) cluster(iso3_id) idcluster(new_iso3): my_cf_boot

    estadd scalar F_first = `f_stat'
    estadd local Controls "Yes"
    estadd local FE "Yes"
    eststo iv_`var'
}

esttab iv_in_degree iv_out_degree iv_closeness iv_eigenvector ///
    using "$TABLES/regression_results_5_centrality_GravityIV_CF_Bootstrap_Adaptation.csv", replace ///
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

global controls gdppc fdi gov_effectiveness disaster_total_impact_pc agri_va rural_pop co2 total_erp_pct trade
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
    global current_iv ln_iv_seaborne_trade_total_ma3

    reghdfe $current_endog $current_iv $controls, absorb(region_id year) vce(cluster iso3_id)
    quietly test $current_iv
    local f_stat = r(F)

    capture drop new_iso3
    set seed 123
    bootstrap, reps($BOOT_REPS) cluster(iso3_id) idcluster(new_iso3): my_cf_boot

    estadd scalar F_first = `f_stat'
    estadd local Controls "Yes"
    estadd local FE "Yes"
    eststo iv_`var'
}

esttab iv_in_degree iv_out_degree iv_closeness iv_eigenvector ///
    using "$TABLES/regression_results_5_centrality_BartikIV_CF_Bootstrap_Adaptation.csv", replace ///
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
