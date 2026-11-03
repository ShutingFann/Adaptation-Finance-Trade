* =========================================================
* 04_bilateral_multilateral.do  Bilateral vs multilateral adaptation finance (gravity IV, Bartik IV)
* Origin : Code/Reg_model.do, sections 3.1-3.2 (lines 884-1040)
* Output : output/tables/
*            regression_results_Bi_Multi_GravityIV_CF.csv
*            regression_results_Bi_Multi_BartikIV_CF.csv
* =========================================================

do "$ANALYSIS/00_prepare.do"

* =========================================================
* 3.1 IV Regression: Gravity IV -- CF Approach -- Bilateral & Multilateral Adaptation
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
        global current_iv "`var'_iv_ma3"

        quietly reghdfe $current_endog $current_iv $controls, absorb(region_id year) vce(cluster iso3_id)
        quietly test $current_iv
		local f_stat = r(F)

        capture drop new_iso3
        set seed 123456
        quietly bootstrap, reps($BOOT_REPS) cluster(iso3_id) idcluster(new_iso3): my_cf_boot

        estadd scalar F_first = `f_stat'
        estadd local Controls "Yes"
        estadd local FE "Yes"
        estadd local DV_Type "`dv_label'"

        eststo `prefix'_`var'
    }
}

esttab bi_in_degree bi_out_degree bi_closeness bi_eigenvector ///
       multi_in_degree multi_out_degree multi_closeness multi_eigenvector ///
    using "$TABLES/regression_results_Bi_Multi_GravityIV_CF.csv", replace ///
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
global controls gdppc fdi gov_effectiveness disaster_total_impact_pc agri_va rural_pop co2 total_erp_pct trade
global current_endog closeness_ma3
global current_iv ln_iv_seaborne_trade_total_ma3

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
        global current_iv ln_iv_seaborne_trade_total_ma3

        quietly reghdfe $current_endog $current_iv $controls, absorb(region_id year) vce(cluster iso3_id)
        quietly test $current_iv
		local f_stat = r(F)

        capture drop new_iso3
        set seed 1234
        quietly bootstrap, reps($BOOT_REPS) cluster(iso3_id) idcluster(new_iso3): my_cf_boot

        estadd scalar F_first = `f_stat'
        estadd local Controls "Yes"
        estadd local FE "Yes"
        estadd local DV_Type "`dv_label'"

        eststo `prefix'_`var'
    }
}

esttab bi_in_degree bi_out_degree bi_closeness bi_eigenvector ///
       multi_in_degree multi_out_degree multi_closeness multi_eigenvector ///
    using "$TABLES/regression_results_Bi_Multi_BartikIV_CF.csv", replace ///
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
* Means of the dependent variables (for relative effects in Fig. 5a-b)
* =========================================================
do "$ANALYSIS/04_dv_means.do"
