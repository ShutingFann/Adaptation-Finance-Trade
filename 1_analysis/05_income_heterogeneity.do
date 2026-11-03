* =========================================================
* 05_income_heterogeneity.do  Heterogeneity by World Bank income group (gravity IV, Bartik IV)
* Origin : Code/Reg_model.do, sections 4.1-4.2 (lines 1044-1202)
* Output : output/tables/
*            regression_results_GravityIV_CF_Bootstrap_by_income.csv
*            regression_results_BartikIV_CF_Bootstrap_by_income.csv
*          output/bootstrap/mod05_{gravity,bartik}_<income group>_<centrality>.dta
*
* Bootstrap replicates are saved for diagnostics only: the income-group samples are
* small and a few replicates can diverge (e.g. Bartik IV, upper middle income,
* out-degree: plain bootstrap SE 1.6e17 in the 2026-10-24 run, before the FDI
* transform). Tables and Fig. 5c-d use the plain bootstrap SE.
* =========================================================

do "$ANALYSIS/00_prepare.do"

if "$BOOTDRAWS" == "" global BOOTDRAWS "$WORKFLOW/output/bootstrap"
capture mkdir "$BOOTDRAWS"

* =========================================================
* 4.1 IV Regression: Gravity IV -- CF Approach -- by Income Group
* =========================================================

global controls gdppc fdi gov_effectiveness disaster_total_impact_pc agri_va rural_pop co2 total_erp_pct trade
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
        global current_iv "`var'_iv_ma3"

        quietly reghdfe $current_endog $current_iv $controls, absorb(region_id year) vce(cluster iso3_id)
        quietly test $current_iv
        local f_stat = r(F)

        capture drop new_iso3
        set seed 123
        quietly bootstrap, reps($BOOT_REPS) cluster(iso3_id) idcluster(new_iso3) ///
            saving("$BOOTDRAWS/mod05_gravity_`clean_g'_`var'.dta", replace): my_cf_boot_group

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
    using "$TABLES/regression_results_GravityIV_CF_Bootstrap_by_income.csv", ///
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

global controls gdppc fdi gov_effectiveness disaster_total_impact_pc agri_va rural_pop co2 total_erp_pct trade
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
        global current_iv "ln_iv_seaborne_trade_total_ma3"

        quietly reghdfe $current_endog $current_iv $controls, absorb(region_id year) vce(cluster iso3_id)
        quietly test $current_iv
        local f_stat = r(F)

        capture drop new_iso3
        set seed 123
        quietly bootstrap, reps($BOOT_REPS) cluster(iso3_id) idcluster(new_iso3) ///
            saving("$BOOTDRAWS/mod05_bartik_`clean_g'_`var'.dta", replace): my_cf_boot_bartik

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
    using "$TABLES/regression_results_BartikIV_CF_Bootstrap_by_income.csv", ///
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
