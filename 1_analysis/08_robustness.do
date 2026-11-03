* =========================================================
* 08_robustness.do  Robustness: MA5, winsorised, lagged controls, excluding COVID years,
*                    contemporaneous (year-t) IV
* Origin : Code/Reg_model.do, sections 6.1-6.2 (lines 2168-2490)
* Output : output/tables/
*            regression_results_IV_CF_robustness_all.csv
*            regression_results_BartikIV_CF_robustness_all.csv
* =========================================================

do "$ANALYSIS/00_prepare.do"

* =========================================================
* 6.1 IV Regression Robustness Tests -- Control Function Approach (CF)
* =========================================================

global controls gdppc fdi gov_effectiveness disaster_total_impact_pc agri_va rural_pop co2 total_erp_pct trade
global centralities in_degree out_degree closeness eigenvector

* *_ma5 and *_iv_ma5 (t-5..t-1) are built in 00_prepare.do

capture xtset iso3_id year
global lag_controls ""
foreach var in $controls {
    capture drop l_`var'
    gen l_`var' = L.`var'
    global lag_controls "$lag_controls l_`var'"
}

local vars_to_winsor "total_adaptation_finance"
foreach var in $centralities {
    local vars_to_winsor "`vars_to_winsor' `var'_ma3 `var'_iv_ma3"
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
	tempvar ppml_d

    quietly reghdfe $current_endog $current_iv $current_controls $if_cond, absorb(region_id year) resid
    quietly predict v_resid_boot if e(sample), resid

    quietly ppmlhdfe $depvar $current_endog $current_controls v_resid_boot $if_cond, absorb(region_id year) d(`ppml_d')
	quietly margins, dydx($current_endog) predict(mu) post
end

* ---------------------------------------------------------

* --- Group 1: Alternative Measure (MA5) ---
foreach var in $centralities {
    display "Group 1 (MA5 IV): `var'"
    global depvar "total_adaptation_finance"
    global current_endog "`var'_ma5"
    global current_iv "`var'_iv_ma5"
    global current_controls "$controls"
    global if_cond ""

    quietly reghdfe $current_endog $current_iv $current_controls, absorb(region_id year) vce(cluster iso3_id)
    quietly test $current_iv
	local f_stat = r(F)

    capture drop new_iso3
    set seed 123456
    quietly bootstrap, reps($BOOT_REPS) cluster(iso3_id) idcluster(new_iso3): my_cf_boot_robust

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
    global current_iv "`var'_iv_ma3_w"
    global current_controls "$controls"
    global if_cond ""

    quietly reghdfe $current_endog $current_iv $current_controls, absorb(region_id year) vce(cluster iso3_id)
	quietly test $current_iv
	local f_stat = r(F)

    capture drop new_iso3
    set seed 123456
    quietly bootstrap, reps($BOOT_REPS) cluster(iso3_id) idcluster(new_iso3): my_cf_boot_robust

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
    global current_iv "`var'_iv_ma3"
    global current_controls "$lag_controls"  //
    global if_cond ""

    quietly reghdfe $current_endog $current_iv $current_controls, absorb(region_id year) vce(cluster iso3_id)
    quietly test $current_iv
	local f_stat = r(F)

    capture drop new_iso3
    set seed 123456
    quietly bootstrap, reps($BOOT_REPS) cluster(iso3_id) idcluster(new_iso3): my_cf_boot_robust

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
    global current_iv "`var'_iv_ma3"
    global current_controls "$controls"
    global if_cond "if year != 2020 & year != 2021" //

    quietly reghdfe $current_endog $current_iv $current_controls $if_cond, absorb(region_id year) vce(cluster iso3_id)
    quietly test $current_iv
	local f_stat = r(F)

    capture drop new_iso3
    set seed 123456
    quietly bootstrap, reps($BOOT_REPS) cluster(iso3_id) idcluster(new_iso3): my_cf_boot_robust

    estadd scalar F_first = `f_stat'
    estadd local Region_FE "Yes"
    estadd local Year_FE "Yes"
    eststo G4_`var'
}

* --- Group 5: Contemporaneous (year-t) IV ---
foreach var in $centralities {
    display "Group 5 (Year-t IV): `var'"
    global depvar "total_adaptation_finance"
    global current_endog "`var'_ma3"
    global current_iv "`var'_iv"
    global current_controls "$controls"
    global if_cond ""

    quietly reghdfe $current_endog $current_iv $current_controls, absorb(region_id year) vce(cluster iso3_id)
    quietly test $current_iv
	local f_stat = r(F)

    capture drop new_iso3
    set seed 123456
    quietly bootstrap, reps($BOOT_REPS) cluster(iso3_id) idcluster(new_iso3): my_cf_boot_robust

    estadd scalar F_first = `f_stat'
    estadd local Region_FE "Yes"
    estadd local Year_FE "Yes"
    eststo G5_`var'
}

* ---------------------------------------------------------
esttab G1_* G2_* G3_* G4_* G5_* ///
    using "$TABLES/regression_results_IV_CF_robustness_all.csv", replace ///
    b(3) se(3) star(* 0.1 ** 0.05 *** 0.01) compress nogaps nodepvars ///
    mtitles("MA5" "MA5" "MA5" "MA5" ///
            "Winsor" "Winsor" "Winsor" "Winsor" ///
            "Lag-Ctrl" "Lag-Ctrl" "Lag-Ctrl" "Lag-Ctrl" ///
            "Ex-Covid" "Ex-Covid" "Ex-Covid" "Ex-Covid" ///
            "IV-t" "IV-t" "IV-t" "IV-t") ///
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

global controls gdppc fdi gov_effectiveness disaster_total_impact_pc agri_va rural_pop co2 total_erp_pct trade
global centralities in_degree out_degree closeness eigenvector

* *_ma5 and *_iv_ma5 (t-5..t-1) are built in 00_prepare.do

capture xtset iso3_id year
global lag_controls ""
foreach var in $controls {
    capture drop l_`var'
    gen l_`var' = L.`var'
    global lag_controls "$lag_controls l_`var'"
}

local vars_to_winsor "total_adaptation_finance ln_iv_seaborne_trade_total_ma3"
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
    tempvar ppml_d

    quietly reghdfe $current_endog $current_iv $current_controls $if_cond, absorb(region_id year) resid
    quietly predict v_resid_boot if e(sample), resid

    quietly ppmlhdfe $depvar $current_endog $current_controls v_resid_boot $if_cond, absorb(region_id year) d(`ppml_d')
	quietly margins, dydx($current_endog) predict(mu) post

end

* ---------------------------------------------------------

* --- Group 1: Alternative Measure (MA5) ---
foreach var in $centralities {
    display "Group 1 (MA5 | Bartik IV): `var'"
    global depvar "total_adaptation_finance"
    global current_endog "`var'_ma5"
    global current_iv "ln_iv_seaborne_trade_total_ma5"
    global current_controls "$controls"
    global if_cond ""

    quietly reghdfe $current_endog $current_iv $current_controls, absorb(region_id year) vce(cluster iso3_id)
    quietly test $current_iv
	local f_stat = r(F)

    capture drop new_iso3
    set seed 123456
    quietly bootstrap, reps($BOOT_REPS) cluster(iso3_id) idcluster(new_iso3): my_cf_boot_bartik_rob

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
    global current_iv "ln_iv_seaborne_trade_total_ma3_w"
    global current_controls "$controls"
    global if_cond ""

    quietly reghdfe $current_endog $current_iv $current_controls, absorb(region_id year) vce(cluster iso3_id)
    quietly test $current_iv
	local f_stat = r(F)

    capture drop new_iso3
    set seed 123456
    quietly bootstrap, reps($BOOT_REPS) cluster(iso3_id) idcluster(new_iso3): my_cf_boot_bartik_rob

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
    global current_iv "ln_iv_seaborne_trade_total_ma3"
    global current_controls "$lag_controls"
    global if_cond ""

    quietly reghdfe $current_endog $current_iv $current_controls, absorb(region_id year) vce(cluster iso3_id)
    quietly test $current_iv
	local f_stat = r(F)

    capture drop new_iso3
    set seed 123456
    quietly bootstrap, reps($BOOT_REPS) cluster(iso3_id) idcluster(new_iso3): my_cf_boot_bartik_rob

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
    global current_iv "ln_iv_seaborne_trade_total_ma3"
    global current_controls "$controls"
    global if_cond "if year != 2020 & year != 2021"

    quietly reghdfe $current_endog $current_iv $current_controls $if_cond, absorb(region_id year) vce(cluster iso3_id)
    quietly test $current_iv
	local f_stat = r(F)

    capture drop new_iso3
    set seed 123456
    quietly bootstrap, reps($BOOT_REPS) cluster(iso3_id) idcluster(new_iso3): my_cf_boot_bartik_rob

    estadd scalar F_first = `f_stat'
    estadd local Region_FE "Yes"
    estadd local Year_FE "Yes"
    eststo G4_`var'
}

* --- Group 5: Contemporaneous (year-t) IV ---
foreach var in $centralities {
    display "Group 5 (Year-t | Bartik IV): `var'"
    global depvar "total_adaptation_finance"
    global current_endog "`var'_ma3"
    global current_iv "ln_iv_seaborne_trade_total"
    global current_controls "$controls"
    global if_cond ""

    quietly reghdfe $current_endog $current_iv $current_controls, absorb(region_id year) vce(cluster iso3_id)
    quietly test $current_iv
	local f_stat = r(F)

    capture drop new_iso3
    set seed 123456
    quietly bootstrap, reps($BOOT_REPS) cluster(iso3_id) idcluster(new_iso3): my_cf_boot_bartik_rob

    estadd scalar F_first = `f_stat'
    estadd local Region_FE "Yes"
    estadd local Year_FE "Yes"
    eststo G5_`var'
}

* ---------------------------------------------------------
esttab G1_* G2_* G3_* G4_* G5_* ///
    using "$TABLES/regression_results_BartikIV_CF_robustness_all.csv", replace ///
    b(3) se(3) star(* 0.1 ** 0.05 *** 0.01) compress nogaps nodepvars ///
    mtitles("MA5" "MA5" "MA5" "MA5" ///
            "Winsor" "Winsor" "Winsor" "Winsor" ///
            "Lag-Ctrl" "Lag-Ctrl" "Lag-Ctrl" "Lag-Ctrl" ///
            "Ex-Covid" "Ex-Covid" "Ex-Covid" "Ex-Covid" ///
            "IV-t" "IV-t" "IV-t" "IV-t") ///
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
