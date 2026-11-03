* =========================================================
* 06_moderation_disaster.do  Moderation by disaster exposure (gravity IV, Bartik IV)
* Origin : Code/Reg_model.do, sections 5.1 and 5.3 (lines 1205-1437, 1686-1919)
* Output : output/tables/
*            Moderation_disaster_Table1_Main_Regression_CF[_Bartik].csv
*            Moderation_disaster_Table2_Marginal_Effects_CF[_Bartik].csv
*            Moderation_disaster_Table3_Marginal_Effects_Differences_CF[_Bartik].csv
*          + *_pctCI.csv versions of Tables 2-3 (percentile bootstrap CIs);
*            replicate draws saved in output/bootstrap/
* Marginal effects at P10/P25/P50/P75/P90/P95/P99 of the moderator (P10/P25 added
* 2026-10-26; the main contrast in Fig. 6a-b is ME(P90) - ME(P10)).
* =========================================================

do "$ANALYSIS/00_prepare.do"
do "$ANALYSIS/00_bootstrap_pctci.do"

* =========================================================
* 5.1 IV Regression Moderation Effects: Gravity IV -- CF Approach -- disaster
* =========================================================

global controls gdppc fdi gov_effectiveness agri_va rural_pop co2 total_erp_pct trade
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
        at(${moderator}=($p10 $p25 $p50 $p75 $p90 $p95 $p99)) post
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
        at(${moderator}=($p10 $p25 $p50 $p75 $p90 $p95 $p99)) post

    return scalar me_p10 = _b[1._at]
    return scalar me_p25 = _b[2._at]
    return scalar me_p50 = _b[3._at]
    return scalar me_p75 = _b[4._at]
    return scalar me_p90 = _b[5._at]
    return scalar me_p95 = _b[6._at]
    return scalar me_p99 = _b[7._at]

    return scalar diff_p90_p10 = _b[5._at] - _b[1._at]
    return scalar diff_p75_p25 = _b[4._at] - _b[2._at]
    return scalar diff_p75_p50 = _b[4._at] - _b[3._at]
    return scalar diff_p90_p50 = _b[5._at] - _b[3._at]
    return scalar diff_p95_p50 = _b[6._at] - _b[3._at]
    return scalar diff_p99_p50 = _b[7._at] - _b[3._at]
    return scalar diff_p90_p75 = _b[5._at] - _b[4._at]
    return scalar diff_p95_p90 = _b[6._at] - _b[5._at]
	return scalar diff_p95_p75 = _b[6._at] - _b[4._at]
    return scalar diff_p99_p95 = _b[7._at] - _b[6._at]
end

local centralities in_degree out_degree closeness eigenvector

* Percentile bootstrap CIs for Tables 2-3 (see 00_bootstrap_pctci.do)
tempfile pct_me pct_diff
pctci_open pf_me "`pct_me'"
pctci_open pf_diff "`pct_diff'"

foreach var in `centralities' {

    global current_endog "`var'_ma3"
    global current_iv "`var'_iv_ma3"

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
    global p95 = r(p95)
    global p99 = r(p99)

    display "P10 = $p10"
    display "P25 = $p25"
    display "P50 = $p50"
    display "P75 = $p75"
    display "P90 = $p90"
    display "P95 = $p95"
    display "P99 = $p99"

    capture drop new_iso3_reg
    set seed 123456

    quietly bootstrap, reps($BOOT_REPS) cluster(iso3_id) idcluster(new_iso3_reg): ///
        my_cf_reg_boot if estimation_sample == 1

    estadd local Controls "Yes"
    estadd local FE "Yes"
    eststo reg_`var'

    capture drop new_iso3_mar
    set seed 123456

    quietly bootstrap, reps($BOOT_REPS) cluster(iso3_id) idcluster(new_iso3_mar) ///
        saving("$BOOTDRAWS/mod06_gravity_margins_`var'.dta", replace): ///
        my_cf_margins_boot if estimation_sample == 1

    eststo margins_`var'
    pctci_post pf_me "$BOOTDRAWS/mod06_gravity_margins_`var'.dta" `var'

    capture drop new_iso3_diff
    set seed 123456

    quietly bootstrap ///
        me_p10=r(me_p10) ///
        me_p25=r(me_p25) ///
        me_p50=r(me_p50) ///
        me_p75=r(me_p75) ///
        me_p90=r(me_p90) ///
        me_p95=r(me_p95) ///
        me_p99=r(me_p99) ///
        diff_p90_p10=r(diff_p90_p10) ///
        diff_p75_p25=r(diff_p75_p25) ///
        diff_p75_p50=r(diff_p75_p50) ///
        diff_p90_p50=r(diff_p90_p50) ///
        diff_p95_p50=r(diff_p95_p50) ///
        diff_p99_p50=r(diff_p99_p50) ///
        diff_p90_p75=r(diff_p90_p75) ///
        diff_p95_p75=r(diff_p95_p75) ///
        diff_p95_p90=r(diff_p95_p90) ///
        diff_p99_p95=r(diff_p99_p95), ///
        reps($BOOT_REPS) cluster(iso3_id) idcluster(new_iso3_diff) ///
        saving("$BOOTDRAWS/mod06_gravity_mediff_`var'.dta", replace): ///
        my_cf_me_diff_boot if estimation_sample == 1

    eststo mediff_`var'
    pctci_post pf_diff "$BOOTDRAWS/mod06_gravity_mediff_`var'.dta" `var'
}


esttab reg_in_degree reg_out_degree reg_closeness reg_eigenvector ///
    using "$TABLES/Moderation_disaster_Table1_Main_Regression_CF.csv", replace ///
    keep(*_ma3 ${moderator} *#* v_resid1_boot v_resid2_boot) ///
    b(3) se(3) star(* 0.1 ** 0.05 *** 0.01) ///
    mtitles("In-Degree" "Out-Degree" "Closeness" "Eigenvector") ///
    scalars("Controls Controls" "FE Region & Year FE") ///
    compress nogaps label ///
    title("Second-Stage PPML Regression Results")


esttab margins_in_degree margins_out_degree margins_closeness margins_eigenvector ///
    using "$TABLES/Moderation_disaster_Table2_Marginal_Effects_CF.csv", replace ///
    b(3) se(3) star(* 0.1 ** 0.05 *** 0.01) ///
    mtitles("In-Degree" "Out-Degree" "Closeness" "Eigenvector") ///
    coeflabels( ///
        1._at "ME at P10" ///
        2._at "ME at P25" ///
        3._at "ME at P50" ///
        4._at "ME at P75" ///
        5._at "ME at P90" ///
        6._at "ME at P95" ///
        7._at "ME at P99" ///
    ) ///
    compress nogaps label ///
    title("Marginal Effects at Different Percentiles")

esttab mediff_in_degree mediff_out_degree mediff_closeness mediff_eigenvector ///
    using "$TABLES/Moderation_disaster_Table3_Marginal_Effects_Differences_CF.csv", replace ///
    keep(diff_p90_p10 diff_p75_p25 diff_p75_p50 diff_p90_p50 diff_p95_p50 diff_p99_p50 ///
         diff_p90_p75 diff_p95_p75 diff_p95_p90 diff_p99_p95) ///
    b(3) se(3) star(* 0.1 ** 0.05 *** 0.01) ///
    mtitles("In-Degree" "Out-Degree" "Closeness" "Eigenvector") ///
    coeflabels( ///
        diff_p90_p10 "ME(P90) - ME(P10)" ///
        diff_p75_p25 "ME(P75) - ME(P25)" ///
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

pctci_export pf_me "`pct_me'" "$TABLES/Moderation_disaster_Table2_Marginal_Effects_CF_pctCI.csv" "P10 P25 P50 P75 P90 P95 P99"
pctci_export pf_diff "`pct_diff'" "$TABLES/Moderation_disaster_Table3_Marginal_Effects_Differences_CF_pctCI.csv" "P10 P25 P50 P75 P90 P95 P99"


* =========================================================
* 5.3 IV Regression Moderation Effects: Bartik IV -- CF Approach -- disaster
* =========================================================

global controls gdppc fdi gov_effectiveness agri_va rural_pop co2 total_erp_pct trade
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
        at(${moderator}=($p10 $p25 $p50 $p75 $p90 $p95 $p99)) post
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
        at(${moderator}=($p10 $p25 $p50 $p75 $p90 $p95 $p99)) post

    return scalar me_p10 = _b[1._at]
    return scalar me_p25 = _b[2._at]
    return scalar me_p50 = _b[3._at]
    return scalar me_p75 = _b[4._at]
    return scalar me_p90 = _b[5._at]
    return scalar me_p95 = _b[6._at]
    return scalar me_p99 = _b[7._at]

    return scalar diff_p90_p10 = _b[5._at] - _b[1._at]
    return scalar diff_p75_p25 = _b[4._at] - _b[2._at]
    return scalar diff_p75_p50 = _b[4._at] - _b[3._at]
    return scalar diff_p90_p50 = _b[5._at] - _b[3._at]
    return scalar diff_p95_p50 = _b[6._at] - _b[3._at]
    return scalar diff_p99_p50 = _b[7._at] - _b[3._at]
    return scalar diff_p90_p75 = _b[5._at] - _b[4._at]
	return scalar diff_p95_p75 = _b[6._at] - _b[4._at]
    return scalar diff_p95_p90 = _b[6._at] - _b[5._at]
    return scalar diff_p99_p95 = _b[7._at] - _b[6._at]
end

local centralities in_degree out_degree closeness eigenvector

* Percentile bootstrap CIs for Tables 2-3 (see 00_bootstrap_pctci.do)
tempfile pct_me pct_diff
pctci_open pf_me "`pct_me'"
pctci_open pf_diff "`pct_diff'"

foreach var in `centralities' {

    global current_endog "`var'_ma3"
    global current_iv "ln_iv_seaborne_trade_total_ma3"

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
    global p95 = r(p95)
    global p99 = r(p99)

    display "P10 = $p10"
    display "P25 = $p25"
    display "P50 = $p50"
    display "P75 = $p75"
    display "P90 = $p90"
    display "P95 = $p95"
    display "P99 = $p99"

    capture drop new_iso3_reg
    set seed 123456

    quietly bootstrap, reps($BOOT_REPS) cluster(iso3_id) idcluster(new_iso3_reg): ///
        my_cf_reg_boot if estimation_sample == 1

    estadd local Controls "Yes"
    estadd local FE "Yes"
    eststo reg_`var'

    capture drop new_iso3_mar
    set seed 123456

    quietly bootstrap, reps($BOOT_REPS) cluster(iso3_id) idcluster(new_iso3_mar) ///
        saving("$BOOTDRAWS/mod06_bartik_margins_`var'.dta", replace): ///
        my_cf_margins_boot if estimation_sample == 1

    eststo margins_`var'
    pctci_post pf_me "$BOOTDRAWS/mod06_bartik_margins_`var'.dta" `var'

    capture drop new_iso3_diff
    set seed 123456

    quietly bootstrap ///
        me_p10=r(me_p10) ///
        me_p25=r(me_p25) ///
        me_p50=r(me_p50) ///
        me_p75=r(me_p75) ///
        me_p90=r(me_p90) ///
        me_p95=r(me_p95) ///
        me_p99=r(me_p99) ///
        diff_p90_p10=r(diff_p90_p10) ///
        diff_p75_p25=r(diff_p75_p25) ///
        diff_p75_p50=r(diff_p75_p50) ///
        diff_p90_p50=r(diff_p90_p50) ///
        diff_p95_p50=r(diff_p95_p50) ///
        diff_p99_p50=r(diff_p99_p50) ///
        diff_p90_p75=r(diff_p90_p75) ///
		diff_p95_p75=r(diff_p95_p75) ///
        diff_p95_p90=r(diff_p95_p90) ///
        diff_p99_p95=r(diff_p99_p95), ///
        reps($BOOT_REPS) cluster(iso3_id) idcluster(new_iso3_diff) ///
        saving("$BOOTDRAWS/mod06_bartik_mediff_`var'.dta", replace): ///
        my_cf_me_diff_boot if estimation_sample == 1

    eststo mediff_`var'
    pctci_post pf_diff "$BOOTDRAWS/mod06_bartik_mediff_`var'.dta" `var'

}


esttab reg_in_degree reg_out_degree reg_closeness reg_eigenvector ///
    using "$TABLES/Moderation_disaster_Table1_Main_Regression_CF_Bartik.csv", replace ///
    keep(*_ma3 ${moderator} *#* v_resid1_boot v_resid2_boot) ///
    b(3) se(3) star(* 0.1 ** 0.05 *** 0.01) ///
    mtitles("In-Degree" "Out-Degree" "Closeness" "Eigenvector") ///
    scalars("Controls Controls" "FE Region & Year FE") ///
    compress nogaps label ///
    title("Second-Stage PPML Regression Results")


esttab margins_in_degree margins_out_degree margins_closeness margins_eigenvector ///
    using "$TABLES/Moderation_disaster_Table2_Marginal_Effects_CF_Bartik.csv", replace ///
    b(3) se(3) star(* 0.1 ** 0.05 *** 0.01) ///
    mtitles("In-Degree" "Out-Degree" "Closeness" "Eigenvector") ///
    coeflabels( ///
        1._at "ME at P10" ///
        2._at "ME at P25" ///
        3._at "ME at P50" ///
        4._at "ME at P75" ///
        5._at "ME at P90" ///
        6._at "ME at P95" ///
        7._at "ME at P99" ///
    ) ///
    compress nogaps label ///
    title("Marginal Effects at Different Percentiles")

esttab mediff_in_degree mediff_out_degree mediff_closeness mediff_eigenvector ///
    using "$TABLES/Moderation_disaster_Table3_Marginal_Effects_Differences_CF_Bartik.csv", replace ///
    keep(diff_p90_p10 diff_p75_p25 diff_p75_p50 diff_p90_p50 diff_p95_p50 diff_p99_p50 ///
         diff_p90_p75 diff_p95_p75 diff_p95_p90 diff_p99_p95) ///
    b(3) se(3) star(* 0.1 ** 0.05 *** 0.01) ///
    mtitles("In-Degree" "Out-Degree" "Closeness" "Eigenvector") ///
    coeflabels( ///
        diff_p90_p10 "ME(P90) - ME(P10)" ///
        diff_p75_p25 "ME(P75) - ME(P25)" ///
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

pctci_export pf_me "`pct_me'" "$TABLES/Moderation_disaster_Table2_Marginal_Effects_CF_Bartik_pctCI.csv" "P10 P25 P50 P75 P90 P95 P99"
pctci_export pf_diff "`pct_diff'" "$TABLES/Moderation_disaster_Table3_Marginal_Effects_Differences_CF_Bartik_pctCI.csv" "P10 P25 P50 P75 P90 P95 P99"
