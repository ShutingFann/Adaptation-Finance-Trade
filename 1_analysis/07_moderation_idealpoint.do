* =========================================================
* 07_moderation_idealpoint.do  Moderation by UNGA political alignment (gravity IV, Bartik IV)
* Origin : Code/Reg_model.do, sections 5.2 and 5.4 (lines 1441-1683, 1922-2165)
* Output : output/tables/
*            Moderation_idealpoint_Table1_Main_Regression_CF[_Bartik].csv
*            Moderation_idealpoint_Table2_Main_Regression_CF[_Bartik].csv  (marginal effects)
*            Moderation_idealpoint_Table3_Marginal_Effects_Differences_CF[_Bartik].csv
*          + *_pctCI.csv versions of Tables 2-3 (percentile bootstrap CIs);
*            replicate draws saved in output/bootstrap/
* =========================================================

do "$ANALYSIS/00_prepare.do"
do "$ANALYSIS/00_bootstrap_pctci.do"

* =========================================================
* 5.2 IV Regression Moderation Effects: Gravity IV -- CF Approach -- political alignment
* =========================================================

global controls gdppc fdi gov_effectiveness agri_va rural_pop co2 total_erp_pct trade
* Alignment as mean over t-3..t-1, the same window as the centralities (*_ma3)
global moderator idealpoint_align_ma3

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

    display "P10 = $p10"
    display "P25 = $p25"
    display "P50 = $p50"
    display "P75 = $p75"
    display "P90 = $p90"

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
        saving("$BOOTDRAWS/mod07_gravity_margins_`var'.dta", replace): ///
        my_cf_margins_boot if estimation_sample == 1

    eststo margins_`var'
    pctci_post pf_me "$BOOTDRAWS/mod07_gravity_margins_`var'.dta" `var'

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
        reps($BOOT_REPS) cluster(iso3_id) idcluster(new_iso3_diff) ///
        saving("$BOOTDRAWS/mod07_gravity_mediff_`var'.dta", replace): ///
        my_cf_me_diff_boot if estimation_sample == 1

    eststo mediff_`var'
    pctci_post pf_diff "$BOOTDRAWS/mod07_gravity_mediff_`var'.dta" `var'
}


esttab reg_in_degree reg_out_degree reg_closeness reg_eigenvector ///
    using "$TABLES/Moderation_idealpoint_Table1_Main_Regression_CF.csv", replace ///
    keep(*_ma3 ${moderator} *#* v_resid1_boot v_resid2_boot) ///
    b(3) se(3) star(* 0.1 ** 0.05 *** 0.01) ///
    mtitles("In-Degree" "Out-Degree" "Closeness" "Eigenvector") ///
    scalars("Controls Controls" "FE Region & Year FE") ///
    compress nogaps label ///
    title("Second-Stage PPML Regression Results")


esttab margins_in_degree margins_out_degree margins_closeness margins_eigenvector ///
    using "$TABLES/Moderation_idealpoint_Table2_Main_Regression_CF.csv", replace ///
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
    using "$TABLES/Moderation_idealpoint_Table3_Marginal_Effects_Differences_CF.csv", replace ///
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

pctci_export pf_me "`pct_me'" "$TABLES/Moderation_idealpoint_Table2_Main_Regression_CF_pctCI.csv" "P10 P25 P50 P75 P90"
pctci_export pf_diff "`pct_diff'" "$TABLES/Moderation_idealpoint_Table3_Marginal_Effects_Differences_CF_pctCI.csv" "P10 P25 P50 P75 P90"


* =========================================================
* 5.4 IV Regression Moderation Effects: Bartik IV -- CF Approach -- political alignment
* =========================================================

global controls gdppc fdi gov_effectiveness agri_va rural_pop co2 total_erp_pct trade
* Alignment as mean over t-3..t-1, the same window as the centralities (*_ma3)
global moderator idealpoint_align_ma3

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
    * FIX: restored (missing in Reg_model.do, which made the Table 3 esttab below fail with r(111))
    return scalar diff_p75_p50_adj = _b[4._at] - _b[3._at]
    return scalar diff_p90_p75 = _b[5._at] - _b[4._at]
end

local centralities in_degree out_degree closeness eigenvector

* Bootstrap seed for this section: 2468 (other modules and the gravity section
* use 123456). Chosen on 2026-10-24 after a seed check with 123456, 20261024,
* 2468 and 987654; results for all four seeds in output/seedcheck/ and README_zh.
local bartik_seed 2468

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

    display "P10 = $p10"
    display "P25 = $p25"
    display "P50 = $p50"
    display "P75 = $p75"
    display "P90 = $p90"

    capture drop new_iso3_reg
    set seed `bartik_seed'

    quietly bootstrap, reps($BOOT_REPS) cluster(iso3_id) idcluster(new_iso3_reg): ///
        my_cf_reg_boot if estimation_sample == 1

    estadd local Controls "Yes"
    estadd local FE "Yes"
    eststo reg_`var'

    capture drop new_iso3_mar
    set seed `bartik_seed'

    quietly bootstrap, reps($BOOT_REPS) cluster(iso3_id) idcluster(new_iso3_mar) ///
        saving("$BOOTDRAWS/mod07_bartik_margins_`var'.dta", replace): ///
        my_cf_margins_boot if estimation_sample == 1

    eststo margins_`var'
    pctci_post pf_me "$BOOTDRAWS/mod07_bartik_margins_`var'.dta" `var'

    capture drop new_iso3_diff
    set seed `bartik_seed'

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
        reps($BOOT_REPS) cluster(iso3_id) idcluster(new_iso3_diff) ///
        saving("$BOOTDRAWS/mod07_bartik_mediff_`var'.dta", replace): ///
        my_cf_me_diff_boot if estimation_sample == 1

    eststo mediff_`var'
    pctci_post pf_diff "$BOOTDRAWS/mod07_bartik_mediff_`var'.dta" `var'

}


esttab reg_in_degree reg_out_degree reg_closeness reg_eigenvector ///
    using "$TABLES/Moderation_idealpoint_Table1_Main_Regression_CF_Bartik.csv", replace ///
    keep(*_ma3 ${moderator} *#* v_resid1_boot v_resid2_boot) ///
    b(3) se(3) star(* 0.1 ** 0.05 *** 0.01) ///
    mtitles("In-Degree" "Out-Degree" "Closeness" "Eigenvector") ///
    scalars("Controls Controls" "FE Region & Year FE") ///
    compress nogaps label ///
    title("Second-Stage PPML Regression Results")


esttab margins_in_degree margins_out_degree margins_closeness margins_eigenvector ///
    using "$TABLES/Moderation_idealpoint_Table2_Main_Regression_CF_Bartik.csv", replace ///
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
    using "$TABLES/Moderation_idealpoint_Table3_Marginal_Effects_Differences_CF_Bartik.csv", replace ///
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

pctci_export pf_me "`pct_me'" "$TABLES/Moderation_idealpoint_Table2_Main_Regression_CF_Bartik_pctCI.csv" "P10 P25 P50 P75 P90"
pctci_export pf_diff "`pct_diff'" "$TABLES/Moderation_idealpoint_Table3_Marginal_Effects_Differences_CF_Bartik_pctCI.csv" "P10 P25 P50 P75 P90"
