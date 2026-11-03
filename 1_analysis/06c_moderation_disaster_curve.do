* =========================================================
* 06c_moderation_disaster_curve.do  Marginal effect of centrality over the
*        range of disaster impact (Fig. 6a-b; gravity IV, Bartik IV)
* Output : output/tables/Moderation_disaster_Curve_CF[_Bartik].csv
*            one row per centrality x grid point: label (P50/P75/P90/P95 or
*            "grid"), moderator value, marginal effect (US$ million per 1 SD
*            of centrality), bootstrap SD; the figure computes the trimmed SE
*          output/tables/Moderation_disaster_Curve_moderator[_Bartik].csv
*            moderator values of the estimation sample (for the histogram)
*          replicate draws: output/bootstrap/mod06c_<iv>_curve_<centrality>.dta
*
* Same model, moderator (disaster_total_impact_pc), estimation sample and seed
* (123456) as 06_moderation_disaster.do, so the replicates draw the same
* countries and the marginal effects at P50/P75/P90/P95 equal module 06 Table 2.
* Range P50-P95 as in Fig. 6a-b. The moderator is right-skewed even after the
* ln(1 + 100 x) transform (P10-P25 nearly equal), so the grid is 14 points evenly
* spaced in value between P50 and P95, plus P50, P75, P90 and P95 themselves.
* A P10-P90 version (annotated ME(P90) - ME(P10)) was tried on 2026-10-26 and
* kept in output/backup_disaster_p10p90_20261026/.
* =========================================================

do "$ANALYSIS/00_prepare.do"

if "$BOOTDRAWS" == "" global BOOTDRAWS "$WORKFLOW/output/bootstrap"
capture mkdir "$BOOTDRAWS"

global controls gdppc fdi gov_effectiveness agri_va rural_pop co2 total_erp_pct trade
global moderator disaster_total_impact_pc

capture xtset, clear
capture tsset, clear

local n_even 14

capture program drop cf_curve_boot
program cf_curve_boot, eclass
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

    quietly margins, dydx(${current_endog}) at(${moderator}=(${grid_vals})) post
end

if "$IVTYPES" == "" global IVTYPES gravity bartik

local centralities in_degree out_degree closeness eigenvector

foreach ivtype of global IVTYPES {

    tempfile res modvals
    capture postclose pf_curve
    postfile pf_curve str8 iv str12 centrality str5 label double(mod_value me se_boot) ///
        long(reps nobs) using "`res'", replace

    foreach var of local centralities {

        global current_endog "`var'_ma3"
        if "`ivtype'" == "gravity" global current_iv "`var'_iv_ma3"
        else                       global current_iv "ln_iv_seaborne_trade_total_ma3"

        capture drop `var'_X_mod `var'_iv_X_mod
        quietly gen `var'_X_mod = ${current_endog} * ${moderator}
        quietly gen `var'_iv_X_mod = ${current_iv} * ${moderator}
        global current_endog_X_mod "`var'_X_mod"
        global current_iv_X_mod "`var'_iv_X_mod"

        * Estimation sample as in module 06
        capture drop estimation_sample
        quietly reghdfe total_adaptation_finance ///
            ${current_endog} ${current_endog_X_mod} ///
            ${current_iv} ${current_iv_X_mod} ///
            ${moderator} $controls, absorb(region_id year)
        gen byte estimation_sample = e(sample)
        quietly count if estimation_sample
        local nobs = r(N)

        * Percentiles as in module 06 (summarize, detail), plus an even grid in value
        quietly summarize $moderator if estimation_sample == 1, detail
        local p50 = r(p50)
        local p75 = r(p75)
        local p90 = r(p90)
        local p95 = r(p95)

        * Grid points sorted by value, each with its label: interior points of
        * the even grid (its endpoints are P50 and P95 themselves) + P50/P75/P90/P95
        tempname G
        matrix `G' = J(`n_even' + 2, 2, .)
        forvalues i = 2/`=`n_even' - 1' {
            matrix `G'[`i' - 1, 1] = `p50' + (`p95' - `p50') * (`i' - 1) / (`n_even' - 1)
            matrix `G'[`i' - 1, 2] = 0
        }
        local k = `n_even' - 2
        foreach p in 50 75 90 95 {
            local ++k
            matrix `G'[`k', 1] = `p`p''
            matrix `G'[`k', 2] = `p'
        }
        preserve
        clear
        quietly svmat double `G', names(g)
        sort g1
        global grid_vals ""
        local grid_labels ""
        forvalues i = 1/`=_N' {
            global grid_vals "${grid_vals} `=g1[`i']'"
            local lab = cond(g2[`i'] > 0, "P" + string(g2[`i']), "grid")
            local grid_labels "`grid_labels' `lab'"
        }
        restore

        display as text _n "== `ivtype' IV, `var': N = `nobs', `: word count ${grid_vals}' grid points"

        local draws "$BOOTDRAWS/mod06c_`ivtype'_curve_`var'.dta"
        capture drop new_iso3_curve
        set seed 123456
        quietly bootstrap, reps($BOOT_REPS) cluster(iso3_id) idcluster(new_iso3_curve) ///
            saving("`draws'", replace): ///
            cf_curve_boot if estimation_sample == 1

        tempname B
        matrix `B' = e(b)
        preserve
        quietly use "`draws'", clear
        local k = 0
        foreach v of varlist _all {
            local ++k
            local lab : word `k' of `grid_labels'
            local val : word `k' of ${grid_vals}
            quietly summarize `v'
            post pf_curve ("`ivtype'") ("`var'") ("`lab'") (`val') (el(`B', 1, `k')) ///
                (r(sd)) (r(N)) (`nobs')
        }
        restore

        * Moderator values of the estimation sample (in_degree sample for the histogram)
        if "`var'" == "in_degree" {
            preserve
            quietly keep if estimation_sample == 1
            keep iso3 year $moderator
            gen iv = "`ivtype'"
            quietly save "`modvals'", replace
            restore
        }
    }

    postclose pf_curve

    * One file per IV, so the two IVs can run in separate Stata processes
    if "`ivtype'" == "gravity" local suffix ""
    else                       local suffix "_Bartik"
    preserve
    quietly use "`res'", clear
    sort centrality mod_value
    export delimited using "$TABLES/Moderation_disaster_Curve_CF`suffix'.csv", replace
    quietly use "`modvals'", clear
    export delimited using "$TABLES/Moderation_disaster_Curve_moderator`suffix'.csv", replace
    restore
}
