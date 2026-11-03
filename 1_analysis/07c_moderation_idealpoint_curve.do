* =========================================================
* 07c_moderation_idealpoint_curve.do  Marginal effect of centrality over the
*        full range of political alignment (Fig. 6c-d plots P10-P75;
*        gravity IV, Bartik IV)
* Output : output/tables/Moderation_idealpoint_Curve_CF[_Bartik].csv
*            one row per centrality x grid point: percentile, alignment
*            value, marginal effect (US$ million per 1 SD of centrality),
*            bootstrap SD; the figure computes the trimmed SE from the draws
*          output/tables/Moderation_idealpoint_Curve_moderator[_Bartik].csv
*            alignment values of the estimation sample (for the histogram)
*          replicate draws: output/bootstrap/mod07c_<iv>_curve_<centrality>.dta
*
* Same model, moderator (idealpoint_align_ma3), estimation sample and seeds as
* 07_moderation_idealpoint.do (gravity 123456, Bartik 2468), so the replicates
* draw the same countries and the marginal effects at P10/P25/P50/P75/P90
* equal module 07 Table 2. Only the grid is denser: P10, P15, ..., P90.
* =========================================================

do "$ANALYSIS/00_prepare.do"

if "$BOOTDRAWS" == "" global BOOTDRAWS "$WORKFLOW/output/bootstrap"
capture mkdir "$BOOTDRAWS"

global controls gdppc fdi gov_effectiveness agri_va rural_pop co2 total_erp_pct trade
global moderator idealpoint_align_ma3

capture xtset, clear
capture tsset, clear

local grid_pcts 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90

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

    if "`ivtype'" == "gravity" local seed 123456
    else                       local seed 2468

    tempfile res modvals
    capture postclose pf_curve
    postfile pf_curve str8 iv str12 centrality int pct double(mod_value me se_boot) ///
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

        * Estimation sample as in module 07
        capture drop estimation_sample
        quietly reghdfe total_adaptation_finance ///
            ${current_endog} ${current_endog_X_mod} ///
            ${current_iv} ${current_iv_X_mod} ///
            ${moderator} $controls, absorb(region_id year)
        gen byte estimation_sample = e(sample)
        quietly count if estimation_sample
        local nobs = r(N)

        * Grid of alignment values at percentiles of the estimation sample
        global grid_vals ""
        foreach p of local grid_pcts {
            _pctile $moderator if estimation_sample == 1, percentiles(`p')
            local v_`p' = r(r1)
            global grid_vals "${grid_vals} `v_`p''"
        }
        display as text _n "== `ivtype' IV, `var': N = `nobs', seed `seed'"

        local draws "$BOOTDRAWS/mod07c_`ivtype'_curve_`var'.dta"
        capture drop new_iso3_curve
        set seed `seed'
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
            local p : word `k' of `grid_pcts'
            quietly summarize `v'
            post pf_curve ("`ivtype'") ("`var'") (`p') (`v_`p'') (el(`B', 1, `k')) ///
                (r(sd)) (r(N)) (`nobs')
        }
        restore

        * Alignment values of the estimation sample (in_degree sample for the histogram)
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
    sort centrality pct
    export delimited using "$TABLES/Moderation_idealpoint_Curve_CF`suffix'.csv", replace
    quietly use "`modvals'", clear
    export delimited using "$TABLES/Moderation_idealpoint_Curve_moderator`suffix'.csv", replace
    restore
}
