* =========================================================
* 04_dv_means.do  Mean of bilateral / multilateral adaptation finance in each
*                 module 04 estimation sample (called at the end of module 04)
* Output : output/tables/regression_results_Bi_Multi_DVmean.csv
*
* With fixed effects, PPML fitted values average to the sample mean of y, so
* AME / mean(y) = PPML coefficient. Fig. 5a-b divides the module 04 AMEs (and
* their SEs) by these means to show effects relative to each channel's size.
* Point estimates only: re-runs the full-sample CF model, no bootstrap.
* =========================================================

do "$ANALYSIS/00_prepare.do"

global controls gdppc fdi gov_effectiveness disaster_total_impact_pc agri_va rural_pop co2 total_erp_pct trade

tempfile means
tempname pf
postfile `pf' str8 iv str12 centrality str12 channel double(dv_mean ppml_b) long N using "`means'", replace

foreach iv in gravity bartik {
    foreach var in in_degree out_degree closeness eigenvector {
        if "`iv'" == "gravity" local z "`var'_iv_ma3"
        else                   local z ln_iv_seaborne_trade_total_ma3

        capture drop v_resid
        quietly reghdfe `var'_ma3 `z' $controls, absorb(region_id year) resid
        quietly predict v_resid, resid

        foreach dv in bilateral_adaptation multilateral_adaptation {
            quietly ppmlhdfe `dv' `var'_ma3 $controls v_resid, absorb(region_id year)
            local b = _b[`var'_ma3]
            quietly summarize `dv' if e(sample)
            local channel = cond("`dv'" == "bilateral_adaptation", "Bilateral", "Multilateral")
            post `pf' ("`iv'") ("`var'") ("`channel'") (r(mean)) (`b') (r(N))
        }
    }
}
postclose `pf'

preserve
use "`means'", clear
list, sepby(iv)
export delimited using "$TABLES/regression_results_Bi_Multi_DVmean.csv", replace
restore
