* =========================================================
* 07b_moderation_idealpoint_robustness.do  Political-alignment moderation:
*        competing interactions and moderator timing (gravity IV, Bartik IV)
* Output : output/tables/
*            Moderation_idealpoint_Robustness_CF_pctCI.csv          (gravity IV)
*            Moderation_idealpoint_Robustness_CF_Bartik_pctCI.csv   (Bartik IV)
*          with suffix _<ROB_TAG> if global ROB_TAG is set (keeps earlier tables);
*          replicate draws saved in output/bootstrap/
*
* Model (CF-PPML, as in 07_moderation_idealpoint.do):
*   y = exp(b1 x + b3 x*m + sum_k c_k x*Z_k + b2 m + Z + controls + v_hat + FE)
*   Endogenous: x, x*m, x*Z_k; excluded instruments: iv, iv*m, iv*Z_k.
*   One first stage (on all instruments) and one CF residual per endogenous term.
*
* Specifications (one common estimation sample per IV, so that differences
* across rows come from the specification, not from the sample); by default
* only base, z_all, ma3 and ma3_zall are run (global SPECS below):
*   base      m = alignment (t),     no Z               (module 07 before 2026-10-24)
*   z_gdp     m = alignment (t),     Z = GDP per capita
*   z_gov     m = alignment (t),     Z = government effectiveness
*   z_umic    m = alignment (t),     Z = upper-middle-income dummy
*   z_all     m = alignment (t),     Z = all three
*   lag1      m = alignment (t-1),   no Z
*   ma3       m = alignment, mean over t-3..t-1 (same window as centrality;
*             moderator of module 07, here on the common sample)
*   ma3_zall  m = alignment MA3,     Z = all three
*
* Reported for each specification and centrality (percentile bootstrap CIs):
*   b_x, b_xm (= b3), b_xz*      PPML coefficients
*   semi_pXX = b1 + b3 m_pXX + sum_k c_k mean(Z_k)   proportional effect
*   me_pXX   = d E[y] / d x at m = pXX (US$ million per 1 SD of centrality)
*   d_pA_pB  = me_pA - me_pB
* Lagged alignment comes from Step 06 (IdealPoint_Alignment_lags_2000_2024.csv).
* =========================================================

do "$ANALYSIS/00_prepare.do"
do "$ANALYSIS/00_bootstrap_pctci.do"

global controls gdppc fdi gov_effectiveness agri_va rural_pop co2 total_erp_pct trade

* Lagged alignment (idealpoint_align_l1, idealpoint_align_ma3) is built in 00_prepare.do

* Upper-middle-income dummy (World Bank group in the panel)
gen byte umic = (incomegroup == "Upper middle income") if incomegroup != ""
label variable umic "Upper middle income"

* Cluster bootstrap with idcluster() requires the panel setting to be cleared
capture xtset, clear
capture tsset, clear

* ---------------------------------------------------------
* Specifications: moderator and competing moderators Z
* ---------------------------------------------------------
* Default: reduced set. For all eight, set beforehand
*   global SPECS base z_gdp z_gov z_umic z_all lag1 ma3 ma3_zall
if "$SPECS" == "" global SPECS base z_all ma3 ma3_zall
local specs $SPECS

local mod_base     idealpoint_align
local mod_z_gdp    idealpoint_align
local mod_z_gov    idealpoint_align
local mod_z_umic   idealpoint_align
local mod_z_all    idealpoint_align
local mod_lag1     idealpoint_align_l1
local mod_ma3      idealpoint_align_ma3
local mod_ma3_zall idealpoint_align_ma3

local z_base
local z_z_gdp    gdppc
local z_z_gov    gov_effectiveness
local z_z_umic   umic
local z_z_all    gdppc gov_effectiveness umic
local z_lag1
local z_ma3
local z_ma3_zall gdppc gov_effectiveness umic

local all_mods idealpoint_align idealpoint_align_l1 idealpoint_align_ma3
local all_z    gdppc gov_effectiveness umic

* ---------------------------------------------------------
* One bootstrap replicate: first stages, CF-PPML, semi-elasticities, margins
* Globals: current_endog, moderator, endog_list, iv_list, zvars, zmain,
*          xz_terms, p10..p90, zbar1..zbar3
* ---------------------------------------------------------
capture program drop cf_mod_rob_boot
program cf_mod_rob_boot, rclass
    version 14.0
    syntax [if] [in]

    marksample touse
    capture drop vr_*

    local k = 0
    foreach e of global endog_list {
        local ++k
        quietly reghdfe `e' $iv_list $moderator $zmain $controls ///
            if `touse', absorb(region_id year) resid
        quietly predict vr_`k' if e(sample), resid
    }

    quietly ppmlhdfe total_adaptation_finance ///
        c.${current_endog}##c.${moderator} $xz_terms ///
        $zmain $controls vr_* ///
        if `touse', absorb(region_id year) d

    * Store everything in locals first: margins below overwrites r()
    local b_x  = _b[${current_endog}]
    local b_xm = _b[c.${current_endog}#c.${moderator}]
    local zshift = 0
    local j = 0
    foreach z of global zvars {
        local ++j
        local b_xz`j' = _b[c.${current_endog}#c.`z']
        local zshift = `zshift' + `b_xz`j'' * ${zbar`j'}
    }
    foreach p in 10 25 50 75 90 {
        local semi_p`p' = `b_x' + `b_xm' * ${p`p'} + `zshift'
    }

    quietly margins, dydx(${current_endog}) ///
        at(${moderator}=($p10 $p25 $p50 $p75 $p90)) post

    return scalar b_x  = `b_x'
    return scalar b_xm = `b_xm'
    forvalues i = 1/`j' {
        return scalar b_xz`i' = `b_xz`i''
    }
    foreach p in 10 25 50 75 90 {
        return scalar semi_p`p' = `semi_p`p''
    }
    return scalar me_p10 = _b[1._at]
    return scalar me_p25 = _b[2._at]
    return scalar me_p50 = _b[3._at]
    return scalar me_p75 = _b[4._at]
    return scalar me_p90 = _b[5._at]
    return scalar d_p10_p50 = _b[1._at] - _b[3._at]
    return scalar d_p25_p50 = _b[2._at] - _b[3._at]
    return scalar d_p75_p50 = _b[4._at] - _b[3._at]
    return scalar d_p90_p50 = _b[5._at] - _b[3._at]
    return scalar d_p75_p25 = _b[4._at] - _b[2._at]
    return scalar d_p90_p10 = _b[5._at] - _b[1._at]
end

* Percentile CIs from a draw file, posted with the specification label
* (same rules as pctci_post in 00_bootstrap_pctci.do, plus the bootstrap SD)
capture program drop rob_post
program rob_post
    args handle drawfile iv spec moderator zlist centrality nobs
    tempname B
    matrix `B' = e(b)
    local names : colfullnames `B'
    preserve
    quietly use "`drawfile'", clear
    local k = 0
    foreach v of varlist _all {
        local ++k
        local term : word `k' of `names'
        quietly summarize `v'
        local reps = r(N)
        local sd = r(sd)
        quietly count if `v' <= 0 & !missing(`v')
        local share = r(N) / `reps'
        local p = 2 * min(`share', 1 - `share')
        _pctile `v', percentiles(0.5 2.5 5 95 97.5 99.5)
        post `handle' ("`iv'") ("`spec'") ("`moderator'") ("`zlist'") ("`centrality'") ///
            ("`term'") (el(`B', 1, `k')) (`sd') ///
            (r(r3)) (r(r4)) (r(r2)) (r(r5)) (r(r1)) (r(r6)) (`p') (`reps') (`nobs')
    }
    restore
end

local centralities in_degree out_degree closeness eigenvector

* IVs to run (set global IVTYPES beforehand to run one IV per Stata process)
if "$IVTYPES" == "" global IVTYPES gravity bartik

foreach ivtype of global IVTYPES {

    tempfile res
    capture postclose pf_rob
    postfile pf_rob str8 iv str10 spec str24 moderator str40 zvars str12 centrality ///
        str16 term double(estimate se_boot ci90_lo ci90_hi ci95_lo ci95_hi ci99_lo ci99_hi p_boot) ///
        long(reps nobs) using "`res'", replace

    foreach var of local centralities {

        global current_endog "`var'_ma3"
        if "`ivtype'" == "gravity" global current_iv "`var'_iv_ma3"
        else                       global current_iv "ln_iv_seaborne_trade_total_ma3"

        * Interaction terms for every moderator and Z used below
        capture drop xm_* ivm_* xz_* ivz_*
        foreach m of local all_mods {
            quietly gen xm_`m'  = ${current_endog} * `m'
            quietly gen ivm_`m' = ${current_iv} * `m'
        }
        foreach z of local all_z {
            quietly gen xz_`z'  = ${current_endog} * `z'
            quietly gen ivz_`z' = ${current_iv} * `z'
        }

        * Common estimation sample for all specifications of this centrality and IV
        capture drop estimation_sample
        quietly reghdfe total_adaptation_finance ${current_endog} ${current_iv} ///
            `all_mods' `all_z' $controls, absorb(region_id year)
        gen byte estimation_sample = e(sample)
        quietly count if estimation_sample
        local nobs = r(N)
        display as text _n "== `ivtype' IV, `var': common sample N = `nobs'"

        foreach s of local specs {

            global moderator "`mod_`s''"
            global zvars "`z_`s''"

            * umic is not among the controls, so its main effect is added here
            global zmain ""
            if strpos(" ${zvars} ", " umic ") global zmain "umic"

            global endog_list "${current_endog} xm_${moderator}"
            global iv_list    "${current_iv} ivm_${moderator}"
            global xz_terms   ""
            local j = 0
            foreach z of global zvars {
                local ++j
                global endog_list "${endog_list} xz_`z'"
                global iv_list    "${iv_list} ivz_`z'"
                global xz_terms   "${xz_terms} c.${current_endog}#c.`z'"
                quietly summarize `z' if estimation_sample
                global zbar`j' = r(mean)
            }

            quietly summarize $moderator if estimation_sample, detail
            foreach p in 10 25 50 75 90 {
                global p`p' = r(p`p')
            }

            local stats b_x=r(b_x) b_xm=r(b_xm)
            forvalues i = 1/`j' {
                local stats `stats' b_xz`i'=r(b_xz`i')
            }
            foreach p in 10 25 50 75 90 {
                local stats `stats' semi_p`p'=r(semi_p`p')
            }
            foreach p in 10 25 50 75 90 {
                local stats `stats' me_p`p'=r(me_p`p')
            }
            local stats `stats' d_p10_p50=r(d_p10_p50) d_p25_p50=r(d_p25_p50) ///
                d_p75_p50=r(d_p75_p50) d_p90_p50=r(d_p90_p50) ///
                d_p75_p25=r(d_p75_p25) d_p90_p10=r(d_p90_p10)

            display as text "   spec `s': m = ${moderator}, Z = ${zvars}"

            local draws "$BOOTDRAWS/mod07b_`ivtype'_`s'_`var'.dta"
            capture drop new_iso3_rob
            set seed 123456
            quietly bootstrap `stats', ///
                reps($BOOT_REPS) cluster(iso3_id) idcluster(new_iso3_rob) ///
                saving("`draws'", replace): ///
                cf_mod_rob_boot if estimation_sample == 1

            rob_post pf_rob "`draws'" "`ivtype'" "`s'" "${moderator}" "${zvars}" "`var'" `nobs'
        }
    }

    postclose pf_rob

    preserve
    quietly use "`res'", clear
    quietly gen str3 stars = ""
    quietly replace stars = "*"   if ci90_lo > 0 | ci90_hi < 0
    quietly replace stars = "**"  if ci95_lo > 0 | ci95_hi < 0
    quietly replace stars = "***" if ci99_lo > 0 | ci99_hi < 0
    order iv spec moderator zvars centrality term estimate ci95_lo ci95_hi stars p_boot se_boot
    if "`ivtype'" == "gravity" local suffix ""
    else                       local suffix "_Bartik"
    if "$ROB_TAG" != "" local suffix "`suffix'_$ROB_TAG"
    export delimited using "$TABLES/Moderation_idealpoint_Robustness_CF`suffix'_pctCI.csv", replace
    restore
}
