* =========================================================
* 00_bootstrap_pctci.do  Percentile bootstrap confidence intervals
* Called by 06_moderation_disaster.do and 07_moderation_idealpoint.do
*
* Marginal effects of the PPML interaction models are highly non-linear in the
* tails of the moderator: a few bootstrap replicates with a large interaction
* coefficient produce extreme values, so the bootstrap SD (normal-approximation
* CI) is dominated by one or two replicates. These modules therefore also report
* percentile CIs (2.5/97.5 etc. of the replicate distribution) and a bootstrap
* p-value, 2 * min(share of replicates <= 0, share > 0).
*
* Usage:
*   pctci_open  <handle> <tempfile>            open a postfile for results
*   bootstrap ..., saving(<drawfile>, replace): ...
*   pctci_post  <handle> <drawfile> <centrality>   after each bootstrap
*   pctci_export <handle> <tempfile> <csv> "<percentile labels>"
* =========================================================

if "$BOOTDRAWS" == "" global BOOTDRAWS "$WORKFLOW/output/bootstrap"
capture mkdir "$BOOTDRAWS"

capture program drop pctci_open
program pctci_open
    args handle file
    capture postclose `handle'
    postfile `handle' str12 centrality str40 term double(estimate ci90_lo ci90_hi ///
        ci95_lo ci95_hi ci99_lo ci99_hi p_boot) long reps using "`file'", replace
end

capture program drop pctci_post
program pctci_post
    args handle drawfile centrality
    * Point estimates of the bootstrap just run (columns match the draw variables)
    tempname B
    matrix `B' = e(b)
    local names : colfullnames `B'
    preserve
    quietly use "`drawfile'", clear
    local k = 0
    foreach v of varlist _all {
        local ++k
        local term : word `k' of `names'
        quietly count if !missing(`v')
        local reps = r(N)
        quietly count if `v' <= 0 & !missing(`v')
        local share = r(N) / `reps'
        local p = 2 * min(`share', 1 - `share')
        _pctile `v', percentiles(0.5 2.5 5 95 97.5 99.5)
        post `handle' ("`centrality'") ("`term'") (el(`B', 1, `k')) ///
            (r(r3)) (r(r4)) (r(r2)) (r(r5)) (r(r1)) (r(r6)) (`p') (`reps')
    }
    restore
end

capture program drop pctci_export
program pctci_export
    args handle file csv pctlist
    postclose `handle'
    preserve
    quietly use "`file'", clear

    * Readable row labels
    quietly gen str40 label = term
    local k = 0
    foreach p of local pctlist {
        local ++k
        quietly replace label = "ME at `p'" if ustrregexm(term, "(^|:)`k'(bn)?\._at$")
    }
    quietly replace label = "ME at " + upper(substr(term, 4, .)) if substr(term, 1, 3) == "me_"
    tempvar t
    quietly gen `t' = subinstr(term, "_", " ", .)
    quietly replace label = "ME(" + upper(word(`t', 2)) + ") - ME(" + upper(word(`t', 3)) + ")" ///
        if substr(term, 1, 5) == "diff_"
    quietly replace label = label + " (adj)" if strpos(term, "_adj")
    drop `t'

    * Stars from the percentile CIs (CI excludes zero)
    quietly gen str3 stars = ""
    quietly replace stars = "*"   if ci90_lo > 0 | ci90_hi < 0
    quietly replace stars = "**"  if ci95_lo > 0 | ci95_hi < 0
    quietly replace stars = "***" if ci99_lo > 0 | ci99_hi < 0

    order centrality term label estimate ci95_lo ci95_hi stars p_boot ci90_lo ci90_hi ci99_lo ci99_hi reps
    export delimited using "`csv'", replace
    restore
end
