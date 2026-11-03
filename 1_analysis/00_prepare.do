* =========================================================
* 00_prepare.do  Load the main sample and build analysis variables
* Origin : Code/Reg_model.do, lines 1-91
* Called at the top of every analysis module (globals set in 00_master.do)
*
* - Sample       : non-high-income countries, 2000-2024
* - Controls     : z-standardised; skewed controls transformed first:
*                  fdi -> asinh(fdi) (% of GDP, can be negative);
*                  disaster_total_impact_pc -> ln(1 + 100 x) (% of population)
* - Centralities : in/out-degree -> ln + std; betweenness/eigenvector ->
*                  ln(x*1e6+1) + std; closeness -> std (same for *_iv)
* - *_ma3 / *_ma5: mean over years t-3..t-1 / t-5..t-1 (lagged moving average),
*                  for the centralities and for both IVs (*_iv_ma3,
*                  ln_iv_seaborne_trade_total_ma3); windows built in Step 12
*                  on the full source panels
* - Alignment    : idealpoint_align (t), idealpoint_align_l1 (t-1) and
*                  idealpoint_align_ma3 (mean t-3..t-1, from Step 06), std
* =========================================================
clear
eststo clear
set more off

import delimited "$DATA/regression_panel_treated_non_high_income.csv", clear

keep if inrange(year, 2000, 2024)

* Skewed controls (2026-10-25): compress the tails before standardisation.
* Raw FDI reaches 103% of GDP (z = 16.5) and caused a diverging bootstrap replicate
* in module 05; disaster impact per capita has skewness 6.3.
replace fdi = asinh(fdi)
replace disaster_total_impact_pc = ln(1 + 100 * disaster_total_impact_pc)

label variable in_degree 			"In-Degree Centrality"
label variable gdppc 				"GDP per capita"
label variable fdi 					"Foreign Direct Investment"
label variable gain_index 			"GAIN Index"
label variable gov_effectiveness 	"Government Effectiveness"
label variable disaster_affected 	"Disaster Affected"
label variable total_erp_pct 		"Total Equity Risk Premium (%)"
label variable out_degree    		"Out-Degree Centrality"
label variable closeness     		"Closeness Centrality"
label variable betweenness   		"Betweenness Centrality"
label variable katz   				"katz Centrality"
label variable eigenvector   		"Eigenvector Centrality"
label variable recap   				"Renewable Capacity Installed"
label variable agri_va        		"Agricultural Value Added (%)"
label variable rural_pop     		"Rural Population (%)"
label variable co2           		"CO2 Emissions per capita"
label variable trade           		"Trade (% of GDP)"

local control_vars gdppc fdi gov_effectiveness recap agri_va rural_pop co2 total_erp_pct trade ///
		disaster_total_impact_pc idealpoint_align idealpoint_align_g7 iv_seaborne_trade_total iv_seaborne_trade_total_ch iv_seaborne_trade_total_gr trade_loaded trade_discharged trade_total ln_iv_seaborne_trade_total

foreach var of local control_vars {
    quietly summarize `var'
    local m_`var' = r(mean)
    local s_`var' = r(sd)
    egen temp_`var' = std(`var')
    replace `var' = temp_`var'
    drop temp_`var'
    local var_label : variable label `var'
    label var `var' "`var_label' (Standardized)"
}

ds, has(type numeric)
local numvars `r(varlist)'

foreach var of local numvars {

    local lower_var = strlower("`var'")

    if inlist("`lower_var'", "in_degree", "out_degree", "in_degree_iv", "out_degree_iv") {
        replace `var' = ln(`var')

        quietly summarize `var'
        local m_`var' = r(mean)
        local s_`var' = r(sd)
        egen tem_`var' = std(`var')
        replace `var' = tem_`var'
        drop tem_`var'
    }

	else if inlist("`lower_var'", "betweenness", "eigenvector", "betweenness_iv", "eigenvector_iv") {
        replace `var' = ln(`var' * 1000000 + 1)

        quietly summarize `var'
        local m_`var' = r(mean)
        local s_`var' = r(sd)
        egen tem_`var' = std(`var')
        replace `var' = tem_`var'
        drop tem_`var'
    }

    else if inlist("`lower_var'", "closeness", "closeness_iv") {
        quietly summarize `var'
        local m_`var' = r(mean)
        local s_`var' = r(sd)
        egen tem_`var' = std(`var')
        replace `var' = tem_`var'
        drop tem_`var'
    }
}

* Lagged alignment from Step 06: t-1 value and mean over t-3..t-1 (same window
* as the centralities), standardised over the loaded panel like idealpoint_align.
* idealpoint_align_ma3 is the moderator of module 07; both are used in 07b.
preserve
import delimited "$DATA/IdealPoint_Alignment_lags_2000_2024.csv", clear
rename iso3c iso3
keep iso3 year idealpoint_alignment_l1 idealpoint_alignment_ma3
tempfile ip_lags
save `ip_lags'
restore

merge m:1 iso3 year using `ip_lags', keep(master match) nogenerate

foreach v in l1 ma3 {
    egen idealpoint_align_`v' = std(idealpoint_alignment_`v')
    drop idealpoint_alignment_`v'
}
label variable idealpoint_align_l1  "IdealPoint_Align, t-1 (Standardized)"
label variable idealpoint_align_ma3 "IdealPoint_Align, mean t-3..t-1 (Standardized)"

encode iso3, gen(iso3_id)
encode m49_region, gen(region_id)
xtset iso3_id year

tab year, gen(yd_)
drop yd_1


* Lagged moving averages: ma{3,5}_<var> come from Step 12 on the transformed
* scale above (before standardisation). Standardising them with the moments of
* the contemporaneous variable gives x_ma3 = mean of standardised x over
* t-3..t-1, the same definition as the former rangestat on the standardised
* variable. Each IV gets the same window as the centrality it instruments.
local ma_vars in_degree out_degree closeness betweenness eigenvector ///
	in_degree_iv out_degree_iv closeness_iv betweenness_iv eigenvector_iv ///
	ln_iv_seaborne_trade_total

foreach var of local ma_vars {
    foreach w in 3 5 {
        gen `var'_ma`w' = (ma`w'_`var' - `m_`var'') / `s_`var''
        label var `var'_ma`w' "`: var label `var''"
        drop ma`w'_`var'
    }
}
