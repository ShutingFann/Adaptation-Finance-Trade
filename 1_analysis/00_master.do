* =========================================================
* 00_master.do  Run the regression analysis
*
* Prerequisite: the data pipeline in Workflow/1_data/ has produced
*               output/data/regression_panel_treated_non_high_income.csv
*
* Every module calls 00_prepare.do first, so each one can be switched
* on/off below and run independently.
* =========================================================
clear all
set more off

* ---------------------------------------------------------
* Paths: WORKFLOW defaults to the parent of the current directory
* (run from Workflow/2_analysis, as run_all.sh does); set it
* beforehand to run from elsewhere
* ---------------------------------------------------------
if "$WORKFLOW" == "" {
    quietly cd ..
    global WORKFLOW "`c(pwd)'"
    quietly cd "$WORKFLOW/2_analysis"
}
global ANALYSIS "$WORKFLOW/2_analysis"
global DATA     "$WORKFLOW/output/data"
global TABLES   "$WORKFLOW/output/tables"

* ---------------------------------------------------------
* Bootstrap settings (paper: 1000 cluster-bootstrap reps,
* 50 Monte Carlo draws per rep in the counterfactual matrix).
* Set lower values beforehand for a quick test run.
* ---------------------------------------------------------
if "$BOOT_REPS" == "" global BOOT_REPS 1000
if "$MC_REPS"   == "" global MC_REPS 50

* ---------------------------------------------------------
* Required packages (install once)
* ---------------------------------------------------------
* ssc install ftools
* ssc install reghdfe
* ssc install ppmlhdfe
* ssc install estout
* ssc install rangestat
* ssc install astile
* ssc install winsor2

* ---------------------------------------------------------
* Modules (1 = run, 0 = skip)                      Paper output
* ---------------------------------------------------------
local run_01 1   // Baseline PPML                       Fig. 4a
local run_01b 1  // Baseline with stepwise controls     SI table
local run_02 1   // IV-CF: gravity & Bartik IV          Fig. 4a
local run_03 1   // Counterfactual 4x4 matrices         Fig. 4b  (slowest)
local run_04 1   // Bilateral vs multilateral           Fig. 5
local run_04b 1  // Bi vs multi, relative (PPML beta)   main text
local run_05 1   // By income group                     Fig. 5
local run_06 1   // Moderation: disasters               Fig. 6a-b
local run_06c 1  // Disaster: ME curve P50-P95          Fig. 6a-b
local run_07 1   // Moderation: political alignment     Fig. 6c-d
local run_07b 1  // Alignment: competing Z, timing      SI table
local run_07c 1  // Alignment: ME curve P10-P90         Fig. 6c-d
local run_08 1   // Robustness                          SI tables

cd "$ANALYSIS"
capture mkdir "$TABLES"

if `run_01' do "$ANALYSIS/01_baseline_ppml.do"
if `run_01b' do "$ANALYSIS/01b_controls_robustness.do"
if `run_02' do "$ANALYSIS/02_iv_cf_main.do"
if `run_03' do "$ANALYSIS/03_counterfactual_matrix.do"
if `run_04' do "$ANALYSIS/04_bilateral_multilateral.do"
if `run_04b' do "$ANALYSIS/04b_bilateral_multilateral_relative.do"
if `run_05' do "$ANALYSIS/05_income_heterogeneity.do"
if `run_06' do "$ANALYSIS/06_moderation_disaster.do"
if `run_06c' do "$ANALYSIS/06c_moderation_disaster_curve.do"
if `run_07' do "$ANALYSIS/07_moderation_idealpoint.do"
if `run_07b' do "$ANALYSIS/07b_moderation_idealpoint_robustness.do"
if `run_07c' do "$ANALYSIS/07c_moderation_idealpoint_curve.do"
if `run_08' do "$ANALYSIS/08_robustness.do"

display "All selected modules finished. Tables in $TABLES"
