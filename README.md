# Global trade networks and international adaptation finance


## Repository structure

```
.
├── README.md
├── 1_analysis/     Stata code: econometric analysis (00_master.do runs all modules)
└── 2_figures/      Jupyter notebooks: main-text and Supplementary figures
```


## Requirements

**Stata 18** (tested with Stata MP) and these community-contributed packages:

```stata
ssc install ftools
ssc install reghdfe
ssc install ppmlhdfe
ssc install estout
ssc install rangestat
ssc install astile
ssc install winsor2
```

**Python 3.11** with Jupyter and:

```
pandas numpy scipy matplotlib seaborn geopandas shapely networkx adjustText
```

## Data availability

The analysis uses a country-year panel built from the sources below. The panel is not included
in this repository because some sources cannot be redistributed. [Describe how to obtain the
processed data, e.g. "available from the corresponding author on request" or a link to a data
repository.]

| Data | Source |
|---|---|
| Bilateral climate finance | OECD Creditor Reporting System project descriptions, classified with ClimateFinanceBERT |
| Multilateral climate finance | [Source] |
| Bilateral trade | IMF International Merchandise Trade Statistics (exports, FOB) |
| Gravity variables (gravity IV) | CEPII Gravity database V202211; CEPII GeoDist |
| Seaborne trade (Bartik IV) | UNCTADstat |
| Disasters | EM-DAT, CRED (not redistributable) |
| Political alignment | UNGA ideal points, Bailey, Strezhnev & Voeten (Harvard Dataverse) |
| Macro and governance controls | World Bank WDI and WGI, ND-GAIN, UNDP HDI, Damodaran country risk premiums, World Bank projects, WTO RTA database |
| Income groups and regions | World Bank country classification; UN M49 |
| Country boundaries (maps) | Natural Earth, 1:110m admin-0 countries |


## How to run

### 1. Econometric analysis

```bash
cd 1_analysis
stata-mp -b do 00_master.do
```

Run `00_master.do` from inside `1_analysis/`. It takes the parent folder as the project folder.
To run it from somewhere else, first set the project folder in Stata:
`global WORKFLOW "/path/to/project"`. You can switch individual modules on or off at the top of
`00_master.do`. Each module first runs `00_prepare.do`, so modules can be run independently.
Random seeds are set inside each module.

The full run takes about 8–9 hours on a laptop. Almost all of that is the cluster bootstrap
(1,000 replications; 50 Monte Carlo draws per replication in the counterfactual matrices). For a
quick test, set fewer replications first:

```stata
global BOOT_REPS 3
global MC_REPS 2
do 00_master.do
```

### 2. Figures

Open the notebooks in `2_figures/` and run them after the analysis has finished. The
regression figures read the tables in `output/tables/`.

## Econometric specification

All models are PPML with region (UN M49) and year fixed effects. Centralities enter as 3-year
lagged moving averages (mean of t-3 to t-1) and are standardised. Results are reported as average
marginal effects in USD million per 1 SD of centrality. In the IV models, the first-stage
residual enters the PPML (control function). Standard errors come from a country-cluster
bootstrap. Each instrument has the same timing as the centrality it instruments.

## Analysis modules (`1_analysis/`)

| Module | Content | Paper |
|---|---|---|
| `00_master.do` | Sets paths and bootstrap settings, runs all modules | – |
| `00_prepare.do` | Loads the panel; transformations, standardisation, moving averages | – |
| `00_bootstrap_pctci.do` | Helper programs for percentile-CI tables (reference only) | – |
| `01_baseline_ppml.do` | Baseline PPML | Fig. 4a |
| `01b_controls_robustness.do` | Baseline with control sets added step by step | SI |
| `02_iv_cf_main.do` | Gravity IV and Bartik IV, control-function PPML | Fig. 4a |
| `03_counterfactual_matrix.do` | Counterfactual reassignment across centrality quartiles (4×4 matrices) | Fig. 4b |
| `04_bilateral_multilateral.do`, `04_dv_means.do` | Bilateral vs multilateral adaptation finance (marginal effects, channel means) | Main text |
| `04b_bilateral_multilateral_relative.do` | Bilateral vs multilateral, effects relative to each channel's mean | Fig. 5a–b |
| `05_income_heterogeneity.do` | Heterogeneity by income group | Fig. 5c–d |
| `06_moderation_disaster.do`, `06c_moderation_disaster_curve.do` | Moderation by disaster exposure | Fig. 6a–b, SI |
| `07_moderation_idealpoint.do`, `07c_moderation_idealpoint_curve.do` | Moderation by political alignment with donors | Fig. 6c–d, SI |
| `07b_moderation_idealpoint_robustness.do` | Political alignment: competing interactions and moderator timing | SI |
| `08_robustness.do` | 5-year moving average, winsorised variables, lagged controls, excluding COVID-19 years, contemporaneous IV | SI |

## Figure notebooks (`2_figures/`)

| Notebook | Figure |
|---|---|
| `Fig1a_climate_finance_trend` | Fig. 1a: trends in international climate finance |
| `Fig1b_climate_finance_sankey` | Fig. 1b: flows of climate finance by region |
| `Fig1cd_climate_finance_map` | Fig. 1c–d: geographical distribution of adaptation finance |
| `Fig1_combined` | Fig. 1, combined layout |
| `Fig2_trade_network_centrality` | Fig. 2: countries' positions in the global trade network |
| `Fig3_finance_vs_centrality` | Fig. 3: adaptation finance and trade-network centrality |
| `Fig4_baseline_and_counterfactual` | Fig. 4: baseline, IV and counterfactual estimates |
| `Fig5_bilateral_multilateral_income` | Fig. 5: bilateral vs multilateral finance; income groups |
| `Fig6ab_moderation_disaster`, `Fig6cd_moderation_idealpoint` | Fig. 6a–b and 6c–d as separate files |
| `Fig6_combined` | Fig. 6: marginal effects across disaster exposure and political alignment |
| `SuppFig_centrality_maps` | SI: maps of standardised centralities, 2000/2010/2020/2024 |
| `SuppFig_centrality_by_income` | SI: centralities by income group, 2024 |
| `SuppFig_moderation_disaster` | SI: marginal effects by disaster-exposure percentile |
| `SuppFig_moderation_idealpoint` | SI: marginal effects by political-alignment percentile |
| `SuppFig_robustness` | SI: robustness checks compared with the main IV estimates |

