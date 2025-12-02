Structural Oil Market Analysis with VAR/SVAR, Copulas and Stress Testing

This repository contains the full, reproducible pipeline developed for my bachelor thesis “Structural Identification of Oil Supply, Global Demand and Shocks: A Factor Analysis with Industrial Applications”.

The project builds a four-variable VAR/SVAR model of the global crude oil market (1990–2024), identified through sign and elasticity restrictions (Kilian–Murphy 2014), and extends it with copula-based dependence modelling, Monte Carlo scenario generation, and an industrial fuel-risk application for airlines.

- **Overview**

The objective is to disentangle the structural drivers of oil price dynamics — flow supply, aggregate demand, and precautionary (inventory) demand shocks — and map them into risk-relevant scenario paths for decision-making in energy-intensive industries (e.g., aviation).

The repository includes:

Fully reproducible data construction pipeline

VAR(12) estimation on stationary transformations

SVAR identification via sign restrictions

Impulse response functions (IRFs)

Forecast error variance decomposition (FEVD)

Historical decomposition

Marginal distribution fitting of structural shocks (Gaussian, t-Student, Logistic)

Bivariate copulas (Gaussian, t, Clayton, Gumbel, Frank)

Tail-dependence analysis

Monte Carlo simulation engine (bootstrap residuals)

Hormuz-type extreme supply-shock stress test

WTI → Jet Fuel pass-through model

Airline hedging application (ITA Airways case)

Every component is implemented in MATLAB, following a transparent, modular structure.

- **Repository Structure**
├── data/                    # Raw and processed data (EIA, FRED, OECD)
├── scripts/
│   ├── build_oil_dataset.m
│   ├── ols_prevar_check.m
│   ├── var_main.m
│   ├── svar_sign_restrictions.m
│   ├── plot_irf_fevd.m
│   ├── fit_marginal_shocks.m
│   ├── copula_bivariate_shocks.m
│   ├── scenario_mc_WTI.m
│   ├── stress_test_hormuz_VAR.m
│   ├── extract_wti_irf_from_var.m
│   ├── estimate_pass_through_jetfuel.m
│   ├── hedging_application_VAR.m
│   ├── ita_hedging_full.m
│   └── utils/
├── results/
│   ├── IRF/
│   ├── FEVD/
│   ├── Copulas/
│   ├── Scenarios/
│   ├── Hedging/
│   └── Figures/
├── thesis/                  # Full LaTeX manuscript
├── README.md
└── LICENSE

- **Methodology**
1. Reduced-form VAR

A VAR(12) is estimated on stationary transformations of:

Δ log U.S. Crude Oil Production

Δ OECD Industrial Activity

Δ log Real WTI

Δ log U.S. Petroleum Inventories

Diagnostics confirm:

no residual autocorrelation

eigenvalues inside unit circle

mild non-normality (heavy tails)

acceptable stability

2. Structural Identification (SVAR)

Identification follows:

sign restrictions

elasticity bounds (Kilian–Murphy, 2014)

Recovering:

Flow Supply Shock

Aggregate Demand Shock

Precautionary/Storage Shock

3. Impulse Responses & FEVD

Key finding:

Global aggregate demand is the dominant driver of medium-run oil price movements,
while supply shocks are short-lived and inventory shocks matter mainly in extremes.

4. Distribution of Structural Shocks

Structural shocks exhibit:

strong leptokurtosis

skewness (demand shocks negative, precautionary positive)

marked tail dependence (especially Demand ↔ Precautionary)

t-Student marginals consistently outperform Gaussian fits.

5. Copula Dependence

Copulas reveal:

strong upper-tail dependence between Demand and Precautionary shocks

weak dependence between Supply and other shocks

moderate dependence between Supply and Demand

6. Scenario Generation & Stress Testing

Two engines:

Monte Carlo (bootstrap residuals)

Extreme scenario (Hormuz-type supply disruption)

WTI paths are transformed to levels and scaled to USD.

7. Jet Fuel Mapping & Hedging

A log–log pass-through model:

log(JF) = α + β log(WTI_real)
β ≈ 1.14


→ Linearised around 60 USD for computational use.

Applied to ITA Airways:

annual fuel demand ≈ 3.5 million bbl

full linear hedge at 200 USD/bbl

**Critical result**:
A full hedge reduces exposure by 30–40% under a Hormuz-like shock.

- **Key Outputs**
Structural Analysis

IRFs for supply, demand, precautionary shocks

FEVD of real WTI

Historical decomposition

Non-Gaussian shock distributions

Copula surfaces and contour plots

Risk Engineering

Monte Carlo WTI distributions

Price paths for extreme supply shocks

Jet fuel price trajectories

ITA Airways monthly/annual fuel bills

Hedging P&L comparison

- **How to Run**

Requirements:

MATLAB R2022+

Statistics Toolbox

Econometrics Toolbox

To rebuild the full pipeline:

run('scripts/build_oil_dataset.m')
run('scripts/var_main.m')
run('scripts/svar_sign_restrictions.m')
run('scripts/fit_marginal_shocks.m')
run('scripts/copula_bivariate_shocks.m')
run('scripts/scenario_mc_WTI.m')
run('scripts/stress_test_hormuz_VAR.m')
run('scripts/estimate_pass_through_jetfuel.m')
run('scripts/ita_hedging_full.m')

- **Citation**

If you use this material, please cite:

Addorisio, G.A.J. (2025). Structural Identification of Oil Supply, Global Demand and Shocks: A Factor Analysis with Industrial Applications. Università Campus Bio-Medico di Roma.

- **License**

MIT License — feel free to use, modify and build upon this work.

- **Contact**

For questions or collaboration:

Gianluca Addorisio
📧 gianluca.addorisio@gmail.com


🔗 https://github.com/gianlucaaddorisio20
