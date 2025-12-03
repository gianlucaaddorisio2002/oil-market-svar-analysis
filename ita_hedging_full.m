%% ============================================================
% ITA_HEDGING_FULL
% Hedging al 100% del consumo carburante ITA Airways
% usando scenari WTI dal VAR + pass-through JetFuel_USGC
% - Supply: shock calibrato (shockSigma_supply_calib, target_increase_supply)
% - Demand: shockSigma_demand (da VAR)
%% ============================================================

clear; clc; close all;

%% 0) Stile grafico uniforme (come negli altri script)
set(groot, ...
    'defaultFigureColor'   , 'w', ...
    'defaultAxesColor'     , 'w', ...
    'defaultAxesXColor'    , 'k', ...
    'defaultAxesYColor'    , 'k', ...
    'defaultAxesGridColor' , [0.85 0.85 0.85], ...
    'defaultAxesFontName'  , 'Helvetica', ...
    'defaultAxesFontSize'  , 10, ...
    'defaultLineLineWidth' , 1.5);

%% 1) Carica scenari WTI dal VAR (stress test calibrato)
if ~exist('stress_hormuz_VAR_results.mat','file')
    error('stress_hormuz_VAR_results.mat mancante. Esegui prima stress_test_hormuz_VAR.m');
end

load('stress_hormuz_VAR_results.mat', ...
     'hvec','P0_USD', ...
     'P_baseline_USD','P_supply_USD','P_demand_USD', ...
     'shockSigma_supply_calib','shockSigma_demand', ...
     'target_increase_supply');

H = numel(hvec);

% Alias per leggibilità
shockSigma_supply = shockSigma_supply_calib;

%% 2) Carica risultati pass-through Jet Fuel
if ~exist('jetfuel_pass_results.mat','file')
    error('jetfuel_pass_results.mat mancante. Esegui prima estimate_pass_through_jetfuel.m');
end

load('jetfuel_pass_results.mat', ...
     'alpha_lin','beta_lin','JF0_hat','WTI0');

%% 3) Consumo ITA Airways (barili)

Q_year_med = 3.5e6;   % stima centrale: 3,5 milioni di barili/anno
Q_year_min = 3.2e6;   % scenario basso (opzionale)
Q_year_max = 4.0e6;   % scenario alto (opzionale)

Q_month_med = Q_year_med / 12;

fprintf('Consumo ITA (centrale): %.1f mln bbl/anno (%.0f bbl/mese)\n', ...
        Q_year_med/1e6, Q_month_med);

%% 4) Trasforma scenari WTI -> Jet Fuel (scala modello)

JF_baseline_mod = alpha_lin + beta_lin * P_baseline_USD;   % vettori h x 1
JF_supply_mod   = alpha_lin + beta_lin * P_supply_USD;
JF_demand_mod   = alpha_lin + beta_lin * P_demand_USD;

% Safety: nessun prezzo negativo (scala modello)
JF_baseline_mod = max(JF_baseline_mod, 0);
JF_supply_mod   = max(JF_supply_mod,   0);
JF_demand_mod   = max(JF_demand_mod,   0);

%% 4b) CALIBRAZIONE A LIVELLO REALISTICO (USD/bbl)

fprintf('\n[DEBUG] Jet Fuel modello a WTI = %.0f: JF0_hat = %.2f (unità modello)\n', ...
        WTI0, JF0_hat);

% Livello realistico di Jet Fuel in USD/bbl a WTI ≈ WTI0 (es. 200 USD/bbl)
JF_target = 200;   

% Fattore di scala per riportare la serie del modello in USD/bbl
scale_JF  = JF_target / JF0_hat;
fprintf('[CALIBRAZIONE] Scala Jet Fuel: scale_JF = %.4f (target %.2f USD/bbl)\n', ...
        scale_JF, JF_target);

% Applica la riscalatura a tutte le traiettorie Jet Fuel
JF_baseline = JF_baseline_mod * scale_JF;
JF_supply   = JF_supply_mod   * scale_JF;
JF_demand   = JF_demand_mod   * scale_JF;

% Safety post-scaling
JF_baseline = max(JF_baseline, 0);
JF_supply   = max(JF_supply,   0);
JF_demand   = max(JF_demand,   0);

%% 5) Prezzo forward di copertura (hedge al 100%)

F0_JF = JF_target;   % prezzo di riferimento per l'hedge
fprintf('Prezzo Jet Fuel di riferimento (WTI = %.0f): %.2f USD/bbl\n', WTI0, F0_JF);

%% 6) Costo carburante mensile (scenario centrale di consumo)

% Senza hedge
Cost_baseline_unhedged = Q_month_med * JF_baseline;   % h x 1
Cost_supply_unhedged   = Q_month_med * JF_supply;
Cost_demand_unhedged   = Q_month_med * JF_demand;

% Con hedge al 100%: costo netto = Q * F0_JF (stesso per tutti gli scenari)
Cost_baseline_hedged = Q_month_med * F0_JF * ones(H,1);
Cost_supply_hedged   = Q_month_med * F0_JF * ones(H,1);
Cost_demand_hedged   = Q_month_med * F0_JF * ones(H,1);

%% 7) Costo totale annuo e risparmio da hedge

Tot_base_unhedged   = sum(Cost_baseline_unhedged);
Tot_supply_unhedged = sum(Cost_supply_unhedged);
Tot_supply_hedged   = sum(Cost_supply_hedged);
Tot_demand_unhedged = sum(Cost_demand_unhedged);
Tot_demand_hedged   = sum(Cost_demand_hedged);

Saving_supply = Tot_supply_unhedged - Tot_supply_hedged;
Saving_demand = Tot_demand_unhedged - Tot_demand_hedged;

% → Risparmi percentuali rispetto al costo carburante SENZA copertura
perc_save_supply_vs_unhedged = 100 * Saving_supply / Tot_supply_unhedged;
perc_save_demand_vs_unhedged = 100 * Saving_demand / Tot_demand_unhedged;

% (opzionale) quanto pesa il risparmio rispetto alla fuel bill di baseline
perc_save_supply_vs_baseline = 100 * Saving_supply / Tot_base_unhedged;
perc_save_demand_vs_baseline = 100 * Saving_demand / Tot_base_unhedged;

%% 8) Stampa risultati (in milioni USD)

mill = 1e6;

fprintf('\n=== COSTO ANNUO CARBURANTE – ITA Airways (scenario centrale Q) ===\n');
fprintf('Baseline (nessuno shock)                         : %8.2f mln USD\n', Tot_base_unhedged/mill);

fprintf('Supply shock calibrato (%.1fσ, target ~+%d%%%% al picco) - senza hedge : %8.2f mln USD\n', ...
        shockSigma_supply, round(100*target_increase_supply), Tot_supply_unhedged/mill);
fprintf('Supply shock calibrato (%.1fσ) - con hedge 100%%%%              : %8.2f mln USD\n', ...
        shockSigma_supply, Tot_supply_hedged/mill);
fprintf('   → Risparmio da hedge (supply)                : %8.2f mln USD (%.1f%%%% del costo non coperto, %.1f%%%% della fuel bill di baseline)\n', ...
        Saving_supply/mill, perc_save_supply_vs_unhedged, perc_save_supply_vs_baseline);

fprintf('\nDomanda forte (%.1fσ) - senza hedge             : %8.2f mln USD\n', ...
        shockSigma_demand, Tot_demand_unhedged/mill);
fprintf('Domanda forte (%.1fσ) - con hedge 100%%%%        : %8.2f mln USD\n', ...
        shockSigma_demand, Tot_demand_hedged/mill);
fprintf('   → Risparmio da hedge (demand)               : %8.2f mln USD (%.1f%%%% del costo non coperto, %.1f%%%% della fuel bill di baseline)\n', ...
        Saving_demand/mill, perc_save_demand_vs_unhedged, perc_save_demand_vs_baseline);


%% 9) Grafici – Prezzo Jet Fuel e costo carburante (stile uniforme)

figure('Color','w');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

colBaseline = [0.2 0.2 0.2];   % grigio scuro
colSupply   = [0.85 0.33 0.10];% arancio
colDemand   = [0 0.45 0.74];   % blu

% (1) Prezzo Jet Fuel – Supply
nexttile;
hJFBaseS = plot(hvec, JF_baseline,'--','LineWidth',1.3,'Color',colBaseline); hold on;
hJFSup   = plot(hvec, JF_supply,'-','LineWidth',1.8,'Color',colSupply);
grid on; box on;
xlabel('Orizzonte (mesi)');
ylabel('Jet Fuel (USD/bbl)');
title(sprintf('Prezzo Jet Fuel – Supply shock (%.1f\\sigma)', shockSigma_supply), 'Color', [0.1 0.1 0.1]);
lg1 = legend([hJFBaseS hJFSup], {'Baseline','Scenario supply'}, 'Location','northwest');
set(lg1,'Box','off','Color','none','TextColor',[0.1 0.1 0.1]);

% (2) Costo carburante – Supply
nexttile;
hCostSupU = plot(hvec, Cost_supply_unhedged/mill,'-','LineWidth',1.8,'Color',colSupply); hold on;
hCostSupH = plot(hvec, Cost_supply_hedged/mill,'--','LineWidth',1.5,'Color',colBaseline);
grid on; box on;
xlabel('Orizzonte (mesi)');
ylabel('Costo carburante (mln USD)');
title(sprintf('ITA – Costo carburante (Supply shock, %.1f\\sigma)', shockSigma_supply), 'Color',[0.1 0.1 0.1]);
lg2 = legend([hCostSupU hCostSupH], {'Senza hedge','Con hedge 100%'}, 'Location','northwest');
set(lg2,'Box','off','Color','none','TextColor',[0.1 0.1 0.1]);

% (3) Prezzo Jet Fuel – Domanda
nexttile;
hJFBaseD = plot(hvec, JF_baseline,'--','LineWidth',1.3,'Color',colBaseline); hold on;
hJFDem   = plot(hvec, JF_demand,'-','LineWidth',1.8,'Color',colDemand);
grid on; box on;
xlabel('Orizzonte (mesi)');
ylabel('Jet Fuel (USD/bbl)');
title(sprintf('Prezzo Jet Fuel – Domanda forte (%.1f\\sigma)', shockSigma_demand), 'Color', [0.1 0.1 0.1]);
lg3 = legend([hJFBaseD hJFDem], {'Baseline','Scenario domanda'}, 'Location','northwest');
set(lg3,'Box','off','Color','none','TextColor',[0.1 0.1 0.1]);

% (4) Costo carburante – Domanda
nexttile;
hCostDemU = plot(hvec, Cost_demand_unhedged/mill,'-','LineWidth',1.8,'Color',colDemand); hold on;
hCostDemH = plot(hvec, Cost_demand_hedged/mill,'--','LineWidth',1.5,'Color',colBaseline);
grid on; box on;
xlabel('Orizzonte (mesi)');
ylabel('Costo carburante (mln USD)');
title(sprintf('ITA – Costo carburante (Domanda forte, %.1f\\sigma)', shockSigma_demand), 'Color',[0.1 0.1 0.1]);
lg4 = legend([hCostDemU hCostDemH], {'Senza hedge','Con hedge 100%'}, 'Location','northwest');
set(lg4,'Box','off','Color','none','TextColor',[0.1 0.1 0.1]);

%% 10) Salva risultati
save('ita_hedging_full_results.mat', ...
     'hvec', ...
     'Q_year_med','Q_year_min','Q_year_max', ...
     'Q_month_med', ...
     'JF_baseline','JF_supply','JF_demand', ...
     'Cost_baseline_unhedged','Cost_supply_unhedged','Cost_demand_unhedged', ...
     'Cost_baseline_hedged','Cost_supply_hedged','Cost_demand_hedged', ...
     'Tot_base_unhedged','Tot_supply_unhedged','Tot_supply_hedged', ...
     'Tot_demand_unhedged','Tot_demand_hedged', ...
     'Saving_supply','Saving_demand', ...
     'F0_JF','scale_JF','JF0_hat','JF_target', ...
     'shockSigma_supply','shockSigma_demand','target_increase_supply');

disp('Risultati salvati in ita_hedging_full_results.mat');