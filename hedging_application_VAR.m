%% ============================================================
% HEDGING APPLICATION basata sullo stress test VAR (calibrato)
% - Usa gli scenari di prezzo WTI generati da stress_test_hormuz_VAR.m
% - Supply: scenario calibrato tipo Hormuz (shockSigma_supply_calib)
% - Demand: shock domanda kσ (da VAR)
%% ============================================================

clear; clc; close all;

%% 0) Stile grafico uniforme (per figure "da tesi")
set(groot, ...
    'defaultFigureColor'   , 'w', ...
    'defaultAxesColor'     , 'w', ...
    'defaultAxesXColor'    , 'k', ...
    'defaultAxesYColor'    , 'k', ...
    'defaultAxesGridColor' , [0.85 0.85 0.85], ...
    'defaultAxesFontName'  , 'Helvetica', ...
    'defaultAxesFontSize'  , 10, ...
    'defaultLineLineWidth' , 1.5);

%% 1) Carica risultati dello stress test
if ~exist('stress_hormuz_VAR_results.mat','file')
    error('stress_hormuz_VAR_results.mat mancante. Esegui prima stress_test_hormuz_VAR.m');
end

load('stress_hormuz_VAR_results.mat', ...
     'hvec', 'P0_USD', ...
     'P_baseline_USD', ...
     'P_supply_USD', 'P_demand_USD', ...
     'shockSigma_supply_calib', 'shockSigma_demand', ...
     'sigma_d', 'target_increase_supply');

H = numel(hvec);           % orizzonte (mesi)

% Alias per leggibilità
shockSigma_supply = shockSigma_supply_calib;

%% 2) Parametri azienda di trasporto
Q_year_barrels  = 3.5e6;           % consumo annuo di carburante (barili)
Q_month_barrels = Q_year_barrels / 12;

% Prezzo del future/swap che l'azienda può bloccare
F0_USD = P0_USD;                   % assumiamo forward ~ prezzo corrente

% Hedge ratio: quota del fabbisogno coperto (0 = nessuna copertura, 1 = full hedge)
h_hedge = 1.0;                     % copertura al 100% 

%% 3) Costo carburante non coperto (baseline + scenari)

% Baseline: prezzo piatto a P0_USD
Cost_baseline = Q_month_barrels * P_baseline_USD;   % H x 1

% Scenario supply (calibrato tipo Hormuz)
Cost_supply_unhedged = Q_month_barrels * P_supply_USD;

% Scenario demand (shockSigma_demand σ)
Cost_demand_unhedged = Q_month_barrels * P_demand_USD;

%% 4) Costo con hedging (futures/swap lineare sul WTI)
% Costo_t = Q [(1-h)*P_t + h*F0]

% Supply scenario (calibrato)
Cost_supply_hedged = Q_month_barrels * ( (1 - h_hedge)*P_supply_USD + h_hedge*F0_USD );

% Demand scenario
Cost_demand_hedged = Q_month_barrels * ( (1 - h_hedge)*P_demand_USD + h_hedge*F0_USD );

%% 5) Costo cumulato (12 mesi)

Tot_baseline        = sum(Cost_baseline);
Tot_supply_unhedged = sum(Cost_supply_unhedged);
Tot_supply_hedged   = sum(Cost_supply_hedged);
Tot_demand_unhedged = sum(Cost_demand_unhedged);
Tot_demand_hedged   = sum(Cost_demand_hedged);

Saving_supply = Tot_supply_unhedged - Tot_supply_hedged;
Saving_demand = Tot_demand_unhedged - Tot_demand_hedged;

fprintf('\n=== RISULTATI ANNUALI (12 mesi) ===\n');
fprintf('Consumo annuo Q = %.0f barili\n', Q_year_barrels);
fprintf('Prezzo iniziale / future F0 = %.2f USD\n\n', F0_USD);

fprintf('Baseline (nessuno shock):                                  %12.2f USD\n', Tot_baseline);

fprintf('Supply shock calibrato (%.1fσ, target ~+%d%%%% al picco) - senza hedge: %12.2f USD\n', ...
        shockSigma_supply, round(100*target_increase_supply), Tot_supply_unhedged);
fprintf('Supply shock calibrato (%.1fσ) - con hedge:                         %12.2f USD\n', ...
        shockSigma_supply, Tot_supply_hedged);
fprintf('   → Risparmio da hedge (supply):                          %12.2f USD\n\n', Saving_supply);

fprintf('Domanda forte %.1fσ - senza hedge:                         %12.2f USD\n', ...
        shockSigma_demand, Tot_demand_unhedged);
fprintf('Domanda forte %.1fσ - con hedge:                           %12.2f USD\n', ...
        shockSigma_demand, Tot_demand_hedged);
fprintf('   → Risparmio da hedge (demand):                          %12.2f USD\n\n', Saving_demand);

%% 6) Grafici: supply scenario (calibrato)

figure('Color','w');
tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

colBaseline = [0.2 0.2 0.2];   % grigio scuro
colSupplyP  = [0.85 0.33 0.10];% arancio per il PREZZO
colSupplyC  = [0 0.45 0.74];   % blu per il COSTO

% 6.1 Prezzo WTI (baseline vs scenario supply)
nexttile;
hBaseS = plot(hvec, P_baseline_USD, '--', ...
    'LineWidth',1.3, 'Color',colBaseline); hold on;
hSupP  = plot(hvec, P_supply_USD, '-', ...
    'LineWidth',1.8, 'Color',colSupplyP);
grid on; box on;
xlabel('Horizon (months)');
ylabel('WTI price (USD)');
title(sprintf('Calibrated supply shock (%.1f\\sigma) – WTI price', shockSigma_supply), ...
      'Color',[0.1 0.1 0.1]);

lgSupTop = legend([hBaseS hSupP], ...
    {'Baseline price','Supply scenario price'}, ...
    'Location','northwest');
set(lgSupTop,'Box','off','Color','none','TextColor',[0.1 0.1 0.1]);

% 6.2 Costo carburante (unhedged vs hedged)
nexttile;
hSupCostU = plot(hvec, Cost_supply_unhedged/1e6, '-', ...
    'LineWidth',1.6, 'Color',colSupplyC); hold on;
hSupCostH = plot(hvec, Cost_supply_hedged/1e6, '--', ...
    'LineWidth',1.6, 'Color',colBaseline);
grid on; box on;
xlabel('Horizon (months)');
ylabel('Fuel cost (USD millions)');
title(sprintf('Fuel cost – supply shock (%.1f\\sigma, h = %.0f%%%%)', ...
      shockSigma_supply, 100*h_hedge), ...
      'Color',[0.1 0.1 0.1]);

lgSupBot = legend([hSupCostU hSupCostH], ...
    {'Unhedged fuel cost','Hedged fuel cost'}, ...
    'Location','northwest');
set(lgSupBot,'Box','off','Color','none','TextColor',[0.1 0.1 0.1]);

%% 7) Grafici: demand scenario

figure('Color','w');
tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

colBaseline = [0.2 0.2 0.2];
colDemand   = [0 0.45 0.74];     % blu per i prezzi
colCostDem  = [0.85 0.33 0.10];  % arancio per i costi

% 7.1 Prezzo WTI
nexttile;
hBaseD = plot(hvec, P_baseline_USD, '--', 'LineWidth',1.3, 'Color',colBaseline); hold on;
hDem   = plot(hvec, P_demand_USD,   '-',  'LineWidth',1.8, 'Color',colDemand);
grid on; box on;
xlabel('Horizon (months)');
ylabel('WTI price (USD)');
title(sprintf('Demand shock (%.1f\\sigma) – WTI price', shockSigma_demand), ...
      'Color',[0.1 0.1 0.1]);

lgTop = legend([hBaseD hDem], {'Baseline price','Demand scenario price'}, ...
               'Location','northwest');
set(lgTop,'Box','off','Color','none','TextColor',[0.1 0.1 0.1]);

% 7.2 Costo carburante con e senza hedge
nexttile;
hCostDemU = plot(hvec, Cost_demand_unhedged/1e6, '-',  'LineWidth',1.6, 'Color',colCostDem); hold on;
hCostDemH = plot(hvec, Cost_demand_hedged/1e6,   '--', 'LineWidth',1.6, 'Color',colBaseline);
grid on; box on;
xlabel('Horizon (months)');
ylabel('Fuel cost (USD millions)');
title(sprintf('Fuel cost – demand shock (%.1f\\sigma, h = %.0f%%%%)', ...
      shockSigma_demand, 100*h_hedge), 'Color',[0.1 0.1 0.1]);

lgBot = legend([hCostDemU hCostDemH], ...
               {'Unhedged fuel cost','Hedged fuel cost'}, ...
               'Location','northwest');
set(lgBot,'Box','off','Color','none','TextColor',[0.1 0.1 0.1]);
