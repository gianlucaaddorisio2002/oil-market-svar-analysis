%% ============================================================
% STRESS TEST HORMUZ basato SOLO sul VAR (Cholesky)
% - Scenario supply shock 3σ
% - Scenario domanda forte 2σ (opzionale)
%% ============================================================

clear; clc; close all;

%% 1) Carica risultati VAR e dataset
if ~exist('var_svar_results.mat','file')
    error('var_svar_results.mat mancante. Esegui prima var_main.m');
end
if ~exist('clean_data.mat','file')
    error('clean_data.mat mancante. Esegui prima build_oil_dataset.m');
end

load('var_svar_results.mat', 'IRF', 'varNames');   % IRF: h x n x n
load('clean_data.mat', 'All');                     % contiene WTI_real

%% 2) Calcola sigma_d dei log-returns WTI (prima dello z-score)

% Ricostruisco i log-returns del WTI reale (come in build_oil_dataset)
wti_levels = All.WTI_real;               % prezzo reale in livello
wti_d      = diff(log(wti_levels));      % log-differenze
sigma_d    = std(wti_d, 'omitnan');      % deviazione standard

fprintf('Deviazione standard log-returns WTI (sigma_d): %.4f\n', sigma_d);

%% 3) Indici delle variabili nel VAR
idxProd = find(strcmp(varNames,'Production_DL'));
idxMacro = find(ismember(varNames, {'OCSE_DL','REA_DL'}));
idxWTI  = find(strcmp(varNames,'WTI_real_DL'));
idxInv  = find(strcmp(varNames,'Inventories_DL'));

if any(isempty([idxProd, idxMacro, idxWTI, idxInv]))
    error('Controlla che varNames contenga Production_DL, OCSE/REA_DL, WTI_real_DL, Inventories_DL');
end

[h, n, ~] = size(IRF);
hScenario = min(12, h-1);     % orizzonte 12 mesi (se possibile)
hvec      = (1:hScenario)';   % partiamo da 1 (mese successivo allo shock)
target_increase_supply = 0.30;   % +30%
%% 4) IRF del WTI in z-score (1 s.d. shock)

% 4.1 Supply: voglio che UNO shock di offerta NEGATIVO faccia SALIRE il WTI
% Calcolo la media della risposta del WTI a uno shock POSITIVO di Production
irf_raw_supply = IRF(1:hScenario, idxWTI, idxProd);
mean12 = mean(irf_raw_supply(1:min(12,hScenario)));

% Se un aumento di Production fa scendere il prezzo (media < 0),
% allora un calo di Production = -shock → segnoSupply = -1 (come prima).
% Se invece media > 0 (VAR "sbagliato" dal punto di vista economico),
% allora definisco direttamente lo shock negativo come +IRF.
if mean12 < 0
    signSupply = -1;   % ribalta: calo Production → prezzo su
else
    signSupply = +1;   % non ribaltare: il VAR dice già che shock negativo fa salire il prezzo
end

irf_wti_supply_z = signSupply * irf_raw_supply;


% 4.2 Demand: aumento della variabile Macro = "domanda forte"
% Qui NON ribaltiamo il segno: uno shock positivo di domanda
% fa salire il WTI se il modello è sensato.
irf_wti_demand_z =  IRF(1:hScenario, idxWTI, idxMacro);

%% 5) Parametri di scenario (statistico vs calibrato)

P0_USD = 60;   % livello iniziale di riferimento

% 5.1 Shock "statistico" (quello vecchio, solo per confronto teorico)
shockSigma_supply_stat = 3;    % 3σ sulla produzione
shockSigma_demand      = 2;    % 2σ sulla domanda (teniamo così)

% 5.2 Shock "calibrato" tipo Hormuz: voglio +30% al picco IRF del WTI
target_increase_supply = 0.30;                       % +30%
target_log_change      = log(1 + target_increase_supply);

% Trovo l'orizzonte h dove la IRF del WTI (supply) è massima
[~, h_peak_supply] = max(irf_wti_supply_z(1:hScenario));
irf_peak_supply    = irf_wti_supply_z(h_peak_supply);

% Fattore k (in sigma) necessario per ottenere quel +30% al picco
shockSigma_supply_calib = target_log_change / (irf_peak_supply * sigma_d);

fprintf('Supply shock statistico: %.1fσ\n', shockSigma_supply_stat);
fprintf('Supply shock calibrato tipo Hormuz: %.1fσ (target +%d% al picco, h=%d)\n', ...
    shockSigma_supply_calib, round(100*target_increase_supply), h_peak_supply);

%% 6) Da IRF (z) a Δlog prezzo reale

% 6.1 Supply "statistico" 3σ (solo per eventuale confronto)
dlogP_supply_stat = irf_wti_supply_z * (sigma_d * shockSigma_supply_stat);   % h x 1

% 6.2 Supply "calibrato" tipo Hormuz (questo è quello che userai per l'hedging)
dlogP_supply = irf_wti_supply_z * (sigma_d * shockSigma_supply_calib);       % h x 1

% 6.3 Demand scenario 2σ (resta come prima)
dlogP_demand = irf_wti_demand_z * (sigma_d * shockSigma_demand);            % h x 1

%% 7) Da Δlog a percorso di prezzo in USD

% Accumulo nel tempo (log-livelli rispetto a t=0)
cumlog_supply_stat = cumsum(dlogP_supply_stat);
cumlog_supply_cal  = cumsum(dlogP_supply);
cumlog_demand      = cumsum(dlogP_demand);

% Baseline piatta
P_baseline_USD = P0_USD * ones(hScenario,1);

% Scenari
P_supply_3sigma_USD = P0_USD * exp(cumlog_supply_stat);  % solo per confronto
P_supply_USD        = P0_USD * exp(cumlog_supply_cal);   % scenario Hormuz calibrato
P_demand_USD        = P0_USD * exp(cumlog_demand);       % domanda forte


%% 8) Plot: confronto baseline vs scenari (versione "paper ready")

figure('Color','w');
tl = tiledlayout(2,1,'Padding','compact','TileSpacing','compact'); %#ok<NASGU>

% Colori "seri"
colBaseline = [0.2 0.2 0.2];     % grigio scuro
colSupply   = [0.85 0.33 0.10];  % arancio
colDemand   = [0 0.45 0.74];     % blu

% ---------- 8.1 Supply shock calibrato ----------
nexttile;
hBase1 = plot(hvec, P_baseline_USD, '--', ...
    'LineWidth',1.3, 'Color',colBaseline); 
hold on;
hSupply = plot(hvec, P_supply_USD, '-', ...
    'LineWidth',1.8, 'Color',colSupply);
grid on; box on;
xlabel('Horizon (months)');
ylabel('WTI price (USD)');
title(sprintf('VAR stress test – calibrated supply shock %.1f\\sigma (≈+%d%% peak)', ...
    shockSigma_supply_calib, round(100*target_increase_supply)), ...
    'Color',[0.1 0.1 0.1]);

% ---------- 8.2 Domanda forte 2σ ----------
nexttile;
hBase2 = plot(hvec, P_baseline_USD, '--', ...
    'LineWidth',1.3, 'Color',colBaseline); 
hold on;
hDemand = plot(hvec, P_demand_USD, '-', ...
    'LineWidth',1.8, 'Color',colDemand);
grid on; box on;
xlabel('Horizon (months)');
ylabel('WTI price (USD)');
title(sprintf('VAR stress test – strong demand %.1f\\sigma (OCSE/REA)', shockSigma_demand), ...
      'Color',[0.1 0.1 0.1]);

% ---------- Legenda unica sopra i pannelli ----------
ax1 = gobjects(1);   % solo per chiarezza (non strettamente necessario)
ax1 = gca;           % current axes è l'ultimo nexttile, ma non importa

lgd = legend([hBase1 hSupply hDemand], ...
    {'Baseline','Supply scenario','Demand scenario'}, ...
    'Orientation','horizontal', ...
    'Box','off', ...
    'FontSize',11);

set(lgd,'Color','none');   % sfondo trasparente
lgd.Layout.Tile = 'north'; % sposta la legenda sopra il tiledlayout

%% 9) Salva risultati per la parte di hedging
save('stress_hormuz_VAR_results.mat', ...
     'hvec','P0_USD', ...
     'P_baseline_USD', ...
     'P_supply_USD', 'P_supply_3sigma_USD', ...
     'P_demand_USD', ...
     'shockSigma_supply_stat', 'shockSigma_supply_calib', ...
     'shockSigma_demand', ...
     'target_increase_supply', 'h_peak_supply', 'sigma_d');

disp('Stress test Hormuz VAR completato e salvato in stress_hormuz_VAR_results.mat');
