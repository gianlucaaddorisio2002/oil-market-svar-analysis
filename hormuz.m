%% STRESS TEST HORMUZ – Prezzo reale WTI in USD (da IRF strutturali calibrate)
clear; clc;

%% 0) Stile grafico uniforme
set(groot, ...
    'defaultFigureColor'   , 'w', ...
    'defaultAxesColor'     , 'w', ...
    'defaultAxesXColor'    , 'k', ...
    'defaultAxesYColor'    , 'k', ...
    'defaultAxesGridColor' , [0.9 0.9 0.9], ...
    'defaultAxesFontName'  , 'Helvetica', ...
    'defaultAxesFontSize'  , 10, ...
    'defaultLineLineWidth' , 1.8);

%% 1) Carica IRF strutturali e dati in livelli
if ~exist('svar_sign_results.mat','file')
    error('svar_sign_results.mat mancante. Esegui prima svar_sign_restrictions.m');
end
if ~exist('clean_data.mat','file')
    error('clean_data.mat mancante. Esegui prima build_oil_dataset.m');
end

load('svar_sign_results.mat', ...
     'IRF_struct','idxProd','idxWTI','horizon');

load('clean_data.mat','All');    % contiene WTI_real in livelli reali (indice)

%% 2) Impostazioni orizzonte
H = min(24, size(IRF_struct,1));   % max 24 mesi o quanto disponibile
t = (0:H-1)';

idxSupply = idxProd;   % shock di offerta (structural supply shock)

%% 3) IRF del WTI in z-score per shock di offerta
irf_WTI_supply_z = squeeze(IRF_struct(1:H, idxWTI, idxSupply));   % H x 1

%% 4) Dal z-score alle variazioni di log-prezzo reale

% Serie del prezzo reale WTI (indice)
wti_real = All.WTI_real;
wti_real = fillmissing(wti_real,'previous');

% Δlog prezzo reale (coerente con All_d)
dlog_wti = diff(log(wti_real));
sigma_d  = std(dlog_wti, 'omitnan');   % deviazione standard del log-return

% === Calibrazione shock Hormuz: voglio ~+30% al picco IRF ===
target_increase_supply = 0.30;             % +30%
target_log_change      = log(1 + target_increase_supply);

% Orizzonte dove la IRF del WTI è massima (risposta a shock di offerta)
[~, h_peak] = max(irf_WTI_supply_z(1:H));
irf_peak    = irf_WTI_supply_z(h_peak);

% Fattore in sigma necessario:
%   ΔlogP_target = sigma_d * shockSigma_calib * irf_peak
shockSigma_calib = target_log_change / (sigma_d * irf_peak);

fprintf('[Stress test Hormuz] Supply shock calibrato: %.1f sigma (target +%d%%%% al picco, h=%d)\n', ...
    shockSigma_calib, round(100*target_increase_supply), h_peak);

% Variazione di log-prezzo generata dallo shock di offerta calibrato:
%   Δlog P_t^shock = sigma_d * (shockSigma_calib * IRF_z_t)
dlog_effect = sigma_d * shockSigma_calib * irf_WTI_supply_z;   % H x 1

% Effetto cumulato sul log-prezzo
cum_dlog_effect = cumsum(dlog_effect);                         % H x 1

%% 5) Dal log-prezzo al prezzo in USD

% Prezzo reale osservato alla fine del campione (indice)
P0_index = wti_real(end);

% Target di partenza "realistico" in USD (es. 60 USD al barile)
P0_USD = 60;
scale  = P0_USD / P0_index;

% Baseline: nessuno shock → prezzo piatto a P0_USD
P_baseline_USD = P0_USD * ones(H,1);

% Scenario Hormuz: applico solo l'effetto dello shock (nessun trend)
P_scenario_index = P0_index * exp(cum_dlog_effect);   % livelli in indice
P_scenario_USD   = P_scenario_index * scale;          % convertito in USD

%% 6) Variazione percentuale dopo 12 mesi
h_12  = min(12, H-1);      % per sicurezza
ret_12 = (P_scenario_USD(h_12+1)/P_baseline_USD(h_12+1) - 1) * 100;

fprintf('Variazione percentuale del WTI dopo %d mesi: %.2f%%%%\n', ...
        h_12, ret_12);

%% ===========================
%   FIGURA – Stress Test Hormuz (full black theme)
% ============================

% --- Parametri modificabili ---
offset_label_x      = 0.6;
offset_label_y      = 0.5;
line_width_base     = 1.6;
line_width_scenario = 2.2;
font_size           = 11;

figure('Name','Stress test Hormuz – Prezzo reale WTI in USD', ...
       'Color','w','Position',[100 100 900 550]);

hold on; grid on; box on;

% --- Imposto lo stile generale tutto nero ---
ax = gca;
ax.XAxis.Color     = 'k';
ax.YAxis.Color     = 'k';
ax.GridColor       = [0.7 0.7 0.7];        % griglia più leggera
ax.MinorGridColor  = [0.7 0.7 0.7];
ax.FontSize        = font_size;
ax.XColor          = 'k';
ax.YColor          = 'k';

% --- Linea baseline (nera tratteggiata) ---
baseline_plot = plot(t, P_baseline_USD, '--', 'Color','k', ...
    'LineWidth', line_width_base);

% --- Linea scenario shock (nera solida) ---
scenario_plot = plot(t, P_scenario_USD, 'k-', ...
    'LineWidth', line_width_scenario);

% --- Annotazione percentuale dopo 12 mesi ---
text(h_12 + offset_label_x, ...
     P_scenario_USD(h_12+1) + offset_label_y, ...
     sprintf('%.1f%% dopo %d mesi', ret_12, h_12), ...
     'FontWeight','bold', ...
     'HorizontalAlignment','left', ...
     'VerticalAlignment','bottom', ...
     'FontSize', font_size, ...
     'Color','k');

% --- Assi ---
xlabel('Orizzonte (mesi)', 'FontSize', font_size, 'Color','k');
ylabel('Prezzo reale WTI (USD)', 'FontSize', font_size, 'Color','k');

% --- Titolo (tutto nero) ---
title({sprintf('Stress test – Calibrated supply shock (%.1f\\sigma)', shockSigma_calib), ...
       '(scenario Hormuz, livelli in USD)'}, ...
       'FontWeight','bold', 'FontSize', font_size+1, 'Color','k');

% --- Legenda nero su nero ---
legend([baseline_plot scenario_plot], ...
    {'Baseline (nessuno shock)', ...
     sprintf('Scenario Hormuz (supply shock %.1f\\sigma)', shockSigma_calib)}, ...
     'Location','best', ...
     'Box','off', ...
     'FontSize', font_size, ...
     'TextColor','k');

% --- Range ---
ylim([min(P_scenario_USD)-1, max(P_scenario_USD)+1]);
xlim([0, H-1]);