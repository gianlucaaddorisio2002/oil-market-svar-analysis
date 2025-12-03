%% HORMUZ_MIX_SCENARIO – Stress test con mix di shock
%   - Scenario 1: solo shock di offerta (3 sigma)
%   - Scenario 2: shock di offerta (3 sigma) + shock precauzionale (2 sigma)
%
%   Usa:
%     - svar_sign_results.mat  (IRF_struct, indici variabili)
%     - clean_data.mat         (WTI_real in livelli reali, indice)

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
     'IRF_struct','idxProd','idxWTI','idxInv','horizon');

load('clean_data.mat','All');    % contiene WTI_real in livelli reali (indice)

%% 2) Impostazioni orizzonte e parametri di scenario
H = min(24, size(IRF_struct,1));   % max 24 mesi o quanto disponibile
t = (0:H-1)';

% Indici shock
idxSupply = idxProd;   % shock di offerta
idxPrec   = idxInv;    % shock precauzionale (storage / inventari)

% Ampiezza shock in sigma
shockSigma_supply = +3;     % supply shock avverso 3 sigma
shockSigma_prec   = +2;     % precautionary shock 2 sigma (puoi cambiare)

fprintf('[Hormuz mix] Supply = %+g sigma, Precautionary = %+g sigma.\n', ...
        shockSigma_supply, shockSigma_prec);

%% 3) IRF del WTI in z-score per i due shock
irf_WTI_supply_z = squeeze(IRF_struct(1:H, idxWTI, idxSupply));   % H x 1
irf_WTI_prec_z   = squeeze(IRF_struct(1:H, idxWTI, idxPrec));     % H x 1

%% 4) Dal z-score alle variazioni di log-prezzo reale

% Serie del prezzo reale WTI (indice)
wti_real = All.WTI_real;
wti_real = fillmissing(wti_real,'previous');

% Δlog prezzo reale (coerente con All_d)
dlog_wti = diff(log(wti_real));
sigma_d  = std(dlog_wti, 'omitnan');   % deviazione standard del ritorno log

% Effetto singolo shock di offerta
dlog_effect_supply = sigma_d * shockSigma_supply * irf_WTI_supply_z;   % H x 1

% Effetto shock precauzionale
dlog_effect_prec   = sigma_d * shockSigma_prec   * irf_WTI_prec_z;     % H x 1

% Effetto combinato (lineare)
dlog_effect_mix    = dlog_effect_supply + dlog_effect_prec;

% Effetti cumulati
cum_dlog_supply = cumsum(dlog_effect_supply);
cum_dlog_mix    = cumsum(dlog_effect_mix);

%% 5) Dal log-prezzo al prezzo in USD

% Prezzo reale osservato alla fine del campione (indice)
P0_index = wti_real(end);

% Target di partenza "realistico" in USD (es. 70 USD al barile)
P0_USD   = 70;
scale    = P0_USD / P0_index;

% Baseline: nessuno shock → prezzo piatto a P0_USD
P_baseline_USD = P0_USD * ones(H,1);

% Scenario 1: solo shock di offerta
P_supply_index = P0_index * exp(cum_dlog_supply);
P_supply_USD   = P_supply_index * scale;

% Scenario 2: mix supply + precautionary
P_mix_index = P0_index * exp(cum_dlog_mix);
P_mix_USD   = P_mix_index * scale;

%% 6) Variazione percentuale dopo 12 mesi per i due scenari
h_12 = min(12, H-1);      % per sicurezza

ret12_supply = (P_supply_USD(h_12+1)/P_baseline_USD(h_12+1) - 1) * 100;
ret12_mix    = (P_mix_USD(h_12+1)/P_baseline_USD(h_12+1)    - 1) * 100;

fprintf('Variazione dopo %d mesi:\n', h_12);
fprintf('  - Solo supply:        %.2f%%%%\n', ret12_supply);
fprintf('  - Supply + precaution:%.2f%%%%\n', ret12_mix);

%% 7) FIGURA – Stress Test Hormuz: solo supply vs mix

offset_label_x = 0.6;
offset_label_y = 0.5;
line_width_base     = 1.4;
line_width_supply   = 1.8;
line_width_mix      = 2.0;
font_size = 11;

figure('Name','Stress test Hormuz – Mix di shock', ...
       'Color','w','Position',[100 100 950 550]);

hold on; grid on; box on;

ax = gca;
ax.XAxis.Color  = 'k';
ax.YAxis.Color  = 'k';
ax.GridColor    = [0.7 0.7 0.7];
ax.MinorGridColor = [0.7 0.7 0.7];
ax.FontSize     = font_size;

% Baseline (tratteggiata)
h_base = plot(t, P_baseline_USD, '--', 'Color','k', ...
    'LineWidth', line_width_base);

% Solo supply (linea nera)
h_supply = plot(t, P_supply_USD, '-', 'Color','k', ...
    'LineWidth', line_width_supply);

% Mix supply + precautionary (linea nera più spessa)
h_mix = plot(t, P_mix_USD, '-', 'Color','k', ...
    'LineWidth', line_width_mix);

% Annotazioni
text(h_12 + offset_label_x, P_supply_USD(h_12+1) + offset_label_y, ...
     sprintf('Solo supply: %.1f%%%% dopo %d mesi', ret12_supply, h_12), ...
     'FontWeight','bold','HorizontalAlignment','left', ...
     'VerticalAlignment','bottom','FontSize', font_size, 'Color','k');

text(h_12 + offset_label_x, P_mix_USD(h_12+1) + 1.0 + offset_label_y, ...
     sprintf('Mix supply+prec: %.1f%%%% dopo %d mesi', ret12_mix, h_12), ...
     'FontWeight','bold','HorizontalAlignment','left', ...
     'VerticalAlignment','bottom','FontSize', font_size, 'Color','k');

% Assi e titolo
xlabel('Orizzonte (mesi)', 'FontSize', font_size, 'Color','k');
ylabel('Prezzo reale WTI (USD)', 'FontSize', font_size, 'Color','k');

title({'Stress test – Scenario Hormuz: solo offerta vs mix', ...
       '(livelli in USD, IRF strutturali)'}, ...
      'FontWeight','bold', 'FontSize', font_size+1, 'Color','k');

% Legenda
legend([h_base h_supply h_mix], ...
    {'Baseline (nessuno shock)', ...
     sprintf('Supply %.0f\\sigma', shockSigma_supply), ...
     sprintf('Supply %.0f\\sigma + Prec %.0f\\sigma', ...
             shockSigma_supply, shockSigma_prec)}, ...
     'Location','best', 'Box','off', ...
     'FontSize', font_size, 'TextColor','k');

% --- Range verticale: prendi min/max su tutte le serie, scalari ---
ymin = min([P_baseline_USD; P_supply_USD; P_mix_USD], [], 'all', 'omitnan');
ymax = max([P_baseline_USD; P_supply_USD; P_mix_USD], [], 'all', 'omitnan');

ylim([ymin - 1, ymax + 1]);
xlim([0, H-1]);


