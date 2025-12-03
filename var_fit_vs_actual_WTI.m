%% VAR – FIT WTI vs PREZZO OSSERVATO - meglio non inserirlo, fa ridere (o piangere)
% Confronto tra prezzo reale WTI osservato e prezzo implicito dal VAR
clear; clc;

%% ---------- 0. STILE GRAFICO UNIFORME ----------
set(groot, ...
    'defaultFigureColor'   , 'w', ...
    'defaultAxesColor'     , 'w', ...
    'defaultAxesXColor'    , 'k', ...
    'defaultAxesYColor'    , 'k', ...
    'defaultAxesGridColor' , [0.9 0.9 0.9], ...
    'defaultAxesFontName'  , 'Helvetica', ...
    'defaultAxesFontSize'  , 10, ...
    'defaultLineLineWidth' , 1.8);

%% ---------- 1. CARICA DATI E VAR ----------
if ~exist('clean_data.mat','file')
    error('clean_data.mat mancante. Esegui prima build_oil_dataset.m');
end
if ~exist('var_svar_results.mat','file')
    error('var_svar_results.mat mancante. Esegui prima var_main.m');
end

load('clean_data.mat','All','All_d');                   % All: livelli, All_d: differenze z-score
load('var_svar_results.mat','Mdl_final','varNames','pChosen');

% Matrice Y (T x n) nello stesso ordine del VAR
TT = All_d(:, varNames);
TT = rmmissing(TT);
Y  = TT{:,:};
[T, n] = size(Y);

fprintf('Osservazioni utili (All_d): %d, numero variabili: %d\n', T, n);
disp('Ordine VAR (All_d):');
disp(varNames(:));

% Indice della variabile WTI (differenze log reali z-score)
idxWTI = find(contains(varNames,'WTI','IgnoreCase',true), 1);
if isempty(idxWTI)
    error('Variabile WTI_real_DL non trovata in varNames.');
end

%% ---------- 2. RESIDUI E VALORI FITTATI DEL VAR ----------
% Residui ridotta forma: U ha T_eff = T - pChosen righe
[U, ~] = infer(Mdl_final, Y);        % T_eff x n
[T_eff, k] = size(U);

if k ~= n
    error('Dimensione di U non coerente con Y.');
end

% Allineo Y agli ultimi T_eff mesi (quelli per cui ho i residui)
Y_eff = Y(end-T_eff+1:end, :);       % T_eff x n

% Valori "previsti" dal VAR (in termini di z-score delle differenze)
Y_hat = Y_eff - U;                   % T_eff x n

% Componente WTI in termini di z-score di Δlog(WTI_real)
z_wti_hat = Y_hat(:, idxWTI);        % T_eff x 1

%% ---------- 3. DA z-score A Δlog(WTI_real) "REALI" ----------
% Serie dei livelli reali del WTI
wti_real = All.WTI_real;             % livello (indice) del prezzo reale
if any(isnan(wti_real))
    wti_real = fillmissing(wti_real,'previous');
end

% Per coerenza temporale: tengo solo l'ultimo segmento coerente con Y_eff
% Y_eff ha T_eff osservazioni; mi servono T_eff+1 livelli per i ritorni
if numel(wti_real) < T_eff + 1
    error('Serie WTI_real troppo corta rispetto a T_eff.');
end

wti_real_eff = wti_real(end-(T_eff+1)+1 : end);   % T_eff+1 livelli
log_wti_eff  = log(wti_real_eff);
wti_d_raw    = diff(log_wti_eff);                 % Δlog(WTI_real), T_eff x 1

mu_wti_raw    = mean(wti_d_raw,'omitnan');
sigma_wti_raw = std(wti_d_raw,'omitnan');

% Inversione dello z-score:
% WTI_real_DL = (Δlog(WTI_real) - mu) / sigma
% ⇒ Δlog_hat = mu + sigma * z_hat
dlog_wti_hat = mu_wti_raw + sigma_wti_raw * z_wti_hat;   % T_eff x 1

%% ---------- 4. RICOSTRUZIONE PREZZO FITTATO ----------
% Base: primo livello del segmento effettivo
P0 = wti_real_eff(1);                      % livello iniziale per la ricostruzione

logP_hat = log(P0) + cumsum(dlog_wti_hat); % T_eff x 1
P_hat    = exp(logP_hat);                  % prezzo "fittato" T_eff x 1

% Prezzo osservato corrispondente (stesse date interne)
P_obs = wti_real_eff(2:end);              % anche questo T_eff x 1

if numel(P_obs) ~= T_eff || numel(P_hat) ~= T_eff
    warning('Lunghezze non allineate. P_obs=%d, P_hat=%d, T_eff=%d', ...
            numel(P_obs), numel(P_hat), T_eff);
end

% Asse temporale: uso le ultime T_eff row times del timetable All_d
time_d_full = All_d.Properties.RowTimes;     % più robusto di .Time
time_d      = time_d_full(end-T_eff+1:end);  % datetime, T_eff x 1

%% ---------- 5. PLOT CONFRONTO PREZZO OSSERVATO VS FITTATO ----------
figure('Name','WTI reale – Prezzo osservato vs prezzo implicito dal VAR', ...
       'Color','w','Position',[100 100 950 550]);

plot(time_d, P_obs, 'k', 'LineWidth', 2); hold on;
plot(time_d, P_hat, 'r--', 'LineWidth', 1.8);

grid on; box on;
xlabel('Tempo','Color','k');
ylabel('Prezzo reale WTI (indice)','Color','k');
title('WTI reale – prezzo osservato vs prezzo "fittato" dal VAR', ...
      'FontWeight','bold','Color','k');

lgd = legend({'Prezzo osservato','Prezzo fittato dal VAR'}, ...
             'Location','best','Box','off');
set(lgd,'TextColor','k','FontWeight','bold');

%% ---------- 6. METRICHE DI ADATTAMENTO (LOG-PREZZI) ----------
logP_obs = log(P_obs);
logP_fit = log(P_hat);

rmse_log = sqrt(mean((logP_fit - logP_obs).^2,'omitnan'));
mae_log  = mean(abs(logP_fit - logP_obs),'omitnan');

fprintf('\n== ADATTAMENTO VAR SU LOG-PREZZO REALE WTI ==\n');
fprintf('RMSE (log-prezzo) = %.4f\n', rmse_log);
fprintf('MAE  (log-prezzo) = %.4f\n', mae_log);

%% ---------- 7. METRICHE DI ADATTAMENTO (LIVELLI) ----------
% Qui usiamo i livelli del WTI reale: P_obs (osservato) e P_hat (fittato)

y     = P_obs;    % WTI reale osservato
y_hat = P_hat;    % WTI reale fittato dal VAR

mask  = ~isnan(y) & ~isnan(y_hat);
y     = y(mask);
y_hat = y_hat(mask);

err   = y_hat - y;

RMSE = sqrt(mean(err.^2, 'omitnan'));          % Root Mean Squared Error
MAE  = mean(abs(err), 'omitnan');              % Mean Absolute Error
MAPE = mean(abs(err ./ y), 'omitnan') * 100;   % Mean Absolute Percentage Error

fprintf('\n== ADATTAMENTO VAR SU PREZZO REALE WTI (LIVELLI) ==\n');
fprintf('RMSE = %.4f\n', RMSE);
fprintf('MAE  = %.4f\n', MAE);
fprintf('MAPE = %.2f %s\n', MAPE, '%%');

save('var_wti_fit_metrics.mat','RMSE','MAE','MAPE');
