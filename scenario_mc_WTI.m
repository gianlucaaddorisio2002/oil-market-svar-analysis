%% SCENARIO MONTE CARLO VAR – PERCORSI FUTURI DEL PREZZO REALE WTI
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

%% ---------- 1. PARAMETRI MONTE CARLO ----------
H     = 12;      % orizzonte in mesi
Nsim  = 5000;    % numero di scenari simulati
alpha = 0.10;    % livello per bande (10%–90%)

fprintf('Simulazioni Monte Carlo VAR: N = %d, orizzonte = %d mesi\n', Nsim, H);

%% ---------- 2. CARICA DATI E VAR FINALE ----------
if ~exist('clean_data.mat','file')
    error('clean_data.mat mancante. Esegui prima build_oil_dataset.m');
end
if ~exist('var_svar_results.mat','file')
    error('var_svar_results.mat mancante. Esegui prima var_main.m');
end

load('clean_data.mat','All','All_d');              % All: livelli, All_d: differenze z-score
load('var_svar_results.mat','Mdl_final','varNames','pChosen');

% Costruisci Y (T x n) nello stesso ordine usato nel VAR
TT = All_d(:, varNames);
TT = rmmissing(TT);
Y  = TT{:,:};
[T, n] = size(Y);

fprintf('Osservazioni utili (All_d): %d, numero variabili: %d\n', T, n);
disp('Ordine VAR (All_d):');
disp(varNames(:));

% Trova indice del WTI nelle variabili del VAR
idxWTI = find(contains(varNames,'WTI','IgnoreCase',true), 1);
if isempty(idxWTI)
    error('Variabile WTI_real_DL non trovata in varNames.');
end

%% ---------- 3. RESIDUI DEL VAR E PARAMETRI PER LA SIMULAZIONE ----------
[U, ~] = infer(Mdl_final, Y);      % T_eff x n, con T_eff = T - pChosen
[T_eff, k] = size(U);

fprintf('Osservazioni utili per residui VAR: %d (T - p), numero variabili: %d\n', T_eff, k);

Sigma_u = cov(U, 1);               %#ok<NASGU> % (se vuoi usarla dopo)

p = pChosen;

if isempty(Mdl_final.Constant)
    Const = zeros(1, k);
else
    Const = Mdl_final.Constant(:)';       % 1 x k
end

ARcoeff = cell(p,1);
for lag = 1:p
    ARcoeff{lag} = Mdl_final.AR{lag};     % ciascuno k x k
end

% Condizioni iniziali: ultimi p valori osservati di Y
if T <= p
    error('Serie troppo corta rispetto al lag p scelto.');
end
Y0 = Y(T-p+1:T, :);     % p x k

%% ---------- 4. MOMENTI GREZZI Δlog(WTI) PER RIENTRARE IN LIVELLI ----------
wti_real = All.WTI_real;   % livello prezzo reale (indice base)
if any(isnan(wti_real))
    wti_real = fillmissing(wti_real,'previous');
end

log_wti   = log(wti_real);
wti_d_raw = diff(log_wti);                % Δlog(WTI_real)

mu_wti_raw    = mean(wti_d_raw,'omitnan');
sigma_wti_raw = std(wti_d_raw,'omitnan');

% Prezzo reale di partenza: ultimo osservato (in unità "indice")
P0_index = wti_real(end);
fprintf('Prezzo reale WTI di partenza (indice) P0 = %.2f\n', P0_index);

% ---- RISCALE A UN PREZZO REALISTICO (ES. 60 USD/BARILE) ----
P_target = 60;                          % scegli tu il livello di partenza in USD
scale_factor = P_target / P0_index;     % riscalatura lineare
fprintf('Riscalo i livelli in modo che il WTI parta da %.2f USD.\n', P_target);

%% ---------- 5. SIMULAZIONI MONTE CARLO ----------
rng(123);   % per riproducibilità

% Contenitori:
Z_WTI_sim = zeros(H, Nsim);   % VAR su z-score della Δlog WTI
P_sim_idx = zeros(H, Nsim);   % livelli di prezzo (in indice base)

fprintf('Inizio simulazioni bootstrap sui residui del VAR...\n');

for s = 1:Nsim

    % 5.1 Inizializza cammino VAR con gli ultimi p valori osservati
    Y_sim = zeros(H + p, k);
    Y_sim(1:p, :) = Y0;

    % 5.2 Simula orizzonte H con residual bootstrap
    for t = p+1 : H + p
        
        % bootstrap residuo ridotta forma
        idx_r = randi(T_eff);          % 1..T_eff
        u_t   = U(idx_r, :);           % 1 x k

        % parte deterministica (VAR)
        mu = Const;
        for lag = 1:p
            mu = mu + Y_sim(t-lag, :) * ARcoeff{lag}';
        end

        Y_sim(t, :) = mu + u_t;
    end

    % 5.3 Estrai solo il tratto futuro (H step) e la componente WTI
    Y_future  = Y_sim(p+1:end, :);           % H x k
    z_wti_sim = Y_future(:, idxWTI);         % H x 1, z-score Δlog(WTI_real)

    Z_WTI_sim(:, s) = z_wti_sim;

    % 5.4 Converti z_t in Δlog(WTI_real) "reali"
    dlog_wti_sim = mu_wti_raw + sigma_wti_raw * z_wti_sim;  % H x 1

    % 5.5 Ricostruisci il livello di prezzo reale WTI (in indice base)
    logP_path = log(P0_index) + cumsum(dlog_wti_sim);
    P_path_idx = exp(logP_path);

    P_sim_idx(:, s) = P_path_idx;
    
    if mod(s, 1000) == 0
        fprintf('  Simulazioni completate: %d / %d\n', s, Nsim);
    end
end

fprintf('Simulazioni Monte Carlo completate.\n');

%% ---------- 6. RISCALE I LIVELLI A USD ----------
% Ora trasformiamo l'indice in "USD equivalenti"
P_sim = P_sim_idx * scale_factor;   % H x Nsim
P0    = P_target;                   % prezzo iniziale in USD

%% ---------- 7. FAN CHART DEL PREZZO REALE WTI ----------
q05 = quantile(P_sim, alpha/2,     2);  % 2.5%
q50 = quantile(P_sim, 0.5,         2);  % mediana
q95 = quantile(P_sim, 1-alpha/2,   2);  % 97.5%

t = (1:H)';

figure('Name','Scenario Monte Carlo – Prezzo reale WTI', ...
       'Color','w', 'Position',[100 100 900 550]);
hold on;

% Banda (fan) 10–90%
fill([t; flipud(t)], [q05; flipud(q95)], ...
     [0.85 0.90 1.0], 'EdgeColor','none', 'FaceAlpha',0.6);
plot(t, q50, 'k-', 'LineWidth', 2);                    % mediana
yline(P0, '--', 'LineWidth',1.3, 'Color',[0.3 0.3 0.3]);  % livello iniziale

grid on; box on;
xlabel('Orizzonte (mesi)', 'Color','k');
ylabel('Prezzo reale WTI (USD)', 'Color','k');
title('Monte Carlo VAR – Scenari futuri del prezzo reale WTI', ...
      'FontWeight','bold','Color','k');

lgd = legend({'Banda 10–90%', 'Mediana', 'Prezzo iniziale'}, ...
             'Location','best','Box','off');
set(lgd,'TextColor','k','FontWeight','bold');

%% ---------- 8. PROBABILITÀ DI EVENTI ESTREMI ----------
% Variazione cumulata a H mesi in termini di prezzo (in USD)
P_H   = P_sim(end, :);
ret_H = P_H / P0 - 1;         % rendimento % rispetto a P0 (invariato dallo scaling)

th_down = -0.30;
th_up   =  0.30;

prob_crash = mean(ret_H <= th_down);
prob_rally = mean(ret_H >= th_up);

fprintf('\n== PROBABILITÀ EVENTI ESTREMI (orizzonte %d mesi) ==\n', H);
fprintf('P[Prezzo WTI scende di almeno 30%%] ≈ %.2f%%%%\n', 100*prob_crash);
fprintf('P[Prezzo WTI sale di almeno 30%%]   ≈ %.2f%%%%\n', 100*prob_rally);

%% ---------- 9. ESPORTA RISULTATI ----------
results_MC = struct();
results_MC.H         = H;
results_MC.Nsim      = Nsim;
results_MC.P0_index  = P0_index;
results_MC.P0_USD    = P0;
results_MC.P_sim_USD = P_sim;
results_MC.q05       = q05;
results_MC.q50       = q50;
results_MC.q95       = q95;
results_MC.ret_H     = ret_H;
results_MC.probCrash = prob_crash;
results_MC.probRally = prob_rally;
results_MC.scale     = scale_factor;

save('scenario_mc_WTI_results.mat','results_MC');

fprintf('\nRisultati Monte Carlo salvati in scenario_mc_WTI_results.mat\n');