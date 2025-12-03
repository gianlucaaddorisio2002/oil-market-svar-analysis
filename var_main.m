%% VAR + SVAR – VAR su serie stazionarie (All_d, OCSE_DL)
clear; clc;

%% 1) Carica dati stazionari
% All_d: timetable con almeno le variabili:
%   Production_DL, OCSE_DL, WTI_real_DL, Inventories_DL
load('clean_data.mat','All_d');

%% 2) Prepara dati per il VAR
% ORDINE ECONOMICO (da replicare nello SVAR):
%   1) Production_DL  -> offerta
%   2) OCSE_DL        -> domanda globale
%   3) WTI_real_DL    -> prezzo reale WTI
%   4) Inventories_DL -> scorte (precautionary/storage)
TT = All_d(:, {
    'Production_DL', ...  % 1) offerta
    'OCSE_DL', ...        % 2) ciclo OCSE (domanda globale)
    'WTI_real_DL', ...    % 3) prezzo reale WTI in differenze log
    'Inventories_DL' ...  % 4) scorte in differenze log
});

TT = rmmissing(TT);   % sicurezza

Y            = TT{:,:};                 % T × 4
[nObs,nVars] = size(Y);
varNames     = TT.Properties.VariableNames;

%% 3) Loop empirico VAR(p) per p = 1…pMax

pMax = 15;

results = table( ...
    'Size'         , [pMax 3], ...
    'VariableTypes', {'double','logical','double'}, ...
    'VariableNames', {'p','Stable','LB_min_pval'} ...
);

for p = 1:pMax

    fprintf('\n==== VAR(%d) su All_d (differenze/zscore) ====\n', p);

    % Stima VAR(p)
    Mdl = varm(nVars, p);
    [EstMdl,~,~,E] = estimate(Mdl, Y);

    % --- Stabilità ---
    ARcell = [{eye(nVars)} EstMdl.AR];
    L      = LagOp(ARcell);
    L      = reflect(L);
    stab   = isStable(L);

    % --- Ljung–Box: p-value minimo sui residui di ciascuna eq. ---
    pvals = zeros(nVars,1);
    for j = 1:nVars
        [~, p_tmp] = lbqtest(E(:,j), 'Lags', 1:12);
        pvals(j) = min(p_tmp);
    end
    lb_min = min(pvals);

    % Salva risultati
    results.p(p)           = p;
    results.Stable(p)      = stab;
    results.LB_min_pval(p) = lb_min;

    fprintf('Stable = %d   min p(LB) = %.4f\n', stab, lb_min);
end

fprintf('\n=== RIEPILOGO EMPIRICO VAR (All_d) ===\n');
disp(results);

%% 4) Scelta dell'ordine p definitivo in modo automatico
idx = find(results.Stable & results.LB_min_pval > 0.05, 1, 'first');
if isempty(idx)
    % nessun modello con LB>0.05 → prendi quello con LB_min_pval massimo
    [~, idx] = max(results.LB_min_pval);
    fprintf('Nessun VAR(p) con LB_min_pval>0.05. Scelgo p con LB_min_pval massimo.\n');
end

pChosen = results.p(idx);
fprintf('\n*** VAR(%d) scelto come modello finale ***\n', pChosen);

%% 5) Stima VAR(pChosen) definitivo + diagnostica

Mdl_final = varm(nVars, pChosen);
[Mdl_final,~,~,E_final] = estimate(Mdl_final, Y);

% --- Stabilità finale ---
ARcell_final = [{eye(nVars)} Mdl_final.AR];
L_final      = LagOp(ARcell_final);
L_final      = reflect(L_final);

if isStable(L_final)
    fprintf('VAR(%d) stabile.\n', pChosen);
else
    warning('VAR(%d) NON stabile. Rivalutare p.', pChosen);
end

% Scatter degli autovalori della matrice companion
p = pChosen;
k = nVars;

% Costruisci matrice companion A
ARmats = Mdl_final.AR;              % cell array {A1,...,Ap}, ciascuna k×k
A1 = [ARmats{:}];                   % riga superiore k × (k*p)

if p > 1
    Ak = [eye(k*(p-1)) zeros(k*(p-1), k)];  % blocco inferiore
    Acomp = [A1; Ak];                       % matrice companion (k*p) × (k*p)
else
    Acomp = ARmats{1};
end

lambda = eig(Acomp);

figure('Name','Autovalori VAR','Color','w');
scatter(real(lambda), imag(lambda), 40, 'filled'); hold on;

% cerchio unitario
theta = linspace(0, 2*pi, 400);
plot(cos(theta), sin(theta), 'k--','LineWidth',1.2);

xlabel('Re(\lambda)','Color','k');
ylabel('Im(\lambda)','Color','k');
title(sprintf('Autovalori matrice companion VAR(%d)', p), 'Color','k');
axis equal; grid on;

% --- Ljung–Box residui finali (tutte le equazioni) ---
pvals_final = zeros(nVars,1);
for j = 1:nVars
    [~, p_tmp] = lbqtest(E_final(:,j),'Lags',1:12);
    pvals_final(j) = min(p_tmp);
end

fprintf('Min p-value Ljung–Box residui VAR(%d): %.4f\n', ...
    pChosen, min(pvals_final));

%% 5a) Equation-by-equation fit statistics (VAR(pChosen))
%   - Number of parameters
%   - Std. error of regression
%   - R^2 e Adjusted R^2
%   → per la tabella "Reduced-form VAR(p) equation fit statistics"

T      = nObs;                 % numero di osservazioni totali
T_eff  = T - pChosen;          % osservazioni effettive dopo i ritardi
k_par  = 1 + pChosen*nVars;    % costante + pChosen*lags*4 variabili

fitStats = table( ...
    'Size'         , [nVars 4], ...
    'VariableTypes', {'double','double','double','double'}, ...
    'VariableNames', {'NumParams','StdError','R2','AdjR2'}, ...
    'RowNames'     , varNames);

for i = 1:nVars

    % Serie dipendente allineata sui ritardi
    yi = Y(pChosen+1:end, i);

    % Residui corrispondenti
    ei = E_final(pChosen+1:end, i);

    % Somma dei quadrati residui e totale
    SSR = sum(ei.^2);
    TSS = sum( (yi - mean(yi)).^2 );

    % Std. error of regression
    stderr_i = sqrt(SSR / (T_eff - k_par));

    % R^2 e Adj. R^2
    R2_i    = 1 - SSR/TSS;
    AdjR2_i = 1 - (1 - R2_i) * (T_eff - 1) / (T_eff - k_par);

    % Scrivi nella tabella
    fitStats.NumParams(i) = k_par;
    fitStats.StdError(i)  = stderr_i;
    fitStats.R2(i)        = R2_i;
    fitStats.AdjR2(i)     = AdjR2_i;
end

fprintf('\n== Equation-by-equation VAR(%d) fit statistics ==\n', pChosen);
disp(fitStats);

% Salva su .mat per usarla in LaTeX
save('var_equation_fit_stats.mat','fitStats');

%% 5bis) Analisi distribuzione residui VAR – equazione WTI_real_DL

% Nuovo ordine (prova): Production_DL, OCSE_DL, WTI_real_DL, Inventories_DL
% → WTI_real_DL è SEMPRE la 3ª variabile (fissa)
e_var = E_final(:,3);              % residui VAR per WTI_real_DL
e_var = e_var(~isnan(e_var));      % pulizia

% Standardizzazione
e_var_std = (e_var - mean(e_var)) / std(e_var);

% Skewness e kurtosis
skew_var      = skewness(e_var_std);
kurt_var      = kurtosis(e_var_std);
excess_kurt_v = kurt_var - 3;

fprintf('\n== MOMENTI DI ORDINE SUPERIORE (RESIDUI VAR – WTI_real_DL) ==\n');
fprintf('Skewness     = %.4f\n', skew_var);
fprintf('Kurtosis     = %.4f (Normale ≈ 3)\n', kurt_var);
fprintf('Excess kurt. = %.4f (Normale = 0)\n', excess_kurt_v);

% Istogramma + normale standard
figure('Name','Residui VAR WTI_real_DL','Color','w');
histogram(e_var_std, 'Normalization','pdf', 'NumBins',30, ...
          'FaceAlpha',0.5, 'EdgeColor',[0.2 0.2 0.2]);
hold on;
xgrid = linspace(min(e_var_std), max(e_var_std), 400);
plot(xgrid, normpdf(xgrid,0,1), 'k','LineWidth',1.8);
xlabel('Residui VAR standardizzati (WTI\_real\_DL)', 'Color', 'k');
ylabel('Densità', 'Color', 'k');
title('Distribuzione residui VAR WTI\_real\_DL vs Normale N(0,1)', 'Color', 'k');
grid on; box on;

% Test di normalità
[p_jb_var, ~] = jbtest(e_var);         % H0: normale
[~, p_ks_var] = kstest(e_var_std);     % H0: N(0,1)

fprintf('\n== TEST NORMALITÀ RESIDUI VAR (WTI_real_DL) ==\n');
fprintf('JbTest: p = %.4f\n', p_jb_var);
fprintf('KsTest: p = %.4f\n', p_ks_var);

% QQ-plot vs Normale standard
figure('Name','QQ-plot residui VAR WTI_real_DL','Color','w');
qqplot(e_var_std);
title('QQ-plot residui VAR WTI\_real\_DL vs Normale N(0,1)','Color','k');
xlabel('Quantili teorici N(0,1)','Color','k');
ylabel('Quantili empirici residui','Color','k');
grid on;

%% ACF + PACF residui VAR (Paper-Ready)

maxLag = 20;

% Calcolo numerico ACF/PACF
[acf_vals,  lags_acf,  bounds_acf]  = autocorr(e_var,  'NumLags', maxLag);
[pacf_vals, lags_pacf, bounds_pacf] = parcorr(e_var, 'NumLags', maxLag);

% Ricava banda simmetrica (95%)
ci_acf  = bounds_acf(2);
ci_pacf = bounds_pacf(2);

%% Figure layout
fig = figure('Name','ACF/PACF Residuals VAR – WTI_real_DL',...
             'Color','w','Units','normalized','Position',[0.1 0.1 0.75 0.55]);

tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

%% ================== ACF ==================
nexttile; hold on;

% Banda di confidenza
patch([-0.5 maxLag+0.5 maxLag+0.5 -0.5],...
      [-ci_acf -ci_acf ci_acf ci_acf],...
      [0.92 0.92 0.92],'EdgeColor','none');

% Zero line
plot([-0.5 maxLag+0.5],[0 0],'k-','LineWidth',1);

% Spike + marker
for i = 1:length(lags_acf)
    plot([lags_acf(i) lags_acf(i)], [0 acf_vals(i)], 'k-', 'LineWidth',1.1);
    plot(lags_acf(i), acf_vals(i), 'ko', 'MarkerSize',4, 'MarkerFaceColor','k');
end

% Scaling dinamico
upper = max(max(acf_vals), ci_acf) * 1.10;
lower = min(min(acf_vals), -ci_acf) * +3.25;

ylim([lower upper]);
xlim([-0.5 maxLag+0.5]);

title('ACF residui VAR – WTI\_real\_DL','Color','k');
xlabel('Lag','Color','k');
ylabel('Autocorrelation','Color','k');

set(gca,'Box','on','FontSize',11,'XColor','k','YColor','k');
grid on;
hold off;

%% ================== PACF ==================
nexttile; hold on;

patch([-0.5 maxLag+0.5 maxLag+0.5 -0.5],...
      [-ci_pacf -ci_pacf ci_pacf ci_pacf],...
      [0.92 0.92 0.92],'EdgeColor','none');

plot([-0.5 maxLag+0.5],[0 0],'k-','LineWidth',1);

for i = 1:length(lags_pacf)
    plot([lags_pacf(i) lags_pacf(i)], [0 pacf_vals(i)], 'k-', 'LineWidth',1.1);
    plot(lags_pacf(i), pacf_vals(i), 'ko', 'MarkerSize',4, 'MarkerFaceColor','k');
end

upper = max(max(pacf_vals), ci_pacf) * 1.10;
lower = min(min(pacf_vals), -ci_pacf) * 3.25;

ylim([lower upper]);
xlim([-0.5 maxLag+0.5]);

title('PACF residui VAR – WTI\_real\_DL','Color','k');
xlabel('Lag','Color','k');
ylabel('Partial autocorrelation','Color','k');

set(gca,'Box','on','FontSize',11,'XColor','k','YColor','k');
grid on;
hold off;

%% 6) IRF + FEVD (Cholesky) – ridotta forma su All_d

horizon = 30;

% IRF ortogonalizzate (Cholesky) dal modello stimato
IRF  = irf(Mdl_final, 'NumObs', horizon);

% FEVD
FEVD = fevd(Mdl_final, 'NumObs', horizon);

% --- Pre-estrai costante e coefficienti AR per velocizzare il bootstrap ---
p = pChosen;
k = nVars;

% Costante come riga 1×k
if isempty(Mdl_final.Constant)
    Const = zeros(1,k);
else
    Const = Mdl_final.Constant(:)';   % 1×k
end

% Coefficienti AR
ARcoeff = cell(p,1);
for lag = 1:p
    ARcoeff{lag} = Mdl_final.AR{lag};    % ciascuno k×k
end

%% 7) Bootstrap delle IRF (residual bootstrap alla Kilian)

numBoot = 50;        % valore variabile
T       = nObs;

% Residui "utili": salto i primi p
E_pool = E_final(p+1:end, :);
T_eff  = size(E_pool,1);

% Condizioni iniziali: primi p valori osservati
Y0 = Y(1:p, :);

% Contenitore IRF bootstrap: (time × risposta × shock × bootstrap)
IRF_boot = zeros(horizon, k, k, numBoot);

fprintf('\n== BOOTSTRAP IRF (residual bootstrap) ==\n');
for b = 1:numBoot

    % 1) ricampiona residui
    idx = randi(T_eff, T_eff, 1);
    E_b = E_pool(idx, :);

    % 2) simula serie VAR
    Y_sim = zeros(T_eff + p, k);
    Y_sim(1:p, :) = Y0;

    for t = p+1 : T_eff + p
        mu = Const;
        for lag = 1:p
            mu = mu + Y_sim(t-lag,:) * ARcoeff{lag}';  % (1×k)*(k×k)
        end
        Y_sim(t,:) = mu + E_b(t-p,:);
    end

    Y_sim_final = Y_sim(p+1:end, :);

    % 3) ristima VAR
    Mdl_b = varm(k, p);
    [Mdl_b,~,~,~] = estimate(Mdl_b, Y_sim_final, 'Display','off');

    % 4) IRF ortogonalizzate
    IRF_b = irf(Mdl_b, 'NumObs', horizon);
    IRF_boot(:,:,:,b) = IRF_b;

    if mod(b,10)==0
        fprintf('  Bootstrap rep %d / %d completata\n', b, numBoot);
    end
end

% 5) Mediana + bande percentile 5–95%
IRF_med     = median(IRF_boot, 4);
IRF_CI_low  = prctile(IRF_boot, 5,  4);
IRF_CI_high = prctile(IRF_boot, 95, 4);

fprintf('Bootstrap IRF completato (%d repliche).\n', numBoot);

%% 8) Salvataggio risultati

save('var_svar_results.mat', ...
     'Mdl_final', ...      % modello VAR stimato
     'IRF', ...            % IRF singola stima (Cholesky)
     'FEVD', ...           % decomposizione varianza
     'varNames', ...       % {'Production_DL','OCSE_DL','WTI_real_DL','Inventories_DL'}
     'pChosen', ...        % lag scelto
     'IRF_boot', ...       % IRF bootstrap complete
     'IRF_med', ...        % mediana bootstrap
     'IRF_CI_low', ...     % banda inferiore 5%%
     'IRF_CI_high');       % banda superiore 95%%

fprintf('\nRISULTATI SALVATI in var_svar_results.mat\n');
