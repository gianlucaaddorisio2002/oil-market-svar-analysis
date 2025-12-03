%% EXTRACT STRUCTURAL SHOCKS (SVAR con restrizioni di segno, OCSE)
clear; clc;

%% ----- Stile grafico globale -----
set(groot, ...
    'defaultFigureColor','w', ...
    'defaultAxesColor','w', ...
    'defaultAxesXColor','k', ...
    'defaultAxesYColor','k', ...
    'defaultAxesGridColor',[0.85 0.85 0.85], ...
    'defaultAxesFontName','Helvetica', ...
    'defaultAxesFontWeight','normal', ...
    'defaultAxesFontSize',12, ...
    'defaultLineLineWidth',2, ...
    'defaultLegendTextColor','k');

%% 1) Carica dati necessari
load('clean_data.mat','All_d');
load('var_svar_results.mat','Mdl_final','varNames');
load('svar_sign_results.mat','acceptedQ');

fprintf('Estraggo gli shock strutturali dal VAR/SVAR...\n');

%% 2) Costruisci Y con stesso ordine del VAR
TT = All_d(:, varNames);
TT = rmmissing(TT);
Y  = TT{:,:};
[T, n] = size(Y);

fprintf('Osservazioni utili: %d, variabili: %d\n', T, n);

%% 3) Residui ridotta forma u_t
[U, ~] = infer(Mdl_final, Y);

%% 4) Covarianza e Cholesky
Sigma_u = cov(U,1);
P = chol(Sigma_u,'lower');

%% 5) Scegli una rotazione accettata
Q_star = acceptedQ(:,:,1);
B = P * Q_star;

%% 6) Shock strutturali eps_t = B^{-1} u_t
Eps = (B \ U')';

%% 7) Nomina shock
shockNames = { ...
    'Supply_shock', ...
    'AggregateDemand_OCSE_shock', ...
    'Precautionary_shock', ...
    'Residual_shock'};

%% 8) Salva output in UNICO file coerente per gli script successivi
shocks_matrix = Eps;
shockLabels   = shockNames;

save('structural_shocks.mat',...
    'Eps','shocks_matrix','shockLabels','varNames',...
    'B','Q_star','Sigma_u');

fprintf('>> Shock strutturali salvati in structural_shocks.mat\n');
%% 9) Descriptive statistics of structural shocks (per Tabella 3.4)

% Eps: T x 4 (Supply, Aggregate Demand, Precautionary, Residual)
% shockLabels: cell array con i nomi

% Se vuoi solo i primi 3 shock "economici":
Eps_use   = Eps(:,1:3);
shockCore = shockLabels(1:3);

stats_shock = table( ...
    'Size',[3 4], ...
    'VariableTypes',{'double','double','double','double'}, ...
    'VariableNames',{'Mean','StdDev','Skewness','Kurtosis'}, ...
    'RowNames',shockCore);

for i = 1:3
    x = Eps_use(:,i);
    x = x(~isnan(x));

    stats_shock.Mean(i)     = mean(x);
    stats_shock.StdDev(i)   = std(x);
    stats_shock.Skewness(i) = skewness(x);
    stats_shock.Kurtosis(i) = kurtosis(x);  % NON excess kurtosis
end

fprintf('\n== Descriptive statistics of structural shocks ==\n');
disp(stats_shock);

save('shock_descriptive_stats.mat','stats_shock','shockCore');
fprintf('>> Statistiche shock salvate in shock_descriptive_stats.mat\n');
