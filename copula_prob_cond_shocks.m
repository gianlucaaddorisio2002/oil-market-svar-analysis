%% COPULA – PROBABILITÀ CONDIZIONATE TRA SHOCK
clear; clc;

%% 1) Carica dati (shock + marginali)
load('structural_shocks.mat', ...
     'shocks_matrix','shockLabels','distFits');

fprintf('Shock disponibili:\n');
disp(shockLabels(:));

% Indici:
% 1 = Supply_shock
% 2 = AggregateDemand_OCSE_shock
% 3 = Precautionary_shock
% 4 = Residual_shock

% Scegli la coppia: qui usiamo AggregateDemand_OCSE vs Precautionary
i = 2;   % X = AggregateDemand_OCSE_shock
j = 3;   % Y = Precautionary_shock

X = shocks_matrix(:,i);
Y = shocks_matrix(:,j);

margX = distFits{i}{2};
margY = distFits{j}{2};

fprintf('\nCoppia scelta: %s (X) vs %s (Y)\n', ...
    shockLabels{i}, shockLabels{j});

%% 2) Trasforma in pseudo-osservazioni U(0,1)
Ux = cdf(margX, X);
Uy = cdf(margY, Y);

epsU = 1e-6;
Ux = min(max(Ux,epsU),1-epsU);
Uy = min(max(Uy,epsU),1-epsU);

U = [Ux Uy];

%% 3) Fitta la t-copula direttamente QUI
[R, nu] = copulafit('t', U);
fprintf('Copula t stimata: nu = %.2f\n', nu);

%% 4) Definisci le soglie (es. 90-esimo percentile = shock "molto alto")
tauX = 0.90;
tauY = 0.90;

x0 = quantile(X, tauX);
y0 = quantile(Y, tauY);

u0 = cdf(margX, x0);
v0 = cdf(margY, y0);

u0 = min(max(u0,epsU),1-epsU);
v0 = min(max(v0,epsU),1-epsU);

uPair = [u0 v0];

%% 5) Valuta la cdf congiunta C(u0,v0)
Cuv = copulacdf('t', uPair, R, nu);

%% 6) Probabilità condizionata
% P(Y > y0 | X > x0) = [1 - u0 - v0 + C(u0,v0)] / (1 - u0)

num  = 1 - u0 - v0 + Cuv;
den  = 1 - u0;
Pcond   = num / den;
Pincond = 1 - v0;

fprintf('\nSoglie: tauX = %.2f, tauY = %.2f\n', tauX, tauY);
fprintf('P(Y > y0) (incond.)       = %.4f\n', Pincond);
fprintf('P(Y > y0 | X > x0)        = %.4f\n', Pcond);
fprintf('Rapporto Pcond / Pincond  = %.2f\n', Pcond/Pincond);
