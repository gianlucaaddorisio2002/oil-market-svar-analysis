%% ============================================================
% PLOT_OLS_DIAGNOSTICS.M
% Figura 2.4 – OLS diagnostics (Ljung-Box, Breusch-Pagan, Jarque-Bera)
%% ============================================================

clear; clc;

if ~exist('clean_data.mat','file')
    error('clean_data.mat non trovato. Esegui prima build_oil_dataset.m');
end

load('clean_data.mat','All');

tbl = timetable2table(All(:,{'WTI_real','Production','OCSE','Inventories'}),...
                      'ConvertRowTimes',false);

Y = tbl.WTI_real;
X = [ones(size(Y)), tbl.Production, tbl.OCSE, tbl.Inventories];
n = size(X,1);

% OLS
beta = X\Y;
e    = Y - X*beta;

% --- Ljung-Box su residui (serial correlation) ---
[~, p_lb] = lbqtest(e,'Lags',12);

% --- Breusch-Pagan per eteroschedasticità ---
% Regressione di e^2 su regressori
e2 = e.^2;
Z  = X;  % stessa matrice regressori
beta_bp = Z\e2;
e2_hat  = Z*beta_bp;
SSR     = sum( (e2 - e2_hat).^2 );
SST     = sum( (e2 - mean(e2)).^2 );
R2      = 1 - SSR/SST;
LM      = n * R2;
df      = size(Z,2)-1;  % escludi intercetta
p_bp    = 1 - chi2cdf(LM, df);

% --- Jarque-Bera per normalità ---
[~, p_jb] = jbtest(e);

pvals = [p_lb; p_bp; p_jb];
tests = {'Ljung-Box (12)','Breusch-Pagan','Jarque-Bera'};

outdir = fullfile('figures','chapter2');
if ~exist(outdir,'dir'); mkdir(outdir); end

figure('Position',[100 100 800 500]);
bar(pvals);
set(gca,'XTick',1:3,'XTickLabel',tests,'XTickLabelRotation',15);
yline(0.05,'r--','LineWidth',1.2);
ylim([0 1]);
ylabel('p-value');
title('Diagnostic tests for static OLS regression of real WTI');

print('-dpdf', fullfile(outdir,'fig_ols_diag.pdf'));
