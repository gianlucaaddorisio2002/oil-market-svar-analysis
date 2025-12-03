%% ============================================================
% PLOT_LAG_SELECTION.M — CLEAN VERSION
%% ============================================================

clear; clc;

if ~exist('clean_data.mat','file')
    error('clean_data.mat non trovato. Esegui prima build_oil_dataset.m');
end

load('clean_data.mat','ALL_VAR');

Y = [ALL_VAR.OCSE_cycle, ALL_VAR.Prod_log, ALL_VAR.WTI_log, ALL_VAR.Inv_log];
Y = Y(~any(isnan(Y),2),:);   
[T, K] = size(Y);

pmax = 12;
AIC  = nan(pmax,1);
BIC  = nan(pmax,1);
HQIC = nan(pmax,1);

for p = 1:pmax
    Mdl = varm(K,p);
    [~,~,logL] = estimate(Mdl, Y, 'Display','off');
    numParams = K^2*p + K;

    [aic_val,bic_val] = aicbic(logL,numParams,T);
    AIC(p)  = aic_val;
    BIC(p)  = bic_val;
    HQIC(p) = -2*logL + 2*numParams*log(log(T));
end

lags = (1:pmax)';

outdir = fullfile('figures','chapter2');
if ~exist(outdir,'dir'); mkdir(outdir); end

figure('Position',[100 100 900 600],'Color','w');
set(gcf,'Color','w');

plot(lags,AIC,'-o','LineWidth',1.8,'MarkerSize',8); hold on;
plot(lags,BIC,'-s','LineWidth',1.8,'MarkerSize',8);
plot(lags,HQIC,'-d','LineWidth',1.8,'MarkerSize',8);
hold off;

% ---- Formatting ----
xlabel('Lag order p','FontSize',14,'Color','k');
ylabel('Information criterion','FontSize',14,'Color','k');
title('VAR lag length selection','FontSize',16,'Color','k');

lgd = legend({'AIC','BIC','HQIC'},'Location','best','FontSize',12);
set(lgd,'Color','w','EdgeColor','k','LineWidth',0.8); % FIX: test legibility

grid on;

print('-dpdf', fullfile(outdir,'fig_lag_selection.pdf'));
