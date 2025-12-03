%% ============================================================
% PLOT_VAR_SERIES.M
% Figura 2.2 – Transformed VAR series (OCSE_cycle, Prod_log, WTI_log, Inv_log)
%% ============================================================

clear; clc;

if ~exist('clean_data.mat','file')
    error('clean_data.mat non trovato. Esegui prima build_oil_dataset.m');
end

load('clean_data.mat','All','ALL_VAR');

t = All.Time;

OCSE_cycle = ALL_VAR.OCSE_cycle;
Prod_log   = ALL_VAR.Prod_log;
WTI_log    = ALL_VAR.WTI_log;
Inv_log    = ALL_VAR.Inv_log;

outdir = fullfile('figures','chapter2');
if ~exist(outdir,'dir'); mkdir(outdir); end

figure('Position',[100 100 1100 900]);

subplot(2,2,1);
plot(t,WTI_log,'k','LineWidth',1.4);
title('Real WTI (log)', 'Color', 'k'); grid on;

subplot(2,2,2);
plot(t,Prod_log,'r','LineWidth',1.4);
title('Global Oil Production (log)', 'Color', 'k'); grid on;

subplot(2,2,3);
plot(t,OCSE_cycle,'b','LineWidth',1.4);
title('OECD Activity (cycle)', 'Color', 'k'); grid on;

subplot(2,2,4);
plot(t,Inv_log,'m','LineWidth',1.4);
title('U.S. Inventories (log)', 'Color', 'k'); grid on;

sgtitle('Transformed VAR series, 1990--2024', 'Color', 'k');

print('-dpdf', fullfile(outdir,'fig_var_series.pdf'));
