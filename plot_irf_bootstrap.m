%% PLOT IRF bootstrap – tutte le risposte e tutti gli shock
clear; clc;

load('var_svar_results.mat');   % carica IRF_med, IRF_CI_low, IRF_CI_high, varNames

horizon = size(IRF_med,1);
t = (1:horizon)';

k = length(varNames);

figure('Color','w','Position',[50 50 1200 800]);
tiledlayout(k,k,"TileSpacing","compact","Padding","compact");

for shock = 1:k
    for resp = 1:k
        
        nexttile;
        
        % estrai IRF bootstrap
        med  = squeeze(IRF_med(:,resp,shock));
        low  = squeeze(IRF_CI_low(:,resp,shock));
        high = squeeze(IRF_CI_high(:,resp,shock));

        % banda (area)
        fill([t; flipud(t)], [low; flipud(high)], ...
             [0.8 0.82 1], 'EdgeColor','none', 'FaceAlpha',0.4);
        hold on;

        % linea mediana
        plot(t, med, 'k', 'LineWidth', 2);

        % linea zero
        yline(0,'k--','LineWidth',1);

        title(sprintf('%s → %s', varNames{shock}, varNames{resp}), ...
            'FontWeight','bold', 'FontSize',10, 'Color', 'k');

        grid on;
    end
end

sgtitle('IRF VAR Bootstrap (5–95% percentile bands)', ...
    'FontWeight','bold','FontSize',14, 'Color', 'k');

%% Verifica numerica significatività IRF bootstrap

H = 8;  % orizzonte breve dove cercare significatività
nVars = length(varNames);

sig = false(H, nVars, nVars);

for shock = 1:nVars
    for resp = 1:nVars
        low  = IRF_CI_low(1:H, resp, shock);
        high = IRF_CI_high(1:H, resp, shock);

        % TRUE quando la banda non contiene lo zero
        sig(:,resp,shock) = (low .* high > 0);
    end
end

fprintf('\n=== MATRICE SIGNIFICATIVITÀ (TRUE = banda non include 0) ===\n');
disp(sig);