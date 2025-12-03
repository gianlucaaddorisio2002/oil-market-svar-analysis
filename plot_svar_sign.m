%% ============================================================
%   PLOT SVAR (Sign Restrictions) + STRESS TEST HORMUZ
%   Variabili: OCSE_DL, Production_DL, WTI_real_DL, Inventories_DL
%% ============================================================

clear; clc;

%% ----- 0. Stile grafico uniforme -----
set(groot, ...
    'defaultFigureColor','w', ...
    'defaultAxesColor','w', ...
    'defaultAxesXColor','k', ...
    'defaultAxesYColor','k', ...
    'defaultAxesGridColor',[0.9 0.9 0.9], ...
    'defaultAxesFontName','Helvetica', ...
    'defaultAxesFontSize',10, ...
    'defaultLineLineWidth',1.8);

%% ----- 1. Carica risultati SVAR -----
if ~exist('svar_sign_results.mat','file')
    error('svar_sign_results.mat mancante. Esegui prima svar_sign_restrictions.m');
end

load('svar_sign_results.mat', ...
     'IRF_struct','IRF_low','IRF_high', ...
     'varNames','horizon','pChosen', ...
     'idxDemand','idxProd','idxWTI','idxInv');

[h, n, ~] = size(IRF_struct);
t = 0:h-1;

% Etichette leggibili
cleanNames = strrep(varNames,'_',' ');

% Indici shock strutturali (ordine: 1 supply, 2 demand, 3 precautionary, 4 residuo)
idxSupply      = idxProd;      % flow supply shock
idxAggDemand   = idxDemand;    % aggregate demand (OCSE)
idxPrecaution  = idxInv;       % precautionary / storage

%% ============================================================
%   2. GRIGLIA COMPLETA DELLE IRF STRUTTURALI
%% ============================================================
figure('Name','SVAR – IRF strutturali (Sign Restrictions)', ...
       'Position',[50 50 1200 800]);

tiledlayout(n,n,"TileSpacing","compact","Padding","compact");

for i = 1:n      % variabile risposta
    for j = 1:n  % shock strutturale
        nexttile;

        med = squeeze(IRF_struct(:,i,j));
        lo  = squeeze(IRF_low(:,i,j));
        hi  = squeeze(IRF_high(:,i,j));

        % Banda di confidenza
        fill([t fliplr(t)], [lo' fliplr(hi')], ...
             [0.85 0.88 1.0], 'EdgeColor','none', 'FaceAlpha',0.5);
        hold on;

        % IRF mediana
        plot(t, med, 'k-', 'LineWidth', 1.8);

        % Asse zero
        yline(0, '--', 'Color', [0.3 0.3 0.3]);

        title(sprintf('%s \x2192 shock %s', ...
              cleanNames{i}, cleanNames{j}), ...
              'FontSize', 9, 'FontWeight', 'bold', 'Color', 'k');

        grid on; box on;
        if i == n
            xlabel('Orizzonte (mesi)','Color','k');
        end
    end
end

sgtitle(sprintf('SVAR (Sign Restrictions) – IRF strutturali, VAR(%d)', pChosen), ...
    'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k');

%% ============================================================
%   3. FOCUS: IRF DEL WTI AI TRE SHOCK PRINCIPALI
%% ============================================================
figure('Name','SVAR – IRF del WTI (3 shock strutturali)', ...
       'Position',[70 70 900 800]);

tiledlayout(3,1,"TileSpacing","compact","Padding","compact");

% Range comune
allResp = [ IRF_struct(:,idxWTI,idxSupply); ...
            IRF_struct(:,idxWTI,idxAggDemand); ...
            IRF_struct(:,idxWTI,idxPrecaution) ];
ymin = min(allResp) - 0.05;
ymax = max(allResp) + 0.05;

% ----- Shock di offerta -----
nexttile;
fill([t fliplr(t)], ...
     [squeeze(IRF_low(:,idxWTI,idxSupply))' ...
      fliplr(squeeze(IRF_high(:,idxWTI,idxSupply))')], ...
     [0.85 0.88 1.0], 'EdgeColor', 'none');
hold on;
plot(t, IRF_struct(:,idxWTI,idxSupply), 'k-', 'LineWidth', 1.8);
yline(0, '--', 'Color', [0.3 0.3 0.3]);
ylim([ymin ymax]); grid on; box on;
title('IRF del WTI a shock di offerta', ...
      'FontSize', 11, 'FontWeight', 'bold', 'Color', 'k');
ylabel('\Delta WTI (z-score)', 'Color', 'k');

% ----- Shock di domanda aggregata (OCSE) -----
nexttile;
fill([t fliplr(t)], ...
     [squeeze(IRF_low(:,idxWTI,idxAggDemand))' ...
      fliplr(squeeze(IRF_high(:,idxWTI,idxAggDemand))')], ...
     [0.85 0.88 1.0], 'EdgeColor', 'none');
hold on;
plot(t, IRF_struct(:,idxWTI,idxAggDemand), 'k-', 'LineWidth', 1.8);
yline(0, '--', 'Color', [0.3 0.3 0.3]);
ylim([ymin ymax]); grid on; box on;
title('IRF del WTI a shock di domanda (OCSE)', ...
      'FontSize', 11, 'FontWeight', 'bold', 'Color', 'k');
ylabel('\Delta WTI (z-score)', 'Color', 'k');

% ----- Shock di domanda precauzionale / storage -----
nexttile;
fill([t fliplr(t)], ...
     [squeeze(IRF_low(:,idxWTI,idxPrecaution))' ...
      fliplr(squeeze(IRF_high(:,idxWTI,idxPrecaution))')], ...
     [0.85 0.88 1.0], 'EdgeColor', 'none');
hold on;
plot(t, IRF_struct(:,idxWTI,idxPrecaution), 'k-', 'LineWidth', 1.8);
yline(0, '--', 'Color', [0.3 0.3 0.3]);
ylim([ymin ymax]); grid on; box on;
title('IRF del WTI a shock precauzionale / storage', ...
      'FontSize', 11, 'FontWeight', 'bold', 'Color', 'k');
ylabel('\Delta WTI (z-score)', 'Color', 'k');
xlabel('Orizzonte (mesi)', 'Color', 'k');

sgtitle('SVAR – IRF del WTI ai tre shock strutturali (restrizioni di segno)', ...
    'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k');

%% ============================================================
%   4. FEVD STRUTTURALE DEL WTI (da IRF strutturali)
%% ============================================================
H = h;
nShocksToPlot = 3;   % offerta, domanda, precauzionale
FEVD_svar = zeros(H, nShocksToPlot);

for jShock = 1:nShocksToPlot
    num = zeros(H,1);
    den = zeros(H,1);
    for tau = 1:H
        num(tau) = sum( IRF_struct(1:tau,idxWTI,jShock).^2 );
        den(tau) = sum( sum( IRF_struct(1:tau,idxWTI,:).^2 ,3) );
    end
    FEVD_svar(:,jShock) = num ./ den;
end

fevdLabels = {'Shock di offerta','Shock di domanda (OCSE)','Shock precauzionale'};

figure('Name','SVAR – FEVD WTI','Position',[100 100 900 700]);
tiledlayout(nShocksToPlot,1,"TileSpacing","compact","Padding","compact");

for jShock = 1:nShocksToPlot
    nexttile;
    fevd_j = FEVD_svar(:,jShock)*100;  % %

    plot(t(2:end), fevd_j(2:end), 'k-', 'LineWidth', 1.8);
    grid on; box on;
    xlim([t(2) t(end)]);

    fevd_final = fevd_j(end);
    yline(fevd_final, '--', 'Color', [0 0 0], 'LineWidth', 1.3);

    text(0.78, 0.85, sprintf('Asintoto: %.2f%%', fevd_final), ...
         'Units','normalized', ...
         'FontSize', 11, 'FontWeight','bold', 'Color', 'k');

    ylabel('% varianza WTI', 'Color', 'k');
    title(fevdLabels{jShock}, 'Color', 'k', 'FontWeight','bold');
end

xlabel('Orizzonte (mesi)', 'Color', 'k');
sgtitle('SVAR – FEVD del WTI (restrizioni di segno)', ...
        'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k');

%% ============================================================
%   5. STRESS TEST HORMUZ: SHOCK DI OFFERTA ESTREMO
%      (chiusura Stretto di Hormuz)
%% ============================================================
% --> Usa IRF del WTI a uno shock di offerta per simulare
%     un percorso in LIVELLI di prezzo (USD) in caso di
%     shock di offerta molto grande (k_sigma std).

% Serve il livello reale del WTI:
if ~exist('clean_data.mat','file')
    warning('clean_data.mat mancante: salto lo stress test Hormuz.');
    return;
end
load('clean_data.mat','All');

if ~ismember('WTI_real', All.Properties.VariableNames)
    warning('Variabile WTI_real non trovata in All. Impossibile fare stress test Hormuz.');
    return;
end

wti_real = All.WTI_real;
if any(isnan(wti_real))
    wti_real = fillmissing(wti_real,'previous');
end

% Ritorni logaritmici reali (per passare da z-score a %)
log_wti   = log(wti_real);
dlog_wti  = diff(log_wti);                     % Δlog(WTI_real)
sigma_raw = std(dlog_wti,'omitnan');           % deviazione std "vera"

% Prezzo di partenza (reale, ultimo dato)
P0 = wti_real(end);

% Orizzonte stress test (max 24 mesi o quanto disponibile dalle IRF)
Hstress = min(24, h);
tau = (0:Hstress-1)';

% IRF del WTI a 1 shock di offerta (in z-score di Δlog)
irf_wti_supply = IRF_struct(1:Hstress, idxWTI, idxSupply);

% Magnitudo scenario Hormuz: k_sigma std dev del supply shock (negativo)
k_sigma = +3;   % es.: -3σ supply shock (molto pesante)
fprintf('\n[Stress test Hormuz] Sto simulando uno shock di offerta di %.1f sigma.\n', k_sigma);

% Δlog(WTI) in termini reali indotti dallo shock
dlog_stress = (k_sigma * sigma_raw) * irf_wti_supply;   % Hstress x 1
dlog_base   = zeros(Hstress,1);                        % baseline: nessuno shock

% Percorsi in livelli
logP0 = log(P0);
logP_stress = logP0 + cumsum(dlog_stress);
logP_base   = logP0 + cumsum(dlog_base);

P_stress = exp(logP_stress);
P_base   = exp(logP_base);

% Drop % dopo 12 mesi (se disponibile)
h12 = min(12, Hstress);
drop12 = (P_stress(h12)/P0 - 1)*100;

figure('Name','Stress test – Chiusura Stretto di Hormuz', ...
       'Position',[120 120 900 550]);

plot(tau, P_base,   '--', 'Color',[0.3 0.3 0.3], 'LineWidth',1.8); hold on;
plot(tau, P_stress, 'k-',  'LineWidth',2.0);

grid on; box on;
xlabel('Orizzonte (mesi)','Color','k');
ylabel('Prezzo reale WTI (USD)','Color','k');
title('Stress test – Shock di offerta estremo (scenario Hormuz)', ...
      'FontWeight','bold','Color','k');

lgd = legend({'Baseline (nessuno shock)','Scenario Hormuz (supply shock estremo)'}, ...
             'Location','best','Box','off');
set(lgd,'TextColor','k','FontWeight','bold');

text(h12, P_stress(h12), sprintf('  %.1f%% dopo %d mesi', drop12, h12), ...
     'FontSize',10,'FontWeight','bold','Color','k', ...
     'HorizontalAlignment','left','VerticalAlignment','bottom');

fprintf('Drop percentuale del WTI dopo %d mesi nello scenario Hormuz: %.2f%%\n', ...
        h12, drop12);