%% OLS e controlli preliminari pre-VAR
clear; clc;
load('clean_data.mat','All_d');

% === Restrizione del campione temporale: 1990–2024  ===
startSub = datetime(1990,1,1);
endSub   = datetime(2024,12,31);

%% SCATTER OLS (dati grezzi + retta OLS)
x = All_d.OCSE_DL;
y = All_d.WTI_real_DL;

figure; hold on; grid on;
scatter(x, y, 20, 'k','filled');

% Regressione OLS
b = regress(y, [ones(length(x),1), x]);
xline_plot = linspace(min(x), max(x), 200);
y_fit = b(1) + b(2).*xline_plot;

plot(xline_plot, y_fit, 'k', 'LineWidth', 2);

xlabel('OCSE\_DL');
ylabel('WTI\_real\_DL');
title('Scatter OLS – relazione statistica tra OCSE e WTI', 'Color', 'k');


%% Scatter matrix – variabili del modello statico OLS
X_ols = [All_d.OCSE_DL, All_d.Production_DL, All_d.WTI_real_DL, All_d.Inventories_DL];
names_ols = {'OCSE\_DL','Prod\_DL','WTI\_DL','Inv\_DL'};

figure('Color','w');
[H,AX,BIGAX] = plotmatrix(X_ols);

% stile nero per tutti gli scatter
for i = 1:numel(H)
    set(H(i), 'Color','k', 'Marker','.', 'MarkerSize',6);
end

% etichette assi
for i = 1:numel(names_ols)
    xlabel(AX(end,i), names_ols{i}, 'Color','k');
    ylabel(AX(i,1),   names_ols{i}, 'Color','k');
end

set(BIGAX,'XColor','k','YColor','k','FontSize',11);
title(BIGAX,'Scatter matrix – variabili OLS statico','Color','k');
grid(BIGAX,'on');


%% OLS e controlli preliminari pre-VAR
clear; clc;
load('clean_data.mat','All','All_d');

% === Restrizione del campione temporale: 1990–2024 ===
startSub = datetime(1990,1,1);
endSub   = datetime(2024,12,31);

All   = All(timerange(startSub,endSub), {'WTI_real','Production','OCSE','Inventories'});
All_d = All_d(timerange(startSub,endSub), ...
              {'WTI_real_DL','Production_DL','OCSE_DL','Inventories_DL'});

%% --- TEST DI STAZIONARIETÀ (ADF) PER LE TABELLE 2.2 E 2.3 ---

% ---- Tabella 2.2: livelli, specifica constant + trend ----
series_lvl = {All.WTI_real, All.Production, All.OCSE, All.Inventories};
names_lvl  = {'WTI_real','Production','OECD','Inventories'};

fprintf('\n================ ADF IN LEVELS (Tabella 2.2) ================\n');
fprintf('%-12s %-18s %12s %12s\n','Series','Specification','ADF stat','p-value');
fprintf('--------------------------------------------------------------\n');

for i = 1:numel(series_lvl)
    y = series_lvl{i};
    y = y(~isnan(y));

    % modello con costante + trend
    [~, pValue, stat] = adftest(y, 'Model','TS');

    fprintf('%-12s %-18s %12.3f %12.4f\n', ...
        names_lvl{i}, 'constant + trend', stat, pValue);
end

% ---- Tabella 2.3: trasformate, specifica solo constant ----
series_tr  = {All_d.WTI_real_DL, All_d.Production_DL, ...
              All_d.OCSE_DL,    All_d.Inventories_DL};
names_tr   = {'WTI_real_DL','Production_DL','OECD_DL','Inventories_DL'};

fprintf('\n=========== ADF SULLLE SERIE TRASFORMATE (Tabella 2.3) =======\n');
fprintf('%-15s %-12s %12s %12s\n','Series','Specification','ADF stat','p-value');
fprintf('--------------------------------------------------------------\n');

for i = 1:numel(series_tr)
    y = series_tr{i};
    y = y(~isnan(y));

    % modello con sola costante
    [~, pValue, stat] = adftest(y, 'Model','ARD');

    fprintf('%-15s %-12s %12.3f %12.4f\n', ...
        names_tr{i}, 'constant', stat, pValue);
end



%% --- Controllo multicollinearità grossa tra regressori ---
X = [All_d.Production_DL, All_d.OCSE_DL, All_d.Inventories_DL];
Rcorr = corrcoef(X);
disp('Matrice di correlazione tra regressori:');
disp(Rcorr);

%% === MODELLO LINEARE CON FITLM ===
T = timetable2table(All_d);
T(:,1) = [];   % rimuove la colonna datetime

mdl = fitlm(T, 'WTI_real_DL ~ Production_DL + OCSE_DL + Inventories_DL');
disp(mdl);

%% === TEST WALD GLOBALE SUI REGRESSORI (con matrice R) ===
% Ordine coefficienti in mdl.Coefficients:
% 1: (Intercept)
% 2: Production_DL
% 3: OCSE_DL
% 4: Inventories_DL

R = [0 1 0 0;   % Production_DL = 0
     0 0 1 0;   % OCSE_DL        = 0
     0 0 0 1];  % Inventories_DL= 0 

r = [0; 0; 0];

% coefTest restituisce p-value
pWald = coefTest(mdl, R, r);

fprintf('\n== TEST WALD GLOBALE SUI REGRESSORI ==\n');
fprintf('H0: Production_DL = OCSE_DL = Inventories_DL = 0\n');
fprintf('p-value (Wald) = %.4f\n', pWald);

if pWald < 0.05
    fprintf('→ Rifiuto H0: almeno un coefficiente è diverso da zero (il modello statico ha un minimo contenuto informativo).\n');
else
    fprintf('→ NON rifiuto H0: congiuntamente i regressori NON sono significativi → modello statico praticamente inutile.\n');
end

%% === ANALISI DEI RESIDUI DEL MODELLO LINEARE ===
% Vettore residui
e = mdl.Residuals.Raw;
e = double(e(:));
e = e(~isnan(e));

% Residui standardizzati (media 0, var 1)
e_std = (e - mean(e))/std(e);

% === Skewness e Kurtosis dei residui ===
skew_e      = skewness(e_std);          % indice di asimmetria
kurt_e      = kurtosis(e_std);          % curtosi "totale" (Normale ≈ 3)
excess_kurt = kurt_e - 3;               % curtosi in eccesso (Normale = 0)

fprintf('\n== MOMENTI DI ORDINE SUPERIORE (RESIDUI STD) ==\n');
fprintf('Skewness     = %.4f\n', skew_e);
fprintf('Kurtosis     = %.4f (Normale ≈ 3)\n', kurt_e);
fprintf('Excess kurt. = %.4f (Normale = 0)\n', excess_kurt);


%% === CONFRONTO EMPIRICO: Normale vs Logistica vs t-Student ===
% Residui standardizzati: e_std

% 1) Stima delle tre distribuzioni
distNorm = fitdist(e_std, 'Normal');            % N(mu, sigma)
distLog  = fitdist(e_std, 'Logistic');          % Logistic(mu, sigma)
distT    = fitdist(e_std, 'tLocationScale');    % t-Student (location-scale)

% 2) Grafico: istogramma + tre densità sovrapposte
figure('Name','Fit Normale vs Logistica vs t-Student', ...
       'Color','w', 'Position',[200 200 900 450]);

histogram(e_std, 'Normalization','pdf', 'NumBins',30, ...
          'FaceColor',[0.8 0.8 0.8], 'FaceAlpha',0.5, 'EdgeColor',[0.3 0.3 0.3]);
hold on;

xgrid = linspace(min(e_std), max(e_std), 400);

plot(xgrid, pdf(distNorm, xgrid), 'k-',  'LineWidth',1.8);  % Normale
plot(xgrid, pdf(distLog,  xgrid), 'b--', 'LineWidth',1.8);  % Logistica
plot(xgrid, pdf(distT,    xgrid), 'm-',  'LineWidth',1.8);  % t-Student

lgd = legend({'Residui','Normale','Logistica','t-Student'}, ...
             'Location','best','Box','off');  % NIENTE 'Color','k' QUI

lgd.FontWeight = 'bold';
lgd.TextColor  = 'k';      % testo nero
lgd.FontSize   = 11;       % opzionale


xlabel('Residui standardizzati');
ylabel('Densità');
title('\bf Confronto empirico: Normale vs Logistica vs t-Student', 'Color', 'k');
grid on; box on;
set(gca,'FontSize',11,'LineWidth',1);


%% 0) Analisi distribuzione residui con dfittool (GUI interattiva)
USE_DFITTOOL = true;  

set(groot,'defaultFigureColor','w');
set(groot,'defaultAxesColor','w');
set(groot,'defaultAxesXColor','k','defaultAxesYColor','k','defaultAxesZColor','k');

figure('Position',[200 200 900 450]);

histogram(e_std, 'Normalization','pdf', 'NumBins',30, ...
          'FaceColor',[0.3 0.6 1], 'FaceAlpha',0.45, 'EdgeColor',[0.2 0.2 0.2]);
hold on;

% Normale standard N(0,1)
xgrid = linspace(min(e_std), max(e_std), 400);
plot(xgrid, normpdf(xgrid,0,1), 'k', 'LineWidth',2);

xlabel('Residui standardizzati','FontSize',12, 'Color', 'k');
ylabel('Densità','FontSize',12, 'Color', 'k');
title('\bf Distribuzione dei residui standardizzati','FontSize',15,'Color','k');

lgd = legend({'Residui','Normale N(0,1)'}, 'Location','best','Box','off');
set(lgd, 'TextColor','k', 'FontWeight','bold', 'FontSize',12);

grid on;
box on;
set(gca,'FontSize',11,'LineWidth',1);

if USE_DFITTOOL
    fprintf('\nApro dfittool sui residui standardizzati...\n');
    dfittool(e_std);
end

%% 1) Serie temporale dei residui
figure;
plot(e,'Color',[1.0 0.5 0.0],'LineWidth',1.2);
hold on;

% *** Linea tratteggiata sullo zero ***
yline(0,'--','Color',[0 0 0],'LineWidth',1.2);

grid on;
title('Residui del modello', 'FontWeight','bold', 'Color',[0.1 0.1 0.1]);
xlabel('Osservazioni (mesi)');
ylabel('e_t');

hold off;


%% 2) Test di normalità + QQ-plot - JBTEST e KSTEST

% Test di normalità classici (sempre vs normale)
[h_jb, p_jb, jbStat] = jbtest(e);   % H0: residui ~ Normale
         % H0: residui ~ Normale
[~, p_ks] = kstest(e_std);     % H0: residui ~ N(0,1)

% === Stima distribuzioni per i QQ-plot ===
% Normale
distN = fitdist(e_std, 'Normal');      % già stimata sopra ma la ristimo qui per chiarezza
muN   = distN.mu;
sigmaN= distN.sigma;

% t-Student (location-scale)
distT = fitdist(e_std, 'tLocationScale');
nu    = distT.nu;
muT   = distT.mu;
sigmaT= distT.sigma;

% === Quantili empirici e teorici ===
e_sorted = sort(e_std);
n        = numel(e_sorted);
p        = (1:n) / (n + 1);           % probabilità per i quantili

% Quantili teorici
qN = muN + sigmaN * norminv(p,0,1);  % Normale stimata
qT = muT + sigmaT * tinv(p,nu);      % t-Student stimata

% === Figura con 2 QQ-plot: Normale vs t-Student ===
figure('Name','QQ-plot: Normale vs t-Student','Color','w');
tiledlayout(1,2,"TileSpacing","compact","Padding","compact");

% --- QQ-plot vs Normale ---
nexttile;
plot(qN, e_sorted, 'b+','MarkerSize',5); hold on;
plot(qN, qN, 'r--','LineWidth',1.2);   % bisettrice
grid on; box on;
title('QQ residui vs Normale','FontWeight','bold', 'Color', 'k');
xlabel('Quantili teorici Normale', 'Color', 'k');
ylabel('Residui ordinati', 'Color', 'k');

% --- QQ-plot vs t-Student ---
nexttile;
plot(qT, e_sorted, 'b+','MarkerSize',5); hold on;
plot(qT, qT, 'r--','LineWidth',1.2);   % bisettrice
grid on; box on;
title(sprintf('QQ residui vs t-Student (\\nu = %.2f)', nu), ...
      'FontWeight','bold', 'Color', 'k');
xlabel('Quantili teorici t-Student', 'Color', 'k');
ylabel('Residui ordinati', 'Color', 'k');

%% 3) Autocorrelazione (ACF e PACF)

% --- ACF ---
figure('Color','w');             % sfondo figura bianco
autocorr(e,'NumLags',20);

ax = gca;
set(ax, 'Color','w', ...         % sfondo assi bianco
        'XColor','k', ...
        'YColor','k', ...
        'FontSize',11, ...
        'LineWidth',1);

title('ACF residui', 'FontWeight','bold', 'Color','k');

lgd = legend;                    % legenda creata automaticamente da autocorr
set(lgd, 'Color','none', ...     % niente riquadro nero
         'Box','on', ...
         'TextColor','k', ...
         'FontWeight','bold', ...
         'FontSize',10, ...
         'AutoUpdate','off');

grid on;
box on;


% --- PACF ---
figure('Color','w');
parcorr(e,'NumLags',20);

ax = gca;
set(ax, 'Color','w', ...
        'XColor','k', ...
        'YColor','k', ...
        'FontSize',11, ...
        'LineWidth',1);

title('PACF residui', 'FontWeight','bold', 'Color','k');

lgd = legend;
set(lgd, 'Color','none', ...
         'Box','on', ...
         'TextColor','k', ...
         'FontWeight','bold', ...
         'FontSize',10, ...
         'AutoUpdate','off');

grid on;
box on;

%% 4) Riepilogo
fprintf('\n== ANALISI RESIDUI ==\n');
fprintf('JbTest (normalità): p = %.4f\n', p_jb);
fprintf('KsTest:             p = %.4f\n', p_ks);


%% === Test formali di autocorrelazione ===
% Durbin-Watson
num = sum(diff(e).^2);
den = sum(e.^2);
DW = num / den;

% Ljung-Box (test su più lag: 1, 6, 12, 24)
lags = [1 6 12 24];
[h_lb, p_lb, stat_lb] = lbqtest(e, "Lags", lags);

fprintf('\n== AUTOCORRELAZIONE RESIDUI ==\n');
fprintf('Durbin-Watson: %.4f (circa 2 = no autocorrelazione)\n', DW);

for i = 1:length(lags)
    fprintf('Ljung-Box lag %2d: stat = %.3f  p = %.4f  (h=%d)\n', ...
        lags(i), stat_lb(i), p_lb(i), h_lb(i));
end

% Per la tabella 2.4 usi SOLO il caso lag 12:
idx12      = (lags == 12);
LB12_stat  = stat_lb(idx12);
LB12_p     = p_lb(idx12);


%% === Test eteroschedasticità: Breusch-Pagan ===
Z = T{:, {'Production_DL', 'OCSE_DL', 'Inventories_DL'}};
Z = [ones(size(Z,1),1), Z];

% Residui al quadrato
s2 = e.^2;

% Aux regression
[~, ~, ~, ~, stats_bp] = regress(s2, Z);
R2_bp = stats_bp(1);
n = length(e);
k = size(Z,2) - 1;

LM = n * R2_bp;
p_BP = 1 - chi2cdf(LM, k);

fprintf('\n== ETEROSCHEDASTICITÀ RESIDUI ==\n');
fprintf('Breusch-Pagan LM=%.4f  p-value=%.4f\n', LM, p_BP);

%% === CONCLUSIONE ===

fprintf('\n=======================\nMODEL VALIDATION SUMMARY\n=======================\n');

if DW > 1.7 && DW < 2.3 && all(h_lb == 0)
    fprintf(' → Nessuna forte evidenza di autocorrelazione residua.\n');
else
    fprintf(' → Evidenza di autocorrelazione residua: modello statico insufficiente.\n');
end

if p_BP > 0.05
    fprintf(' → Nessuna evidenza di eteroschedasticità.\n');
else
    fprintf(' → Eteroschedasticità rilevata: OLS non ottimale.\n');
end

if p_jb > 0.05 && p_ks > 0.05
    fprintf(' → Residui compatibili con normalità.\n');
else
    fprintf(' → Residui non normali.\n');
end

fprintf('--------------------------------\nDecisione: ');

if any(h_lb == 1) || p_BP < 0.05
    fprintf('Passare a VAR è giustificato.\n');
else
    fprintf('Modello OLS accettabile, ma VAR può migliorare dinamica.\n');
end

fprintf('=======================\n\n');

fprintf('\n== RIEPILOGO PER TABELLA 2.4 ==\n');
fprintf('Ljung-Box (12 lags): stat = %.3f, p-value = %.4f\n', LB12_stat, LB12_p);
fprintf('Breusch-Pagan:       LM   = %.3f, p-value = %.4f\n', LM, p_BP);
fprintf('Jarque-Bera:         JB   = %.3f, p-value = %.4f\n', jbStat, p_jb);

[h_jb, p_jb, jbStat] = jbtest(e);   % l'hai già sopra

fprintf('Jarque–Bera:         JB   = %.3f, p-value = %.4f\n', jbStat, p_jb);
