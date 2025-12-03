%% ============================================================
% PLOT_VAR_CHOLESKY.M
% Plot IRF & FEVD del VAR ridotto (Cholesky) + focus sul WTI
% Coerente con var_main.m
%
% Ordine economico atteso nel VAR (var_main):
%   1) Production_DL   -> offerta (supply)
%   2) OCSE_DL         -> domanda aggregata globale
%   3) WTI_real_DL     -> prezzo reale del WTI
%   4) Inventories_DL  -> scorte (precautionary/storage)
%
% ATTENZIONE: - aggiornato perchè prima c'era incoerenza 
% - Gli indici per supply/demand/precautionary sono ricavati dai NOMI,
%   non hard-coded, per evitare incoerenze.
%% ============================================================

clear; clc;

%% ----- Stile generale figure -----
set(groot, ...
    'defaultFigureColor'  , 'w', ...
    'defaultAxesColor'    , 'w', ...
    'defaultAxesXColor'   , 'k', ...
    'defaultAxesYColor'   , 'k', ...
    'defaultAxesGridColor', [0.9 0.9 0.9], ...
    'defaultAxesFontName' , 'Helvetica', ...
    'defaultAxesFontSize' , 10, ...
    'defaultLineLineWidth', 2);

%% ----- Carica risultati VAR ridotto + Cholesky -----
if ~exist('var_svar_results.mat','file')
    error('File var_svar_results.mat mancante. Esegui prima var_main.m');
end

load('var_svar_results.mat', ...
     'Mdl_final','IRF','FEVD','varNames','pChosen');

fprintf('Plotting IRF and FEVD for reduced-form VAR(%d) – Cholesky.\n', pChosen);

% Controllo dimensioni IRF
[h, n, n2] = size(IRF);
if n ~= n2
    error('IRF ha dimensioni inattese: non è (h × n × n).');
end
t = 0:h-1;

disp('Nomi variabili nel VAR (ordine come stimato in var_main):');
disp(varNames(:));

%% ----- Ricava indici per supply/demand/preca/WTI in modo robusto -----
idxProd = find(strcmp(varNames,'Production_DL'));
idxOCSE = find(strcmp(varNames,'OCSE_DL'));
idxWTI  = find(strcmp(varNames,'WTI_real_DL'));
idxInv  = find(strcmp(varNames,'Inventories_DL'));

if isempty(idxProd) || isempty(idxOCSE) || isempty(idxWTI) || isempty(idxInv)
    error(['Uno o più nomi variabili non trovati in varNames. ', ...
           'Controlla coerenza tra var_main.m e var_svar_results.mat.']);
end

fprintf('\nIndici trovati:\n');
fprintf('  Production_DL   (supply)       -> %d\n', idxProd);
fprintf('  OCSE_DL        (agg. demand)  -> %d\n', idxOCSE);
fprintf('  WTI_real_DL    (price)        -> %d\n', idxWTI);
fprintf('  Inventories_DL (precaution.)  -> %d\n\n', idxInv);

%% =====================================================
% 1) IRF – griglia completa (diagnostica)
%    IRF(i,j): risposta di variabile i a shock in variabile j
%% =====================================================

figure('Name','VAR – IRF (Cholesky)');
tiledlayout(n, n, "TileSpacing","compact","Padding","compact");

for i = 1:n          % variabile risposta (righe)
    for j = 1:n      % shock nella variabile j (colonne)
        nexttile;

        resp_ij = squeeze(IRF(:, i, j));
        plot(t, resp_ij, 'Color',[0 0.2 0.6]); hold on;
        yline(0,'Color',[0.5 0.5 0.5],'LineWidth',1);

        grid on; box on;

        title(sprintf('%s: risposta a shock in %s', ...
            strrep(varNames{i},'_',' '), ...
            strrep(varNames{j},'_',' ')), ...
            'FontSize',8,'FontWeight','bold','Color','k');

        if i == n
            xlabel('Mesi','Color','k');
        end
    end
end

sgtitle(sprintf('IRF (VAR ridotto, Cholesky) – VAR(%d)', pChosen), ...
    'FontSize',13,'FontWeight','bold','Color','k');


%% =====================================================
% 2) FEVD – griglia completa (diagnostica)
%    FEVD(i,j,h): quota varianza di variabile i spiegata dallo shock j
%% =====================================================

[h_f, n_f, n2_f] = size(FEVD);
if h_f ~= h || n_f ~= n || n2_f ~= n
    warning('Dimensioni FEVD inattese rispetto a IRF (atteso: h×n×n).');
end

figure('Name','VAR – FEVD (Cholesky)');
tiledlayout(n, n, "TileSpacing","compact","Padding","compact");

for i = 1:n
    for j = 1:n
        nexttile;

        fevd_ij = squeeze(FEVD(:, i, j)) * 100;  % in percentuale
        plot(1:h_f, fevd_ij, 'Color',[0 0.2 0.6]); hold on;
        ylim([0 100]);

        yline(0,'Color',[0.5 0.5 0.5],'LineWidth',1);
        grid on; box on;

        title(sprintf('%s: quota spiegata da shock in %s', ...
            strrep(varNames{i},'_',' '), ...
            strrep(varNames{j},'_',' ')), ...
            'FontSize',8,'FontWeight','bold','Color','k');

        xlabel('Orizzonte (mesi)','Color','k');
        ylabel('%','Color','k');
    end
end

sgtitle(sprintf('FEVD (VAR ridotto, Cholesky) – VAR(%d)', pChosen), ...
    'FontSize',13,'FontWeight','bold','Color','k');


%% =====================================================
% 3) IRF specifiche per il WTI (innovazioni VAR)
%    Focus: risposta del WTI a shock di:
%      - supply (Production_DL)
%      - aggregate demand (OCSE_DL)
%      - precautionary (Inventories_DL)
%
%    Convenzione: rappresentiamo shock NEGATIVO di supply/demand
%                 (↓ Production, ↓ OCSE) → cambiamo segno alle IRF.
%% =====================================================

IRF_prod_WTI = squeeze(IRF(:, idxWTI, idxProd));   % shock in Production (supply)
IRF_ocse_WTI = squeeze(IRF(:, idxWTI, idxOCSE));   % shock in OCSE (agg. demand)
IRF_inv_WTI  = squeeze(IRF(:, idxWTI, idxInv));    % shock in Inventories (precautionary)

% Convenzione: consideriamo shock NEGATIVO di supply/demand
IRF_prod_WTI = -IRF_prod_WTI;   % calo della produzione = shock di offerta negativo
IRF_ocse_WTI = -IRF_ocse_WTI;   % calo della domanda globale = shock di domanda negativo
% Precautionary lo lasciamo con segno "naturale": aumento delle scorte.

figure('Name','VAR – IRF WTI (Cholesky)');
tiledlayout(3,1,"TileSpacing","compact","Padding","compact");

allResp = [IRF_prod_WTI; IRF_ocse_WTI; IRF_inv_WTI];
ymin = min(allResp); ymax = max(allResp);
if ymax - ymin < 1e-8
    marg = 0.1;
else
    marg = 0.05*(ymax-ymin);
end
ymin = ymin - marg; ymax = ymax + marg;

% 3.1 Shock di offerta (Production ↓)
nexttile;
plot(t, IRF_prod_WTI,'LineWidth',2,'Color',[0 0.2 0.6]); hold on;
yline(0,'Color',[0.5 0.5 0.5],'LineWidth',1);
ylim([ymin ymax]);
grid on; box on;
title('VAR – risposta del WTI a un calo della Production (shock di offerta negativo)', ...
      'FontSize',11,'FontWeight','bold','Color','k');
ylabel('Risposta WTI','Color','k');

% 3.2 Shock di domanda aggregata (OCSE ↓)
nexttile;
plot(t, IRF_ocse_WTI,'LineWidth',2,'Color',[0.8 0.2 0.1]); hold on;
yline(0,'Color',[0.5 0.5 0.5],'LineWidth',1);
ylim([ymin ymax]);
grid on; box on;
title('VAR – risposta del WTI a un calo dell''OCSE IPI (shock di domanda aggregata negativo)', ...
      'FontSize',11,'FontWeight','bold','Color','k');
ylabel('Risposta WTI','Color','k');

% 3.3 Shock precauzionale (Inventories ↑)
nexttile;
plot(t, IRF_inv_WTI,'LineWidth',2,'Color',[0 0.6 0.2]); hold on;
yline(0,'Color',[0.5 0.5 0.5],'LineWidth',1);
ylim([ymin ymax]);
grid on; box on;
title('VAR – risposta del WTI a un aumento delle Inventories (shock precauzionale)', ...
      'FontSize',11,'FontWeight','bold','Color','k');
xlabel('Orizzonte (mesi)','Color','k');
ylabel('Risposta WTI','Color','k');


%% =====================================================
% 4) FEVD WTI – 3 subplot con asintoto e valore
%    Mostra la quota di varianza del WTI spiegata:
%      - da supply (Production)
%      - da aggregate demand (OCSE)
%      - da precautionary (Inventories)
%% =====================================================

t_fevd = 1:h_f;

fevd_supply_WTI = squeeze(FEVD(:, idxWTI, idxProd)) * 100;
fevd_demand_WTI = squeeze(FEVD(:, idxWTI, idxOCSE)) * 100;
fevd_prec_WTI   = squeeze(FEVD(:, idxWTI, idxInv )) * 100;

figure('Name','VAR – FEVD WTI per shock (Cholesky)');
tiledlayout(3,1,"TileSpacing","compact","Padding","compact");

%% 4.1 Supply (Production_DL)
nexttile;
plot(t_fevd, fevd_supply_WTI, 'LineWidth',2,'Color',[0 0.2 0.6]); hold on;
asymp = fevd_supply_WTI(end);
yline(asymp,'--','Color',[0.3 0.3 0.3],'LineWidth',1);
if max(fevd_supply_WTI) > 0
    label_offset = 0.07 * max(fevd_supply_WTI);
else
    label_offset = 1;
end
text(h_f*0.85, asymp - label_offset, sprintf(' %.2f%%', asymp), ...
     'Color',[0.2 0.2 0.2],'FontSize',12,'FontWeight','bold');
ylim([0 max(fevd_supply_WTI)*1.05 + 1e-6]);
yline(0,'Color',[0.5 0.5 0.5],'LineWidth',1);
grid on; box on;
xlabel('Orizzonte (mesi)','Color','k');
ylabel('% varianza spiegata','Color','k');
title(sprintf('WTI – quota spiegata da %s (supply)', ...
      strrep(varNames{idxProd},'_',' ')), ...
      'FontSize',12,'FontWeight','bold','Color','k');

%% 4.2 Aggregate Demand (OCSE_DL)
nexttile;
plot(t_fevd, fevd_demand_WTI, 'LineWidth',2,'Color',[0 0.2 0.6]); hold on;
asymp = fevd_demand_WTI(end);
yline(asymp,'--','Color',[0.3 0.3 0.3],'LineWidth',1);
if max(fevd_demand_WTI) > 0
    label_offset = 0.07 * max(fevd_demand_WTI);
else
    label_offset = 1;
end
text(h_f*0.85, asymp - label_offset, sprintf(' %.2f%%', asymp), ...
     'Color',[0.2 0.2 0.2],'FontSize',12,'FontWeight','bold');
ylim([0 max(fevd_demand_WTI)*1.05 + 1e-6]);
yline(0,'Color',[0.5 0.5 0.5],'LineWidth',1);
grid on; box on;
xlabel('Orizzonte (mesi)','Color','k');
ylabel('% varianza spiegata','Color','k');
title(sprintf('WTI – quota spiegata da %s (aggregate demand)', ...
      strrep(varNames{idxOCSE},'_',' ')), ...
      'FontSize',12,'FontWeight','bold','Color','k');

%% 4.3 Precautionary (Inventories_DL)
nexttile;
plot(t_fevd, fevd_prec_WTI, 'LineWidth',2,'Color',[0 0.2 0.6]); hold on;
asymp = fevd_prec_WTI(end);
yline(asymp,'--','Color',[0.3 0.3 0.3],'LineWidth',1);
if max(fevd_prec_WTI) > 0
    label_offset = 0.07 * max(fevd_prec_WTI);
else
    label_offset = 1;
end
text(h_f*0.85, asymp - label_offset, sprintf(' %.2f%%', asymp), ...
     'Color',[0.2 0.2 0.2],'FontSize',12,'FontWeight','bold');
ylim([0 max(fevd_prec_WTI)*1.05 + 1e-6]);
yline(0,'Color',[0.5 0.5 0.5],'LineWidth',1);
grid on; box on;
xlabel('Orizzonte (mesi)','Color','k');
ylabel('% varianza spiegata','Color','k');
title(sprintf('WTI – quota spiegata da %s (precautionary)', ...
      strrep(varNames{idxInv},'_',' ')), ...
      'FontSize',12,'FontWeight','bold','Color','k');

sgtitle('VAR – FEVD del WTI per i tre shock principali (Cholesky)', ...
    'FontSize',13,'FontWeight','bold','Color','k');
