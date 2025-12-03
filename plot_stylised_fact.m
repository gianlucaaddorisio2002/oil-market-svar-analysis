%% ============================================================
% PLOT_STYLISED_FACT.M
% Figura 1 — Stylised Facts of the Oil Market (1990–2024)
% Usa il dataset costruito da build_oil_dataset.m
% Output: fig_stylised_oil_market.pdf
%% ============================================================

clear; clc;

% 1) Costruisci/aggiorna il dataset (crea All, All_d, ALL_VAR)
run('build_oil_dataset.m');   % deve essere nello stesso folder o in path

% 2) Estrai le serie dalla timetable All
if ~exist('All','var')
    error('La variabile "All" non esiste dopo build_oil_dataset.m');
end

t           = All.Time;
WTI_real    = All.WTI_real;
OCSE        = All.OCSE;
Production  = All.Production;
Inventories = All.Inventories;

% 3) Rimuovi eventuali osservazioni con NaN
valid = ~isnan(WTI_real) & ~isnan(OCSE) & ...
        ~isnan(Production) & ~isnan(Inventories);

t           = t(valid);
WTI_real    = WTI_real(valid);
OCSE        = OCSE(valid);
Production  = Production(valid);
Inventories = Inventories(valid);

% 4) Normalizza tutte le serie a Base 100 (prima osservazione = 100)
base100 = @(x) 100 * x ./ x(1);

WTI100  = base100(WTI_real);
OCSE100 = base100(OCSE);
Prod100 = base100(Production);
Inv100  = base100(Inventories);

% 5) Grafico 2x2
figure('Position',[100 100 1100 900]);

subplot(2,2,1);
plot(t,WTI100,'k','LineWidth',1.4);
title('Real WTI Price (Base 100)', 'Color', 'k');
grid on;

subplot(2,2,2);
plot(t,Prod100,'r','LineWidth',1.4);
title('Global Oil Production (Base 100)', 'Color', 'k');
grid on;

subplot(2,2,3);
plot(t,OCSE100,'b','LineWidth',1.4);
title('OECD Activity (Base 100)', 'Color', 'k');
grid on;

subplot(2,2,4);
plot(t,Inv100,'m','LineWidth',1.4);
title('U.S.A Inventories (Base 100)', 'Color', 'k');
grid on;

sgtitle('Stylised Facts of the Oil Market, 1990--2024', 'Color', 'k');

% 6) Esporta in PDF
print(gcf,'-dpdf','fig_stylised_oil_market.pdf');
