%% Estrai IRF del WTI a shock chiave dal VAR Cholesky
clear; clc;

load('var_svar_results.mat', 'IRF', 'varNames');

% Indici delle variabili (devono riflettere l'ordine nel VAR)
idxProd = find(strcmp(varNames,'Production_DL'));
idxMacro = find(ismember(varNames, {'OCSE_DL','REA_DL'}));
idxWTI  = find(strcmp(varNames,'WTI_real_DL'));
idxInv  = find(strcmp(varNames,'Inventories_DL'));

if any(isempty([idxProd, idxMacro, idxWTI, idxInv]))
    error('Controlla i nomi in varNames.');
end

% IRF del WTI a 1 s.d. shocks
irf_wti_supply  = -IRF(:, idxWTI, idxProd);   % calo Production = shock di offerta negativo
irf_wti_demand  = -IRF(:, idxWTI, idxMacro);  % calo OCSE/REA = shock di domanda negativa
irf_wti_storage =  IRF(:, idxWTI, idxInv);    % aumento Inventories = shock precauzionale/storage

horizon = (0:size(IRF,1)-1)';

save('irf_wti_shocks.mat', ...
     'horizon','irf_wti_supply','irf_wti_demand','irf_wti_storage');

disp('IRF del WTI per supply/demand/storage salvate in irf_wti_shocks.mat');