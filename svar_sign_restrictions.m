%% SVAR a restrizioni di segno – versione Kilian-style (4 variabili)
% Ordine VAR (obbligatorio):
%   1) Production_DL      -> offerta
%   2) OCSE_DL            -> domanda globale (attività economica)
%   3) WTI_real_DL        -> prezzo reale WTI
%   4) Inventories_DL     -> scorte (precautionary/storage)
%
% Shock strutturali:
%   1 = Supply shock (flow supply)
%   2 = Aggregate demand shock
%   3 = Precautionary / storage demand shock
%   4 = Residual (nessun vincolo)

clear; clc;

%% 0) Carica risultati VAR (ridotta forma)
if ~exist('var_svar_results.mat','file')
    error('var_svar_results.mat mancante. Esegui prima var_main.m');
end

load('var_svar_results.mat','Mdl_final','varNames','pChosen');

fprintf('SVAR con restrizioni di segno su VAR(%d)\n', pChosen);
disp('Ordine variabili (VAR ridotta forma):');
disp(varNames(:));

%% Controllo ordine variabili
idxProd = find(strcmpi(varNames,'Production_DL'));
idxOCSE = find(strcmpi(varNames,'OCSE_DL'));
idxWTI  = find(strcmpi(varNames,'WTI_real_DL'));
idxInv  = find(strcmpi(varNames,'Inventories_DL'));

if any(isempty([idxProd idxOCSE idxWTI idxInv]))
    error('Le variabili non sono nell''ordine Production, OCSE, WTI, Inventories. Controlla var_main.m');
end

%% 1) IRF ridotta forma con Cholesky
horizon  = 24;   % mesi
IRF_chol = irf(Mdl_final,'NumObs',horizon);   % h x n x n
[h, n, ~] = size(IRF_chol);

if n ~= 4
    error('Mi aspetto esattamente 4 variabili nel VAR.');
end

t = (0:h-1)';

%% 2) Definizione restrizioni di segno (Kilian-style)

% signR(iVar, jShock, tau) = segno richiesto (+1, -1, 0)
%
% VAR order:
%   1: Production_DL
%   2: OCSE_DL
%   3: WTI_real_DL
%   4: Inventories_DL

H_restrict = 1;         % SOLO impatto (t = 0) - si può variare per test
signR      = zeros(n, n, H_restrict);

tau = 1;   % t = 0 in Matlab (prima riga delle IRF)

% ---------- SHOCK 1: FLOW SUPPLY ----------
% Idee di base ( secondo Kilian 2009):
%   - Produzione ↓ (shock negativo di offerta)
%   - Prezzo WTI ↑ (meno offerta → prezzo sale)
%   - Inventories ↓ o 0 (meno disponibilità fisica)
%
signR(1,1,tau) = -1;   % Production_DL ↓
signR(3,1,tau) = +1;   % WTI_real_DL ↑
% signR(4,1,tau) = -1; % opzionale: scorte ↓ (puoi attivarlo se vuoi più forza)

% ---------- SHOCK 2: AGGREGATE DEMAND ----------
% Idee di base:
%   - OCSE_DL ↑ (più attività globale)
%   - WTI_real_DL ↑ (più domanda di petrolio)
%   - Production_DL ↑ (le imprese reagiscono)
%   - Inventories_DL ↑ (più uso capacità / accumulo scorte)
%
signR(2,2,tau) = +1;   % OCSE_DL ↑
signR(3,2,tau) = +1;   % WTI_real_DL ↑
% segni più deboli sulle altre:
signR(1,2,tau) = +1;   % Production_DL ↑
signR(4,2,tau) = +1;   % Inventories_DL ↑

% ---------- SHOCK 3: PRECAUTIONARY / STORAGE DEMAND ----------
% Idee di base:
%   - WTI_real_DL ↑ (timore futuro, prezzo sale)
%   - Inventories_DL ↑ (si accumulano scorte)
%   - Nessun effetto contemporaneo su produzione e attività reale (0)
%
signR(3,3,tau) = +1;   % WTI_real_DL ↑
signR(4,3,tau) = +1;   % Inventories_DL ↑
% Production_DL e OCSE_DL = 0 (nessun vincolo)

% ---------- SHOCK 4: RESIDUAL ----------
% Nessun vincolo
% (colonna 4 di signR tutta a zero)

%% 3) Ricerca rotazioni ortogonali

nDraws      = 100000;
maxAccepted = 1000;

acceptedIRF = zeros(h, n, n, maxAccepted);
acceptedQ   = zeros(n, n, maxAccepted);
nAcc        = 0;

tol = 1e-4;

fprintf('\nInizio ricerca di rotazioni (restrizioni di segno)...\n');

for draw = 1:nDraws
    
    % 3.1 Matrice ortogonale Q ~ Haar
    [Q, ~] = qr(randn(n));
    if det(Q) < 0
        Q(:,1) = -Q(:,1);
    end

    % 3.2 Ruota IRF ridotta forma
    IRF_Q = zeros(h, n, n);
    for tt = 1:h
        IRF_Q(tt,:,:) = squeeze(IRF_chol(tt,:,:)) * Q;
    end

    % 3.3 Controllo segni a impatto
    ok = true;
    for jShock = 1:n
        for iVar = 1:n
            sgn = signR(iVar, jShock, tau);
            if sgn == 0
                continue;
            end
            val = IRF_Q(tau, iVar, jShock);
            if (sgn > 0 && val <  tol) || (sgn < 0 && val > -tol)
                ok = false;
                break;
            end
        end
        if ~ok, break; end
    end

    if ~ok
        continue;
    end

    % 3.4 Salva rotazione accettata
    nAcc = nAcc + 1;
    acceptedIRF(:,:,:,nAcc) = IRF_Q;
    acceptedQ(:,:,nAcc)     = Q;

    if mod(nAcc,50) == 0
        fprintf('Rotazioni accettate: %d\n', nAcc);
    end

    if nAcc >= maxAccepted
        break;
    end
end

if nAcc == 0
    error('Nessuna rotazione soddisfa le restrizioni di segno. Il problema è il VAR o le restrizioni troppo forti.');
end

fprintf('\nTotale rotazioni accettate: %d\n', nAcc);

%% 4) IRF strutturali finali (mediana + banda 16–84%)

IRF_struct = median(acceptedIRF(:,:,:,1:nAcc), 4);
IRF_low    = prctile(acceptedIRF(:,:,:,1:nAcc), 16, 4);
IRF_high   = prctile(acceptedIRF(:,:,:,1:nAcc), 84, 4);

%% 5) Salva risultati SVAR

idxDemand = idxOCSE;
idxInv    = idxInv;   % solo per chiarezza

save('svar_sign_results.mat', ...
     'IRF_struct','IRF_low','IRF_high', ...
     'acceptedQ','varNames','horizon','pChosen', ...
     'idxDemand','idxProd','idxWTI','idxInv');

fprintf('\nRisultati SVAR con nuove restrizioni salvati in svar_sign_results.mat\n');