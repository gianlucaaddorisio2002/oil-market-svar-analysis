%% ============================================================
% JETFUEL_PASS_THROUGH
% Stima del legame di lungo periodo tra WTI reale e Jet fuel
% log(JF_bbl) = alpha + beta * log(WTI_real) + eps
% + linearizzazione locale intorno a WTI = 60 USD/bbl
%% ============================================================

clear; clc;

if ~exist('clean_data.mat','file')
    error('clean_data.mat mancante. Esegui prima build_oil_dataset.m');
end

load('clean_data.mat','All');   % timetable mensile

%% 1) Controllo variabili
if ~all(ismember({'WTI_real','JetFuel'}, All.Properties.VariableNames))
    error('In All devono esserci WTI_real e JetFuel (controlla build_oil_dataset).');
end

% Opzionale: restringi il campione a regime "moderno"
sampleStart = datetime(2000,1,1);   % puoi cambiare in 1991 se vuoi tutto
T = All(timerange(sampleStart, All.Time(end), 'closed'), ...
        {'WTI_real','JetFuel'});

% Rimuovi eventuali NaN
T = rmmissing(T);

wti  = T.WTI_real;      % USD/bbl (reale)
jf   = T.JetFuel;       % USD/gallon

%% 2) Conversione Jet fuel in USD per barile
gallon_per_barrel = 42;
jf_bbl = jf * gallon_per_barrel;   % USD/bbl

%% 3) Regressione in log-livelli (elasticità di lungo periodo)

y = log(jf_bbl);        % log prezzo jet fuel (bbl)
x = log(wti);           % log WTI reale

X = [ones(size(x)) x];
[b,~,e,~,stats] = regress(y, X);

alpha = b(1);
beta  = b(2);
R2    = stats(1);
sigma_eps = std(e);

fprintf('\n=== Pass-through WTI → Jet fuel (log-log) ===\n');
fprintf('log(JF_bbl) = %.4f + %.4f * log(WTI_real)\n', alpha, beta);
fprintf('R^2 = %.3f, sigma_eps = %.4f\n\n', R2, sigma_eps);

%% 4) Linearizzazione locale intorno a WTI0 = 60 USD/bbl

WTI0 = 60;                       % livello di riferimento (coerente con stress test)
JF0_hat = exp(alpha) * WTI0^beta;  % prezzo jet fuel stimato a WTI=60

% derivata dJF/dWTI in log-log: beta * JF / WTI
beta_lin  = beta * (JF0_hat / WTI0);
alpha_lin = JF0_hat - beta_lin * WTI0;

fprintf('Linearizzazione intorno a WTI = %.2f\n', WTI0);
fprintf('JF_hat(WTI) ≈ alpha_lin + beta_lin * WTI\n');
fprintf('alpha_lin = %.4f, beta_lin = %.4f (USD/bbl)\n\n', alpha_lin, beta_lin);

%% 5) Salva risultati per il modello di hedging

save('jetfuel_pass_results.mat', ...
     'alpha','beta','R2','sigma_eps', ...
     'alpha_lin','beta_lin', ...
     'WTI0','JF0_hat', ...
     'sampleStart');

disp('Risultati salvati in jetfuel_pass_results.mat');
