%% ============================================================
% ESTIMATE_PASS_THROUGH_JETFUEL
% Stima del legame di lungo periodo tra WTI reale e Jet fuel USGC
% log(JF_bbl) = alpha + beta * log(WTI_real) + eps
% + linearizzazione locale intorno a WTI = 60 USD/bbl
%% ============================================================

clear; clc;

if ~exist('clean_data.mat','file')
    error('clean_data.mat mancante. Esegui prima build_oil_dataset.m');
end

load('clean_data.mat','All');   % timetable mensile

%% 1) Controllo variabili
if ~all(ismember({'WTI_real','JetFuel_USGC'}, All.Properties.VariableNames))
    error('In All devono esserci WTI_real e JetFuel_USGC (controlla build_oil_dataset).');
end

% Restringi il campione a un regime moderno (puoi allargare a 1991 se vuoi)
sampleStart = datetime(2000,1,1);
T = All(timerange(sampleStart, All.Time(end), 'closed'), ...
        {'WTI_real','JetFuel_USGC'});

T = rmmissing(T);

wti  = T.WTI_real;        % USD/bbl (reale)
jf_gal = T.JetFuel_USGC;  % USD per gallon

%% 2) Conversione Jet fuel in USD per barile
gallon_per_barrel = 42;
jf_bbl = jf_gal * gallon_per_barrel;   % USD/bbl

%% 3) Regressione in log-livelli

y = log(jf_bbl);          % log prezzo jet fuel (bbl)
x = log(wti);             % log WTI reale

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

WTI0 = 60;
JF0_hat = exp(alpha) * WTI0^beta;   % prezzo jet fuel stimato a WTI=60

% derivata locale dJF/dWTI in log-log: beta * JF / WTI
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

%% 6) Fitted values and scatter plot WTI–Jet Fuel

% fitted in log-levels and in USD/bbl
JF_hat_log = X*b;              % log(JF_hat)
JF_hat_bbl = exp(JF_hat_log);  % JF_hat in USD/bbl

figure;
scatter(wti, jf_bbl, 15, 'filled'); hold on;
plot(wti, JF_hat_bbl, 'LineWidth', 1.5);
grid on;
xlabel('WTI real price (USD/bbl)');
ylabel('Jet Fuel USGC (USD/bbl)');
title('Observed vs fitted Jet Fuel prices', 'Color','k');

print('-dpng','figures/chapter4/fig_pass_through_scatter.png');

%% 7) Salva anche le serie per eventuali usi futuri
save('jetfuel_pass_results.mat', ...
     'alpha','beta','R2','sigma_eps', ...
     'alpha_lin','beta_lin', ...
     'WTI0','JF0_hat', ...
     'sampleStart', ...
     'wti','jf_bbl','JF_hat_bbl');


