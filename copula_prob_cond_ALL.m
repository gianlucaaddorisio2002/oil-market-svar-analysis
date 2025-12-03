%% COPULA – Probabilità condizionate per tutte le coppie (τ = 0.95)
clear; clc;

%% 1) Carica dati e marginali
load('structural_shocks.mat','shocks_matrix','shockLabels','distFits');

fprintf('Shock disponibili:\n');
disp(shockLabels(:));

% coppie da testare (solo quelle economiche sensate)
pairs = {
    1,2;  % Supply -> AggregateDemand_OCSE
    1,3;  % Supply -> Precautionary
    2,3;  % AggregateDemand_OCSE -> Precautionary
};

tau = 0.95;  % soglia tail

results = [];

for k = 1:size(pairs,1)

    i = pairs{k,1};   % X
    j = pairs{k,2};   % Y

    X = shocks_matrix(:,i);
    Y = shocks_matrix(:,j);

    margX = distFits{i}{2};
    margY = distFits{j}{2};

    % pseudo osservazioni
    Ux = cdf(margX, X);
    Uy = cdf(margY, Y);

    epsU = 1e-6;
    Ux = min(max(Ux,epsU),1-epsU);
    Uy = min(max(Uy,epsU),1-epsU);

    U = [Ux Uy];

    %% Fit copula migliore (Gaussian, t, Clayton, Gumbel, Frank)
    families = {'Gaussian','t','Clayton','Gumbel','Frank'};
    AIC = inf;   % temporaneo
    bestFam = '';
    bestParams = [];

    for f = 1:length(families)
        fam = families{f};
        try
            switch fam
                case 'Gaussian'
                    theta = copulafit('Gaussian',U);
                    ll = sum(log(copulapdf('Gaussian',U,theta)));
                    kpar = 1;

                case 't'
                    [R,nu] = copulafit('t',U);
                    theta = {R,nu};
                    ll = sum(log(copulapdf('t',U,R,nu)));
                    kpar = 2;

                case 'Clayton'
                    theta = copulafit('Clayton',U);
                    ll = sum(log(copulapdf('Clayton',U,theta)));
                    kpar = 1;

                case 'Gumbel'
                    theta = copulafit('Gumbel',U);
                    ll = sum(log(copulapdf('Gumbel',U,theta)));
                    kpar = 1;

                case 'Frank'
                    theta = copulafit('Frank',U);
                    ll = sum(log(copulapdf('Frank',U,theta)));
                    kpar = 1;
            end

            AICtest = -2*ll + 2*kpar;

            if AICtest < AIC
                AIC = AICtest;
                bestFam = fam;
                bestParams = theta;
            end
        catch
        end
    end

    %% Probabilità condizionata

    x0 = quantile(X,tau);
    y0 = quantile(Y,tau);

    u0 = cdf(margX,x0);
    v0 = cdf(margY,y0);
    u0 = min(max(u0,epsU),1-epsU);
    v0 = min(max(v0,epsU),1-epsU);

    switch bestFam
        case 't'
            R = bestParams{1};
            nu = bestParams{2};
            Cuv = copulacdf('t',[u0 v0],R,nu);
        case 'Gaussian'
            R = bestParams;
            Cuv = copulacdf('Gaussian',[u0 v0],R);
        otherwise
            theta = bestParams;
            Cuv = copulacdf(bestFam,[u0 v0],theta);
    end

    Pincond = 1 - v0;
    Pcond = (1 - u0 - v0 + Cuv) / (1 - u0);
    Ratio = Pcond / Pincond;

    results = [results; {shockLabels{i},shockLabels{j},bestFam,Pcond,Pincond,Ratio,AIC}];

end

%% Mostra tabella finale
T = cell2table(results,...
    'VariableNames',{'Shock_X','Shock_Y','Best_Copula','P_cond',...
                     'P_incond','Ratio','AIC'});

disp('=== RISULTATI FINALI (τ = 0.95) ===')
disp(T)

writetable(T,'copula_conditional_results_tau95.csv');
fprintf('\nTabella salvata in: copula_conditional_results_tau95.csv\n');