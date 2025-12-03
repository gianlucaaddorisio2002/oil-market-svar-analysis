%% COPULA BIVARIATE TRA SHOCK STRUTTURALI
clear; clc;

%% ----- Stile grafico globale -----
set(groot, ...
    'defaultFigureColor','w', ...
    'defaultAxesColor','w', ...
    'defaultAxesXColor','k', ...
    'defaultAxesYColor','k', ...
    'defaultAxesGridColor',[0.85 0.85 0.85], ...
    'defaultAxesFontName','Helvetica', ...
    'defaultAxesFontWeight','normal', ...
    'defaultAxesFontSize',12, ...
    'defaultLineLineWidth',2, ...
    'defaultLegendTextColor','k');

%% 1) Carica shock + marginali
load('structural_shocks.mat',...
     'shocks_matrix','shockLabels','distFits');

[T,nShocks] = size(shocks_matrix);

fprintf('Osservazioni: %d\n',T);
disp('Shock disponibili:');
disp(shockLabels(:));

pairs = [1 2;
         1 3;
         2 3];

allResults = cell(size(pairs,1),1);

%% 2) Loop sulle coppie
for p = 1:size(pairs,1)

    i = pairs(p,1);
    j = pairs(p,2);

    Xi = shocks_matrix(:,i);
    Xj = shocks_matrix(:,j);

    marg_i = distFits{i}{2};   % risulta: t-Student è la migliore
    marg_j = distFits{j}{2};

    Ui = cdf(marg_i,Xi);
    Uj = cdf(marg_j,Xj);

    epsU = 1e-6;
    Ui = min(max(Ui,epsU),1-epsU);
    Uj = min(max(Uj,epsU),1-epsU);

    U = [Ui Uj];

    % Famiglie copula
    fams = {'Gaussian','t','Clayton','Gumbel','Frank'};
    AIC = zeros(5,1);
    BIC = zeros(5,1);
    LL  = zeros(5,1);
    param = cell(5,1);

    for f = 1:5
        fam = fams{f};
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

        LL(f) = ll;
        AIC(f) = -2*ll + 2*kpar;
        BIC(f) = -2*ll + kpar*log(T);
        param{f} = theta;
    end

    % Tabella risultati
    tab = table(fams',AIC,BIC,LL,'VariableNames',...
         {'Family','AIC','BIC','LogLik'});
    disp(tab);

    % Copula migliore
    [~,idx] = min(AIC);
    bestFam = fams{idx};
    bestParam = param{idx};

    allResults{p} = struct('pair',pairs(p,:),...
                           'bestFamily',bestFam,...
                           'results',tab,...
                           'params',bestParam);

    % ------------ PLOT ------------
    figure('Color','w','Position',[100 100 1200 400]);
    tiledlayout(1,3);

    % Scatter originali
    nexttile;
    scatter(Xi,Xj,15,'k','filled','MarkerFaceAlpha',0.4);
    xlabel(strrep(shockLabels{i},'_',' '),'Color','k');
    ylabel(strrep(shockLabels{j},'_',' '),'Color','k');
    title('Scatter shock','Color','k','FontWeight','bold');
    grid on; box on;

    % Scatter uniformi
    nexttile;
    scatter(Ui,Uj,15,'k','filled','MarkerFaceAlpha',0.4);
    xlabel('U1','Color','k'); ylabel('U2','Color','k');
    title('Scatter Uniform(0,1)','Color','k','FontWeight','bold');
    xlim([0 1]); ylim([0 1]); grid on; box on;

    % Contour copula
    nexttile;
    u = linspace(0.001,0.999,50);
    [U1,U2] = meshgrid(u,u);
    switch bestFam
        case 'Gaussian'
            Z = copulapdf('Gaussian',[U1(:) U2(:)],bestParam);
        case 't'
            Z = copulapdf('t',[U1(:) U2(:)],bestParam{1},bestParam{2});
        case 'Clayton'
            Z = copulapdf('Clayton',[U1(:) U2(:)],bestParam);
        case 'Gumbel'
            Z = copulapdf('Gumbel',[U1(:) U2(:)],bestParam);
        case 'Frank'
            Z = copulapdf('Frank',[U1(:) U2(:)],bestParam);
    end
    Z = reshape(Z,size(U1));
    contour(U1,U2,Z,10,'k');
    xlim([0 1]); ylim([0 1]);
    xlabel('U1','Color','k'); ylabel('U2','Color','k');
    title(['Copula ' bestFam],'Color','k','FontWeight','bold');
    grid on; box on;

    sgtitle(sprintf('Copula %s – %s vs %s',bestFam,...
            strrep(shockLabels{i},'_',' '),strrep(shockLabels{j},'_',' ')),...
            'FontWeight','bold','Color','k');
end

save('copula_results_shocks.mat','allResults','pairs','shockLabels');

fprintf('\n>> Copule bivariate salvate.\n');