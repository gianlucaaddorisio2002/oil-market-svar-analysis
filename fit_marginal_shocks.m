%% FIT MARGINAL DISTRIBUTIONS FOR STRUCTURAL SHOCKS
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

if ~exist('structural_shocks.mat','file')
    error('Esegui prima extract_structural_shocks.m');
end

load('structural_shocks.mat','Eps','shockLabels');

[T, n] = size(Eps);
fprintf('Numero osservazioni: %d, shock: %d\n', T, n);

distFits = cell(n,1);  % <- qui salviamo le marginali

for j = 1:3
    x = Eps(:,j);

    fprintf('\n==============================\n');
    fprintf('Shock %d: %s\n', j, shockLabels{j});
    fprintf('==============================\n');

    % Statistiche
    fprintf('Skewness = %.4f\n', skewness(x));
    fprintf('Kurtosis = %.4f\n', kurtosis(x));

    % Fit distribuzioni
    distNorm     = fitdist(x,'Normal');
    distT        = fitdist(x,'tLocationScale');
    distLogistic = fitdist(x,'Logistic');

    distFits{j} = {distNorm, distT, distLogistic};

    % ---------------- PLOT marginali ----------------
    figure('Color','w','Position',[100 100 900 500]);

    histogram(x,'Normalization','pdf',...
              'FaceColor',[0.8 0.8 0.8],...
              'EdgeColor',[0.3 0.3 0.3]);
    hold on;

    xg = linspace(min(x),max(x),400);
    plot(xg,pdf(distNorm,xg),'k','LineWidth',2);
    plot(xg,pdf(distT,xg),'b','LineWidth',2);
    plot(xg,pdf(distLogistic,xg),'r-.','LineWidth',2);

    xlabel('Shock value','Color','k');
    ylabel('Density','Color','k');
    title(['Marginal distribution – ' strrep(shockLabels{j},'_',' ')],...
          'FontWeight','bold','Color','k');

    lgd = legend({'Data','Normal','t-Student','Logistic'},...
                 'Location','best','Box','off');
    set(lgd,'TextColor','k','FontWeight','bold');

    % QQ-plot
    figure('Color','w','Position',[200 200 600 500]);
    qqplot(x);
    title(['QQ-plot – ' strrep(shockLabels{j},'_',' ')],...
          'FontWeight','bold','Color','k');
    grid on; box on;
end

% Salva anche distFits
save('structural_shocks.mat','distFits','-append');

fprintf('\n>> Stima marginali completata.\n');