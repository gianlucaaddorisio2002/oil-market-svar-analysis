%% COPULA & ANALISI DENSITA' (PROB.) CONGIUNTA
% Analisi di dipendenza non lineare tra WTI e variabili macro
% (OCSE, Production, Inventories) tramite copule archimedeane
% + stima della densità congiunta empirica.

clear; clc;
load('clean_data.mat','All_d');

%% ================== SCELTA SOTTOINTERVALLO TEMPORALE ====================
% In questo caso sono inutili perchè non ho selezionato un sotto-campione,
% ma utile per testare robustezza modello
startSub = datetime(1990,1,1);
endSub   = datetime(2024,12,31);

All_sub = All_d(timerange(startSub, endSub, 'closed'), :);

fprintf('Intervallo usato per copule: %s – %s\n', ...
        datestr(startSub,'yyyy-mm'), datestr(endSub,'yyyy-mm'));

%% ================== SCELTA DELLE COPPIE DA ANALIZZARE ===================
% Usiamo variazioni log-differenziate (già in All_d)
pairs = {
    'WTI_real_DL', 'OCSE_DL';          % WTI vs domanda globale (OCSE)
    'WTI_real_DL', 'Production_DL';   % WTI vs produzione
    'WTI_real_DL', 'Inventories_DL'   % WTI vs scorte
};

pairNames = {
    'WTI vs OCSE';
    'WTI vs Production';
    'WTI vs Inventories'
};

% Per ogni coppia stimo copule + densità congiunta
for k = 1:size(pairs,1)

    varY = pairs{k,1};   % sempre WTI_real_DL
    varX = pairs{k,2};   % OCSE / Production / Inventories
    niceName = pairNames{k};

    fprintf('\n=============================\n');
    fprintf('Analisi copula: %s\n', niceName);
    fprintf('=============================\n');

    % Estrazione dati grezzi (scala originale delle differenze log)
    Y = All_sub.(varY);
    X = All_sub.(varX);

    % Rimozione NaN
    valid = ~isnan(X) & ~isnan(Y);
    X = X(valid);
    Y = Y(valid);

    n = length(X);
    fprintf('Numero osservazioni utili: %d\n', n);

    %% ===== STEP 1: trasformazione in uniformi (empirical CDF) ===========
    % Le copule lavorano su U,V ~ Uniform(0,1)
    u = tiedrank(X) / (n + 1);
    v = tiedrank(Y) / (n + 1);

    %% ===== STEP 2: stima delle copule archimedeane =======================

    % Clayton: dipendenza in coda inferiore
    theta_clay = copulafit('clayton', [u v]);

    % Gumbel: dipendenza in coda superiore
    theta_gum  = copulafit('gumbel',  [u v]);

    % Frank: dipendenza simmetrica (no prevalenza di coda)
    theta_fr   = copulafit('frank',   [u v]);

    % Log-likelihood (utile per confronto informale tra copule)
    ll_clay = sum( log( copulapdf('clayton',[u v],theta_clay) ) );
    ll_gum  = sum( log( copulapdf('gumbel', [u v],theta_gum ) ) );
    ll_fr   = sum( log( copulapdf('frank',  [u v],theta_fr  ) ) );

    fprintf('Log-likelihood copula Clayton: %.2f\n', ll_clay);
    fprintf('Log-likelihood copula Gumbel : %.2f\n', ll_gum);
    fprintf('Log-likelihood copula Frank  : %.2f\n', ll_fr);

    %% ===== STEP 3: grafico scatter nello spazio uniforme (U,V) ==========
    figure('Name',['Scatter U,V - ' niceName], 'Color','w');
    scatter(u, v, 35, 'filled', 'MarkerFaceAlpha',0.5);
    grid on; box on;
    xlabel('u = F_X(X)', 'FontSize',12);
    ylabel('v = F_Y(Y)', 'FontSize',12);
    title(['Scatter nello spazio uniforme - ' niceName], ...
          'FontWeight','bold','Color', 'k');
    xlim([0 1]); ylim([0 1]);

    %% ===== STEP 4: superfici 3D delle densità di copula =================
    Ugrid = linspace(0.001, 0.999, 50);
    [U,V] = meshgrid(Ugrid, Ugrid);

    C_clay = copulapdf('clayton', [U(:) V(:)], theta_clay);
    C_clay = reshape(C_clay, size(U));

    C_gum  = copulapdf('gumbel',  [U(:) V(:)], theta_gum);
    C_gum  = reshape(C_gum, size(U));

    C_fr   = copulapdf('frank',   [U(:) V(:)], theta_fr);
    C_fr   = reshape(C_fr, size(U));

    % --- Clayton
    figure('Name',['Copula Clayton - ' niceName], 'Color','w');
    surf(U, V, C_clay, 'EdgeColor','none');
    xlabel('u'); ylabel('v'); zlabel('c(u,v)');
    title(['Copula Clayton - ' niceName], 'FontWeight','bold', 'Color', 'k');
    colorbar; view(135,30); grid on;

    % --- Gumbel
    figure('Name',['Copula Gumbel - ' niceName], 'Color','w');
    surf(U, V, C_gum, 'EdgeColor','none');
    xlabel('u'); ylabel('v'); zlabel('c(u,v)');
    title(['Copula Gumbel - ' niceName], 'FontWeight','bold', 'Color', 'k');
    colorbar; view(135,30); grid on;

    % --- Frank
    figure('Name',['Copula Frank - ' niceName], 'Color','w');
    surf(U, V, C_fr, 'EdgeColor','none');
    xlabel('u'); ylabel('v'); zlabel('c(u,v)');
    title(['Copula Frank - ' niceName], 'FontWeight','bold', 'Color', 'k');
    colorbar; view(135,30); grid on;

    %% ===== STEP 5: densità congiunta empirica (X,Y) =====================
    % Kernel density 2D sulla scala originale (Δlog)
    % ATTENZIONE: pura descrizione NON parametrica.
    gridSize = 60;
    xmin = prctile(X,2);  xmax = prctile(X,98);
    ymin = prctile(Y,2);  ymax = prctile(Y,98);

    xgrid = linspace(xmin, xmax, gridSize);
    ygrid = linspace(ymin, ymax, gridSize);
    [Xg,Yg] = meshgrid(xgrid, ygrid);

    % ksdensity bivariata
    [fXY, Xi] = ksdensity([X Y], [Xg(:) Yg(:)]);
    fXY = reshape(fXY, size(Xg));

    figure('Name',['Densità congiunta empirica - ' niceName], 'Color','w');
    surf(Xg, Yg, fXY, 'EdgeColor','none');
    xlabel([varX ' (Δlog)']);
    ylabel([varY ' (Δlog)']);
    zlabel('f_{X,Y}(x,y)');
    title(['Densità congiunta empirica - ' niceName], 'FontWeight','bold');
    colorbar; view(135,30); grid on;

    % Contour plot 2D (più leggibile in tesi)
    figure('Name',['Contour densità congiunta - ' niceName], 'Color','w');
    contourf(Xg, Yg, fXY, 12, 'LineColor','none');
    xlabel([varX ' (Δlog)']);
    ylabel([varY ' (Δlog)']);
    title(['Contour densità congiunta - ' niceName], 'FontWeight','bold');
    colorbar; grid on;

end