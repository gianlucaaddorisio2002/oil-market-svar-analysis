%% OLS a finestre temporali – analisi di stabilità
clear; clc;
load('clean_data.mat','All_d');

% All_d è un timetable con rowtimes = dates_d
allStart = All_d.Properties.RowTimes(1);
allEnd   = All_d.Properties.RowTimes(end);


fprintf('Campione completo: %s – %s (%d osservazioni)\n', ...
    datestr(allStart), datestr(allEnd), height(All_d));

%% 1) Definisci le finestre temporali
windows = struct( ...
    'name',  {}, ...
    'start', {}, ...
    'end',   {} );

windows(1).name  = '1990–1999';
windows(1).start = datetime(1990,1,1);
windows(1).end   = datetime(1999,12,31);

windows(2).name  = '2000–2007';
windows(2).start = datetime(2000,1,1);
windows(2).end   = datetime(2007,12,31);

windows(3).name  = '2008–2014';
windows(3).start = datetime(2008,1,1);
windows(3).end   = datetime(2014,12,31);

windows(4).name  = '2015–2024';
windows(4).start = datetime(2015,1,1);
windows(4).end   = datetime(2024,12,31);

nW = numel(windows);

%% 2) Prealloca tabella risultati
results = table( ...
    'Size',[nW 11], ...
    'VariableTypes', {'string','datetime','datetime','double', ...
                      'double','double','double','double', ...
                      'double','double','double'}, ...
    'VariableNames', {'Window','StartDate','EndDate','N', ...
                      'R2','AdjR2','pWald','JB_p', ...
                      'KS_p','DW','LB_min_p'} );

%% 3) Matrice R per il test di Wald (same as ols_prevar_check)
% Coefficienti in mdl.Coefficients:
% 1: Intercetta
% 2: Production_DL
% 3: OCSE_DL
% 4: Inventories_DL
R = [0 1 0 0;   % Production_DL = 0
     0 0 1 0;   % OCSE_DL        = 0
     0 0 0 1];  % Inventories_DL= 0
r = [0; 0; 0];

%% 4) Loop sulle finestre
for k = 1:nW
    wname  = windows(k).name;
    wstart = windows(k).start;
    wend   = windows(k).end;

    % Sotto-campione della finestra
    Sub = All_d(timerange(wstart, wend, 'closed'), :);

    fprintf('\n===== Finestra %s (%s – %s) =====\n', ...
        wname, datestr(wstart), datestr(wend));

    if height(Sub) < 30
        fprintf('Finestra troppo corta (%d osservazioni). Salto.\n', height(Sub));
        results.Window(k)    = wname;
        results.StartDate(k) = wstart;
        results.EndDate(k)   = wend;
        results.N(k)         = height(Sub);
        continue;
    end

    % Converte in table per fitlm
    Tsub = timetable2table(Sub);
    Tsub(:,1) = []; % rimuove colonna date

    % Modello OLS nella finestra
    mdl = fitlm(Tsub, 'WTI_real_DL ~ Production_DL + OCSE_DL + Inventories_DL');
    disp(mdl);

    % Residui
    e = mdl.Residuals.Raw;
    e = double(e(:));
    e = e(~isnan(e));
    e_std = (e - mean(e))/std(e);

    % R^2
    R2    = mdl.Rsquared.Ordinary;
    R2adj = mdl.Rsquared.Adjusted;

    % Test WALD globale
    pWald = coefTest(mdl, R, r);

    % Normalità: JB + KS
    [p_jb, ~] = jbtest(e);
    [~, p_ks] = kstest(e_std);  % vs N(0,1)

    % Autocorrelazione: Durbin–Watson + Ljung–Box (1:12)
    num = sum(diff(e).^2);
    den = sum(e.^2);
    DW  = num / den;

    [~, p_lb] = lbqtest(e, "Lags", 1:12);
    LB_min = min(p_lb);

    % Riempie tabella risultati
    results.Window(k)    = wname;
    results.StartDate(k) = wstart;
    results.EndDate(k)   = wend;
    results.N(k)         = height(Sub);
    results.R2(k)        = R2;
    results.AdjR2(k)     = R2adj;
    results.pWald(k)     = pWald;
    results.JB_p(k)      = p_jb;
    results.KS_p(k)      = p_ks;
    results.DW(k)        = DW;
    results.LB_min_p(k)  = LB_min;

    fprintf('N = %d  |  R2 = %.4f  |  AdjR2 = %.4f  |  pWald = %.4f  |  LB_min_p = %.4f\n', ...
        height(Sub), R2, R2adj, pWald, LB_min);
end

%% 5) Mostra riepilogo
fprintf('\n===== RIEPILOGO FINESTRE OLS =====\n');
disp(results);

% (Opzionale) salva su file .mat e/o .csv
save('ols_windows_results.mat','results');
writetable(results,'ols_windows_results.csv');

fprintf('Risultati salvati in ols_windows_results.mat e .csv\n');