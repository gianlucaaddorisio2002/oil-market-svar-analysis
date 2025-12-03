%% ============================================================
% BUILD OIL DATASET – VERSIONE corretta CON OCSE (Aggregato)
%% ============================================================
clear; clc

startD = datetime(1990,1,1);
endD   = datetime(2024,12,31);

%% ============================================================
% 1) LETTURA FILE GREZZI
%% ============================================================

%% --- WTI (FRED ma salvato come semicolon) ---
opts_wti = detectImportOptions('wti.csv','Delimiter',';');
opts_wti.SelectedVariableNames = {'Date','Value'};
wti = readtable('wti.csv', opts_wti);

%% --- CPI (FRED: separato da virgola) ---
cpi = readtable('cpi.csv','Range','A3:B100000','ReadVariableNames',false);
cpi.Properties.VariableNames = {'Date','Value'};

%% --- JET FUEL (EIA) --- 
% jetfuel1.xls: UNA sola colonna di testo del tipo "1991-01,0.741"
Traw = readtable('jetfuel1.xls', ...
                 'ReadVariableNames', false, ...
                 'TextType','string');

col = Traw.Var1;                % unica colonna

% togli eventuali celle vuote
col = col(~ismissing(col));

% se la prima riga è "Datetime,Value", eliminala
if contains(col(1),"Datetime","IgnoreCase",true)
    col = col(2:end);
end

% spezza "1991-01,0.741" -> ["1991-01","0.741"]
parts = split(col, ',' );       % N×2 string array

jet = table();
jet.Date  = parts(:,1);         % string date "1991-01"
jet.Value = parts(:,2);         % string "0.741"



%% --- OCSE (sostituisce REA; serie ocse.csv) ---
% NIENTE hard-code di 'DateTime': leggiamo i nomi e prendiamo le prime 2 colonne
opts_ocse = detectImportOptions('ocse.csv', ...
    'TextType','string', ...
    'VariableNamingRule','preserve');

% Se vuoi vederli a schermo:
% disp('Nomi colonne OCSE:');
% disp(opts_ocse.VariableNames);

% Seleziona semplicemente le prime 2 colonne (data + valore)
opts_ocse.DataLines = 2;                         % salta solo l'header
opts_ocse.SelectedVariableNames = opts_ocse.VariableNames(1:2);

ocse = readtable('ocse.csv', opts_ocse);
ocse.Properties.VariableNames = {'Date','Value'};

%% --- FFR (FRED) ---
ffr = readtable('ffr.csv','Range','A3:B100000','ReadVariableNames',false);
ffr.Properties.VariableNames = {'Date','Value'};

%% --- Inventories (Excel) ---
inv = readtable('inventories.xlsx','Range','A3:B1048576','ReadVariableNames',false);
inv.Properties.VariableNames = {'Date','Value'};

%% --- Production (semicolon + mesi italiani es. "gen-1990") ---
opts_prd = detectImportOptions('production.csv','Delimiter',';','TextType','string');
opts_prd.SelectedVariableNames = opts_prd.VariableNames(1:2);
prd = readtable('production.csv', opts_prd);
prd.Properties.VariableNames = {'Date','Value'};

%% ============================================================
% 2) PULIZIA NUMERICA (numero con virgole → number)
%% ============================================================

fixnum = @(s) str2double( regexprep(string(s), '[,]', '' ) );
jet.Value = fixnum(jet.Value);

wti.Value  = fixnum(wti.Value);
cpi.Value  = fixnum(cpi.Value);
ocse.Value = fixnum(ocse.Value);
ffr.Value  = fixnum(ffr.Value);
inv.Value  = fixnum(inv.Value);
prd.Value  = fixnum(prd.Value);

%% ============================================================
% 3) PARSING DATE (corretto per ogni file)
%% ============================================================

wti.Date  = datetime(string(wti.Date),'InputFormat','yyyy-MM-dd');
cpi.Date  = datetime(string(cpi.Date),'InputFormat','yyyy-MM-dd');
ocse.Date = datetime(string(ocse.Date),'InputFormat','yyyy-MM-dd');
ffr.Date  = datetime(string(ffr.Date),'InputFormat','yyyy-MM-dd');
inv.Date  = datetime(string(inv.Date),'InputFormat','MM/yyyy');
   
jet.Date = datetime(string(jet.Date),'InputFormat','yyyy-MM','Locale','it_IT');


% produzione con mesi italiani (gen, feb, mar, ...)
try
    prd.Date = datetime(string(prd.Date),'InputFormat','MMM-yyyy','Locale','it_IT');
catch
    prd.Date = datetime(string(prd.Date),'InputFormat','MMM-yyyy','Locale','en_US');
end

%% ============================================================
% 4) FORZA PRIMO GIORNO DEL MESE
%% ============================================================

startMonth = @(d) dateshift(d,'start','month');

wti.Date  = startMonth(wti.Date);
cpi.Date  = startMonth(cpi.Date);
ocse.Date = startMonth(ocse.Date);
ffr.Date  = startMonth(ffr.Date);
inv.Date  = startMonth(inv.Date);
prd.Date  = startMonth(prd.Date);
jet.Date = startMonth(jet.Date);


%% ============================================================
% 5) COSTRUZIONE TIMETABLES
%% ============================================================

WTI  = timetable(wti.Date,  wti.Value,  'VariableNames', {'WTI'});
CPI  = timetable(cpi.Date,  cpi.Value,  'VariableNames', {'CPI'});
OCSE = timetable(ocse.Date, ocse.Value, 'VariableNames', {'OCSE'});
FFR  = timetable(ffr.Date,  ffr.Value,  'VariableNames', {'FFR'});
INV  = timetable(inv.Date,  inv.Value,  'VariableNames', {'Inventories'});
PRD  = timetable(prd.Date,  prd.Value,  'VariableNames', {'Production'});

JET = timetable(jet.Date, jet.Value, 'VariableNames', {'JetFuel_USGC'});

% ordina
WTI  = sortrows(WTI);
CPI  = sortrows(CPI);
OCSE = sortrows(OCSE);
FFR  = sortrows(FFR);
INV  = sortrows(INV);
PRD  = sortrows(PRD);
JET = sortrows(JET);

JET = JET(timerange(startD,endD,'closed'),:);

% taglia periodo
WTI  = WTI(timerange(startD,endD,'closed'),:);
CPI  = CPI(timerange(startD,endD,'closed'),:);
OCSE = OCSE(timerange(startD,endD,'closed'),:);
FFR  = FFR(timerange(startD,endD,'closed'),:);
INV  = INV(timerange(startD,endD,'closed'),:);
PRD  = PRD(timerange(startD,endD,'closed'),:);

%% ============================================================
% 6) SINCRONIZZA SU DATE COMUNI
%% ============================================================

All = synchronize(WTI, CPI, OCSE, INV, PRD, FFR, JET, 'intersection');

%% ============================================================
% 7) COSTRUZIONE WTI REALE
%% ============================================================

All.WTI_real = (All.WTI ./ All.CPI) * 100;
All.WTI = [];   % rimuovo la nominale

%% ============================================================
% 8) DATASET STAZIONARIO: All_d (differenze + zscore)
%% ============================================================

wti_d   = diff(log(All.WTI_real));
prod_d  = diff(log(All.Production));
ocse_d  = diff(All.OCSE);           % OCSE è un ciclo → niente log
inv_d   = diff(log(All.Inventories));
jet_d = diff(log(All.JetFuel_USGC));


dates_d = All.Time(2:end);

All_d = timetable(dates_d, wti_d, prod_d, ocse_d, inv_d, jet_d, ...
    'VariableNames', {'WTI_real_DL','Production_DL','OCSE_DL','Inventories_DL', 'JeTFuel_DL'});

All_d{:,:} = zscore(All_d{:,:});
All_d = rmmissing(All_d);

%% ============================================================
% 9) DATASET VAR STRUTTURALE (ALL_VAR)
%    log-livelli + ciclo OCSE
%% ============================================================

All.WTI_log  = log(All.WTI_real);
All.Prod_log = log(All.Production);
All.Inv_log  = log(All.Inventories);

ALL_VAR = All(:, {
    'OCSE', ...        % ciclo OCSE
    'Prod_log', ...
    'WTI_log', ...
    'Inv_log'
});

ALL_VAR.Properties.VariableNames = {
    'OCSE_cycle', ...
    'Prod_log', ...
    'WTI_log', ...
    'Inv_log'
};

%% ============================================================
% 10) SALVA
%% ============================================================

save('clean_data.mat','All','All_d','ALL_VAR','-v7');
disp('Dati salvati in clean_data.mat (All, All_d, ALL_VAR)');
