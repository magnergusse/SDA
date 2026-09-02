%% Comparação de Dias: FBG vs Sensor Comercial (Sincronizado)
% Este script usa EXATAMENTE a mesma lógica de processamento do
% script "fbg_respiratory_analysis.m" para garantir consistência nos BPMs.

%clc; clear; close all;

%% 1. Seleção de Ficheiros FBG
fprintf('=== PASSO 1: Selecionar ficheiros FBG (Os seus dados) ===\n');
[files_fbg, path_fbg] = uigetfile('*.txt', 'Selecione ficheiros FBG (Múltiplos)', 'MultiSelect', 'on');

if isequal(files_fbg, 0), disp('Cancelado.'); return; end
if ischar(files_fbg), files_fbg = {files_fbg}; end 

num_dias = length(files_fbg);
fprintf('Selecionados %d ficheiros FBG.\n', num_dias);

%% 2. Seleção de Ficheiros Comerciais (Opcional)
tem_comercial = false;
files_bio = {};
path_bio = '';

resp = questdlg('Tem ficheiros do sensor Comercial para estes dias?', 'Comercial', 'Sim', 'Não', 'Não');

if strcmp(resp, 'Sim')
    fprintf('\n=== PASSO 2: Selecionar ficheiros COMERCIAIS correspondentes ===\n');
    [files_b, path_b] = uigetfile('*.txt', 'Selecione ficheiros OpenSignals', 'MultiSelect', 'on');
    
    if ~isequal(files_b, 0)
        if ischar(files_b), files_b = {files_b}; end
        if length(files_b) ~= num_dias
            warning('Aviso: Número de ficheiros comerciais diferente dos FBG.');
        end
        files_bio = files_b;
        path_bio = path_b;
        tem_comercial = true;
    end
end

% Estrutura de Resultados
resultados = struct(); 

% Parâmetros de Filtro (Iguais à Análise Completa)
f_low_target = 0.1; 
f_high_target = 2.0;

%% 3. Processamento (Dia a Dia)
for i = 1:num_dias
    % --- A. Processar Nome e Ler Ficheiro ---
    fname_fbg = files_fbg{i};
    [~, raw_name] = fileparts(fname_fbg);
    token = regexp(raw_name, '^(.*?)[\s\-_]*day[\s\-_]*(\d+)', 'tokens', 'ignorecase');
    if ~isempty(token)
        label_dia = sprintf('%s\n(Dia %s)', strtrim(token{1}{1}), token{1}{2});
    else
        label_dia = strrep(raw_name, '_', ' '); 
    end
    
    full_fbg = fullfile(path_fbg, fname_fbg);
    fprintf('\n--- Processando: %s ---\n', label_dia);
    
    opts = detectImportOptions(full_fbg); opts.Delimiter = '\t'; opts.VariableNamingRule = 'preserve';
    try, opts.CommentStyle = '#'; catch, end
    dad_f = readtable(full_fbg, opts);
    
    % --- B. LÓGICA DE SENSORES (Sincronizada com Análise Completa) ---
    col_inicio_wl = 3; 
    dados_wl = dad_f{:, col_inicio_wl:end};
    colunas_validas = ~all(isnan(dados_wl), 1); 
    indices_ativos = find(colunas_validas);
    n_sens = length(indices_ativos);
    
    if n_sens == 3
        % Lógica Específica: 1=Esq, 2=Temp, 3=Dir
        idx_esq = indices_ativos(1); 
        idx_temp = indices_ativos(2); 
        idx_dir = indices_ativos(3);
        
        s_esq = dados_wl(:, idx_esq);
        s_temp = dados_wl(:, idx_temp);
        s_dir = dados_wl(:, idx_dir);
        
        % Compensação Individual e Fusão
        s_esq_comb = (s_esq - mean(s_esq)) - (s_temp - mean(s_temp));
        s_dir_comb = (s_dir - mean(s_dir)) - (s_temp - mean(s_temp));
        
        raw_sig = (s_esq_comb + s_dir_comb) / 2;
        
    elseif n_sens == 2
        % Lógica Específica: 1=Resp, 2=Temp
        idx_resp = indices_ativos(1); 
        idx_temp = indices_ativos(2);
        s_resp = dados_wl(:, idx_resp);
        s_temp = dados_wl(:, idx_temp);
        
        raw_sig = (s_resp - mean(s_resp)) - (s_temp - mean(s_temp));
    else
        % 1 Canal
        raw_sig = dados_wl(:, indices_ativos(1));
        raw_sig = raw_sig - mean(raw_sig);
    end
    
    % --- C. Tempo e Filtro ---
    col_t = dad_f{:, 2};
    if iscell(col_t) || isstring(col_t), col_t = str2double(strrep(string(col_t), ',', '.')); end
    if any(isnan(col_t)), col_t = (1:height(dad_f))'; end
    dur = col_t(end)-col_t(1); if dur<=0, dur=height(dad_f); end
    Fs_f = height(dad_f)/dur;
    t_f = (0:height(dad_f)-1)/Fs_f;
    
    % Correção Nyquist
    f_high_f = f_high_target;
    if (Fs_f / 2) <= f_high_f, f_high_f = (Fs_f / 2) * 0.95; end
    
    Wn = [f_low_target, f_high_f]/(Fs_f/2); 
    [b,a]=butter(2,Wn,'bandpass');
    sig_f = filtfilt(b,a, raw_sig);
    
    % --- D. DETEÇÃO DE PICOS (Sincronizada - Estratégia Avançada) ---
    min_dist_sec = 1.0; 
    dist_min = Fs_f * min_dist_sec;
    prom_min = std(sig_f) * 0.5;
    altura_min = max(sig_f) * 0.15; % Adicionado: Filtra picos pequenos
    
    [pks_f, locs_f] = findpeaks(sig_f, ...
        'MinPeakDistance', dist_min, ...
        'MinPeakProminence', prom_min, ...
        'MinPeakHeight', altura_min);
        
    rr_f = 0; if ~isempty(locs_f), rr_f = mean(60./diff(t_f(locs_f))); end
    
    
    % --- E. PROCESSAR COMERCIAL ---
    rr_b = 0; sig_b = []; t_b = []; has_bio = false;
    
    if tem_comercial && i <= length(files_bio)
        try
            fname_bio = files_bio{i};
            full_bio = fullfile(path_bio, fname_bio);
            
            % Ler Sampling Rate
            fid=fopen(full_bio); Fs_b=1000; found=false;
            while ~feof(fid)
                l=fgetl(fid); if startsWith(l,'#') && contains(l,'sampling rate')
                   t=regexp(l,'"sampling rate":\s*(\d+)','tokens');
                   if ~isempty(t), Fs_b=str2double(t{1}{1}); found=true; end
                elseif ~startsWith(l,'#'), break; end
            end
            fclose(fid);
            
            optsB = detectImportOptions(full_bio); optsB.CommentStyle='#'; optsB.VariableNamingRule='preserve';
            dad_b = readtable(full_bio, optsB);
            raw_b = dad_b{:, end}; 
            
            t_b = (0:length(raw_b)-1)/Fs_b;
            
            f_high_b = f_high_target;
            if (Fs_b/2) <= f_high_b, f_high_b = (Fs_b/2)*0.95; end
            
            Wn_b = [f_low_target, f_high_b]/(Fs_b/2); 
            [bb,ab]=butter(2,Wn_b,'bandpass');
            sig_b = filtfilt(bb,ab, raw_b - mean(raw_b));
            
            % Picos Comercial (Mesma estratégia)
            dist_min_b = Fs_b * min_dist_sec;
            prom_min_b = std(sig_b) * 0.5;
            altura_min_b = max(sig_b) * 0.15;
            
            [~, locs_b] = findpeaks(sig_b, ...
                'MinPeakDistance', dist_min_b, ...
                'MinPeakProminence', prom_min_b, ...
                'MinPeakHeight', altura_min_b);
                
            if ~isempty(locs_b), rr_b = mean(60./diff(t_b(locs_b))); end
            has_bio = true;
        catch
            warning('Erro ao processar comercial do dia %d', i);
        end
    end
    
    % Guardar
    resultados(i).dia = label_dia;
    resultados(i).RR_fbg = rr_f;
    resultados(i).RR_bio = rr_b;
    resultados(i).sig_f = sig_f;
    resultados(i).t_f = t_f;
    resultados(i).sig_b = sig_b;
    resultados(i).t_b = t_b;
    resultados(i).has_bio = has_bio;
end

%% 4. Visualização
figure('Color', 'white', 'Position', [50, 50, 1200, 800], 'Name', 'Comparacao Sincronizada');

% A. Barras RR
subplot(2, 2, [1, 3]); 
vals_f = [resultados.RR_fbg];
vals_b = [resultados.RR_bio];
labels = {resultados.dia};

if tem_comercial
    b = bar([vals_f', vals_b']);
    legend({'FBG', 'Comercial'}, 'Location', 'northwest');
    b(1).FaceColor = [0 0.4470 0.7410]; 
    b(2).FaceColor = [0.8500 0.3250 0.0980];
else
    bar(vals_f); legend('FBG');
end
title('Comparação de RR (Algoritmo Unificado)');
ylabel('BPM'); xticklabels(labels); xtickangle(45); grid on;

for k = 1:num_dias
    text(k-0.15, vals_f(k), sprintf('%.1f', vals_f(k)), 'Vert','bottom','Horiz','center','FontWeight','bold');
    if tem_comercial && resultados(k).has_bio
        text(k+0.15, vals_b(k), sprintf('%.1f', vals_b(k)), 'Vert','bottom','Horiz','center','FontWeight','bold');
    end
end

% B. Sinais
n_plots = min(4, num_dias);
for k = 1:n_plots
    subplot(4, 2, k*2);
    tf = resultados(k).t_f; sf = resultados(k).sig_f;
    sf_n = (sf-mean(sf)); if std(sf_n)~=0, sf_n=sf_n/std(sf_n); end
    
    plot(tf, sf_n, 'b', 'LineWidth', 1.5); hold on;
    
    if resultados(k).has_bio
        tb = resultados(k).t_b; sb = resultados(k).sig_b;
        sb_n = (sb-mean(sb)); if std(sb_n)~=0, sb_n=sb_n/std(sb_n); end
        plot(tb, sb_n, 'Color', [0.85 0.33 0.1], 'LineWidth', 1.5);
        if k==1, legend({'FBG', 'Comercial'}, 'FontSize', 8, 'Location', 'best'); end
        xlim([0, min(30, min(tf(end), tb(end)))]);
    else
        if k==1, legend({'FBG'}, 'FontSize', 8, 'Location', 'best'); end
        xlim([0, min(30, tf(end))]);
    end
    
    title(strrep(resultados(k).dia, sprintf('\n'), ' '));
    grid on; axis tight;
    if k==n_plots, xlabel('Tempo (s)'); else, set(gca,'XTickLabel',[]); end
end
fprintf('Concluído.\n');