%% Comparação de Dias: Cálculo Individualizado (Algoritmo SmartVest)
% Este script calcula o RR para cada dia de forma completamente isolada,
% usando EXATAMENTE o mesmo algoritmo do ficheiro 'smartvest.m'.
% Não há partilha de variáveis entre iterações para evitar bias.

clc; clear; close all;

%% 1. SELEÇÃO DE FICHEIROS
fprintf('=== 1. SELEÇÃO DE DADOS ===\n');
fprintf('Selecione os ficheiros FBG (SmartVest) para os vários dias.\n');
[files_fbg, path_fbg] = uigetfile('*.txt', 'Selecione Ficheiros FBG', 'MultiSelect', 'on');

if isequal(files_fbg, 0), return; end
if ischar(files_fbg), files_fbg = {files_fbg}; end
num_dias = length(files_fbg);

% Pergunta sobre Comercial
files_bio = {}; path_bio = '';
resp = questdlg('Tem dados do Sensor Comercial para comparar?', 'Comercial', 'Sim', 'Não', 'Não');
tem_comercial = strcmp(resp, 'Sim');

if tem_comercial
    fprintf('Selecione os ficheiros COMERCIAIS na mesma ordem.\n');
    [files_b, path_b] = uigetfile('*.txt', 'Selecione Ficheiros OpenSignals', 'MultiSelect', 'on');
    if ischar(files_b), files_b = {files_b}; end
    if length(files_b) ~= num_dias
        warning('Atenção: Número de ficheiros comerciais diferente dos FBG.'); 
    end
    files_bio = files_b; path_bio = path_b;
end

% Estrutura para guardar resultados finais
ResultadosFinais = struct();

%% 2. LOOP DE PROCESSAMENTO (DIA A DIA)
for i = 1:num_dias
    % --- INÍCIO DO CÁLCULO ISOLADO PARA O DIA 'i' ---
    
    % A. Preparar Nome
    fName = files_fbg{i};
    [~, raw_name] = fileparts(fName);
    tok = regexp(raw_name, '^(.*?)[\s\-_]*day[\s\-_]*(\d+)', 'tokens', 'ignorecase');
    if ~isempty(tok)
        label_dia = sprintf('%s (Dia %s)', strtrim(tok{1}{1}), tok{1}{2});
    else
        label_dia = strrep(raw_name, '_', ' ');
    end
    fprintf('\n>>> A processar Dia %d: %s\n', i, label_dia);
    
    % B. Algoritmo SmartVest (Cópia Exata da Lógica)
    full_path_fbg = fullfile(path_fbg, fName);
    
    % B1. Ler Ficheiro
    opts = detectImportOptions(full_path_fbg); 
    opts.Delimiter = '\t'; 
    opts.VariableNamingRule = 'preserve';
    try, opts.CommentStyle = '#'; catch, end
    dados = readtable(full_path_fbg, opts);
    
    % B2. Identificar Sensores (Lógica SmartVest)
    col_inicio = 3;
    d_wl = dados{:, col_inicio:end};
    valid_cols = ~all(isnan(d_wl), 1);
    idx_active = find(valid_cols);
    n_sens = length(idx_active);
    
    sinal_fbg_proc = []; % Sinal final para este dia
    
    if n_sens == 3
        % Lógica 3 Canais (Esq, Temp, Dir)
        i_esq = idx_active(1); i_temp = idx_active(2); i_dir = idx_active(3);
        s_esq = d_wl(:, i_esq); s_temp = d_wl(:, i_temp); s_dir = d_wl(:, i_dir);
        
        % Compensação (SmartVest: Subtrair Temp de cada lado)
        s_esq_c = (s_esq - mean(s_esq)) - (s_temp - mean(s_temp));
        s_dir_c = (s_dir - mean(s_dir)) - (s_temp - mean(s_temp));
        
        % Fusão (Média dos dois)
        sinal_fbg_proc = (s_esq_c + s_dir_c) / 2;
        
    elseif n_sens == 2
        % Lógica 2 Canais (Resp, Temp)
        i_resp = idx_active(1); i_temp = idx_active(2);
        s_resp = d_wl(:, i_resp); s_temp = d_wl(:, i_temp);
        sinal_fbg_proc = (s_resp - mean(s_resp)) - (s_temp - mean(s_temp));
        
    else
        % 1 Canal
        sinal_fbg_proc = d_wl(:, idx_active(1));
        sinal_fbg_proc = sinal_fbg_proc - mean(sinal_fbg_proc);
    end
    
    % B3. Tempo e Filtro (Lógica SmartVest)
    col_t = dados{:, 2};
    if iscell(col_t) || isstring(col_t)
        col_t = str2double(strrep(string(col_t), ',', '.'));
    end
    if any(isnan(col_t)), col_t = (1:height(dados))'; end
    
    duracao = col_t(end) - col_t(1);
    if duracao <= 0, duracao = height(dados); end
    Fs_fbg = height(dados) / duracao;
    t_fbg = (0:height(dados)-1) / Fs_fbg;
    
    % Filtro Butterworth (0.1 - 2.0 Hz com Proteção Nyquist)
    f_low = 0.1; f_high = 2.0;
    if (Fs_fbg/2) <= f_high, f_high = (Fs_fbg/2) * 0.9; end
    
    % Proteção extra para Fs muito baixa
    if f_low >= f_high, f_low = f_high * 0.1; end
    
    Wn = [f_low, f_high] / (Fs_fbg/2);
    [b, a] = butter(2, Wn, 'bandpass');
    
    sinal_fbg_filt = filtfilt(b, a, sinal_fbg_proc);
    
    % B4. Deteção de Picos (Parâmetros SmartVest)
    dist_min = Fs_fbg * 1.5; % 1.5 seg
    prom_min = std(sinal_fbg_filt) * 0.5;
    
    % Nota: Usamos MinPeakHeight também para garantir consistência
    altura_min = max(sinal_fbg_filt) * 0.15;
    
    [pks, locs] = findpeaks(sinal_fbg_filt, ...
        'MinPeakDistance', dist_min, ...
        'MinPeakProminence', prom_min, ...
        'MinPeakHeight', altura_min);
    
    rr_fbg = 0;
    if ~isempty(locs)
        rr_fbg = mean(60 ./ diff(t_fbg(locs)));
    end
    
    % --- CÁLCULO COMERCIAL (Se existir para este dia) ---
    rr_bio = 0; t_bio = []; s_bio_filt = []; has_bio = false;
    
    if tem_comercial && i <= length(files_bio)
        try
            fNameBio = files_bio{i};
            full_path_bio = fullfile(path_bio, fNameBio);
            
            % Ler Sampling Rate
            fid = fopen(full_path_bio); Fs_bio=1000; found=false;
            while ~feof(fid)
                l = fgetl(fid); 
                if startsWith(l,'#') && contains(l,'sampling rate')
                    tk = regexp(l,'"sampling rate":\s*(\d+)','tokens');
                    if ~isempty(tk), Fs_bio=str2double(tk{1}{1}); found=true; end
                elseif ~startsWith(l,'#'), break; end
            end
            fclose(fid);
            
            optsB = detectImportOptions(full_path_bio); 
            optsB.CommentStyle='#'; optsB.VariableNamingRule='preserve';
            d_bio = readtable(full_path_bio, optsB);
            raw_bio = d_bio{:, end}; % Assumir última coluna
            
            t_bio = (0:length(raw_bio)-1)/Fs_bio;
            
            % Filtro Igual
            f_h_b = 2.0; if (Fs_bio/2)<=f_h_b, f_h_b=(Fs_bio/2)*0.9; end
            Wn_b = [0.1, f_h_b]/(Fs_bio/2);
            [bb, ab] = butter(2, Wn_b, 'bandpass');
            s_bio_filt = filtfilt(bb, ab, raw_bio - mean(raw_bio));
            
            % Picos Igual
            d_min_b = Fs_bio * 1.5;
            p_min_b = std(s_bio_filt) * 0.5;
            alt_min_b = max(s_bio_filt) * 0.15;
            
            [~, l_b] = findpeaks(s_bio_filt, ...
                'MinPeakDistance', d_min_b, ...
                'MinPeakProminence', p_min_b, ...
                'MinPeakHeight', alt_min_b);
            
            if ~isempty(l_b)
                rr_bio = mean(60 ./ diff(t_bio(l_b)));
            end
            has_bio = true;
            
        catch ME
            warning('Erro no comercial dia %d: %s', i, ME.message);
        end
    end
    
    % --- GUARDAR RESULTADOS DESTE DIA ---
    ResultadosFinais(i).DiaLabel = label_dia;
    ResultadosFinais(i).RR_FBG = rr_fbg;
    ResultadosFinais(i).RR_BIO = rr_bio;
    ResultadosFinais(i).SinalFBG = sinal_fbg_filt;
    ResultadosFinais(i).TempoFBG = t_fbg;
    ResultadosFinais(i).SinalBIO = s_bio_filt;
    ResultadosFinais(i).TempoBIO = t_bio;
    ResultadosFinais(i).TemBio = has_bio;
    
    fprintf('   -> FBG: %.2f BPM | BIO: %.2f BPM\n', rr_fbg, rr_bio);
end

%% 3. VISUALIZAÇÃO COMPARATIVA
if isempty(fieldnames(ResultadosFinais)), return; end

figure('Color', 'white', 'Position', [50, 50, 1200, 800], 'Name', 'Comparacao Dias Final');

% A. Gráfico de Barras (RR)
subplot(2, 2, [1, 3]);
vals_f = [ResultadosFinais.RR_FBG];
vals_b = [ResultadosFinais.RR_BIO];

% --- EDIÇÃO DE ETIQUETAS DO EIXO X ---
% Por defeito, usa os nomes extraídos dos ficheiros.
% Se quiser editar manualmente, descomente e altere a linha abaixo:
% lbls = {'Dia 1 (Joao)', 'Dia 2 (Ana)', 'Dia 3 (Teste)'}; 
lbls = {ResultadosFinais.DiaLabel};

if tem_comercial
    b = bar([vals_f', vals_b']);
    legend({'SmartVest (FBG)', 'Comercial'}, 'Location', 'northwest');
    b(1).FaceColor = [0 0.4470 0.7410]; 
    b(2).FaceColor = [0.8500 0.3250 0.0980];
else
    bar(vals_f); legend('SmartVest (FBG)');
end

title('Comparação de Frequência Respiratória (BPM) por Dia');
ylabel('BPM'); 
xticklabels(lbls); % Aplica as etiquetas definidas acima
xtickangle(45); grid on;

% Etiquetas nas barras
for k = 1:num_dias
    text(k-0.15, vals_f(k), sprintf('%.1f', vals_f(k)), 'Vert','bottom','Horiz','center','FontWeight','bold');
    if tem_comercial && ResultadosFinais(k).TemBio
        text(k+0.15, vals_b(k), sprintf('%.1f', vals_b(k)), 'Vert','bottom','Horiz','center','FontWeight','bold');
    end
end

% B. Mini-Gráficos dos Sinais (Primeiros 30s)
n_show = min(4, num_dias);
for k = 1:n_show
    subplot(4, 2, k*2);
    
    % FBG
    tf = ResultadosFinais(k).TempoFBG; 
    sf = ResultadosFinais(k).SinalFBG;
    % Normalizar (Z-score) para ver forma
    sf = (sf - mean(sf)); if std(sf)~=0, sf=sf/std(sf); end
    
    plot(tf, sf, 'b', 'LineWidth', 1.5); hold on;
    
    % BIO
    if ResultadosFinais(k).TemBio
        tb = ResultadosFinais(k).TempoBIO;
        sb = ResultadosFinais(k).SinalBIO;
        sb = (sb - mean(sb)); if std(sb)~=0, sb=sb/std(sb); end
        
        plot(tb, sb, 'Color', [0.85 0.33 0.1], 'LineWidth', 1.5);
        if k==1, legend({'FBG', 'Comercial'}, 'FontSize', 8); end
        
        lim_x = min(30, min(tf(end), tb(end)));
        xlim([0, lim_x]);
    else
        if k==1, legend({'FBG'}, 'FontSize', 8); end
        xlim([0, min(30, tf(end))]);
    end
    
    title(strrep(ResultadosFinais(k).DiaLabel, sprintf('\n'), ' '));
    grid on; axis tight;
    % Define "Tempo (s)" como a etiqueta do eixo X para todos os mini-gráficos
    xlabel('Tempo (s)'); 
    
    if k~=n_show
        set(gca,'XTickLabel',[]); % Esconde números do eixo X exceto no último
    end
end

fprintf('Análise Completa.\n');