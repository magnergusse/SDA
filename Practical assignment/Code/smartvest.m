%% Análise de Sensores Respiratórios FBG + Comparação Biosignalsplux
% Este script processa dados FBG (Esquerda/Direita/Temp) e compara com
% um sensor comercial (Biosignalsplux).

clc; clear; close all;

%% ========================================================================
%% PARTE 1: ANÁLISE DOS DADOS FBG (Série Temporal)
%% ========================================================================
fprintf('--- PARTE 1: DADOS FBG ---\n');
[ficheiro, caminho] = uigetfile('*.txt', '1. Selecione o ficheiro FBG (Ex: Mediçao3FBG.txt)');
if isequal(ficheiro, 0)
    disp('Cancelado.');
    return;
end
caminho_completo = fullfile(caminho, ficheiro);
fprintf('A carregar FBG: %s ...\n', ficheiro);

% Ler ficheiro
opcoes = detectImportOptions(caminho_completo);
opcoes.Delimiter = '\t'; 
opcoes.VariableNamingRule = 'preserve'; 
try opcoes.CommentStyle = '#'; catch, end 
dados = readtable(caminho_completo, opcoes);

% Identificação dos Sensores (Ignora primeiras 2 colunas)
col_inicio_wl = 3; 
dados_wl = dados{:, col_inicio_wl:end};
colunas_validas = ~all(isnan(dados_wl), 1); 
indices_ativos = find(colunas_validas);
num_sensores = length(indices_ativos);

fprintf('FBG: Detetados %d sensores ativos.\n', num_sensores);

% Variáveis para guardar os sinais individuais
sinal_filt_esq = [];
sinal_filt_dir = [];

if num_sensores == 3
    % --- CONFIGURAÇÃO DE 3 CANAIS ---
    % Canal 1: Esq, Canal 2: Temp, Canal 3: Dir
    idx_esq  = indices_ativos(1); 
    idx_temp = indices_ativos(2); 
    idx_dir  = indices_ativos(3); 
    
    sinal_esq  = dados_wl(:, idx_esq);
    sinal_temp = dados_wl(:, idx_temp);
    sinal_dir  = dados_wl(:, idx_dir);
    
    nomes = dados.Properties.VariableNames(col_inicio_wl + indices_ativos - 1);
    fprintf('Config: Esq=%s, Temp=%s, Dir=%s\n', nomes{1}, nomes{2}, nomes{3});
    
    % Compensação (Subtrair Temperatura)
    sinal_esq_cent = sinal_esq - mean(sinal_esq);
    sinal_dir_cent = sinal_dir - mean(sinal_dir);
    sinal_temp_cent = sinal_temp - mean(sinal_temp);
    
    % Sinais Individuais Compensados
    sinal_esq_comb = sinal_esq_cent - sinal_temp_cent;
    sinal_dir_comb = sinal_dir_cent - sinal_temp_cent;
    
    % Sinal Combinado (Média dos dois lados)
    sinal_bruto_combinado = (sinal_esq_comb + sinal_dir_comb) / 2;
    
elseif num_sensores == 2
    % 2 Canais (Resp, Temp)
    idx_resp = indices_ativos(1); idx_temp = indices_ativos(2);
    sinal_resp = dados_wl(:, idx_resp); sinal_temp = dados_wl(:, idx_temp);
    
    sinal_bruto_combinado = (sinal_resp - mean(sinal_resp)) - (sinal_temp - mean(sinal_temp));
    fprintf('Config: 2 Canais (Resp + Temp)\n');
else
    % 1 Canal
    idx_ativo = indices_ativos(1);
    sinal_bruto_combinado = dados_wl(:, idx_ativo);
    fprintf('Config: 1 Canal (Sem compensação)\n');
end

% Processamento Temporal
coluna_tempo = dados{:, 2};
if iscell(coluna_tempo) || isstring(coluna_tempo)
    coluna_tempo = strrep(string(coluna_tempo), ',', '.');
    tempo_bruto_fbg = str2double(coluna_tempo);
else
    tempo_bruto_fbg = double(coluna_tempo);
end

% Correção de tempo se houver NaNs
if any(isnan(tempo_bruto_fbg)), tempo_bruto_fbg = (1:height(dados))'; end

duracao = tempo_bruto_fbg(end) - tempo_bruto_fbg(1);
if duracao <= 0, duracao = height(dados); end

Fs_fbg = height(dados) / duracao;
t_fbg = (0:height(dados)-1) / Fs_fbg;
fprintf('FBG Fs estimada: %.2f Hz\n', Fs_fbg);

% Filtros (0.1 Hz a 2.0 Hz conforme seu pedido)
freq_baixa = 0.1;
freq_alta = 2.0;
if Fs_fbg/2 <= freq_alta, freq_alta = (Fs_fbg/2) * 0.9; end

Wn = [freq_baixa, freq_alta] / (Fs_fbg/2);
[b, a] = butter(2, Wn, 'bandpass');

% Aplicar filtros
sinal_filtrado = filtfilt(b, a, sinal_bruto_combinado); % Média/Fusão

if num_sensores == 3
    sinal_filt_esq = filtfilt(b, a, sinal_esq_comb);
    sinal_filt_dir = filtfilt(b, a, sinal_dir_comb);
end

% --- ESTRATÉGIA AVANÇADA DE DETEÇÃO DE PICOS (FBG) ---
min_dist_sec = 3; 
dist_min = Fs_fbg * min_dist_sec;
prom_min = std(sinal_filtrado) * 0.5;
altura_min = max(sinal_filtrado) * 0.15;

[picos, locs] = findpeaks(sinal_filtrado, ...
    'MinPeakDistance', dist_min, ...
    'MinPeakProminence', prom_min, ...
    'MinPeakHeight', altura_min);

if isempty(locs), rr_media = 0; else, rr_media = mean(60 ./ diff(t_fbg(locs))); end
amp_media = mean(picos);

fprintf('FBG RR: %.1f BPM | Amp: %.4f nm\n', rr_media, amp_media);


%% ========================================================================
%% PARTE 2: DADOS BIOSIGNALSPLUX (Comercial)
%% ========================================================================
fprintf('\n--- PARTE 2: SENSOR COMERCIAL ---\n');
resp = questdlg('Carregar dados Biosignalsplux?', 'Comparar', 'Sim', 'Não', 'Sim');

tem_comercial = false;
t_bio = []; sinal_bio_filt = []; rr_bio = 0; nome_bio = '';

if strcmp(resp, 'Sim')
    [fBio, pBio] = uigetfile('*.txt', '2. Selecione ficheiro OpenSignals');
    if ~isequal(fBio, 0)
        file_bio = fullfile(pBio, fBio);
        fprintf('A ler Biosignals: %s ...\n', fBio);
        
        % 1. Ler Sampling Rate (Parser Manual)
        fid = fopen(file_bio, 'rt');
        Fs_bio = 1000; found_Fs = false;
        while true
            line = fgetl(fid);
            if ~ischar(line), break; end
            if startsWith(line, '#') && contains(line, 'sampling rate')
                tok = regexp(line, '"sampling rate":\s*(\d+)', 'tokens');
                if ~isempty(tok), Fs_bio = str2double(tok{1}{1}); found_Fs = true; end
            end
            if ~startsWith(line, '#'), break; end
        end
        fclose(fid);
        
        if ~found_Fs
            inp = inputdlg('Frequência (Hz):', 'Fs', 1, {'1000'});
            if ~isempty(inp), Fs_bio = str2double(inp{1}); end
        end
        
        % 2. Ler Dados
        optsBio = detectImportOptions(file_bio);
        optsBio.CommentStyle = '#';
        optsBio.VariableNamingRule = 'preserve';
        dados_bio = readtable(file_bio, optsBio);
        
        % Assumir última coluna como sinal
        sinal_bio_bruto = dados_bio{:, end};
        nome_bio = dados_bio.Properties.VariableNames{end};
        
        t_bio = (0:length(sinal_bio_bruto)-1) / Fs_bio;
        
        % 3. Processar
        Wn_bio = [freq_baixa, freq_alta] / (Fs_bio/2);
        [b_bio, a_bio] = butter(2, Wn_bio, 'bandpass');
        sinal_bio_filt = filtfilt(b_bio, a_bio, sinal_bio_bruto - mean(sinal_bio_bruto));
        
        % --- ESTRATÉGIA AVANÇADA DE DETEÇÃO DE PICOS (COMERCIAL) ---
        dist_min_bio = Fs_bio * min_dist_sec;
        prom_min_bio = std(sinal_bio_filt) * 0.5;
        altura_min_bio = max(sinal_bio_filt) * 0.15;

        [picos_bio, locs_bio] = findpeaks(sinal_bio_filt, ...
            'MinPeakDistance', dist_min_bio, ...
            'MinPeakProminence', prom_min_bio, ...
            'MinPeakHeight', altura_min_bio);
        
        if isempty(locs_bio), rr_bio = 0; else, rr_bio = mean(60 ./ diff(t_bio(locs_bio))); end
        tem_comercial = true;
        fprintf('Biosignals RR: %.1f BPM\n', rr_bio);
    end
end


%% ========================================================================
%% PARTE 3: VISUALIZAÇÃO UNIFICADA (TEMPO)
%% ========================================================================
figure('Color', 'white', 'Position', [50, 50, 1000, 900], 'Name', 'Analise FBG vs Comercial');

% --- SUBPLOT 1: DETALHE FBG (Esq vs Dir vs Média) ---
ax1 = subplot(3,1,1);
hold on;
if num_sensores == 3
    plot(t_fbg, sinal_filt_esq, 'Color', [0 0.5 1], 'LineWidth', 1); % Azul (Esq)
    plot(t_fbg, sinal_filt_dir, 'Color', [0 0.8 0], 'LineWidth', 1); % Verde (Dir)
    plot(t_fbg, sinal_filtrado, 'k-', 'LineWidth', 2); % Preto (Média)
    legend('Esquerda', 'Direita', 'Média (Fusão)');
else
    plot(t_fbg, sinal_filtrado, 'k-', 'LineWidth', 1.5);
    legend('Sinal FBG');
end
plot(t_fbg(locs), picos, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5,'DisplayName','Picos');
ylabel('Amplitude (nm)');
title(sprintf('1. FBG: Comparação Lados | RR: %.1f BPM', rr_media));
grid on; axis tight;

% --- SUBPLOT 2: SINAL COMERCIAL ---
ax2 = subplot(3,1,2);
if tem_comercial
    plot(t_bio, sinal_bio_filt, 'Color', [0.85 0.33 0.1], 'LineWidth', 1.5); hold on;
    plot(t_bio(locs_bio), picos_bio, 'mo', 'MarkerFaceColor', 'm', 'MarkerSize', 5);
    ylabel('Amplitude (u.a.)');
    title(sprintf('2. SENSOR COMERCIAL (%s) | RR: %.1f BPM', nome_bio, rr_bio));
    legend('Biosignalsplux', 'Picos');
else
    text(0.5, 0.5, 'Nenhum dado comercial carregado', 'HorizontalAlignment', 'center');
    title('2. Sensor Comercial (Dados em falta)');
end
xlabel('Tempo (s)');
grid on; axis tight;

% --- SUBPLOT 3: SOBREPOSIÇÃO (COMPARACAO DIRETA) ---
ax3 = subplot(3,1,3);
if tem_comercial
    fbg_norm = (sinal_filtrado - mean(sinal_filtrado));
    if std(fbg_norm) ~= 0, fbg_norm = fbg_norm / std(fbg_norm); end
    bio_norm = (sinal_bio_filt - mean(sinal_bio_filt));
    if std(bio_norm) ~= 0, bio_norm = bio_norm / std(bio_norm); end
    
    plot(t_fbg, fbg_norm, 'b-', 'LineWidth', 1.5); hold on;
    plot(t_bio, bio_norm, 'Color', [0.85 0.33 0.1], 'LineWidth', 1.5); 
    
    ylabel('Amp (Norm)'); xlabel('Tempo (s)');
    title('3. Sobreposição Temporal: FBG (Azul) vs Comercial (Laranja)');
    legend('FBG (Média)', 'Sensor Comercial');
    grid on; axis tight;
else
    text(0.5, 0.5, 'Dados comerciais não disponíveis', 'HorizontalAlignment', 'center');
    title('3. Sobreposição FBG vs Comercial');
end

if tem_comercial
    linkaxes([ax1, ax2, ax3], 'x');
    max_t = min(t_fbg(end), t_bio(end));
    xlim([0, max_t]);
else
    xlim([0, t_fbg(end)]);
end


%% ========================================================================
%% PARTE 4: ANÁLISE ESPECTRAL (FREQ. RESPIRATÓRIA - FFT)
%% ========================================================================
% Cria uma NOVA janela para comparar os espectros de frequência
figure('Color', 'white', 'Position', [150, 150, 1000, 800], 'Name', 'Comparacao Espectral (FFT)');

% Definir função simples para FFT (Potência)
calc_psd = @(x, fs) abs(fft(x - mean(x))).^2 / length(x);
% Eixo de frequência
get_freq_axis = @(x, fs) (0:length(x)-1)*(fs/length(x));

% --- SUBPLOT 1: Comparação Esquerda vs Direita vs Média (FBG) ---
subplot(2,1,1);
hold on;

% 1. Média
N_fbg = length(sinal_filtrado);
f_axis_fbg = get_freq_axis(sinal_filtrado, Fs_fbg);
P_fbg = calc_psd(sinal_filtrado, Fs_fbg);
% Cortar para 0-2 Hz (range respiratório)
idx_lim_fbg = f_axis_fbg <= 2.0;

if num_sensores == 3
    % Calcular para Esq e Dir
    P_esq = calc_psd(sinal_filt_esq, Fs_fbg);
    P_dir = calc_psd(sinal_filt_dir, Fs_fbg);
    
    plot(f_axis_fbg(idx_lim_fbg), P_esq(idx_lim_fbg), 'Color', [0 0.5 1], 'LineWidth', 1);
    plot(f_axis_fbg(idx_lim_fbg), P_dir(idx_lim_fbg), 'Color', [0 0.8 0], 'LineWidth', 1);
    plot(f_axis_fbg(idx_lim_fbg), P_fbg(idx_lim_fbg), 'k', 'LineWidth', 2);
    legend('Esquerda', 'Direita', 'Média (Fusão)');
else
    plot(f_axis_fbg(idx_lim_fbg), P_fbg(idx_lim_fbg), 'k', 'LineWidth', 2);
    legend('Sinal FBG');
end
title('Espectro de Frequência Respiratória: Comparação FBG');
xlabel('Frequência (Hz)'); ylabel('Potência (u.a.)');
grid on;

% --- SUBPLOT 2: Comparação Comercial vs Média FBG ---
subplot(2,1,2);
hold on;

% Plotar FBG (Normalizado para facilitar comparação visual de picos)
P_fbg_norm = P_fbg / max(P_fbg(idx_lim_fbg)); % Normalizar pelo pico na banda de interesse
plot(f_axis_fbg(idx_lim_fbg), P_fbg_norm(idx_lim_fbg), 'b', 'LineWidth', 1.5);

if tem_comercial
    N_bio = length(sinal_bio_filt);
    f_axis_bio = get_freq_axis(sinal_bio_filt, Fs_bio);
    P_bio = calc_psd(sinal_bio_filt, Fs_bio);
    
    idx_lim_bio = f_axis_bio <= 2.0;
    
    % Normalizar
    P_bio_norm = P_bio / max(P_bio(idx_lim_bio));
    
    plot(f_axis_bio(idx_lim_bio), P_bio_norm(idx_lim_bio), 'Color', [0.85 0.33 0.1], 'LineWidth', 1.5);
    legend('FBG (Média)', 'Sensor Comercial');
    title('Comparação Espectral: FBG vs Comercial (Normalizado)');
else
    title('Comparação Espectral (Sem dados comerciais)');
end
xlabel('Frequência (Hz)'); ylabel('Potência Normalizada');
grid on;