%% Visualização Detalhada FBG: Esquerda vs Direita vs Média vs Referência
% Este script foca-se na comparação da deformação mecânica captada
% em diferentes zonas do tórax, incluindo o sensor de referência (Temperatura)
% para provar o isolamento mecânico.
% INCLUI: Comparação Raw vs Filtered e Seleção Simplificada de Canais.

clc; clear; close all;

%% 1. Carregar Ficheiro FBG
fprintf('--- ANÁLISE DETALHADA DE DEFORMAÇÃO ---\n');
[ficheiro, caminho] = uigetfile('*.txt', 'Selecione o ficheiro FBG (Ex: Mediçao3FBG.txt)');

if isequal(ficheiro, 0), disp('Cancelado.'); return; end
caminho_completo = fullfile(caminho, ficheiro);
fprintf('A carregar: %s ...\n', ficheiro);

% Nome para o gráfico
nome_input = inputdlg('Nome do Paciente/Ensaio (para o título):', 'Input', 1, {'Ensaio 1'});
if isempty(nome_input), nome_ensaio = 'Dados FBG'; else, nome_ensaio = nome_input{1}; end

%% 2. Processamento de Dados
opts = detectImportOptions(caminho_completo);
opts.Delimiter = '\t'; opts.VariableNamingRule = 'preserve';
try, opts.CommentStyle = '#'; catch, end
dados = readtable(caminho_completo, opts);

% Identificar Sensores
col_inicio = 3; 
dados_wl = dados{:, col_inicio:end};
cols_validas = ~all(isnan(dados_wl), 1);
idx_ativos = find(cols_validas);

if length(idx_ativos) < 3
    errordlg('Este script requer 3 sensores ativos (Esq, Temp, Dir) para funcionar corretamente.', 'Erro');
    return;
end

% --- SELEÇÃO DE CANAIS (Alternar Ordem) ---
modo_sel = questdlg('Qual a ordem dos sensores na fibra?', ...
    'Configuração de Canais', 'Padrão (Esq-Temp-Dir)', 'Alternativa (Esq-Dir-Temp)', 'Padrão (Esq-Temp-Dir)');

if strcmp(modo_sel, 'Alternativa (Esq-Dir-Temp)')
    % Ordem Alternativa: O sensor do MEIO é o da DIREITA, o ÚLTIMO é TEMP
    idx_esq = idx_ativos(1);
    idx_dir = idx_ativos(2); % Troca
    idx_temp = idx_ativos(3); % Troca
    fprintf('Configuração: Esq (Col %d), Dir (Col %d), Temp (Col %d)\n', idx_esq, idx_dir, idx_temp);
else
    % Padrão: 1=Esq, 2=Temp, 3=Dir
    idx_esq = idx_ativos(1);
    idx_temp = idx_ativos(2);
    idx_dir = idx_ativos(3);
    fprintf('Configuração Padrão: Esq (Col %d), Temp (Col %d), Dir (Col %d)\n', idx_esq, idx_temp, idx_dir);
end

sinal_esq_raw = dados_wl(:, idx_esq);
sinal_temp_raw = dados_wl(:, idx_temp);
sinal_dir_raw = dados_wl(:, idx_dir);


% --- COMPENSAÇÃO DE TEMPERATURA E FUSÃO ---
% 1. Centralizar (remover média DC) para ver apenas a variação AC
s_esq_c = sinal_esq_raw - mean(sinal_esq_raw);
s_dir_c = sinal_dir_raw - mean(sinal_dir_raw);
s_temp_c = sinal_temp_raw - mean(sinal_temp_raw); % Sinal do meio puramente centrado

% 2. Subtrair o ruído térmico de cada lado
deformacao_esq = s_esq_c - s_temp_c;
deformacao_dir = s_dir_c - s_temp_c;

% 3. Calcular a Média (Sinal Fusão) - ESTE É O SINAL BRUTO (RAW)
deformacao_media_raw = (deformacao_esq + deformacao_dir) / 2;

% --- TEMPO E FILTRAGEM ---
col_t = dados{:, 2};
if iscell(col_t) || isstring(col_t)
    col_t = str2double(strrep(string(col_t), ',', '.'));
end
if any(isnan(col_t)), col_t = (1:height(dados))'; end
duracao = col_t(end) - col_t(1); if duracao <= 0, duracao = height(dados); end
Fs = height(dados) / duracao;
t = (0:height(dados)-1) / Fs;

% Filtro Passa-Banda (0.1 - 2.0 Hz)
f_low = 0.1; f_high = 2.0;
if (Fs/2) <= f_high, f_high = (Fs/2)*0.95; end
Wn = [f_low, f_high] / (Fs/2);
[b, a] = butter(2, Wn, 'bandpass');

% Aplicar filtro aos sinais de interesse
filt_esq = filtfilt(b, a, deformacao_esq);
filt_dir = filtfilt(b, a, deformacao_dir);
filt_media = filtfilt(b, a, deformacao_media_raw); % ESTE É O SINAL FILTRADO
filt_temp = filtfilt(b, a, s_temp_c); 

% Calcular Amplitudes P2P (CORRIGIDO: max - min)
amp_esq = max(filt_esq) - min(filt_esq);
amp_dir = max(filt_dir) - min(filt_dir);
amp_temp = max(filt_temp) - min(filt_temp);
amp_media = max(filt_media) - min(filt_media);

fprintf('Dados processados. A gerar gráficos...\n');

%% 3. GRÁFICO 1: COMPARAÇÃO TEMPORAL (DEFORMAÇÃO)
figure('Color', 'white', 'Position', [50, 50, 1200, 600], 'Name', 'Deformacao Tempo');

% Plotar as linhas
plot(t, deformacao_esq, 'Color', [0 0.4470 0.7410], 'LineWidth', 1.2); hold on; % Azul (Esq)
plot(t, deformacao_dir, 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 1.2); % Verde (Dir)
plot(t, s_temp_c, 'r--', 'LineWidth', 1.0); % Vermelho Tracejado (Referência)
plot(t, deformacao_media_raw, 'k-', 'LineWidth', 2.5); % Preto Espesso (Média)

% Estética
grid on;
title(sprintf('Validação da Deformação Torácica: %s', nome_ensaio), 'FontSize', 14);
ylabel('Amplitude de Deformação (nm)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Tempo (segundos)', 'FontSize', 12, 'FontWeight', 'bold');

% Legenda Inteligente com Valores P2P
leg_esq = sprintf('Lado Esquerdo ');
leg_dir = sprintf('Lado Direito ');
leg_temp = sprintf('Ref. Temperatura ');
leg_media = sprintf('Média Global ');

legend({leg_esq, leg_dir, leg_temp, leg_media}, 'Location', 'best', 'FontSize', 10);

if t(end) > 30, xlim([0, 30]); else, xlim([0, t(end)]); end
axis tight; 

% txt_info = 'Nota: O sinal de referência (vermelho) deve ter amplitude muito inferior aos sinais respiratórios.';
% text(t(1), min([filt_esq; filt_dir])*1.1, txt_info, 'FontSize', 9, 'Color', [0.3 0.3 0.3]);


%% 4. GRÁFICO 2: COMPARAÇÃO ESPECTRAL (FFT)
figure('Color', 'white', 'Position', [100, 100, 1000, 600], 'Name', 'Comparacao Espectral');

calc_psd = @(x) abs(fft(x)).^2 / length(x);
freq_axis = (0:length(filt_media)-1) * (Fs / length(filt_media));

psd_esq = calc_psd(filt_esq);
psd_dir = calc_psd(filt_dir);
psd_media = calc_psd(filt_media);
psd_temp = calc_psd(filt_temp); 

idx_lim = freq_axis <= 1.5;

plot(freq_axis(idx_lim), psd_esq(idx_lim), 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5); hold on;
plot(freq_axis(idx_lim), psd_dir(idx_lim), 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 1.5);
plot(freq_axis(idx_lim), psd_temp(idx_lim), 'r--', 'LineWidth', 1.5); 
plot(freq_axis(idx_lim), psd_media(idx_lim), 'k-', 'LineWidth', 2.5);

grid on;
title(sprintf('Validação Espectral (FFT): %s', nome_ensaio), 'FontSize', 14);
ylabel('Densidade de Potência (u.a.)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Frequência (Hz)', 'FontSize', 12, 'FontWeight', 'bold');
legend({'Esquerda', 'Direita', 'Referência (Temp)', 'Respiração (Média)'}, ...
    'Location', 'northeast', 'FontSize', 11);

[max_p, idx_max] = max(psd_media(idx_lim));
freq_pico = freq_axis(idx_max);
bpm = freq_pico * 60;

text(freq_pico, max_p, sprintf('  %.1f BPM', bpm), ...
    'VerticalAlignment', 'bottom', 'FontWeight', 'bold', 'FontSize', 12);


%% 5. GRÁFICO 3: EFEITO DO FILTRO (RAW VS FILTERED)
figure('Color', 'white', 'Position', [150, 150, 1000, 800], 'Name', 'Efeito Filtro');

% --- A. Comparação no Tempo ---
subplot(2,1,1);
% Centralizar o Raw para comparar visualmente com o Filtered (que já não tem DC)
raw_centered = deformacao_media_raw - mean(deformacao_media_raw);

plot(t, raw_centered, 'Color', [0.7 0.7 0.7], 'LineWidth', 1.0); hold on; % Cinzento
plot(t, filt_media, 'b-', 'LineWidth', 1.5); % Azul

grid on;
title(sprintf('Efeito do Filtro no Domínio do Tempo (%s)', nome_ensaio), 'FontSize', 12);
ylabel('Amplitude (nm)');
legend('Sinal Bruto (Com Ruído)', 'Sinal Filtrado (Limpo)');
if t(end) > 30, xlim([0, 30]); else, xlim([0, t(end)]); end

% --- B. Comparação na Frequência (FFT) ---
subplot(2,1,2);
psd_raw = calc_psd(raw_centered);

% Mostrar uma gama maior (até 5 Hz) para ver o ruído de alta frequência
idx_noise = freq_axis <= 5.0; 

plot(freq_axis(idx_noise), psd_raw(idx_noise), 'Color', [0.7 0.7 0.7], 'LineWidth', 1.0); hold on;
plot(freq_axis(idx_noise), psd_media(idx_noise), 'b-', 'LineWidth', 2.0);

grid on;
title('Efeito do Filtro no Domínio da Frequência', 'FontSize', 12);
xlabel('Frequência (Hz)'); ylabel('Potência');
legend('Espectro Bruto (Ruído Visível)', 'Espectro Filtrado');

fprintf('Gráficos gerados.\n');
