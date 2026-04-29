%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ntn_dl_sim.m
%   NR NTN 下行链路仿真主程序（基于 3GPP TR 38.811）
%
%   仿真输出：
%     (1) BER vs SNR（不同 Rician K 因子）
%     (2) BLER vs SNR（含 NR LDPC 编码）
%     (3) 多普勒频偏对 BER 的影响
%     (4) 信道估计 NMSE 对比（LS_INTERP / MMSE / 2D_WIENER）
%
%   所有仿真参数集中在 §1 参数配置区，无需修改其他函数。
%
%   依赖文件（需位于同一目录或 MATLAB path 中）：
%     ntn_channel_model.m
%     nr_ldpc_encode.m  nr_ldpc_decode.m  nr_ldpc_bg_tables.m
%     nr_chest_dmrs_standalone.m
%
%   作者: wangzl  日期: 2026-04-28
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;

%% ========================================================================
%  §1  仿真参数配置（全部集中于此）
%  ========================================================================

% ---- NR 帧结构 ----
cfg.mu          = 1;          % Numerology: 0→15kHz, 1→30kHz, 2→60kHz
cfg.n_rb        = 25;         % RB 数（25RB≈10MHz @ 30kHz SCS）
cfg.mod_order   = 4;          % 调制阶数：2=BPSK, 4=QPSK, 16=16QAM, 64=64QAM
cfg.coding      = true;       % true=使用 LDPC 编码；false=无编码（BER验证）
cfg.ldpc_rate   = 1/3;        % LDPC 码率
cfg.ldpc_iter   = 20;         % BP 解码最大迭代次数

% ---- 卫星 / NTN 信道参数 ----
cfg.altitude_m      = 600e3;  % 轨道高度（m），LEO-600
cfg.fc_hz           = 2e9;    % 载波频率（Hz），S 频段
cfg.elevation_deg   = 45;     % 仰角（度）
cfg.mobile_speed    = 3;      % 地面用户速度（m/s）
cfg.enable_shadow   = true;   % 是否叠加对数正态阴影衰落
cfg.shadow_std_db   = 2.0;    % 阴影衰落标准差（dB）

% ---- 仿真实验配置 ----
cfg.snr_dB_vec  = -10:2:20;   % SNR 扫描范围（dB）
cfg.n_frames    = 200;        % 每个 SNR 点仿真帧数（越大结果越准，速度越慢）
cfg.min_errors  = 100;        % 最小误块数（达到后提前停止）

% ---- 实验一：不同 K 因子的 BER（Rician 衰落分析）----
cfg.K_factor_list_db = [-Inf, 0, 7, 15]; % dB；-Inf = 瑞利衰落
cfg.K_labels         = {'Rayleigh', 'K=0dB', 'K=7dB', 'K=15dB'};

% ---- 实验二：不同多普勒频偏的 BER ----
cfg.doppler_list_hz  = [0, 500, 2000, 7000]; % 额外施加的频偏（Hz）
cfg.doppler_labels   = {'0 Hz', '500 Hz', '2 kHz', '7 kHz'};

% ---- 实验三：信道估计方法对比（NMSE）----
cfg.chest_methods    = {'LS_INTERP', 'MMSE', '2D_WIENER'};
cfg.K_for_chest_db   = 7;      % 信道估计对比时使用的固定 K 因子

%% ========================================================================
%  §2  NR 参数推导（由 cfg 自动计算，勿手动修改）
%  ========================================================================
p = nr_numerology(cfg.mu);            % 子载波间隔、CP、符号数等

N_sc    = cfg.n_rb * 12;              % 总子载波数
N_sym   = p.sym_per_slot;            % 每时隙 OFDM 符号数（14）
N_fft   = p.N_fft;                   % FFT 点数
fs_hz   = p.scs_hz * N_fft;          % 采样率（Hz）
cp_vec  = p.cp_lengths;              % 每个符号的 CP 长度（采样点数）[1×14]
N_slot_samples = sum(cp_vec) + N_sym * N_fft; % 每时隙总采样点数

% DMRS 配置：符号2、7、11（Type1 梳状，每2子载波1导频）
dmrs_sym = [3, 8, 12];               % 1-based DMRS 符号位置
dmrs_re  = 1:2:N_sc;                 % 梳状-2 子载波位置（奇数）
N_dmrs   = length(dmrs_re);

% 数据 RE（排除 DMRS 符号的全部 RE 用于信道估计，数据符号的 RE 用于传输）
data_sym = setdiff(1:N_sym, dmrs_sym);
N_data_re_per_slot = length(data_sym) * N_sc; % 数据 RE 总数

% LDPC 每码块信息比特数（取数据 RE 承载比特数的一半，留速率匹配余量）
bits_per_re  = log2(cfg.mod_order);
total_data_bits = N_data_re_per_slot * bits_per_re;
cfg.ldpc_K   = min(256, floor(total_data_bits * cfg.ldpc_rate));

fprintf('NR 参数: mu=%d, SCS=%.0fkHz, N_FFT=%d, N_sc=%d, fs=%.2fMHz\n', ...
    cfg.mu, p.scs_hz/1e3, N_fft, N_sc, fs_hz/1e6);
fprintf('每时隙采样数: %d, 数据RE数: %d, LDPC K=%d bits\n', ...
    N_slot_samples, N_data_re_per_slot, cfg.ldpc_K);

%% ========================================================================
%  §3  DMRS 序列生成（Gold 序列，TS 38.211 §7.4.1）
%  ========================================================================
tx_dmrs = gen_dmrs_seq(N_dmrs, 0);   % 单位功率导频序列 [N_dmrs×1]

%% ========================================================================
%  §4  实验一：BER/BLER vs SNR，不同 K 因子
%  ========================================================================
fprintf('\n=== 实验一：BER/BLER vs SNR（不同 K 因子）===\n');

n_K     = length(cfg.K_factor_list_db);
n_snr   = length(cfg.snr_dB_vec);

ber_K   = zeros(n_K, n_snr);
bler_K  = zeros(n_K, n_snr);

for ki = 1:n_K
    K_dB = cfg.K_factor_list_db(ki);
    ntn_params = make_ntn_params(cfg, K_dB, 0);   % 0 = 无额外频偏
    fprintf('  K=%s dB ...\n', num2str(K_dB));

    for si = 1:n_snr
        snr_dB = cfg.snr_dB_vec(si);
        [ber_K(ki,si), bler_K(ki,si)] = run_link(cfg, p, ntn_params, snr_dB, ...
            tx_dmrs, dmrs_sym, dmrs_re, 'LS_INTERP');
    end
end

%% ========================================================================
%  §5  实验二：多普勒频偏对 BER 的影响
%  ========================================================================
fprintf('\n=== 实验二：多普勒频偏对 BER 的影响 ===\n');

n_dop   = length(cfg.doppler_list_hz);
ber_dop = zeros(n_dop, n_snr);

% 使用固定 K 因子（K=7dB，中等 LOS）
K_fixed_dB = 7;

for di = 1:n_dop
    fd_hz = cfg.doppler_list_hz(di);
    ntn_params = make_ntn_params(cfg, K_fixed_dB, fd_hz);
    fprintf('  fd=%.0f Hz ...\n', fd_hz);

    for si = 1:n_snr
        snr_dB = cfg.snr_dB_vec(si);
        [ber_dop(di,si), ~] = run_link(cfg, p, ntn_params, snr_dB, ...
            tx_dmrs, dmrs_sym, dmrs_re, 'LS_INTERP');
    end
end

%% ========================================================================
%  §6  实验三：信道估计方法 NMSE 对比
%  ========================================================================
fprintf('\n=== 实验三：信道估计 NMSE 对比 ===\n');

n_chest  = length(cfg.chest_methods);
nmse_mat = zeros(n_chest, n_snr);

ntn_params_chest = make_ntn_params(cfg, cfg.K_for_chest_db, 0);

for ci = 1:n_chest
    fprintf('  方法: %s ...\n', cfg.chest_methods{ci});
    for si = 1:n_snr
        snr_dB = cfg.snr_dB_vec(si);
        nmse_mat(ci,si) = run_chest_nmse(cfg, p, ntn_params_chest, snr_dB, ...
            tx_dmrs, dmrs_sym, dmrs_re, cfg.chest_methods{ci});
    end
end

%% ========================================================================
%  §7  绘图与保存
%  ========================================================================
plot_results(cfg, ber_K, bler_K, ber_dop, nmse_mat);
save('ntn_sim_results.mat', 'cfg', 'ber_K', 'bler_K', 'ber_dop', 'nmse_mat');
fprintf('\n仿真完成，结果已保存至 ntn_sim_results.mat\n');

%% ========================================================================
%  §8  内部函数
%  ========================================================================

%--------------------------------------------------------------------------
function [ber, bler] = run_link(cfg, p, ntn_p, snr_dB, tx_dmrs, dsym, dre, chest_method)
% 运行一次完整的 NR NTN 下行链路仿真，返回 BER 和 BLER

N_sc      = cfg.n_rb * 12;
N_sym     = p.sym_per_slot;
N_fft     = p.N_fft;
fs_hz     = p.scs_hz * N_fft;
cp_vec    = p.cp_lengths;
data_sym  = setdiff(1:N_sym, dsym);
N_dmrs    = length(dre);

n_bit_err = 0;
n_blk_err = 0;
n_bit_tot = 0;
n_blk_tot = 0;

for frm = 1:cfg.n_frames
    %-- 发送比特 ---
    if cfg.coding
        info_bits = randi([0 1], cfg.ldpc_K, 1);
        [coded_bits, enc_p] = nr_ldpc_encode(info_bits, cfg.ldpc_rate, 'auto');
        tx_bits   = coded_bits;
    else
        n_tx_bits = length(data_sym) * N_sc * log2(cfg.mod_order);
        tx_bits   = randi([0 1], n_tx_bits, 1);
        info_bits = tx_bits;
        enc_p     = [];
    end

    %-- 调制 ---
    tx_syms = qam_mod(tx_bits, cfg.mod_order);   % 复数符号

    %-- 映射到资源网格 ---
    tx_grid = zeros(N_sc, N_sym);
    % 数据符号
    re_avail = length(data_sym) * N_sc;
    n_fill   = min(length(tx_syms), re_avail);
    data_vec = zeros(re_avail, 1);
    data_vec(1:n_fill) = tx_syms(1:n_fill);
    col_idx = 0;
    for l = data_sym
        tx_grid(:, l) = data_vec(col_idx*N_sc+1 : (col_idx+1)*N_sc);
        col_idx = col_idx + 1;
    end
    % DMRS 符号
    for l = dsym
        tx_grid(dre, l) = tx_dmrs;
        tx_grid(setdiff(1:N_sc, dre), l) = 0;  % 非导频位置置零（保护）
    end

    %-- OFDM 调制 ---
    tx_wave = ofdm_mod(tx_grid, N_fft, cp_vec, N_sc);

    %-- NTN 信道（Rician + 多普勒 + 时延）---
    n_samp = size(tx_wave, 2);
    [rx_wave_ch, ~] = ntn_channel_model(tx_wave, ntn_p, frm, n_samp, fs_hz);

    %-- 多普勒频偏（额外施加，用于实验二）---
    if ntn_p.extra_fd_hz ~= 0
        t_vec = (0:n_samp-1) / fs_hz;
        rx_wave_ch = rx_wave_ch .* exp(1j * 2*pi * ntn_p.extra_fd_hz * t_vec);
    end

    %-- 加 AWGN ---
    snr_lin  = 10^(snr_dB/10);
    sig_pwr  = mean(abs(rx_wave_ch(:)).^2);
    noise_var = sig_pwr / snr_lin;
    noise    = sqrt(noise_var/2) * (randn(size(rx_wave_ch)) + 1j*randn(size(rx_wave_ch)));
    rx_wave  = rx_wave_ch + noise;

    %-- OFDM 解调 ---
    rx_grid = ofdm_demod(rx_wave, N_fft, cp_vec, N_sc, N_sym);

    %-- 信道估计与均衡 ---
    dmrs_cfg_s.re_positions  = dre - 1;    % 0-based
    dmrs_cfg_s.sym_positions = dsym - 1;   % 0-based
    nr_cfg_s.n_rb_max        = cfg.n_rb;
    nr_cfg_s.T_slot_s        = 1 / (p.scs_hz * p.sym_per_slot);
    nr_cfg_s.sym_per_slot    = p.sym_per_slot;
    chest_p.method           = chest_method;
    chest_p.snr_db           = snr_dB;
    chest_p.noise_var        = noise_var / sig_pwr;

    [H_est, ~] = nr_chest_dmrs_standalone(rx_grid, tx_dmrs, dmrs_cfg_s, nr_cfg_s, chest_p);

    % ZF 均衡：每个数据 RE 除以信道估计
    rx_eq = zeros(N_sc, N_sym);
    for l = 1:N_sym
        h_col = H_est(:, l);
        h_col(abs(h_col) < 1e-6) = 1e-6;   % 防止除零
        rx_eq(:, l) = rx_grid(:, l) ./ h_col;
    end

    %-- 提取数据符号并解调 ---
    rx_data = zeros(re_avail, 1);
    col_idx = 0;
    for l = data_sym
        rx_data(col_idx*N_sc+1 : (col_idx+1)*N_sc) = rx_eq(:, l);
        col_idx = col_idx + 1;
    end
    rx_data = rx_data(1:n_fill);

    %-- 软解映射（LLR）---
    llr = qam_demod_soft(rx_data, cfg.mod_order, noise_var / sig_pwr);

    %-- 解码 ---
    if cfg.coding
        [dec_bits, crc_pass, ~] = nr_ldpc_decode(llr, enc_p, cfg.ldpc_iter, 'SMS');
        dec_bits = dec_bits(1:cfg.ldpc_K);
        n_bit_err = n_bit_err + sum(dec_bits ~= info_bits);
        n_bit_tot = n_bit_tot + cfg.ldpc_K;
        n_blk_err = n_blk_err + (~crc_pass);
        n_blk_tot = n_blk_tot + 1;
    else
        hard_bits = double(llr < 0);
        hard_bits = hard_bits(1:length(info_bits));
        n_bit_err = n_bit_err + sum(hard_bits ~= info_bits);
        n_bit_tot = n_bit_tot + length(info_bits);
        n_blk_err = n_blk_err + (sum(hard_bits ~= info_bits) > 0);
        n_blk_tot = n_blk_tot + 1;
    end

    % 提前停止
    if n_blk_err >= cfg.min_errors && frm >= 20
        break;
    end
end

ber  = n_bit_err / max(n_bit_tot, 1);
bler = n_blk_err / max(n_blk_tot, 1);
end


%--------------------------------------------------------------------------
function nmse = run_chest_nmse(cfg, p, ntn_p, snr_dB, tx_dmrs, dsym, dre, method)
% 计算信道估计 NMSE，需要已知真实信道

N_sc    = cfg.n_rb * 12;
N_sym   = p.sym_per_slot;
N_fft   = p.N_fft;
fs_hz   = p.scs_hz * N_fft;
cp_vec  = p.cp_lengths;

nmse_acc = 0;
n_trials = min(cfg.n_frames, 50);   % NMSE 仿真帧数（50帧足够）

for frm = 1:n_trials
    % 发送全1导频时隙（便于提取真实信道）
    tx_grid = zeros(N_sc, N_sym);
    for l = dsym
        tx_grid(dre, l) = tx_dmrs;
    end
    % 数据位置发白噪声（模拟实际发送）
    for l = setdiff(1:N_sym, dsym)
        tx_grid(:, l) = (randn(N_sc,1) + 1j*randn(N_sc,1)) / sqrt(2);
    end

    tx_wave = ofdm_mod(tx_grid, N_fft, cp_vec, N_sc);
    n_samp  = size(tx_wave, 2);

    % 记录信道输入/输出提取真实 H
    rx_wave_ch_clean = ofdm_mod(tx_grid, N_fft, cp_vec, N_sc);  % 无噪参考
    [rx_wave_ch, ~]  = ntn_channel_model(tx_wave, ntn_p, frm, n_samp, fs_hz);

    % AWGN
    snr_lin   = 10^(snr_dB/10);
    sig_pwr   = mean(abs(rx_wave_ch(:)).^2);
    noise_var = sig_pwr / snr_lin;
    rx_wave   = rx_wave_ch + sqrt(noise_var/2) * ...
                (randn(size(rx_wave_ch)) + 1j*randn(size(rx_wave_ch)));

    rx_grid     = ofdm_demod(rx_wave,    N_fft, cp_vec, N_sc, N_sym);
    rx_grid_ref = ofdm_demod(rx_wave_ch, N_fft, cp_vec, N_sc, N_sym);

    % 真实信道：导频位置的信道系数（LS 无噪估计，作为 ground truth）
    H_true = zeros(N_sc, N_sym);
    for l = dsym
        H_true(dre, l) = rx_grid_ref(dre, l) ./ tx_dmrs;
    end
    % 用线性插值填充全频带（ground truth）
    for l = dsym
        x_p = dre(:);
        h_p = H_true(dre, l);
        H_true(:, l) = interp1(x_p, h_p, (1:N_sc)', 'linear', 'extrap');
    end
    % 时域插值 ground truth
    H_true_full = zeros(N_sc, N_sym);
    for sc = 1:N_sc
        h_row = real(H_true(sc, dsym)).' + 1j * imag(H_true(sc, dsym)).';
        if length(dsym) >= 2
            H_true_full(sc,:) = interp1(dsym(:), h_row, (1:N_sym)', 'linear', 'extrap');
        else
            H_true_full(sc,:) = h_row(1);
        end
    end

    % 待评估的信道估计
    dmrs_cfg_s.re_positions  = dre - 1;
    dmrs_cfg_s.sym_positions = dsym - 1;
    nr_cfg_s.n_rb_max        = cfg.n_rb;
    nr_cfg_s.T_slot_s        = 1 / (p.scs_hz * p.sym_per_slot);
    nr_cfg_s.sym_per_slot    = p.sym_per_slot;
    chest_p.method           = method;
    chest_p.snr_db           = snr_dB;
    chest_p.noise_var        = noise_var / sig_pwr;

    [H_est, ~] = nr_chest_dmrs_standalone(rx_grid, tx_dmrs, dmrs_cfg_s, nr_cfg_s, chest_p);

    % NMSE = ||H_est - H_true||^2 / ||H_true||^2
    err    = H_est - H_true_full;
    nmse_acc = nmse_acc + sum(abs(err(:)).^2) / max(sum(abs(H_true_full(:)).^2), 1e-12);
end

nmse = nmse_acc / n_trials;
end


%--------------------------------------------------------------------------
function ntn_p = make_ntn_params(cfg, K_dB, extra_fd_hz)
% 构造 ntn_channel_model 所需的参数结构体

ntn_p.fc_hz          = cfg.fc_hz;
ntn_p.altitude_m     = cfg.altitude_m;
ntn_p.elevation_deg  = cfg.elevation_deg;
ntn_p.mobile_speed   = cfg.mobile_speed;
ntn_p.K_factor_db    = K_dB;
ntn_p.enable_pathloss = false;  % 路径损耗已包含在 SNR 中
ntn_p.enable_shadow  = cfg.enable_shadow;
ntn_p.shadow_std_db  = cfg.shadow_std_db;
ntn_p.enable_rician  = true;
ntn_p.delay_s        = 0;       % 0 = 由轨道高度自动计算
ntn_p.ta_precomp_en  = false;   % 独立链路仿真不做 TA 预补偿（时延已被归一化）
ntn_p.extra_fd_hz    = extra_fd_hz;
end


%--------------------------------------------------------------------------
function p = nr_numerology(mu)
% 根据 Numerology 参数 mu 计算 NR 帧结构参数

scs_hz_list = [15e3, 30e3, 60e3, 120e3, 240e3];
p.mu          = mu;
p.scs_hz      = scs_hz_list(mu + 1);
p.sym_per_slot = 14;   % normal CP

% FFT 点数（依据带宽和 SCS 选择，此处按 10MHz 典型值）
switch mu
    case 0;  p.N_fft = 1024;
    case 1;  p.N_fft = 512;
    case 2;  p.N_fft = 256;
    case 3;  p.N_fft = 128;
    otherwise; p.N_fft = 512;
end

% CP 长度（TS 38.211 Table 5.3.1-1，normal CP）
% 第1个符号 CP 较长，其余相同
Nfft = p.N_fft;
cp_normal = round(Nfft * 144/2048);   % 标准比例
cp_first  = round(Nfft * 160/2048);   % 第0个符号（每半帧首符号）
p.cp_lengths = [cp_first, repmat(cp_normal, 1, 13)];

p.T_slot_s = 1e-3 / (2^mu);   % 每时隙时长（s）
end


%--------------------------------------------------------------------------
function tx_wave = ofdm_mod(grid, N_fft, cp_vec, N_sc)
% 自包含 OFDM 调制（单天线）
% grid: [N_sc × N_sym] 频域资源网格
% 输出: [1 × N_samples] 时域信号

[~, N_sym] = size(grid);
dc = N_fft/2 + 1;
n_lower = floor(N_sc/2);
n_upper = N_sc - n_lower;

N_samp = sum(cp_vec) + N_sym * N_fft;
tx_wave = zeros(1, N_samp);
ptr = 1;

for l = 1:N_sym
    freq_buf = zeros(N_fft, 1);
    freq_buf(dc-n_lower : dc-1)   = grid(1:n_lower, l);
    freq_buf(dc+1       : dc+n_upper) = grid(n_lower+1:end, l);
    td = ifft(ifftshift(freq_buf), N_fft) * sqrt(N_fft);
    cp = cp_vec(l);
    sym = [td(end-cp+1:end); td];
    tx_wave(1, ptr:ptr+length(sym)-1) = sym.';
    ptr = ptr + length(sym);
end
end


%--------------------------------------------------------------------------
function rx_grid = ofdm_demod(rx_wave, N_fft, cp_vec, N_sc, N_sym)
% 自包含 OFDM 解调（单天线）
% rx_wave: [1 × N_samples]
% 输出: [N_sc × N_sym] 频域资源网格

dc = N_fft/2 + 1;
n_lower = floor(N_sc/2);
n_upper = N_sc - n_lower;
rx_grid = zeros(N_sc, N_sym);
ptr = 1;

for l = 1:N_sym
    cp = cp_vec(l);
    sym_sig = rx_wave(1, ptr : ptr + cp + N_fft - 1);
    ptr = ptr + cp + N_fft;
    td  = sym_sig(cp+1 : end).';   % 去 CP
    fd  = fftshift(fft(td, N_fft)) / sqrt(N_fft);
    lower = fd(dc-n_lower : dc-1);
    upper = fd(dc+1       : dc+n_upper);
    rx_grid(:, l) = [lower; upper];
end
end


%--------------------------------------------------------------------------
function syms = qam_mod(bits, M)
% QAM 调制（自然映射）
% bits: [N_bits×1]  M: 调制阶数（2/4/16/64）

bits = bits(:);
k    = log2(M);
n_sym = floor(length(bits) / k);
bits  = bits(1 : n_sym*k);
% 二进制转十进制（不依赖 bi2de）
bit_mat = reshape(bits, k, n_sym)';   % [n_sym × k]
pw      = 2.^(k-1:-1:0);
idx     = bit_mat * pw(:);            % 自然映射索引 [n_sym×1]，值域 0~M-1

switch M
    case 2   % BPSK
        syms = 1 - 2*bits(1:k:end);
        syms = complex(syms, 0);
    case 4   % QPSK
        const = (1/sqrt(2)) * [1+1j, -1+1j, 1-1j, -1-1j];
        syms  = const(idx + 1).';
    case 16  % 16QAM
        const16 = qammod_const(16);
        syms = const16(idx + 1).';
    case 64  % 64QAM
        const64 = qammod_const(64);
        syms = const64(idx + 1).';
    otherwise
        error('不支持的调制阶数: %d', M);
end
syms = syms(:);
end


%--------------------------------------------------------------------------
function llr = qam_demod_soft(rx_syms, M, noise_var)
% 近似软解映射（独立 bit 的 LLR）
% rx_syms: [N_sym×1]  M: 调制阶数  noise_var: 噪声方差

rx_syms = rx_syms(:);
k       = log2(M);
N_sym   = length(rx_syms);
llr     = zeros(N_sym * k, 1);
sigma2  = max(noise_var, 1e-10);

switch M
    case 2   % BPSK
        llr = 2 * real(rx_syms) / sigma2;
    case 4   % QPSK（独立 I/Q 解映射）
        llr(1:2:end) = 2 * sqrt(2) * real(rx_syms) / sigma2;
        llr(2:2:end) = 2 * sqrt(2) * imag(rx_syms) / sigma2;
    otherwise
        % 近似 Max-Log-MAP：对每个 bit，找最近的0类和1类星座点
        const = qammod_const(M);
        const_norm = const / sqrt(mean(abs(const).^2));
        for si = 1:N_sym
            r = rx_syms(si);
            for b = 1:k
                % 构建 bit b 为0/1 的星座子集
                all_idx = 0:M-1;
                bit_b   = floor(all_idx / 2^(k-b));
                mask0   = mod(bit_b, 2) == 0;
                mask1   = ~mask0;
                d0 = min(abs(r - const_norm(mask0)).^2);
                d1 = min(abs(r - const_norm(mask1)).^2);
                llr((si-1)*k + b) = (d1 - d0) / sigma2;
            end
        end
end
end


%--------------------------------------------------------------------------
function const = qammod_const(M)
% 生成归一化 QAM 星座点（自然映射顺序）
sqM   = sqrt(M);
pts   = (-(sqM-1) : 2 : (sqM-1));
[re_g, im_g] = meshgrid(pts, fliplr(pts));
const = (re_g(:) + 1j*im_g(:)).' ;
const = const / sqrt(mean(abs(const).^2));  % 归一化单位平均功率
end


%--------------------------------------------------------------------------
function seq = gen_dmrs_seq(N, n_id)
% 生成长度 N 的 DMRS 序列（Gold 序列，TS 38.211 §5.2.1 简化实现）
% 不依赖 Communications Toolbox
seq = gold_seq(N, n_id, 0) + 1j * gold_seq(N, n_id, 1);
seq = seq / sqrt(mean(abs(seq).^2));   % 归一化单位功率
end

function c = gold_seq(N, n_id, offset)
% 生成 Gold 序列的 I 或 Q 分量
c_init = mod(n_id + offset, 2^31);
% 31位移位寄存器初始化（内置转换，不用 de2bi）
x1 = bitget(c_init + 2^31, 31:-1:1);   % MSB first
x2 = bitget(1 + offset,    31:-1:1);

c = zeros(N, 1);
for n = 1:N
    new_x1 = mod(x1(3) + x1(end), 2);
    new_x2 = mod(x2(3) + x2(1) + x2(end-1) + x2(end), 2);
    c(n)   = mod(x1(end) + x2(end), 2);
    x1 = [x1(2:end), new_x1];
    x2 = [x2(2:end), new_x2];
end
c = (1 - 2*c) / sqrt(2);   % BPSK 映射 ±1/√2
end


%--------------------------------------------------------------------------
function plot_results(cfg, ber_K, bler_K, ber_dop, nmse_mat)

colors   = lines(max(length(cfg.K_factor_list_db), length(cfg.doppler_list_hz)));
markers  = {'o-','s-','d-','^-','v-'};
snr      = cfg.snr_dB_vec;

figure('Name','NR NTN 下行链路仿真结果','NumberTitle','off','Position',[100 100 1300 900]);

% --- 子图1：BER vs SNR（不同K因子）---
subplot(2,2,1);
for ki = 1:length(cfg.K_factor_list_db)
    semilogy(snr, max(ber_K(ki,:), 1e-5), markers{mod(ki-1,5)+1}, ...
        'Color', colors(ki,:), 'LineWidth', 1.5, 'MarkerSize', 6); hold on;
end
grid on; ylim([1e-5 1]);
xlabel('SNR (dB)'); ylabel('BER');
title('BER vs SNR（不同 Rician K 因子）');
legend(cfg.K_labels, 'Location', 'southwest');

% --- 子图2：BLER vs SNR（不同K因子）---
subplot(2,2,2);
for ki = 1:length(cfg.K_factor_list_db)
    semilogy(snr, max(bler_K(ki,:), 1e-3), markers{mod(ki-1,5)+1}, ...
        'Color', colors(ki,:), 'LineWidth', 1.5, 'MarkerSize', 6); hold on;
end
grid on; ylim([1e-3 1]);
xlabel('SNR (dB)'); ylabel('BLER');
title('BLER vs SNR（LDPC 编码，不同 K 因子）');
legend(cfg.K_labels, 'Location', 'southwest');

% --- 子图3：BER vs SNR（不同多普勒频偏）---
subplot(2,2,3);
for di = 1:length(cfg.doppler_list_hz)
    semilogy(snr, max(ber_dop(di,:), 1e-5), markers{mod(di-1,5)+1}, ...
        'Color', colors(di,:), 'LineWidth', 1.5, 'MarkerSize', 6); hold on;
end
grid on; ylim([1e-5 1]);
xlabel('SNR (dB)'); ylabel('BER');
title('BER vs SNR（不同多普勒频偏，K=7dB）');
legend(cfg.doppler_labels, 'Location', 'southwest');

% --- 子图4：信道估计 NMSE vs SNR ---
subplot(2,2,4);
for ci = 1:length(cfg.chest_methods)
    semilogy(snr, max(nmse_mat(ci,:), 1e-4), markers{mod(ci-1,5)+1}, ...
        'Color', colors(ci,:), 'LineWidth', 1.5, 'MarkerSize', 6); hold on;
end
grid on;
xlabel('SNR (dB)'); ylabel('NMSE');
title('信道估计 NMSE 对比（K=7dB）');
legend(cfg.chest_methods, 'Location', 'northeast');

sgtitle('NR NTN 下行链路仿真（3GPP TR 38.811, LEO-600, S频段）');

saveas(gcf, 'ntn_sim_results.png');
fprintf('图片已保存至 ntn_sim_results.png\n');
end
