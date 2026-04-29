%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   filename:    nr_chest_dmrs_standalone.m
%   description: NR PDSCH 下行信道估计，基于 DMRS 导频。
%
%   信道估计算法：
%       LS（最小二乘）：H_LS(k) = Y(k) / X(k)，导频位置直接除法
%       MMSE（最小均方误差）：利用信道相关性平滑，假设时不变或慢变
%       线性插值：LS 估计后在频域线性内插到全子载波
%       2D Wiener 滤波：同时利用时频相关性（高精度但计算量大，NTN高速场景）
%
%   NTN 场景特殊处理：
%       - 高多普勒扩展使时间相关性快速衰减 → 需要密集 DMRS（已在 nr_dmrs_config 配置）
%       - 大时延 → 信道在 CP 范围内有效（若时延超过 CP 需截断或额外处理）
%       - 支持多 DMRS 符号联合估计（时域平均后插值）
%
%   input:
%       rx_grid     [N_sc × N_sym] complex  接收到的时隙资源网格（频域）
%       tx_dmrs     [N_dmrs × 1]  complex   发送的 DMRS 序列（来自 nr_dmrs_config）
%       dmrs_cfg    struct                   DMRS 配置（来自 nr_dmrs_config）
%       nr_cfg      struct                   NR 帧结构（来自 nr_numerology_config）
%       chest_params struct                  信道估计参数：
%           .method     string  'LS' | 'LS_INTERP' | 'MMSE' | '2D_WIENER'
%           .noise_var  scalar  噪声方差（MMSE用，可从导频残差估计）
%           .snr_db     scalar  估计SNR（dB，MMSE初始化用）
%           .interp_method string '线性插值' | '三次样条'（LS_INTERP）
%
%   output:
%       H_est       [N_sc × N_sym] complex  估计的信道矩阵（每RE的复增益）
%       noise_est   scalar                   估计的噪声功率
%
%   note:
%       rx_grid 和 H_est 的维度均为 [N_sc × N_sym]，N_sc 为 n_rb×12，
%       N_sym 为时隙内 OFDM 符号数（14 for normal CP）。
%
%   update note:
%       2026-04-15  created by wangzl  (NR-NTN feature extension)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [H_est, noise_est] = nr_chest_dmrs_standalone(rx_grid, tx_dmrs, dmrs_cfg, nr_cfg, chest_params)

% ---- 参数默认值 ----
if ~isfield(chest_params, 'method'),       chest_params.method       = 'LS_INTERP'; end
if ~isfield(chest_params, 'noise_var'),    chest_params.noise_var    = 1e-3;        end
if ~isfield(chest_params, 'interp_method'), chest_params.interp_method = 'linear';  end
if ~isfield(chest_params, 'snr_db'),       chest_params.snr_db       = 20;          end

[N_sc, N_sym] = size(rx_grid);
re_pos  = dmrs_cfg.re_positions + 1;   % 0-based → 1-based
sym_pos = dmrs_cfg.sym_positions + 1;  % 0-based → 1-based

% ---- Step 1：导频位置 LS 估计 ----
% 对每个 DMRS 符号位置进行 LS
H_dmrs_all = [];  % [N_dmrs × N_dmrs_sym] 存储各符号的 LS 估计
sym_pos_list = [];

for l_idx = 1:length(sym_pos)
    l = sym_pos(l_idx);
    if l < 1 || l > N_sym, continue; end

    % 取导频位置接收信号
    y_pilot = rx_grid(re_pos, l);    % [N_dmrs × 1]

    % LS 估计：H_LS = Y / X（导频已知）
    H_ls = y_pilot ./ tx_dmrs;       % [N_dmrs × 1]

    H_dmrs_all   = [H_dmrs_all, H_ls];
    sym_pos_list = [sym_pos_list, l];
end

if isempty(H_dmrs_all)
    warning('nr_chest_dmrs: 未找到有效 DMRS 符号，返回零信道');
    H_est = zeros(N_sc, N_sym);
    noise_est = chest_params.noise_var;
    return;
end

% ---- Step 2：噪声功率估计（导频残差法）----
% 若有多个 DMRS 符号，利用相邻符号差分估计噪声
if size(H_dmrs_all, 2) >= 2
    diff_H = diff(H_dmrs_all, 1, 2);
    noise_est = mean(abs(diff_H(:)).^2) / 2;
    noise_est = max(noise_est, 1e-10);  % 防止为0
else
    noise_est = chest_params.noise_var;
end

% ---- Step 3：频域插值/平滑 → 全子载波信道估计 ----
% 首先对多个 DMRS 符号取平均（时域方向减噪）
H_dmrs_avg = mean(H_dmrs_all, 2);   % [N_dmrs × 1]，多符号平均

switch upper(chest_params.method)

    case 'LS'
        % ---- 纯 LS：仅导频位置有效，其余置零 ----
        H_est = zeros(N_sc, N_sym);
        for l_idx = 1:length(sym_pos_list)
            l = sym_pos_list(l_idx);
            H_est(re_pos, l) = H_dmrs_all(:, l_idx);
        end

    case 'LS_INTERP'
        % ---- LS + 频域插值：从 DMRS 子载波插值到全子载波 ----
        H_full_one = interp_to_full_sc(re_pos, H_dmrs_avg, N_sc, ...
                                        chest_params.interp_method);

        % 时域插值：将单符号估计复制/插值到所有 OFDM 符号
        H_est = zeros(N_sc, N_sym);
        if length(sym_pos_list) == 1
            % 仅一个 DMRS 符号：直接复制到所有符号（静态信道近似）
            for l = 1:N_sym
                H_est(:, l) = H_full_one;
            end
        else
            % 多个 DMRS 符号：各符号分别插值，时域线性内插
            H_sym_full = zeros(N_sc, length(sym_pos_list));
            for l_idx = 1:length(sym_pos_list)
                H_sym_full(:, l_idx) = interp_to_full_sc(re_pos, ...
                    H_dmrs_all(:, l_idx), N_sc, chest_params.interp_method);
            end
            % 时域线性插值到 N_sym 个符号
            H_est = interp_to_full_sym(sym_pos_list, H_sym_full, N_sym);
        end

    case 'MMSE'
        % ---- MMSE 频域平滑（假设宽平稳频率相关信道）----
        % H_mmse = R_hh * (R_hh + sigma^2 * I)^{-1} * H_LS
        % 简化：假设频率相关性为 sinc 函数（均匀多径）
        sigma2 = noise_est;
        snr_lin = 10^(chest_params.snr_db / 10);

        % 频域相关矩阵（简化为对角，等效于频率平滑）
        % 完整 MMSE 需要完整 R_hh（可用 SCM/CDL 模型参数化）
        alpha  = snr_lin / (snr_lin + 1);  % 平滑系数

        % 先做频域插值
        H_mmse_one = interp_to_full_sc(re_pos, H_dmrs_avg, N_sc, 'linear');

        % 频域平滑（移动平均，等效简单 MMSE 滤波）
        win_len = max(3, round(N_sc / (2 * nr_cfg.n_rb_max)));
        H_smooth = movmean(H_mmse_one, win_len);

        % 加权合并：alpha 由 SNR 控制
        H_combined = alpha * H_smooth + (1 - alpha) * H_mmse_one;

        H_est = zeros(N_sc, N_sym);
        for l = 1:N_sym
            H_est(:, l) = H_combined;
        end

    case '2D_WIENER'
        % ---- 2D Wiener 滤波（时频联合，高精度，NTN高速推荐）----
        % 简化实现：先频域 MMSE，再时域 Wiener 时变滤波
        % 完整实现需要 Jakes 模型的 Doppler 谱参数

        sigma2   = noise_est;
        fd_norm  = 0.1;   % 归一化多普勒（fd*Ts_slot），可由 ntn_ch_params.fd 传入

        % 频域 LS 估计（各 DMRS 符号）
        H_sym_full = zeros(N_sc, length(sym_pos_list));
        for l_idx = 1:length(sym_pos_list)
            H_sym_full(:, l_idx) = interp_to_full_sc(re_pos, ...
                H_dmrs_all(:, l_idx), N_sc, 'linear');
        end

        % 时域 Wiener 系数（基于 Jakes 谱，J0(2*pi*fd*n*T_sym)）
        T_sym_s   = nr_cfg.T_slot_s / nr_cfg.sym_per_slot;
        n_past    = min(length(sym_pos_list), 4);

        % 构建相关向量（Jakes）
        r_td = zeros(n_past, 1);
        for n = 1:n_past
            dt_sym = (sym_pos_list(end) - sym_pos_list(max(1,end-n+1)));
            r_td(n) = besselj(0, 2*pi*fd_norm * dt_sym);
        end
        % Wiener 权重（近似解）
        R_auto = toeplitz(r_td);
        w_wien = (R_auto + sigma2*eye(n_past)) \ r_td;

        % 对每个子载波做时域 Wiener 滤波预测（最新1个符号的估计）
        H_est = zeros(N_sc, N_sym);
        if size(H_sym_full, 2) >= n_past
            H_filtered = H_sym_full(:, end-n_past+1:end) * w_wien;
        else
            H_filtered = H_sym_full(:, end);
        end
        for l = 1:N_sym
            H_est(:, l) = H_filtered;
        end

    otherwise
        error('未知信道估计方法: %s', chest_params.method);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   内部函数：频域插值（从 DMRS 子载波插值到全部子载波）
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function H_full = interp_to_full_sc(re_pos, H_dmrs, N_sc, method)
% re_pos: [N_dmrs×1] DMRS 子载波位置（1-based）
% H_dmrs: [N_dmrs×1] DMRS LS 估计
% N_sc:   全子载波数
% method: 'linear' | 'spline' | 'nearest'

x_pilot  = re_pos(:);
H_pilot  = H_dmrs(:);
x_full   = (1:N_sc)';

% 若 DMRS 位置不覆盖边缘，需要外插（使用最近端点值）
x_pilot_real = x_pilot;
H_pilot_real = real(H_pilot);
H_pilot_imag = imag(H_pilot);

if length(x_pilot) < 2
    % 单点：直接填充
    H_full = H_pilot(1) * ones(N_sc, 1);
    return;
end

% 实部插值
Hr_full = interp1(x_pilot_real, H_pilot_real, x_full, method, 'extrap');
Hi_full = interp1(x_pilot_real, H_pilot_imag, x_full, method, 'extrap');
H_full  = Hr_full + 1j * Hi_full;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   内部函数：时域插值（从 DMRS 符号位置插值到全部 OFDM 符号）
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function H_all = interp_to_full_sym(sym_list, H_sym, N_sym)
% sym_list: [1×K] DMRS 符号位置（1-based）
% H_sym:   [N_sc×K] 各 DMRS 符号处的信道估计
% N_sym:   时隙总符号数

[N_sc, K] = size(H_sym);
H_all     = zeros(N_sc, N_sym);
sym_full  = (1:N_sym)';

for k = 1:N_sc
    hr = real(H_sym(k, :))';
    hi = imag(H_sym(k, :))';
    hr_full = interp1(sym_list(:), hr, sym_full, 'linear', 'extrap');
    hi_full = interp1(sym_list(:), hi, sym_full, 'linear', 'extrap');
    H_all(k, :) = (hr_full + 1j*hi_full).';   % 用非共轭转置 .'，否则虚部符号会被反转
end

end
