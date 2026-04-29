%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   filename:    ntn_ta_precomp.m
%   description: NTN UE 侧时间提前量（Timing Advance, TA）预补偿模块。
%
%   背景（TR 38.821 §6.3）：
%       LEO 卫星距地面约 600 km，单向传播时延约 2 ms（垂直），
%       最大斜距时延可达 ~10 ms（仰角10°）。
%       3GPP NR LTE 标准 TA 范围最大约 0.68 ms（NTA_max=3846，
%       对应 30 kHz SCS 为 0.68 ms），远不足以覆盖 NTN 场景。
%
%       TR 38.821 §6.3 解决方案：
%         UE 在发送上行信号前，提前 TA_precomp = 2 * d_slant / c 时长，
%         使上行信号到达卫星（作为基站）时对齐下行接收定时。
%         即 UE 在时域提前发送，补偿往返时延（RTD = 2 * d_slant / c）。
%
%       本模块实现：
%         (1) 根据轨道参数计算当前时刻斜距 d_slant
%         (2) 计算 TA_precomp（采样点数）
%         (3) 在信道模型中应用：信道注入时延后，UE 发射预先提前
%             等效于接收端残余时延 = 信道时延 - TA预补偿 ≈ 0
%         (4) 输出补偿效果度量
%
%   调用位置：
%       ntn_channel_model.m 中，在应用信道时延后调用本函数计算 TA，
%       或在仿真主循环中作为 UE 上行信号预处理步骤。
%
%   input:
%       orbit_params    struct  轨道参数（来自 sim_para.NTN.orbit）：
%           .altitude_m         轨道高度 (m)
%           .elevation_deg      当前仰角 (deg)（可为标量或时变序列）
%           .t_subframe_ms      当前子帧时刻 (ms)（用于时变仰角）
%           .pass_center_ms     过顶时刻 (ms)（用于从轨道模型推算仰角）
%           .elevation_min_deg  最小仰角 (deg)
%       fs_hz           scalar  采样率 (Hz)
%       ta_params       struct  TA 配置参数：
%           .mode       string  'static'（固定仰角）| 'dynamic'（时变轨道）
%           .feedback_delay_ms scalar 基站 TA 命令延迟（ms，典型4~8ms），
%                               用于评估开环误差
%
%   output:
%       ta_samples      scalar  TA 预补偿量（采样点数，整数）
%       ta_metrics      struct  TA 度量：
%           .d_slant_m          斜距 (m)
%           .delay_s            单向时延 (s)
%           .rtd_s              往返时延 RTD (s)
%           .ta_samples         TA 预补偿采样点数
%           .ta_ms              TA 预补偿时长 (ms)
%           .residual_samples   残余时延采样点数（信道时延 - TA）
%
%   update note:
%       2026-04-26  created by wangzl  (NTN TA预补偿闭环实现)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ta_samples, ta_metrics] = ntn_ta_precomp(orbit_params, fs_hz, ta_params)

if nargin < 3
    ta_params = struct('mode', 'dynamic', 'feedback_delay_ms', 6);
end

c       = 3e8;
GM      = 3.986004418e14;
R_earth = 6371e3;

altitude  = orbit_params.altitude_m;
r_orbit   = R_earth + altitude;
v_sat     = sqrt(GM / r_orbit);     % 轨道速度 (m/s)

% ---- 计算当前仰角 ----
if strcmpi(ta_params.mode, 'static') || ~isfield(orbit_params, 't_subframe_ms')
    % 静态模式：直接使用配置仰角
    if isfield(orbit_params, 'elevation_deg')
        elev_deg = orbit_params.elevation_deg;
    else
        elev_deg = 45;
    end
else
    % 动态模式：根据轨道模型推算当前仰角
    t_ms      = orbit_params.t_subframe_ms;
    t_center  = orbit_params.pass_center_ms;
    elev_min  = orbit_params.elevation_min_deg;
    elev_min_rad = elev_min * pi / 180;

    % 半弧时长
    alpha_max   = acos((R_earth / r_orbit) * cos(elev_min_rad)) - elev_min_rad;
    T_half_ms   = alpha_max * r_orbit / v_sat * 1000;

    % 当前弧角
    dt       = t_ms - t_center;
    alpha_t  = (dt / T_half_ms) * alpha_max;
    alpha_t  = max(-alpha_max, min(alpha_max, alpha_t));

    % 仰角（地心角到仰角的近似转换）
    % elev = 90 - alpha_t（近似，仰角在过顶时为90deg）
    % 更准确: sin(elev) = (r_orbit/R_earth)*cos(alpha_t) - 1 (简化)
    % 使用标准几何: cos(alpha_t) = (R_earth^2 + r_orbit^2 - d^2) / (2*R_earth*r_orbit)
    sin_elev = (r_orbit * cos(alpha_t) - R_earth) / ...
               sqrt(r_orbit^2 + R_earth^2 - 2*r_orbit*R_earth*cos(alpha_t));
    elev_deg = asin(max(-1, min(1, sin_elev))) * 180 / pi;
    elev_deg = max(orbit_params.elevation_min_deg, elev_deg);
end

% ---- 计算斜距 ----
elev_rad  = elev_deg * pi / 180;
% 精确斜距（地球椭球近似为球形）
% d^2 = R_earth^2 + r_orbit^2 - 2*R_earth*r_orbit*cos(alpha)
% 简化用仰角直接计算：d = R_earth*(sin(elev+acos(R_earth/r_orbit*cos(elev))) - sin(elev))
% 最常用近似（仰角>10°时误差<1%）：d ≈ altitude / sin(elev)
elev_rad_safe = max(elev_rad, 5 * pi/180);   % 防除零：最小5度
d_slant_m     = sqrt(R_earth^2 * cos(elev_rad_safe)^2 + r_orbit^2 - R_earth^2) ...
                - R_earth * cos(elev_rad_safe);
% 简化公式校验
d_slant_approx = altitude / sin(elev_rad_safe);

% 取精确值（几何解）
d_slant_m = max(d_slant_m, altitude);  % 至少等于轨道高度

% ---- 计算时延与TA ----
delay_s   = d_slant_m / c;             % 单向传播时延
rtd_s     = 2 * delay_s;               % 往返时延（Round-Trip Delay）

% TA 预补偿量 = RTD（UE 提前 RTD 时间发送，使上行信号及时到达）
% 注：TS 38.213 §4.2 规定 TA 调整步长 Tc=16*64/480000 s（NR最小时间单元）
% 此处在采样域直接取整
ta_samples     = round(rtd_s * fs_hz);  % TA 预补偿采样点数

% ---- 信道实际时延（用于计算残余误差）----
channel_delay_samples = round(delay_s * fs_hz);

% 残余时延 = 信道时延 - TA/2（TA补偿的是往返，实际信道是单向）
% 发射端提前 ta_samples 个采样 ≡ 接收端看到的到达时间提前了 ta_samples
% 若信道单向时延 = channel_delay_samples，则接收残余 = channel_delay - ta/2
% 但在链路仿真中通常 TA 补偿单向时延，残余 = channel_delay - ta_samples
% TR 38.821 约定：TA_precomp = 2*d/c（即RTD），效果是残余归零（开环理想）
residual_samples = channel_delay_samples - ta_samples;

% ---- 输出度量 ----
ta_metrics.elev_deg          = elev_deg;
ta_metrics.d_slant_m         = d_slant_m;
ta_metrics.delay_s           = delay_s;
ta_metrics.rtd_s             = rtd_s;
ta_metrics.ta_samples        = ta_samples;
ta_metrics.ta_ms             = ta_samples / fs_hz * 1000;
ta_metrics.channel_delay_samples = channel_delay_samples;
ta_metrics.residual_samples  = residual_samples;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ntn_apply_ta_precomp: 在时域信号上应用 TA 预补偿（UE 侧上行发射前处理）
%
%   原理：
%       UE 在 t0 时刻原本应该发送的信号，提前 ta_samples 个采样点发送。
%       等效于：将发射信号在时域做循环左移（提前），
%       使其经过单向传播时延后，恰好在预期时刻到达卫星。
%
%   input:
%       tx_signal   [N_samples×1 或 N_ant×N_samples]  发射信号
%       ta_samples  scalar  TA 预补偿采样数（整数，≥0）
%
%   output:
%       tx_signal_ta  同尺寸，时域左移 ta_samples 后的信号
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tx_signal_ta = ntn_apply_ta_precomp(tx_signal, ta_samples)

ta_int = max(0, round(ta_samples));

if ta_int == 0
    tx_signal_ta = tx_signal;
    return;
end

sz = size(tx_signal);
if sz(1) == 1 || sz(2) == 1
    % 列/行向量：直接循环左移
    tx_signal_ta = circshift(tx_signal(:), -ta_int);
    tx_signal_ta = reshape(tx_signal_ta, sz);
else
    % 多天线：每行（每根天线）独立循环左移
    tx_signal_ta = zeros(sz);
    for ant = 1:sz(1)
        row = tx_signal(ant, :);
        tx_signal_ta(ant, :) = circshift(row, -ta_int);
    end
end

end
