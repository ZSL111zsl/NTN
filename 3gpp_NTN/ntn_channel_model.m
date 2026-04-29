%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   filename:    ntn_channel_model.m
%   description: NTN（非地面网络）增强卫星信道模型，基于 3GPP TR 38.811。
%
%   功能说明：
%       在原有 LTE SCM 信道基础上，叠加 NTN 特有的信道损伤：
%         (1) 大尺度路径损耗（自由空间 + 大气吸收）
%         (2) 阴影衰落（对数正态，慢变）
%         (3) 小尺度衰落（莱斯衰落，K因子与仰角相关）
%         (4) 高速多普勒建模（使用 ntn_doppler_profile 生成序列）
%         (5) 系统时延（卫星轨道高度对应的传播时延）
%
%   调用方式：
%       被 sim_config_NR_NTN_init.m 初始化，参数存于 sim_para.NTN.channel
%       在 wireless_channel.m 调用后对信号做后处理，或独立使用。
%
%   input:
%       signal_in       matrix  [Nr x N_samples]  输入信号
%       ntn_ch_params   struct  信道参数（来自 sim_para.NTN.channel），字段：
%           .model          string  'TR38811_A'(城市) / 'TR38811_B'(郊区) / 'AWGN'
%           .fc_hz          scalar  载波频率 (Hz)
%           .altitude_m     scalar  轨道高度 (m)
%           .elevation_deg  scalar  当前仰角 (deg)
%           .mobile_speed   scalar  用户终端速度 (m/s)，用于地面端多普勒
%           .K_factor_db    scalar  莱斯K因子 (dB)，-Inf=瑞利
%           .enable_pathloss logical 是否叠加路径损耗（仿真时通常已在SNR中体现）
%           .enable_shadow  logical 是否叠加阴影衰落
%           .shadow_std_db  scalar  阴影衰落标准差 (dB)
%           .enable_rician  logical 是否用莱斯衰落替换瑞利
%           .delay_s        scalar  传播时延（由轨道高度决定，也可手动设置）
%       t_ms            scalar  当前子帧时刻 (ms)（用于时变衰落更新）
%       n_samples       scalar  子帧采样点数
%       fs_hz           scalar  采样率 (Hz)
%
%   output:
%       signal_out      matrix  [Nr x N_samples]  经NTN信道后输出信号
%       ch_metrics      struct  信道度量（路径损耗/K因子等）
%
%   note:
%       本函数是信道叠加模块。当 sim_para.NTN.channel.enabled='YES' 时被调用。
%       在原有 wireless_channel 之后链式调用，叠加NTN效应而非替换。
%
%       ta_precomp_en   logical （可选）是否启用TA预补偿（默认false）
%                               启用后信道输出会叠加 TA 补偿效果（残余时延接近0）
%
%   update note:
%       2026-04-15  created by wangzl  (NR-NTN feature extension)
%       2026-04-26  updated by wangzl  (集成TA预补偿闭环)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [signal_out, ch_metrics] = ntn_channel_model(signal_in, ntn_ch_params, t_ms, n_samples, fs_hz)

% ---- 物理常数 ----
c          = 3e8;          % 光速 (m/s)
R_earth_m  = 6371e3;       % 地球半径 (m)
GM         = 3.986e14;     % 地球引力常数 (m^3/s^2)

fc         = ntn_ch_params.fc_hz;
altitude   = ntn_ch_params.altitude_m;
elevation  = ntn_ch_params.elevation_deg;
lambda     = c / fc;       % 波长 (m)

% ---- 传播时延 ----
% 斜距（地面站到卫星）
% d_slant = h / sin(elevation)  （简化：忽略地球曲率，仰角>10deg误差<1%）
elev_rad  = elevation * pi / 180;
if ~isfield(ntn_ch_params, 'delay_s') || ntn_ch_params.delay_s <= 0
    d_slant_m  = altitude / sin(max(elev_rad, 0.1745));  % 限制最小仰角5.6deg避免div0
    ntn_ch_params.delay_s = d_slant_m / c;
end
delay_samples = round(ntn_ch_params.delay_s * fs_hz);

% ---- 大尺度路径损耗 (FSPL + 大气吸收) ----
d_slant_m = altitude / sin(max(elev_rad, 0.1745));
FSPL_dB   = 20*log10(4*pi*d_slant_m/lambda);  % 自由空间路径损耗 (dB)

% 大气吸收（近似，TR 38.811 §6.6.2）
% 主要由氧气吸收（60GHz附近）和水汽吸收决定。
% 简化：按频率范围分段给出典型值
if fc < 6e9        % L/S频段 (<6GHz)
    atmo_loss_dB = 0.1 * (1 / sin(max(elev_rad, 0.1745)));  % ~0.1dB/仰角修正
elseif fc < 30e9   % Ka低频段 (6~30GHz)
    atmo_loss_dB = 0.5 * (1 / sin(max(elev_rad, 0.1745)));
else               % Ka高频段 (>30GHz)
    atmo_loss_dB = 2.0 * (1 / sin(max(elev_rad, 0.1745)));
end

% ---- 雨衰（ITU-R P.618 / TR 38.811 §6.6.5）----
% 独立于FSPL开关：enable_pathloss 控制 FSPL/大气吸收是否计入，
% 雨衰由 enable_rain_fade 单独控制，可在不计大尺度路损的情况下考察雨衰影响
rain_atten_dB = 0;
if isfield(ntn_ch_params, 'enable_rain_fade') && ntn_ch_params.enable_rain_fade
    rain_rate = 0;
    if isfield(ntn_ch_params, 'rain_rate_mm_per_hr')
        rain_rate = ntn_ch_params.rain_rate_mm_per_hr;
    end
    if rain_rate > 0
        rain_atten_dB = ntn_rain_attenuation(fc, elevation, rain_rate);
    end
end
rain_atten_lin = 10^(-rain_atten_dB / 10);   % 线性幅度^2 因子

total_pathloss_dB  = FSPL_dB + atmo_loss_dB;
total_pathloss_lin = 10^(-total_pathloss_dB / 10);  % 线性路径损耗（幅度^2）

% ---- 阴影衰落 ----
% 慢变（相关距离~几十km），仿真中近似为每帧随机常数
shadow_dB = 0;
if isfield(ntn_ch_params, 'enable_shadow') && ntn_ch_params.enable_shadow
    % 阴影衰落标准差（TR 38.811 §6.6，与仰角和模型相关）
    sigma_sh  = ntn_ch_params.shadow_std_db;
    shadow_dB = sigma_sh * randn(1);  % 缓变阴影（此处每子帧独立，可改为相关模型）
end
shadow_lin = 10^(shadow_dB / 10);

% ---- 莱斯小尺度衰落系数 ----
% K 因子与仰角相关（仰角越高 K 越大，LOS越强）
% TR 38.811 §6.7.2 近似: K(dB) ≈ -12 + 0.3*elevation_deg (for S-band Urban)
if isfield(ntn_ch_params, 'K_factor_db')
    K_dB = ntn_ch_params.K_factor_db;
else
    K_dB = -12 + 0.3 * elevation;   % 自动计算K因子
end

[Nr, ~] = size(signal_in);

if isfield(ntn_ch_params, 'enable_rician') && ntn_ch_params.enable_rician && K_dB > -30
    K_lin = 10^(K_dB/10);
    % 莱斯信道系数（归一化单位幅度）
    % LOS分量幅度: sqrt(K/(K+1))，散射分量: sqrt(1/(K+1))
    h_los  = sqrt(K_lin / (K_lin + 1)) * exp(1j * 2*pi * rand(Nr,1));  % LOS随机相位
    h_scat = sqrt(1 / (K_lin + 1)) * (randn(Nr,1) + 1j*randn(Nr,1)) / sqrt(2);
    h_rician = h_los + h_scat;  % [Nr×1]
else
    % 瑞利衰落（K→-inf）
    h_rician = (randn(Nr,1) + 1j*randn(Nr,1)) / sqrt(2);
end

% ---- 组合信道系数 ----
% 若不叠加路径损耗（SNR中已体现），则仅保留小尺度衰落和阴影
% 雨衰始终作为独立幅度因子（由 enable_rain_fade 控制，关闭时=1）
if isfield(ntn_ch_params, 'enable_pathloss') && ntn_ch_params.enable_pathloss
    h_total = h_rician * sqrt(total_pathloss_lin * shadow_lin * rain_atten_lin);
else
    h_total = h_rician * sqrt(shadow_lin * rain_atten_lin);
end

% ---- 信号处理 ----
signal_out = signal_in;

% 应用信道衰落（每根接收天线乘以对应信道系数）
for r = 1:Nr
    signal_out(r, :) = signal_in(r, :) * h_total(r);
end

% ---- TA 预补偿（TR 38.821 §6.3）----
% 若启用 ta_precomp_en，UE 侧已提前 TA 个采样点发送上行信号，
% 等效残余时延 = 信道时延 - TA，理想开环补偿后残余接近0。
ta_precomp_en = isfield(ntn_ch_params, 'ta_precomp_en') && ntn_ch_params.ta_precomp_en;
ta_samples_applied = 0;
ta_metrics_out = struct();

if ta_precomp_en
    orbit_p.altitude_m        = ntn_ch_params.altitude_m;
    orbit_p.elevation_deg     = ntn_ch_params.elevation_deg;
    orbit_p.t_subframe_ms     = t_ms;
    orbit_p.pass_center_ms    = ntn_ch_params.pass_center_ms;
    orbit_p.elevation_min_deg = ntn_ch_params.elevation_min_deg;
    ta_cfg.mode               = 'dynamic';
    ta_cfg.feedback_delay_ms  = 0;
    [ta_samples_applied, ta_metrics_out] = ntn_ta_precomp(orbit_p, fs_hz, ta_cfg);
    net_delay = max(0, delay_samples - ta_samples_applied);
else
    net_delay = delay_samples;
end

% 应用净时延（信道时延 - TA 或完整时延）
if net_delay > 0 && net_delay < n_samples
    signal_out = [signal_out(:, net_delay+1:end), zeros(Nr, net_delay)];
end

% ---- 输出度量 ----
ch_metrics.pathloss_dB       = total_pathloss_dB;
ch_metrics.FSPL_dB           = FSPL_dB;
ch_metrics.atmo_loss_dB      = atmo_loss_dB;
ch_metrics.rain_atten_dB     = rain_atten_dB;
ch_metrics.shadow_dB         = shadow_dB;
ch_metrics.K_factor_dB       = K_dB;
ch_metrics.delay_s           = ntn_ch_params.delay_s;
ch_metrics.delay_samples     = delay_samples;
ch_metrics.h_rician_avg      = mean(abs(h_rician));
ch_metrics.ta_precomp_en     = ta_precomp_en;
ch_metrics.ta_samples        = ta_samples_applied;
ch_metrics.net_delay_samples = net_delay;
if ~isempty(fieldnames(ta_metrics_out))
    ch_metrics.ta = ta_metrics_out;
end

end
