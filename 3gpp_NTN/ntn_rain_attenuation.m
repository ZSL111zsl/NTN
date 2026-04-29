%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   filename:    ntn_rain_attenuation.m
%   description: 卫星链路雨衰计算（基于 ITU-R P.618-13 简化模型）
%
%   背景（3GPP TR 38.811 §6.6.5 引用 ITU-R P.618）：
%       卫星链路在雨天会产生额外路径衰减，模型包含三步：
%         (1) 比衰减系数  γ_R = k * R^α   (dB/km)     (ITU-R P.838-3)
%         (2) 有效路径长度 L_eff = Ls * r_{0.01}       (km)
%         (3) 超 0.01% 时间雨衰 A_001 = γ_R * L_eff    (dB)
%
%       k / α 与频率、极化、温度有关，本函数采用水平极化的简化表（P.838-3
%       典型查表值）并按对数插值，覆盖 1~40 GHz。
%
%   input:
%       fc_hz           scalar  载波频率 (Hz)
%       elevation_deg   scalar  仰角 (deg)
%       rain_rate_mm_hr scalar  雨率 (mm/h)，典型值：
%                               0=无雨, 5=小雨, 25=中雨, 50=大雨, 100=暴雨
%       h_s_km          scalar  地面站海拔高度 (km)，默认 0
%
%   output:
%       A_rain_dB       scalar  雨衰 (dB)，≥0
%       info            struct  详细中间量（gamma_R, L_eff 等），可选调试
%
%   note:
%       - 对 S 频段 (2 GHz)，即使 50 mm/h 暴雨，A_rain < 0.5 dB；
%         对 Ka 频段 (20 GHz)，50 mm/h 可达 20 dB 以上，
%         这也是 NTN 高频段部署的主要挑战之一。
%       - 本实现用于链路级仿真的链路预算调整，不考虑雨滴尺寸分布细节。
%
%   update note:
%       2026-04-29  created by wangzl  (NR-NTN 毕设仿真雨衰建模)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A_rain_dB, info] = ntn_rain_attenuation(fc_hz, elevation_deg, rain_rate_mm_hr, h_s_km)

if nargin < 4
    h_s_km = 0;      % 地面站默认海平面
end

% 输入保护
rain_rate_mm_hr = max(rain_rate_mm_hr, 0);
elevation_deg   = max(elevation_deg, 5);     % 最小仰角限制避免发散

% ---- 1) 查 ITU-R P.838-3 表，得到比衰减系数 k 和 α ----
% 水平极化（H极化）近似典型值（常用 LEO 卫星S/X/Ku/Ka）
%   频率 GHz   k            α
fc_GHz = fc_hz / 1e9;
freq_table   = [1,     2,     4,     6,     8,     10,    12,    15,    20,    25,    30,    35,    40    ];
k_h_table    = [2.59e-5, 8.47e-5, 6.50e-4, 1.75e-3, 4.54e-3, 1.01e-2, 1.88e-2, 3.67e-2, 7.51e-2, 1.24e-1, 1.87e-1, 2.63e-1, 3.50e-1];
alpha_h_table= [0.969, 1.066, 1.121, 1.308, 1.327, 1.276, 1.217, 1.154, 1.099, 1.061, 1.021, 0.979, 0.939];

% 对 log(k) 做线性插值（k 本身跨几个数量级），α 直接线性插值
if fc_GHz < freq_table(1)
    k_R    = k_h_table(1);
    alpha_R= alpha_h_table(1);
elseif fc_GHz > freq_table(end)
    k_R    = k_h_table(end);
    alpha_R= alpha_h_table(end);
else
    logk_R = interp1(freq_table, log(k_h_table),    fc_GHz, 'linear');
    k_R    = exp(logk_R);
    alpha_R= interp1(freq_table, alpha_h_table,      fc_GHz, 'linear');
end

% ---- 2) 比衰减系数 γ_R ----
if rain_rate_mm_hr <= 0
    gamma_R = 0;
else
    gamma_R = k_R * rain_rate_mm_hr.^alpha_R;   % dB/km
end

% ---- 3) 雨层顶高度 h_R (ITU-R P.839 近似) ----
% 对中纬度地面站典型 h_R ≈ 4 km（0°C 等温线高度），赤道更高、极地更低
h_R_km = 4.0;                                   % 简化取固定值
if h_R_km <= h_s_km
    A_rain_dB = 0;
    info = pack_info(k_R, alpha_R, gamma_R, 0, 0, 0, h_R_km);
    return;
end

% ---- 4) 斜距中雨层长度 Ls ----
elev_rad = elevation_deg * pi / 180;
if elevation_deg >= 5
    Ls_km = (h_R_km - h_s_km) / sin(elev_rad);                     % km
else
    % 小仰角近似（ITU-R P.618 公式 8b）
    R_e   = 8500;   % km, 等效地球半径（含折射）
    Ls_km = 2*(h_R_km - h_s_km) / ...
            ( sqrt( sin(elev_rad).^2 + 2*(h_R_km - h_s_km)/R_e ) + sin(elev_rad) );
end

% ---- 5) 水平投影 LG 与缩减因子 r_{0.01} ----
LG_km = Ls_km * cos(elev_rad);                                  % 水平投影长度
if gamma_R <= 0
    r_001 = 1;
else
    r_001 = 1 ./ ( 1 + 0.78*sqrt(LG_km.*gamma_R./max(fc_GHz,1)) ...
                    - 0.38*(1 - exp(-2*LG_km)) );
    r_001 = max(min(r_001, 1.0), 0.1);       % 限幅
end

% ---- 6) 有效路径长度与总雨衰 ----
L_eff_km = Ls_km * r_001;                   % 有效斜距雨层长度
A_rain_dB= gamma_R * L_eff_km;              % 超 0.01% 时间雨衰（dB）

A_rain_dB = max(A_rain_dB, 0);

info = pack_info(k_R, alpha_R, gamma_R, Ls_km, LG_km, r_001, h_R_km);
info.A_rain_dB = A_rain_dB;
info.L_eff_km  = L_eff_km;

end

%--------------------------------------------------------------------------
function info = pack_info(k_R, alpha_R, gamma_R, Ls_km, LG_km, r_001, h_R_km)
info.k_R       = k_R;
info.alpha_R   = alpha_R;
info.gamma_R   = gamma_R;       % dB/km
info.Ls_km     = Ls_km;
info.LG_km     = LG_km;
info.r_001     = r_001;
info.h_R_km    = h_R_km;
end
