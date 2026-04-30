%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ntn_sim_main.m
%   NR NTN 下行链路仿真（本科毕业设计中期答辩专用主函数）
%
%   一、设计目标
%       参考 3GPP TR 38.811（NR-NTN 信道模型）和 TR 38.821（NR-NTN 解决方案），
%       搭建基于 LEO-600 卫星的 NR 下行单链路仿真，完成 6 项核心分析：
%         M1 ideal        理想 AWGN 下 BER（仿真 vs 理论，用于校准基线）
%         M2 kfactor      不同 Rician K 因子下的 BER / BLER
%         M3 doppler      LEO 多普勒频偏对 BER 的影响
%         M4 delay        残余传播时延对 BER 的影响，并对比 TA 预补偿
%         M5 pathloss     不同仰角下的 FSPL + 大气吸收 + 雨衰
%         M6 rain         S 频段 vs Ka 频段雨衰对 BER 的影响
%         M7 compensation GNSS 频偏补偿 & 多种信道估计算法补偿效果
%
%   二、使用方法（一键运行 / 按模块运行）
%       >> ntn_sim_main                                % 按 cfg.modules 运行
%       >> cfg.modules={'all'};     ntn_sim_main       % 全部 7 个模块
%       >> cfg.modules={'ideal'};   ntn_sim_main       % 仅跑校准基线
%       >> cfg.modules={'rain','compensation'}; ntn_sim_main
%
%   三、参数修改
%       打开本文件 §1 参数区即可，所有可调项集中于此；
%       其它函数一律不需要修改。
%
%   四、输出
%       ./results/*.png   每个模块对应 1~2 张仿真图（答辩可直接使用）
%       ./results/ntn_sim_main_results.mat  所有数据矩阵
%
%   五、依赖（必须位于同一目录）
%       ntn_channel_model.m, ntn_rain_attenuation.m, ntn_ta_precomp.m
%       nr_ldpc_encode.m, nr_ldpc_decode.m, nr_ldpc_bg_tables.m
%       nr_chest_dmrs_standalone.m
%
%   作者: wangzl  日期: 2026-04-29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ntn_sim_main(cfg_override)

%% ========================================================================
%  §1  参数配置区（按需修改）
%  ========================================================================
cfg = local_default_cfg();
if nargin >= 1 && isstruct(cfg_override)
    % 允许外部传入部分字段覆盖默认值
    fn = fieldnames(cfg_override);
    for k = 1:numel(fn)
        cfg.(fn{k}) = cfg_override.(fn{k});
    end
end

% 兼容从脚本调用：若 base 工作区已存在 cfg 结构体，则以其字段覆盖默认值
if nargin < 1
    try
        base_cfg = evalin('base', 'cfg');
        if isstruct(base_cfg)
            fn = fieldnames(base_cfg);
            for k = 1:numel(fn)
                cfg.(fn{k}) = base_cfg.(fn{k});
            end
        end
    catch
    end
end

%% ========================================================================
%  §2  模块展开与 NR 参数推导
%  ========================================================================
all_modules = {'ideal','kfactor','doppler','delay','pathloss','rain','compensation'};
if any(strcmpi(cfg.modules, 'all'))
    cfg.modules = all_modules;
end

% 固定全局随机种子：保证每次仿真结果完全可复现
if isfield(cfg, 'rng_seed') && ~isempty(cfg.rng_seed)
    rng(cfg.rng_seed, 'twister');
end

[p, cfg] = local_expand_nr_params(cfg);
fprintf('\n============================================================\n');
fprintf('  NR-NTN 下行链路仿真（TR 38.811）\n');
if isfield(cfg,'fast_verify') && cfg.fast_verify
    fprintf('  [模式] FAST_VERIFY —— 3~5 分钟出图，用于验证代码方向性\n');
    fprintf('         正式仿真请在调用前设： cfg.fast_verify = false\n');
else
    fprintf('  [模式] FULL 正式仿真 —— 预计 30~60 分钟\n');
end
fprintf('  载波: %.2f GHz | SCS: %.0f kHz | N_FFT: %d | N_RB: %d\n', ...
        cfg.fc_hz/1e9, p.scs_hz/1e3, p.N_fft, cfg.n_rb);
fprintf('  SNR: %.1f..%.1f dB (%d 点) | 每点最大 %d 帧\n', ...
        cfg.snr_dB_vec(1), cfg.snr_dB_vec(end), ...
        length(cfg.snr_dB_vec), cfg.n_frames);
fprintf('  启用模块: %s\n', strjoin(cfg.modules, ', '));
if isfield(cfg, 'rng_seed') && ~isempty(cfg.rng_seed)
    fprintf('  随机种子: %d (结果完全可复现)\n', cfg.rng_seed);
end
fprintf('============================================================\n\n');

results_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
if exist(results_dir, 'dir') == 0
    mkdir(results_dir);
end

R = struct();   % 结果容器

%% ========================================================================
%  §3  逐模块执行
%  ========================================================================
t0 = tic;

% 每个模块跑完立即保存，防止长时间仿真中途崩溃丢数据
mat_path = fullfile(results_dir, 'ntn_sim_main_results.mat');
n_total_mods = sum(ismember({'ideal','kfactor','doppler','delay','pathloss','rain','compensation'}, cfg.modules));
mod_done = 0;

if ismember('ideal',        cfg.modules), tm=tic; R.ideal        = exp_ideal(cfg, p, results_dir);        mod_done=mod_done+1; log_mod_done('M1 ideal',        tm, mod_done, n_total_mods, t0); save_checkpoint(mat_path, R, cfg, p, t0); end
if ismember('kfactor',      cfg.modules), tm=tic; R.kfactor      = exp_kfactor(cfg, p, results_dir);      mod_done=mod_done+1; log_mod_done('M2 kfactor',      tm, mod_done, n_total_mods, t0); save_checkpoint(mat_path, R, cfg, p, t0); end
if ismember('doppler',      cfg.modules), tm=tic; R.doppler      = exp_doppler(cfg, p, results_dir);      mod_done=mod_done+1; log_mod_done('M3 doppler',      tm, mod_done, n_total_mods, t0); save_checkpoint(mat_path, R, cfg, p, t0); end
if ismember('delay',        cfg.modules), tm=tic; R.delay        = exp_delay(cfg, p, results_dir);        mod_done=mod_done+1; log_mod_done('M4 delay',        tm, mod_done, n_total_mods, t0); save_checkpoint(mat_path, R, cfg, p, t0); end
if ismember('pathloss',     cfg.modules), tm=tic; R.pathloss     = exp_pathloss(cfg, p, results_dir);     mod_done=mod_done+1; log_mod_done('M5 pathloss',     tm, mod_done, n_total_mods, t0); save_checkpoint(mat_path, R, cfg, p, t0); end
if ismember('rain',         cfg.modules), tm=tic; R.rain         = exp_rain(cfg, p, results_dir);         mod_done=mod_done+1; log_mod_done('M6 rain',         tm, mod_done, n_total_mods, t0); save_checkpoint(mat_path, R, cfg, p, t0); end
if ismember('compensation', cfg.modules), tm=tic; R.compensation = exp_compensation(cfg, p, results_dir); mod_done=mod_done+1; log_mod_done('M7 compensation', tm, mod_done, n_total_mods, t0); save_checkpoint(mat_path, R, cfg, p, t0); end

R.cfg      = cfg;
R.p_nr     = p;
R.runtime  = toc(t0);

save(mat_path, '-struct', 'R');
fprintf('\n=== 全部仿真完成，用时 %.1f 分钟 ===\n', R.runtime/60);
fprintf('    数据: %s\n', fullfile(results_dir, 'ntn_sim_main_results.mat'));
fprintf('    图像: %s\\*.png\n\n', results_dir);

end

%% ========================================================================
%  §4  默认参数
%  ========================================================================
function log_mod_done(mod_name, tm, mod_done, n_total, t0)
% 打印单模块耗时 + 总体进度 + ETA
mod_dur = toc(tm);
elapsed = toc(t0);
if mod_done > 0 && mod_done < n_total
    avg_per_mod = elapsed / mod_done;
    eta_min = avg_per_mod * (n_total - mod_done) / 60;
    fprintf('    [done] %-18s | 本模块 %5.1fs | 累计 %5.1f 分钟 | 剩余 %d/%d 模块 ≈ %.1f 分钟\n', ...
            mod_name, mod_dur, elapsed/60, n_total - mod_done, n_total, eta_min);
else
    fprintf('    [done] %-18s | 本模块 %5.1fs | 累计 %5.1f 分钟 | 全部完成\n', ...
            mod_name, mod_dur, elapsed/60);
end
end

function save_checkpoint(mat_path, R, cfg, p, t0)
% 每个模块跑完后立即写 .mat 增量保存（防崩溃）
R.cfg     = cfg;
R.p_nr    = p;
R.runtime = toc(t0);
save(mat_path, '-struct', 'R');
fprintf('    [checkpoint] 已保存 %s （累计 %.1f 分钟）\n', ...
        mat_path, R.runtime/60);
end

function cfg = local_default_cfg()

% ---- 运行模块（可用 'all' 或子集） ----
cfg.modules = {'all'};

% ---- 快速验证模式：true=秒级出 7 张图验证代码正确性；false=真实仿真 ----
cfg.fast_verify = true;

% ---- 随机种子：固定便于复现（答辩可用于 "结果可重复" 说明）----
cfg.rng_seed = 20260429;

% ---- NR 帧结构 ----
cfg.mu          = 1;          % Numerology: 0=15kHz,1=30kHz,2=60kHz
cfg.n_rb        = 25;         % RB 数（25RB≈10MHz @ 30kHz SCS）
cfg.mod_order   = 4;          % 2=BPSK 4=QPSK 16=16QAM 64=64QAM
cfg.coding      = true;       % 是否使用 NR LDPC
cfg.ldpc_rate   = 1/3;
cfg.ldpc_iter   = 10;         % SMS 解码迭代数（快速模式）

% ---- LEO 轨道 / NTN 信道 ----
cfg.altitude_m      = 600e3;  % LEO-600
cfg.fc_hz           = 2e9;    % S 频段
cfg.elevation_deg   = 45;     % 默认仰角
cfg.mobile_speed    = 3;      % 地面用户速度 (m/s)
cfg.enable_shadow   = false;  % 答辩对比图关闭阴影（保证可复现）
cfg.shadow_std_db   = 2.0;

% ---- 仿真实验默认扫描 ----
cfg.snr_dB_vec  = -5:3:20;    % 9 个 SNR 点
% 统计量：按 10·(1/BER) 最低样本经验值设定；Monte-Carlo 相对误差 ~7%
cfg.n_frames    = 600;        % 每点最多帧数（从 60 提至 600）
cfg.min_errors  = 200;        % 达到后提前停止（从 50 提至 200）
cfg.min_frames  = 30;         % 统计有效的最少帧数（避免低 SNR 过早停）

% ---- 各模块独立参数 ----
cfg.K_factor_list_db   = [-Inf, 0, 7, 15];
cfg.K_labels           = {'Rayleigh','K=0 dB','K=7 dB','K=15 dB'};

% 多普勒扫描：覆盖 TR 38.811 LEO-600 @ 2 GHz 真实最大多普勒 (~±48 kHz)
% 子载波间距 30 kHz，48 kHz ≈ 1.6×SCS，ICI 完全破坏正交性
cfg.doppler_list_hz    = [0, 500, 2000, 7000, 24000, 48000];
cfg.doppler_labels     = {'0 Hz','500 Hz','2 kHz','7 kHz','24 kHz','48 kHz (LEO 峰值)'};

% 时延（μs）：扫描残余定时偏差，对比 TA 补偿效果
cfg.delay_list_us      = [0, 1, 5, 10, 30];
cfg.delay_labels       = {'0 μs','1 μs','5 μs','10 μs','30 μs'};

% 雨衰：两个频段、各四个雨率
cfg.rain_bands_hz      = [2e9, 20e9];
cfg.rain_band_labels   = {'S频段 2 GHz','Ka频段 20 GHz'};
cfg.rain_rate_list     = [0, 5, 25, 50];
cfg.rain_labels        = {'无雨','小雨 5 mm/h','中雨 25 mm/h','大雨 50 mm/h'};

% 补偿对比（M7）
cfg.gnss_err_hz_list   = [0, 1e4];    % GNSS 补偿残余误差（0=完美补偿）
cfg.gnss_labels        = {'GNSS 补偿后','无补偿'};
cfg.chest_methods      = {'LS_INTERP','MMSE','2D_WIENER'};

% 路径损耗（M5）
cfg.pathloss_elev_deg  = 10:5:90;

% ---- fast_verify 覆盖：中等规模验证代码正确性 ----
% 目标：3~5 分钟跑完 7 图，BER 曲线形状物理上可信（K 因子方向对、BER ≤ 0.5）
% 若仍要更快，手动把 n_frames 降到 40 即可
if cfg.fast_verify
    cfg.snr_dB_vec        = [-5, 0, 5, 10, 15];      % 5 SNR 点（加密 0 dB / 10 dB）
    cfg.n_frames          = 80;                      % 每点最多 80 帧（15→80 消除方向反转）
    cfg.min_errors        = 25;                      % 25 个块错即停（5→25）
    cfg.min_frames        = 15;                      % 至少 15 帧保证 K 因子差异可辨
    cfg.K_factor_list_db  = [-Inf, 7];               % 2 条 K 曲线
    cfg.K_labels          = {'Rayleigh','K=7 dB'};
    cfg.doppler_list_hz   = [0, 7000, 48000];        % 3 条 fd 曲线
    cfg.doppler_labels    = {'0 Hz','7 kHz','48 kHz'};
    cfg.delay_list_us     = [0, 5, 30];              % 3 条时延
    cfg.delay_labels      = {'0 μs','5 μs','30 μs'};
    cfg.rain_rate_list    = [0, 50];                 % 无雨 vs 大雨
    cfg.rain_labels       = {'无雨','大雨 50 mm/h'};
    cfg.chest_methods     = {'LS_INTERP','MMSE'};    % 2 种方法
    cfg.pathloss_elev_deg = [10, 30, 60, 90];        % 4 个仰角点
end

end

%% ========================================================================
%  §5  NR 参数推导
%  ========================================================================
function [p, cfg] = local_expand_nr_params(cfg)
scs_hz_list = [15e3, 30e3, 60e3, 120e3];
p.mu          = cfg.mu;
p.scs_hz      = scs_hz_list(cfg.mu + 1);
p.sym_per_slot= 14;
switch cfg.mu
    case 0;  p.N_fft = 1024;
    case 1;  p.N_fft = 512;
    case 2;  p.N_fft = 256;
    otherwise; p.N_fft = 512;
end
p.cp_normal    = round(p.N_fft * 144/2048);
p.cp_first     = round(p.N_fft * 160/2048);
p.cp_lengths   = [p.cp_first, repmat(p.cp_normal, 1, 13)];
p.T_slot_s     = 1e-3 / (2^cfg.mu);
p.fs_hz        = p.scs_hz * p.N_fft;
p.N_sc         = cfg.n_rb * 12;
p.N_slot_samps = sum(p.cp_lengths) + p.sym_per_slot * p.N_fft;

% DMRS: 符号 3/8/12（1-based），梳状-2
p.dmrs_sym = [3, 8, 12];
p.dmrs_re  = 1:2:p.N_sc;
p.data_sym = setdiff(1:p.sym_per_slot, p.dmrs_sym);
p.N_data_re= length(p.data_sym) * p.N_sc;

bits_per_re = log2(cfg.mod_order);
cfg.ldpc_K  = min(256, floor(p.N_data_re * bits_per_re * cfg.ldpc_rate));

% 发射端 DMRS 序列（复用同一套，所有实验共用）
p.tx_dmrs = gen_dmrs_seq(length(p.dmrs_re), 0);

end

%% ========================================================================
%  §6  各实验模块
%  ========================================================================

% --- M1 : 理想 AWGN BER（仿真 + 理论）------------------------------------
function out = exp_ideal(cfg, p, results_dir)
fprintf('[M1] 理想 AWGN BER（未编码 BER vs 理论 Q 函数）...\n');
n_snr = length(cfg.snr_dB_vec);
ber_sim = zeros(1, n_snr);

% 临时覆盖：理论参考曲线对应未编码 BER，这里强制关闭 LDPC
cfg_m1 = cfg;
cfg_m1.coding = false;

ntn_p = make_ntn_params(cfg_m1, 100, 0);   % 占位值（awgn_only=true 时不使用）
for si = 1:n_snr
    snr_dB = cfg_m1.snr_dB_vec(si);
    opt = make_default_opt();
    opt.awgn_only  = true;               % 关闭所有衰落，走纯 AWGN 参考
    [ber_sim(si), ~] = run_link_ext(cfg_m1, p, ntn_p, snr_dB, opt);
    fprintf('    SNR=%4.1f dB -> BER=%.3e\n', snr_dB, ber_sim(si));
end

% 理论 BER（AWGN，Gray 编码 QPSK/BPSK 对数域 Es/N0 已等于 SNR 线性）
ber_theory = theory_awgn_ber(cfg.snr_dB_vec, cfg.mod_order);

out.snr_dB     = cfg.snr_dB_vec;
out.ber_sim    = ber_sim;
out.ber_theory = ber_theory;
plot_ideal(out, cfg, results_dir);
end

% --- M2 : K 因子（BER/BLER）----------------------------------------------
function out = exp_kfactor(cfg, p, results_dir)
fprintf('[M2] K 因子对 BER/BLER 影响 ...\n');
n_K = length(cfg.K_factor_list_db);
n_snr = length(cfg.snr_dB_vec);
ber = zeros(n_K, n_snr);
bler= zeros(n_K, n_snr);
for ki = 1:n_K
    K_dB = cfg.K_factor_list_db(ki);
    fprintf('    K=%s dB ...\n', num2str(K_dB));
    ntn_p = make_ntn_params(cfg, K_dB, 0);
    for si = 1:n_snr
        opt = make_default_opt();
        [ber(ki,si), bler(ki,si)] = run_link_ext(cfg, p, ntn_p, cfg.snr_dB_vec(si), opt);
    end
end
out.snr_dB = cfg.snr_dB_vec;
out.ber    = ber;
out.bler   = bler;
plot_kfactor(out, cfg, results_dir);
end

% --- M3 : 多普勒频偏 -----------------------------------------------------
function out = exp_doppler(cfg, p, results_dir)
fprintf('[M3] 多普勒频偏对 BER 影响 ...\n');
n_fd  = length(cfg.doppler_list_hz);
n_snr = length(cfg.snr_dB_vec);
ber = zeros(n_fd, n_snr);
for di = 1:n_fd
    fd = cfg.doppler_list_hz(di);
    fprintf('    fd=%.0f Hz ...\n', fd);
    ntn_p = make_ntn_params(cfg, 7, fd);
    for si = 1:n_snr
        opt = make_default_opt();
        [ber(di,si), ~] = run_link_ext(cfg, p, ntn_p, cfg.snr_dB_vec(si), opt);
    end
end
out.snr_dB = cfg.snr_dB_vec;
out.ber    = ber;
plot_doppler(out, cfg, results_dir);
end

% --- M4 : 残余时延 + TA 补偿对比 ----------------------------------------
% 两子图：
%   (a) 残余时延扫描 0~30 μs，揭示 CP 长度阈值（~4.7 μs @ 30 kHz SCS）
%   (b) LEO-600 真实传播时延 ~2.8 ms：未补偿 vs 启用 TA 预补偿
function out = exp_delay(cfg, p, results_dir)
fprintf('[M4] 残余时延 & TA 补偿对比 ...\n');
n_d   = length(cfg.delay_list_us);
n_snr = length(cfg.snr_dB_vec);
ber_noTA = zeros(n_d, n_snr);

% ---- (a) 残余时延扫描（无 TA 补偿）----
for di = 1:n_d
    d_us = cfg.delay_list_us(di);
    fprintf('    残余时延=%g μs （无 TA 补偿）...\n', d_us);
    ntn_p = make_ntn_params(cfg, 7, 0);
    for si = 1:n_snr
        opt = make_default_opt();
        opt.ta_en       = false;
        opt.residual_us = d_us;
        [ber_noTA(di,si), ~] = run_link_ext(cfg, p, ntn_p, cfg.snr_dB_vec(si), opt);
    end
end

% ---- (b) LEO 真实传播时延场景 ----
% 斜距（仰角 45°，h=600 km）≈ 848 km → τ ≈ 2.83 ms，远超 CP 长度
% 无 TA：残余时延 = τ（mod slot 采样数等效为大的 FFT 窗错位） → BER=0.5 平台
% 启用 TA：历元对齐残余 ~1 μs，接近理想基线
c_light  = 3e8;
elev_rad = cfg.elevation_deg * pi/180;
d_slant  = slant_range_m(cfg.altitude_m, max(elev_rad, 5*pi/180));
tau_prop_us = (d_slant / c_light) * 1e6;
fprintf('    LEO 真实传播时延 = %.2f ms （未补偿 vs 启用 TA）...\n', tau_prop_us/1e3);

ber_leoNoTA = zeros(1, n_snr);
ber_leoTA   = zeros(1, n_snr);
for si = 1:n_snr
    ntn_p = make_ntn_params(cfg, 7, 0);
    ntn_p.delay_s = tau_prop_us * 1e-6;   % 注入真实传播时延，交给 ntn_channel_model
    % 无 TA：信道时延完整作用于信号 → FFT 窗错位
    opt = make_default_opt();
    opt.ta_en          = false;
    opt.use_real_delay = true;
    [ber_leoNoTA(si), ~] = run_link_ext(cfg, p, ntn_p, cfg.snr_dB_vec(si), opt);
    % 启用 TA 闭环：ntn_ta_precomp 动态计算预补偿量，残余 ≈ 0
    opt = make_default_opt();
    opt.ta_en          = true;
    opt.use_real_delay = true;
    [ber_leoTA(si), ~] = run_link_ext(cfg, p, ntn_p, cfg.snr_dB_vec(si), opt);
end

out.snr_dB        = cfg.snr_dB_vec;
out.ber_noTA      = ber_noTA;
out.ber_leoNoTA   = ber_leoNoTA;
out.ber_leoTA     = ber_leoTA;
out.tau_prop_us   = tau_prop_us;
% CP 长度（μs）：用第 2..14 个符号（普通 CP）计算
out.cp_len_us     = (p.cp_normal / p.fs_hz) * 1e6;
plot_delay(out, cfg, results_dir);
end

% --- M5 : 路径损耗与仰角 -------------------------------------------------
function out = exp_pathloss(cfg, p, results_dir)
fprintf('[M5] 仰角 vs 路径损耗（FSPL + 大气 + 雨衰）...\n');
elev_vec = cfg.pathloss_elev_deg;
n_e = length(elev_vec);

fspl_dB = zeros(1, n_e);
atmo_dB = zeros(1, n_e);
rain_dB_S  = zeros(1, n_e);
rain_dB_Ka = zeros(1, n_e);

c = 3e8;
lambda_S  = c / 2e9;
lambda_Ka = c / 20e9;

for ei = 1:n_e
    elev = elev_vec(ei);
    elev_rad = max(elev*pi/180, 5*pi/180);
    d_slant  = slant_range_m(cfg.altitude_m, elev_rad);
    fspl_dB(ei) = 20*log10(4*pi*d_slant/lambda_S);   % 以 S 频段 FSPL 代表
    atmo_dB(ei) = 0.1 / sin(elev_rad);               % S 频段大气吸收近似
    rain_dB_S(ei)  = ntn_rain_attenuation(2e9,  elev, 25);
    rain_dB_Ka(ei) = ntn_rain_attenuation(20e9, elev, 25);
end

% 另做固定 SNR 下的 BER vs 仰角（K 由 TR 38.811 Table 6.7.2-2 查表）
% Suburban S-band LOS 场景下 Rician K 因子（均值 μ，dB），每 10° 一个采样点
%   来源：3GPP TR 38.811 v15.4.0, Table 6.7.2-2a (Suburban, S-band)
% 同时开启 enable_pathloss 让仰角通过路损真正影响链路（归一化到仰角 90° 为 0 dB）
fprintf('    固定 SNR=10dB 下 BER vs 仰角（含路损 + K 查表）...\n');
elev_ref_deg = [10, 20, 30, 40, 50, 60, 70, 80, 90];
K_ref_dB     = [ 2.1, 5.6, 6.4, 7.4, 8.0, 8.5, 8.7, 9.1, 9.2];
% 以 90° 为基准的相对路损（自由空间 + 大气吸收）
PL90_dB = 20*log10(4*pi*cfg.altitude_m/lambda_S) + 0.1;
ber_fix_snr = zeros(1, n_e);
K_used_dB   = zeros(1, n_e);
rel_PL_dB   = zeros(1, n_e);
for ei = 1:n_e
    elev = elev_vec(ei);
    K_dB = interp1(elev_ref_deg, K_ref_dB, elev, 'linear', 'extrap');
    K_used_dB(ei) = K_dB;
    ntn_p = make_ntn_params(cfg, K_dB, 0);
    ntn_p.elevation_deg = elev;
    % 把仰角对路损的影响折算为等效 SNR 损失（避免 FSPL 绝对值压垮链路）
    rel_PL_dB(ei) = fspl_dB(ei) + atmo_dB(ei) - PL90_dB;
    snr_eff = 10 - rel_PL_dB(ei);
    opt = make_default_opt();
    [ber_fix_snr(ei), ~] = run_link_ext(cfg, p, ntn_p, snr_eff, opt);
end

out.elev_deg     = elev_vec;
out.fspl_S_dB    = fspl_dB;
out.atmo_dB      = atmo_dB;
out.rain_S_dB_25 = rain_dB_S;
out.rain_Ka_dB_25= rain_dB_Ka;
out.ber_fixSNR   = ber_fix_snr;
out.K_used_dB    = K_used_dB;
out.rel_PL_dB    = rel_PL_dB;
plot_pathloss(out, cfg, results_dir);
end

% --- M6 : 雨衰 BER（Ka 频段 × 不同雨率，S 频段仅作基线）-----------------
% S 频段 2 GHz 时 ITU-R P.838 的 k_H ~ 3e-5，各雨率下雨衰 <0.01 dB 可忽略，
% 故 S 频段只跑 1 条无雨曲线作为基线，Ka 频段扫描 4 个雨率。
function out = exp_rain(cfg, p, results_dir)
fprintf('[M6] 雨衰对 BER 影响（S 频段基线 vs Ka 频段扫描）...\n');
n_snr = length(cfg.snr_dB_vec);

% S 频段雨衰查表（仅用于副标题说明）
rain_dB_S = zeros(1, length(cfg.rain_rate_list));
for ri = 1:length(cfg.rain_rate_list)
    rain_dB_S(ri) = ntn_rain_attenuation(2e9, cfg.elevation_deg, cfg.rain_rate_list(ri));
end

% S 频段基线：无雨（其它雨率在 S 频段结果等价，不重复跑）
fprintf('    S 频段 2 GHz 基线（雨衰<%.3f dB 可忽略）...\n', max(rain_dB_S)+1e-3);
ber_S = zeros(1, n_snr);
ntn_p = make_ntn_params(cfg, 7, 0);
ntn_p.fc_hz = 2e9;
for si = 1:n_snr
    opt = make_default_opt();
    [ber_S(si), ~] = run_link_ext(cfg, p, ntn_p, cfg.snr_dB_vec(si), opt);
end

% Ka 频段：扫描四个雨率
rains_Ka = cfg.rain_rate_list;
ber_Ka   = zeros(length(rains_Ka), n_snr);
rain_dB_Ka = zeros(1, length(rains_Ka));
for ri = 1:length(rains_Ka)
    rate = rains_Ka(ri);
    rain_dB_Ka(ri) = ntn_rain_attenuation(20e9, cfg.elevation_deg, rate);
    fprintf('    Ka 频段 20 GHz rain=%g mm/h (雨衰=%.2f dB) ...\n', rate, rain_dB_Ka(ri));
    ntn_p = make_ntn_params(cfg, 7, 0);
    ntn_p.fc_hz               = 20e9;
    ntn_p.enable_rain_fade    = (rate > 0);
    ntn_p.rain_rate_mm_per_hr = rate;
    for si = 1:n_snr
        opt = make_default_opt();
        [ber_Ka(ri, si), ~] = run_link_ext(cfg, p, ntn_p, cfg.snr_dB_vec(si), opt);
    end
end

out.snr_dB     = cfg.snr_dB_vec;
out.ber_S      = ber_S;
out.ber_Ka     = ber_Ka;
out.rains_Ka   = rains_Ka;
out.rain_dB_S  = rain_dB_S;
out.rain_dB_Ka = rain_dB_Ka;
plot_rain(out, cfg, results_dir);
end

% --- M7 : 补偿前后对比（GNSS 频偏补偿 & 信道估计算法）-------------------
function out = exp_compensation(cfg, p, results_dir)
fprintf('[M7] GNSS 频偏补偿 + 信道估计算法对比 ...\n');
fd_true = cfg.doppler_list_hz(end);     % 用最大多普勒（最难场景）
n_snr   = length(cfg.snr_dB_vec);

% (a) GNSS 频偏补偿前/后 BER
ber_nocomp = zeros(1, n_snr);
ber_comp   = zeros(1, n_snr);
fprintf('    GNSS 补偿对比（fd=%g Hz, 残余误差=10 Hz）...\n', fd_true);
for si = 1:n_snr
    ntn_p = make_ntn_params(cfg, 7, fd_true);
    opt = make_default_opt();
    opt.compensate_fd = false;
    [ber_nocomp(si), ~] = run_link_ext(cfg, p, ntn_p, cfg.snr_dB_vec(si), opt);
    opt.compensate_fd = true;
    opt.gnss_err_hz   = 10;             % GNSS 历元对齐残余 ~10 Hz
    [ber_comp(si), ~] = run_link_ext(cfg, p, ntn_p, cfg.snr_dB_vec(si), opt);
end

% (b) 不同信道估计方法 BER（固定 K=7dB, 启用 GNSS 补偿）
methods = cfg.chest_methods;
ber_chest = zeros(length(methods), n_snr);
for mi = 1:length(methods)
    fprintf('    信道估计方法: %s ...\n', methods{mi});
    ntn_p = make_ntn_params(cfg, 7, fd_true);
    for si = 1:n_snr
        opt = make_default_opt();
        opt.compensate_fd = true;
        opt.gnss_err_hz   = 10;
        opt.chest_method  = methods{mi};
        [ber_chest(mi,si), ~] = run_link_ext(cfg, p, ntn_p, cfg.snr_dB_vec(si), opt);
    end
end

out.snr_dB      = cfg.snr_dB_vec;
out.ber_nocomp  = ber_nocomp;
out.ber_comp    = ber_comp;
out.methods     = methods;
out.ber_chest   = ber_chest;
plot_compensation(out, cfg, results_dir);
end

%% ========================================================================
%  §7  核心仿真：扩展版 run_link（支持 GNSS 补偿 / 残余时延 / 雨衰）
%  ========================================================================
function opt = make_default_opt()
opt.compensate_fd = false;
opt.gnss_err_hz   = 0;
opt.chest_method  = 'LS_INTERP';
opt.ta_en         = true;         % 默认启用 TA（残余时延≈0）
opt.residual_us   = 0;            % μs，强制残余时延（供 M4a 使用）
opt.awgn_only     = false;        % true => 纯 AWGN，无衰落（供 M1 使用）
opt.use_real_delay= false;        % true => 信道内保留真实传播时延（供 M4b 测试 TA 闭环）
end

function [ber, bler] = run_link_ext(cfg, p, ntn_p, snr_dB, opt)

N_sc     = p.N_sc;
N_sym    = p.sym_per_slot;
N_fft    = p.N_fft;
fs_hz    = p.fs_hz;
cp_vec   = p.cp_lengths;
data_sym = p.data_sym;
tx_dmrs  = p.tx_dmrs;
dsym     = p.dmrs_sym;
dre      = p.dmrs_re;

n_bit_err=0; n_bit_tot=0; n_blk_err=0; n_blk_tot=0;

for frm = 1:cfg.n_frames
    %-- 发送比特
    if cfg.coding
        info_bits = randi([0 1], cfg.ldpc_K, 1);
        [coded_bits, enc_p] = nr_ldpc_encode(info_bits, cfg.ldpc_rate, 'auto');
        tx_bits   = coded_bits;
    else
        n_tx = length(data_sym) * N_sc * log2(cfg.mod_order);
        tx_bits   = randi([0 1], n_tx, 1);
        info_bits = tx_bits;
        enc_p     = [];
    end

    %-- 调制
    tx_syms = qam_mod(tx_bits, cfg.mod_order);

    %-- 映射到资源网格
    tx_grid = zeros(N_sc, N_sym);
    re_avail = length(data_sym) * N_sc;
    n_fill   = min(length(tx_syms), re_avail);
    data_vec = zeros(re_avail, 1);
    data_vec(1:n_fill) = tx_syms(1:n_fill);
    col_idx = 0;
    for l = data_sym
        tx_grid(:, l) = data_vec(col_idx*N_sc+1 : (col_idx+1)*N_sc);
        col_idx = col_idx + 1;
    end
    for l = dsym
        tx_grid(dre, l) = tx_dmrs;
    end

    %-- OFDM 调制
    tx_wave = ofdm_mod(tx_grid, N_fft, cp_vec, N_sc);

    %-- NTN 信道（复用 ntn_channel_model）
    n_samp = size(tx_wave, 2);
    if opt.awgn_only
        rx_wave_ch = tx_wave;   % 纯 AWGN 参考路径，跳过衰落/阴影
    else
        ntn_p_use = ntn_p;
        ntn_p_use.ta_precomp_en = opt.ta_en;        % 由 opt 控制
        if opt.use_real_delay
            % M4b 专用：保留信道内真实传播时延 + TA 闭环
            % delay_s 由 ntn_channel_model 按轨道自动计算（若未传入）
            if ~isfield(ntn_p_use, 'delay_s') || ntn_p_use.delay_s <= 0
                ntn_p_use.delay_s = 0;   % 让 channel_model 自动算
            end
            % 时间基准用实际 slot 时间（ms），让 TA 动态模式有意义
            t_ms = frm * (1e3 / (2^cfg.mu));   % 1 slot = 0.5ms @ μ=1
        else
            ntn_p_use.delay_s       = 1e-12;   % 其它实验：忽略大时延，只保留信道本身
            t_ms = frm;
        end
        [rx_wave_ch, ~] = ntn_channel_model(tx_wave, ntn_p_use, t_ms, n_samp, fs_hz);
    end

    %-- 强制残余时延（M4 专用）
    % 物理正确做法：零填充 + 截尾（接收开头一段没收到，末尾丢失）
    % 对超过 1 个 slot 的大时延，信号整段被推出 FFT 窗 → 接收到的都是 0 → BER≈0.5
    if opt.residual_us > 0
        n_shift = round(opt.residual_us * 1e-6 * fs_hz);
        if n_shift >= n_samp
            rx_wave_ch = zeros(size(rx_wave_ch));   % 完全错过 FFT 窗
        else
            rx_wave_ch = [zeros(size(rx_wave_ch,1), n_shift), ...
                          rx_wave_ch(:, 1:end-n_shift)];
        end
    end

    %-- 施加多普勒频偏（真实值）
    if isfield(ntn_p, 'extra_fd_hz') && ntn_p.extra_fd_hz ~= 0
        t_vec = (0:n_samp-1) / fs_hz;
        rx_wave_ch = rx_wave_ch .* exp(1j * 2*pi * ntn_p.extra_fd_hz * t_vec);
    end

    %-- 加 AWGN：统一 per-RE Es/N0 = SNR 定义
    % 频域每 RE 发送符号能量 Es=1，时域每采样噪声方差 = N0 = 1/SNR_lin
    % （DFT 保持方差，所以时域 = 频域噪声方差）
    noise_var = 1 / 10^(snr_dB/10);
    rx_wave   = rx_wave_ch + sqrt(noise_var/2) * ...
                (randn(size(rx_wave_ch)) + 1j*randn(size(rx_wave_ch)));

    %-- GNSS 频偏补偿
    if opt.compensate_fd && isfield(ntn_p, 'extra_fd_hz') && ntn_p.extra_fd_hz ~= 0
        t_vec  = (0:n_samp-1) / fs_hz;
        fd_est = ntn_p.extra_fd_hz - opt.gnss_err_hz;   % GNSS 残差
        rx_wave = rx_wave .* exp(-1j * 2*pi * fd_est * t_vec);
    end

    %-- OFDM 解调
    rx_grid = ofdm_demod(rx_wave, N_fft, cp_vec, N_sc, N_sym);

    %-- 信道估计与均衡
    if opt.awgn_only
        rx_eq      = rx_grid;                          % H=1, 跳过估计
        sigma2_grid= noise_var * ones(N_sc, N_sym);    % per-RE 噪声方差
    else
        dmrs_cfg_s.re_positions  = dre - 1;
        dmrs_cfg_s.sym_positions = dsym - 1;
        nr_cfg_s.n_rb_max        = cfg.n_rb;
        nr_cfg_s.T_slot_s        = p.T_slot_s;
        nr_cfg_s.sym_per_slot    = p.sym_per_slot;
        chest_p.method           = opt.chest_method;
        chest_p.snr_db           = snr_dB;
        chest_p.noise_var        = noise_var;
        [H_est, ~] = nr_chest_dmrs_standalone(rx_grid, tx_dmrs, dmrs_cfg_s, nr_cfg_s, chest_p);
        rx_eq = zeros(N_sc, N_sym);
        H_pwr = abs(H_est).^2;
        H_pwr = max(H_pwr, 1e-6);           % 防除零
        for l = 1:N_sym
            h_col = H_est(:, l);
            h_col(abs(h_col) < 1e-6) = 1e-6;
            rx_eq(:, l) = rx_grid(:, l) ./ h_col;
        end
        % per-RE 均衡后噪声方差：N0 / |H(k,l)|²（深衰落 RE 会被正确"去权"）
        sigma2_grid = noise_var ./ H_pwr;
    end

    %-- 抽取数据符号（信号 + per-RE 噪声方差）
    rx_data    = zeros(re_avail, 1);
    sigma2_vec = zeros(re_avail, 1);
    col_idx = 0;
    for l = data_sym
        idx = col_idx*N_sc+1 : (col_idx+1)*N_sc;
        rx_data(idx)    = rx_eq(:, l);
        sigma2_vec(idx) = sigma2_grid(:, l);
        col_idx = col_idx + 1;
    end
    rx_data    = rx_data(1:n_fill);
    sigma2_vec = sigma2_vec(1:n_fill);

    %-- 软解映射（per-RE sigma²，深衰落 RE 的 LLR 自动降权）
    llr = qam_demod_soft(rx_data, cfg.mod_order, sigma2_vec);

    %-- 解码
    if cfg.coding
        [dec_bits, ~, ~] = nr_ldpc_decode(llr, enc_p, cfg.ldpc_iter, 'SMS');
        dec_bits = dec_bits(1:cfg.ldpc_K);
        info_err = sum(dec_bits ~= info_bits);
        n_bit_err = n_bit_err + info_err;
        n_bit_tot = n_bit_tot + cfg.ldpc_K;
        % BLER = 信息位块错误率（业界标准，不依赖 BP 校验位收敛）
        n_blk_err = n_blk_err + (info_err > 0);
        n_blk_tot = n_blk_tot + 1;
    else
        hard = double(llr < 0);
        hard = hard(1:length(info_bits));
        n_bit_err = n_bit_err + sum(hard ~= info_bits);
        n_bit_tot = n_bit_tot + length(info_bits);
        n_blk_err = n_blk_err + (sum(hard ~= info_bits) > 0);
        n_blk_tot = n_blk_tot + 1;
    end

    % 自适应停止：达到最小错误数且跑够最少帧数（保证统计有效性）
    min_frm = 10;
    if isfield(cfg, 'min_frames'), min_frm = cfg.min_frames; end
    if n_blk_err >= cfg.min_errors && frm >= min_frm
        break;
    end
end

ber  = n_bit_err / max(n_bit_tot, 1);
bler = n_blk_err / max(n_blk_tot, 1);

% BER 物理上限 = 0.5（超过即可反转比特判决得到 ≤0.5），
% 大于 0.5 一定是统计涨落 → 提醒用户增加 n_frames
if ber > 0.5 && n_blk_tot < 100
    fprintf('        [warn] BER=%.3f>0.5 (samples=%d)，样本量不足，建议增加 n_frames\n', ...
            ber, n_blk_tot);
end
end

%% ========================================================================
%  §8  NTN 信道参数构造
%  ========================================================================
function ntn_p = make_ntn_params(cfg, K_dB, extra_fd_hz)
ntn_p.fc_hz          = cfg.fc_hz;
ntn_p.altitude_m     = cfg.altitude_m;
ntn_p.elevation_deg  = cfg.elevation_deg;
ntn_p.mobile_speed   = cfg.mobile_speed;
ntn_p.K_factor_db    = K_dB;
ntn_p.enable_pathloss= false;    % 路径损耗单独在 M5 展示
ntn_p.enable_shadow  = cfg.enable_shadow;
ntn_p.shadow_std_db  = cfg.shadow_std_db;
% Rician 使能分支（写清楚三种情形）
%   K = -Inf  → 纯 Rayleigh，enable_rician=true，Rician 分支内部 K_lin=0 退化为散射
%   K = +Inf  → 纯 AWGN baseline，enable_rician=false 直接设 H=1
%   其他数值  → Rician 衰落
if isinf(K_dB) && K_dB > 0
    ntn_p.enable_rician = false;
else
    ntn_p.enable_rician = true;
end

ntn_p.delay_s        = 0;
ntn_p.ta_precomp_en  = false;    % 默认由 opt 接管
ntn_p.pass_center_ms = 0;
ntn_p.elevation_min_deg = 10;
ntn_p.enable_rain_fade    = false;
ntn_p.rain_rate_mm_per_hr = 0;
ntn_p.extra_fd_hz    = extra_fd_hz;
end

%% ========================================================================
%  §9  绘图函数（每个模块独立一张）
%  ========================================================================
function plot_ideal(r, cfg, results_dir)
fh = figure('Name','M1 理想 AWGN BER','NumberTitle','off','Position',[120 120 720 480]);
y_sim = r.ber_sim;    y_sim(y_sim<=0) = NaN;   % 0 错误帧断线
y_th  = r.ber_theory; y_th(y_th<=0)   = NaN;
semilogy(r.snr_dB, y_sim, 'bs-', 'LineWidth', 1.8, 'MarkerSize', 7); hold on;
semilogy(r.snr_dB, y_th,  'k--', 'LineWidth', 1.8);
grid on; ylim([1e-5 1]);
xlabel('per-RE E_s/N_0 (dB)'); ylabel('BER');
title(sprintf('M1: 理想 AWGN BER（%s, 仿真 vs 理论）', mod_name(cfg.mod_order)));
legend({'仿真','理论'}, 'Location','southwest');
drawnow; saveas(fh, fullfile(results_dir, 'M1_ideal_AWGN_BER.png'));
end

function plot_kfactor(r, cfg, results_dir)
fh = figure('Name','M2 K 因子','NumberTitle','off','Position',[120 120 1200 480]);
subplot(1,2,1);
plot_curve_set(r.snr_dB, r.ber, cfg.K_labels, 1e-5);
xlabel('SNR (dB)'); ylabel('BER');
title('M2a: BER vs SNR（不同 Rician K 因子）');

subplot(1,2,2);
plot_curve_set(r.snr_dB, r.bler, cfg.K_labels, 1e-3);
xlabel('SNR (dB)'); ylabel('BLER');
title('M2b: BLER vs SNR（LDPC 编码）');

sgtitle('LEO-NTN 小尺度衰落分析（TR 38.811 §6.7.2）');
drawnow; saveas(fh, fullfile(results_dir, 'M2_Kfactor_BER_BLER.png'));
end

function plot_doppler(r, cfg, results_dir)
fh = figure('Name','M3 多普勒','NumberTitle','off','Position',[120 120 720 480]);
plot_curve_set(r.snr_dB, r.ber, cfg.doppler_labels, 1e-5);
hold on;
% BER=0.5 物理上限参考线（CFO 未补偿最坏情况）
yline(0.5, 'k:', '随机判决上限 BER=0.5', 'LineWidth', 1, ...
      'LabelVerticalAlignment','bottom','LabelHorizontalAlignment','left');
xlabel('SNR (dB)'); ylabel('BER');
title('M3: 多普勒频偏对 BER 影响（LEO-600, K=7 dB）');
drawnow; saveas(fh, fullfile(results_dir, 'M3_doppler_BER.png'));
end

function plot_delay(r, cfg, results_dir)
fh = figure('Name','M4 时延 & TA 补偿','NumberTitle','off','Position',[120 120 1200 480]);

% (a) 残余时延扫描：未补偿场景，观察 CP 长度阈值
subplot(1,2,1);
plot_curve_set(r.snr_dB, r.ber_noTA, cfg.delay_labels, 1e-5);
xlabel('SNR (dB)'); ylabel('BER');
title(sprintf('M4a: 残余时延对 BER 影响（CP ≈ %.2f μs）', r.cp_len_us));

% (b) LEO 真实传播时延 vs TA 补偿效果
subplot(1,2,2);
y1 = r.ber_leoNoTA; y1(y1<=0) = NaN;
y2 = r.ber_leoTA;   y2(y2<=0) = NaN;
semilogy(r.snr_dB, y1, 'r-s', 'LineWidth', 2, 'MarkerSize', 7); hold on;
semilogy(r.snr_dB, y2, 'b-o', 'LineWidth', 2, 'MarkerSize', 7);
grid on; ylim([1e-5 1]);
xlabel('SNR (dB)'); ylabel('BER');
title(sprintf('M4b: LEO 真实传播 %.2f ms：未补偿 vs TA 补偿', r.tau_prop_us/1e3));
legend({sprintf('无 TA（残余 %.2f ms，超 CP）', r.tau_prop_us/1e3), ...
        '启用 TA（残余 1 μs）'}, 'Location','southwest');

sgtitle('TR 38.821：LEO 传播时延与 TA 预补偿效果');
drawnow; saveas(fh, fullfile(results_dir, 'M4_delay_TA.png'));
end

function plot_pathloss(r, cfg, results_dir)
fh = figure('Name','M5 路径损耗','NumberTitle','off','Position',[120 120 1400 520]);

% (a) S 频段路损分量重合（物理事实：2 GHz 雨衰<0.01 dB）
%     只画 S-FSPL/ S-FSPL+大气 两条主干；Ka 用独立右轴对比显示
% M5a：改用单 Y 轴 + 统一 dB 坐标，S/Ka 同图直接对比（避免 yyaxis 标签被裁）
ax1 = subplot(1,2,1);
plot(r.elev_deg, r.fspl_S_dB, 'b-o', 'LineWidth', 1.6, 'MarkerSize',5); hold on;
plot(r.elev_deg, r.fspl_S_dB + r.atmo_dB, 'c--s', 'LineWidth', 1.4, 'MarkerSize',5);
plot(r.elev_deg, r.fspl_S_dB + r.atmo_dB + r.rain_Ka_dB_25 + ...
     20*log10(20e9/2e9), 'm-d', 'LineWidth', 1.8, 'MarkerSize',5);
grid on; xlabel('仰角 (deg)'); ylabel('路径损耗 (dB)');
legend({'S: FSPL','S: FSPL + 大气','Ka: FSPL + 大气 + 雨衰(25mm/h)'}, ...
       'Location','northeast','FontSize',9);
title({'M5a: 仰角 vs 路径损耗', ...
       sprintf('S 雨衰 %.3f dB 可忽略；Ka 雨衰 %.1f dB 主导', ...
               max(r.rain_S_dB_25), max(r.rain_Ka_dB_25))});
set(ax1,'Position', get(ax1,'Position') + [0.01 0 -0.01 0]);  % 微调边距

subplot(1,2,2);
y = r.ber_fixSNR; y(y<=0) = NaN;
semilogy(r.elev_deg, y, 'r-o', 'LineWidth', 1.8, 'MarkerSize', 6);
grid on; xlabel('仰角 (deg)'); ylabel('BER');
ylim([1e-5 1]);
title({'M5b: BER vs 仰角（SNR_{90°}=10 dB，等效路损 + K 查表）', ...
       'K: TR 38.811 Table 6.7.2-2a (Suburban S-band)'});

drawnow; saveas(fh, fullfile(results_dir, 'M5_pathloss_elevation.png'));
end

function plot_rain(r, cfg, results_dir)
fh = figure('Name','M6 雨衰影响','NumberTitle','off','Position',[120 120 1200 480]);

% S 子图：单条基线 + 注释
subplot(1,2,1);
y = r.ber_S; y(y<=0) = NaN;
semilogy(r.snr_dB, y, 'b-o', 'LineWidth', 1.8, 'MarkerSize', 6);
grid on; ylim([1e-5 1]);
xlabel('SNR (dB)'); ylabel('BER');
title({'M6-S 频段 2 GHz（基线）', ...
       sprintf('ITU-R P.838：5/25/50 mm/h 雨衰 <%.3f dB，与无雨无差异', ...
               max(r.rain_dB_S)+1e-3)});
legend({'全雨率下 BER 曲线重合'}, 'Location','southwest');

% Ka 子图：四雨率扫描
subplot(1,2,2);
lbls = cell(1, length(r.rains_Ka));
for ri = 1:length(r.rains_Ka)
    lbls{ri} = sprintf('%s (%.1f dB)', cfg.rain_labels{ri}, r.rain_dB_Ka(ri));
end
plot_curve_set(r.snr_dB, r.ber_Ka, lbls, 1e-5);
xlabel('SNR (dB)'); ylabel('BER');
title('M6-Ka 频段 20 GHz：BER vs SNR（不同雨率）');

sgtitle('雨衰对 NR-NTN 下行链路影响（ITU-R P.618 / TR 38.811 §6.6.5）');
drawnow; saveas(fh, fullfile(results_dir, 'M6_rain_BER.png'));
end

function plot_compensation(r, cfg, results_dir)
fh = figure('Name','M7 补偿前后对比','NumberTitle','off','Position',[120 120 1200 480]);
subplot(1,2,1);
y1 = r.ber_nocomp; y1(y1<=0) = NaN;
y2 = r.ber_comp;   y2(y2<=0) = NaN;
semilogy(r.snr_dB, y1, 'r-s', 'LineWidth', 1.8, 'MarkerSize', 6); hold on;
semilogy(r.snr_dB, y2, 'b-o', 'LineWidth', 1.8, 'MarkerSize', 6);
yline(0.5, 'k:', 'BER=0.5', 'LineWidth', 1);
grid on; ylim([1e-5 1]);
xlabel('SNR (dB)'); ylabel('BER');
title('M7a: GNSS 频偏补偿前后 BER 对比（fd=7 kHz）');
legend({'无补偿','GNSS 补偿（残余 10 Hz）'}, 'Location','southwest');

subplot(1,2,2);
plot_curve_set(r.snr_dB, r.ber_chest, r.methods, 1e-5);
xlabel('SNR (dB)'); ylabel('BER');
title('M7b: 不同信道估计方法 BER（启用 GNSS 补偿）');

sgtitle('NR-NTN 补偿算法效果对比');
drawnow; saveas(fh, fullfile(results_dir, 'M7_compensation.png'));
end

% --- 通用簇曲线绘制（legend 由调用方统一处理）
% NaN 化：y<=0（0 错误帧）→ NaN，避免被 y_min 地板误导成"瀑布反弹"
function plot_curve_set(x, y_mat, labels, y_min)
markers = {'-o','-s','-d','-^','-v','-x','-*'};
colors  = lines(max(size(y_mat,1), 5));
for i = 1:size(y_mat,1)
    y = y_mat(i,:);
    y(y <= 0) = NaN;
    semilogy(x, y, markers{mod(i-1,length(markers))+1}, ...
        'Color', colors(i,:), 'LineWidth', 1.6, 'MarkerSize', 6); hold on;
end
grid on; ylim([y_min 1]);
if ~isempty(labels)
    legend(labels, 'Location','southwest');
end
end

function s = mod_name(M)
switch M
    case 2;  s = 'BPSK';
    case 4;  s = 'QPSK';
    case 16; s = '16QAM';
    case 64; s = '64QAM';
    otherwise; s = sprintf('%d-QAM', M);
end
end

function ber = theory_awgn_ber(snr_dB_vec, M)
% AWGN 理论 BER（Gray 编码）
snr_lin = 10.^(snr_dB_vec/10);   % Es/N0 线性
switch M
    case 2      % BPSK: BER = Q(sqrt(2 Eb/N0)), Eb/N0 = SNR
        ber = qfunc_compat(sqrt(2*snr_lin));
    case 4      % QPSK: BER = Q(sqrt(Es/N0)), Eb/N0 = SNR/2
        ber = qfunc_compat(sqrt(snr_lin));
    case 16     % 16QAM 近似 (Gray): BER ≈ (3/4) Q(sqrt(Es/(5 N0)))
        ber = (3/4) * qfunc_compat(sqrt(snr_lin/5));
    case 64     % 64QAM 近似 (Gray): BER ≈ (7/12) Q(sqrt(Es/(21 N0)))
        ber = (7/12) * qfunc_compat(sqrt(snr_lin/21));
    otherwise
        ber = qfunc_compat(sqrt(snr_lin));
end
end

function v = qfunc_compat(x)
% 不依赖 communications toolbox 的 qfunc 替代
v = 0.5 * erfc(x/sqrt(2));
end

%% ========================================================================
%  §10  OFDM / QAM / DMRS 辅助函数（从 ntn_dl_sim.m 搬入，自包含）
%  ========================================================================
function tx_wave = ofdm_mod(grid, N_fft, cp_vec, N_sc)
[~, N_sym] = size(grid);
dc = N_fft/2 + 1;
n_lower = floor(N_sc/2); n_upper = N_sc - n_lower;
N_samp  = sum(cp_vec) + N_sym * N_fft;
tx_wave = zeros(1, N_samp);
ptr = 1;
for l = 1:N_sym
    freq_buf = zeros(N_fft, 1);
    freq_buf(dc-n_lower : dc-1)       = grid(1:n_lower, l);
    freq_buf(dc+1       : dc+n_upper) = grid(n_lower+1:end, l);
    td = ifft(ifftshift(freq_buf), N_fft) * sqrt(N_fft);
    cp = cp_vec(l);
    sym = [td(end-cp+1:end); td];
    tx_wave(1, ptr:ptr+length(sym)-1) = sym.';
    ptr = ptr + length(sym);
end
end

function rx_grid = ofdm_demod(rx_wave, N_fft, cp_vec, N_sc, N_sym)
dc = N_fft/2 + 1;
n_lower = floor(N_sc/2); n_upper = N_sc - n_lower;
rx_grid = zeros(N_sc, N_sym);
ptr = 1;
for l = 1:N_sym
    cp = cp_vec(l);
    sym_sig = rx_wave(1, ptr : ptr + cp + N_fft - 1);
    ptr = ptr + cp + N_fft;
    td  = sym_sig(cp+1 : end).';
    fd  = fftshift(fft(td, N_fft)) / sqrt(N_fft);
    lower = fd(dc-n_lower : dc-1);
    upper = fd(dc+1       : dc+n_upper);
    rx_grid(:, l) = [lower; upper];
end
end

function syms = qam_mod(bits, M)
bits = bits(:);
k    = log2(M);
n_sym= floor(length(bits)/k);
bits = bits(1:n_sym*k);
bit_mat = reshape(bits, k, n_sym)';
pw      = 2.^(k-1:-1:0);
idx     = bit_mat * pw(:);
switch M
    case 2
        syms = 1 - 2*bits(1:k:end);
        syms = complex(syms, 0);
    case 4
        % QPSK: bit0 控制 I, bit1 控制 Q（bit=0 -> +1, bit=1 -> -1）
        % 与 qam_demod_soft 中 LLR = 2√2·real(y)/σ² 的解映射保持一致
        I = 1 - 2*bits(1:2:end);
        Q = 1 - 2*bits(2:2:end);
        syms = (I + 1j*Q) / sqrt(2);
    case {16, 64}
        c = qammod_const(M);
        syms = c(idx+1).';
    otherwise
        error('qam_mod: 不支持的 M=%d', M);
end
syms = syms(:);
end

function llr = qam_demod_soft(rx_syms, M, noise_var)
% noise_var 可为标量（全部 RE 共用）或向量 [N_sym × 1]（per-RE）
rx_syms = rx_syms(:);
k       = log2(M);
N_sym   = length(rx_syms);
llr     = zeros(N_sym*k, 1);
if isscalar(noise_var)
    sigma2 = max(noise_var, 1e-10) * ones(N_sym, 1);
else
    sigma2 = max(noise_var(:), 1e-10);
end
switch M
    case 2
        llr = 2 * real(rx_syms) ./ sigma2;
    case 4
        % QPSK: syms = (bI + j·bQ)/sqrt(2)，每维幅度 1/√2
        % Max-Log LLR = 2·(1/√2)·real(y)/σ² = √2·real(y)/σ²
        llr(1:2:end) = sqrt(2)*real(rx_syms)./sigma2;
        llr(2:2:end) = sqrt(2)*imag(rx_syms)./sigma2;
    otherwise
        c = qammod_const(M);
        c = c / sqrt(mean(abs(c).^2));
        for si = 1:N_sym
            r = rx_syms(si);
            s2 = sigma2(si);
            for b = 1:k
                all_idx = 0:M-1;
                bit_b   = floor(all_idx / 2^(k-b));
                mask0 = mod(bit_b,2)==0; mask1 = ~mask0;
                d0 = min(abs(r - c(mask0)).^2);
                d1 = min(abs(r - c(mask1)).^2);
                llr((si-1)*k + b) = (d1 - d0) / s2;
            end
        end
end
end

function const = qammod_const(M)
sqM = sqrt(M);
pts = -(sqM-1) : 2 : (sqM-1);
[re_g, im_g] = meshgrid(pts, fliplr(pts));
const = (re_g(:) + 1j*im_g(:)).';
const = const / sqrt(mean(abs(const).^2));
end

function d = slant_range_m(alt_m, elev_rad)
% 含地球曲率的精确斜距（球形地球近似）
% d = -R·sinE + √(R²sin²E + 2R·h + h²)
R_earth = 6371e3;
h = alt_m;
sinE = sin(elev_rad);
d = -R_earth*sinE + sqrt(R_earth^2 * sinE^2 + 2*R_earth*h + h^2);
end

function seq = gen_dmrs_seq(N, n_id)
seq = gold_seq(N, n_id, 0) + 1j*gold_seq(N, n_id, 1);
seq = seq / sqrt(mean(abs(seq).^2));
end

function c = gold_seq(N, n_id, offset)
c_init = mod(n_id + offset, 2^31);
x1 = bitget(c_init + 2^31, 31:-1:1);
x2 = bitget(1 + offset,    31:-1:1);
c = zeros(N, 1);
for n = 1:N
    new_x1 = mod(x1(3) + x1(end), 2);
    new_x2 = mod(x2(3) + x2(1) + x2(end-1) + x2(end), 2);
    c(n)   = mod(x1(end) + x2(end), 2);
    x1 = [x1(2:end), new_x1];
    x2 = [x2(2:end), new_x2];
end
c = (1 - 2*c) / sqrt(2);
end
