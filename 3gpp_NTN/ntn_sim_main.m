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

[p, cfg] = local_expand_nr_params(cfg);
fprintf('\n============================================================\n');
fprintf('  NR-NTN 下行链路仿真（TR 38.811）\n');
fprintf('  载波: %.2f GHz | SCS: %.0f kHz | N_FFT: %d | N_RB: %d\n', ...
        cfg.fc_hz/1e9, p.scs_hz/1e3, p.N_fft, cfg.n_rb);
fprintf('  SNR: %.1f..%.1f dB (%d 点) | 每点最大 %d 帧\n', ...
        cfg.snr_dB_vec(1), cfg.snr_dB_vec(end), ...
        length(cfg.snr_dB_vec), cfg.n_frames);
fprintf('  启用模块: %s\n', strjoin(cfg.modules, ', '));
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

if ismember('ideal',        cfg.modules), R.ideal        = exp_ideal(cfg, p, results_dir);        end
if ismember('kfactor',      cfg.modules), R.kfactor      = exp_kfactor(cfg, p, results_dir);      end
if ismember('doppler',      cfg.modules), R.doppler      = exp_doppler(cfg, p, results_dir);      end
if ismember('delay',        cfg.modules), R.delay        = exp_delay(cfg, p, results_dir);        end
if ismember('pathloss',     cfg.modules), R.pathloss     = exp_pathloss(cfg, p, results_dir);     end
if ismember('rain',         cfg.modules), R.rain         = exp_rain(cfg, p, results_dir);         end
if ismember('compensation', cfg.modules), R.compensation = exp_compensation(cfg, p, results_dir); end

R.cfg      = cfg;
R.p_nr     = p;
R.runtime  = toc(t0);

save(fullfile(results_dir, 'ntn_sim_main_results.mat'), '-struct', 'R');
fprintf('\n=== 全部仿真完成，用时 %.1f 分钟 ===\n', R.runtime/60);
fprintf('    数据: %s\n', fullfile(results_dir, 'ntn_sim_main_results.mat'));
fprintf('    图像: %s\\*.png\n\n', results_dir);

end

%% ========================================================================
%  §4  默认参数
%  ========================================================================
function cfg = local_default_cfg()

% ---- 运行模块（可用 'all' 或子集） ----
cfg.modules = {'all'};

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
cfg.n_frames    = 60;         % 快速模式
cfg.min_errors  = 50;         % 达到后提前停止

% ---- 各模块独立参数 ----
cfg.K_factor_list_db   = [-Inf, 0, 7, 15];
cfg.K_labels           = {'Rayleigh','K=0 dB','K=7 dB','K=15 dB'};

cfg.doppler_list_hz    = [0, 500, 2000, 7000];
cfg.doppler_labels     = {'0 Hz','500 Hz','2 kHz','7 kHz'};

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
function out = exp_delay(cfg, p, results_dir)
fprintf('[M4] 残余时延 & TA 补偿对比 ...\n');
n_d   = length(cfg.delay_list_us);
n_snr = length(cfg.snr_dB_vec);
ber_noTA = zeros(n_d, n_snr);
ber_TA   = zeros(1, n_snr);     % 启用 TA 后残余≈0，一条曲线

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

fprintf('    启用 TA 预补偿 （残余≈0） ...\n');
ntn_p = make_ntn_params(cfg, 7, 0);
for si = 1:n_snr
    opt = make_default_opt();
    opt.ta_en       = true;
    opt.residual_us = 0;
    [ber_TA(si), ~] = run_link_ext(cfg, p, ntn_p, cfg.snr_dB_vec(si), opt);
end

out.snr_dB   = cfg.snr_dB_vec;
out.ber_noTA = ber_noTA;
out.ber_TA   = ber_TA;
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
    d_slant  = cfg.altitude_m / sin(elev_rad);
    fspl_dB(ei) = 20*log10(4*pi*d_slant/lambda_S);   % 以 S 频段 FSPL 代表
    atmo_dB(ei) = 0.1 / sin(elev_rad);               % S 频段大气吸收近似
    rain_dB_S(ei)  = ntn_rain_attenuation(2e9,  elev, 25);
    rain_dB_Ka(ei) = ntn_rain_attenuation(20e9, elev, 25);
end

% 另做固定 SNR 下的 BER vs 仰角（K 由仰角自动）
fprintf('    固定 SNR=10dB 下 BER vs 仰角（K 随仰角）...\n');
ber_fix_snr = zeros(1, n_e);
for ei = 1:n_e
    elev = elev_vec(ei);
    K_dB = max(-10, -12 + 0.3*elev);
    ntn_p = make_ntn_params(cfg, K_dB, 0);
    ntn_p.elevation_deg = elev;
    opt = make_default_opt();
    [ber_fix_snr(ei), ~] = run_link_ext(cfg, p, ntn_p, 10, opt);
end

out.elev_deg     = elev_vec;
out.fspl_S_dB    = fspl_dB;
out.atmo_dB      = atmo_dB;
out.rain_S_dB_25 = rain_dB_S;
out.rain_Ka_dB_25= rain_dB_Ka;
out.ber_fixSNR   = ber_fix_snr;
plot_pathloss(out, cfg, results_dir);
end

% --- M6 : 雨衰 BER（S/Ka 双频段 × 不同雨率）----------------------------
function out = exp_rain(cfg, p, results_dir)
fprintf('[M6] 雨衰对 BER 影响（S 频段 vs Ka 频段）...\n');
bands = cfg.rain_bands_hz;
rains = cfg.rain_rate_list;
n_snr = length(cfg.snr_dB_vec);

ber_all = zeros(length(bands), length(rains), n_snr);
rain_dB_table = zeros(length(bands), length(rains));
for bi = 1:length(bands)
    for ri = 1:length(rains)
        fc   = bands(bi);
        rate = rains(ri);
        rain_dB_table(bi, ri) = ntn_rain_attenuation(fc, cfg.elevation_deg, rate);
        fprintf('    band=%.0fGHz rain=%g mm/h (雨衰=%.2f dB) ...\n', ...
                fc/1e9, rate, rain_dB_table(bi, ri));
        ntn_p = make_ntn_params(cfg, 7, 0);
        ntn_p.fc_hz                = fc;
        ntn_p.enable_rain_fade     = (rate > 0);
        ntn_p.rain_rate_mm_per_hr  = rate;
        for si = 1:n_snr
            opt = make_default_opt();
            [ber_all(bi,ri,si), ~] = run_link_ext(cfg, p, ntn_p, cfg.snr_dB_vec(si), opt);
        end
    end
end
out.snr_dB     = cfg.snr_dB_vec;
out.bands_hz   = bands;
out.rains      = rains;
out.ber_all    = ber_all;
out.rain_dB    = rain_dB_table;
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
opt.residual_us   = 0;            % μs，强制残余时延（供 M4 使用）
opt.awgn_only     = false;        % true => 纯 AWGN，无衰落（供 M1 使用）
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

    %-- NTN 信道（复用 ntn_channel_model，不做传播时延，避免与残余时延冲突）
    n_samp = size(tx_wave, 2);
    if opt.awgn_only
        rx_wave_ch = tx_wave;   % 纯 AWGN 参考路径，跳过衰落/阴影
    else
        ntn_p_use = ntn_p;
        ntn_p_use.ta_precomp_en = opt.ta_en;        % 由 opt 控制
        ntn_p_use.delay_s       = 1e-12;            % 忽略大时延，只保留信道本身
        [rx_wave_ch, ~] = ntn_channel_model(tx_wave, ntn_p_use, frm, n_samp, fs_hz);
    end

    %-- 强制残余时延（M4 专用），循环右移等效于接收端延迟采样
    if opt.residual_us > 0
        n_shift = round(opt.residual_us * 1e-6 * fs_hz);
        n_shift = min(n_shift, n_samp-1);
        rx_wave_ch = [zeros(size(rx_wave_ch,1), n_shift), rx_wave_ch(:, 1:end-n_shift)];
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
        rx_eq     = rx_grid;           % H=1, 跳过估计
        sigma2_fd = noise_var;         % rx_grid per-RE 噪声方差
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
        for l = 1:N_sym
            h_col = H_est(:, l);
            h_col(abs(h_col) < 1e-6) = 1e-6;
            rx_eq(:, l) = rx_grid(:, l) ./ h_col;
        end
        % 均衡后等效噪声方差：N0 / |H|²（用均值近似）
        H_pwr_avg = max(mean(abs(H_est(:)).^2), 1e-6);
        sigma2_fd = noise_var / H_pwr_avg;
    end

    %-- 抽取数据符号
    rx_data = zeros(re_avail, 1);
    col_idx = 0;
    for l = data_sym
        rx_data(col_idx*N_sc+1 : (col_idx+1)*N_sc) = rx_eq(:, l);
        col_idx = col_idx + 1;
    end
    rx_data = rx_data(1:n_fill);

    %-- 软解映射（sigma² = 均衡后 per-RE 噪声方差）
    llr = qam_demod_soft(rx_data, cfg.mod_order, sigma2_fd);

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

    if n_blk_err >= cfg.min_errors && frm >= 10
        break;
    end
end

ber  = n_bit_err / max(n_bit_tot, 1);
bler = n_blk_err / max(n_blk_tot, 1);
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
ntn_p.enable_rician  = ~isinf(K_dB) || (K_dB > 0);   % Inf->AWGN 情形
if isinf(K_dB) && K_dB > 0
    ntn_p.enable_rician = false;  % K=+Inf => AWGN baseline
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
semilogy(r.snr_dB, max(r.ber_sim, 1e-6), 'bs-', 'LineWidth', 1.8, 'MarkerSize', 7); hold on;
semilogy(r.snr_dB, max(r.ber_theory, 1e-6), 'k--', 'LineWidth', 1.8);
grid on; ylim([1e-5 1]);
xlabel('SNR = E_s/N_0 (dB)'); ylabel('BER');
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
xlabel('SNR (dB)'); ylabel('BER');
title('M3: 多普勒频偏对 BER 影响（LEO-600, K=7 dB）');
drawnow; saveas(fh, fullfile(results_dir, 'M3_doppler_BER.png'));
end

function plot_delay(r, cfg, results_dir)
fh = figure('Name','M4 时延 & TA 补偿','NumberTitle','off','Position',[120 120 720 480]);
plot_curve_set(r.snr_dB, r.ber_noTA, {}, 1e-5);   % 不先画 legend
semilogy(r.snr_dB, max(r.ber_TA, 1e-5), 'k--*', 'LineWidth', 2, 'MarkerSize', 8);
lg = [cfg.delay_labels, {'启用 TA 补偿'}];
xlabel('SNR (dB)'); ylabel('BER');
title('M4: 残余时延对 BER 影响 & TA 预补偿效果（TR 38.821）');
legend(lg, 'Location','southwest');
drawnow; saveas(fh, fullfile(results_dir, 'M4_delay_TA.png'));
end

function plot_pathloss(r, cfg, results_dir)
fh = figure('Name','M5 路径损耗','NumberTitle','off','Position',[120 120 1200 480]);
subplot(1,2,1);
plot(r.elev_deg, r.fspl_S_dB, 'b-o', 'LineWidth', 1.6, 'MarkerSize',5); hold on;
plot(r.elev_deg, r.fspl_S_dB + r.atmo_dB, 'r-s', 'LineWidth', 1.6);
plot(r.elev_deg, r.fspl_S_dB + r.atmo_dB + r.rain_S_dB_25, 'g-^', 'LineWidth', 1.6);
plot(r.elev_deg, r.fspl_S_dB + r.atmo_dB + r.rain_Ka_dB_25, 'm--d', 'LineWidth', 1.6);
grid on; xlabel('仰角 (deg)'); ylabel('总路径损耗 (dB)');
legend({'FSPL (S频段)','+大气吸收','+雨衰 25mm/h (S)','+雨衰 25mm/h (Ka)'}, 'Location','northeast');
title('M5a: 仰角 vs 路径损耗分量（LEO-600, 雨率 25 mm/h）');

subplot(1,2,2);
semilogy(r.elev_deg, max(r.ber_fixSNR, 1e-5), 'r-o', 'LineWidth', 1.8, 'MarkerSize', 6);
grid on; xlabel('仰角 (deg)'); ylabel('BER');
ylim([1e-5 1]);
title('M5b: 固定 SNR=10 dB 下 BER 随仰角变化（K 自动由仰角决定）');

drawnow; saveas(fh, fullfile(results_dir, 'M5_pathloss_elevation.png'));
end

function plot_rain(r, cfg, results_dir)
fh = figure('Name','M6 雨衰影响','NumberTitle','off','Position',[120 120 1200 480]);
for bi = 1:length(r.bands_hz)
    subplot(1,2,bi);
    ber_mat = squeeze(r.ber_all(bi,:,:));
    lbls = cell(1, length(r.rains));
    for ri = 1:length(r.rains)
        lbls{ri} = sprintf('%s (%.1f dB)', cfg.rain_labels{ri}, r.rain_dB(bi,ri));
    end
    plot_curve_set(r.snr_dB, ber_mat, lbls, 1e-5);
    xlabel('SNR (dB)'); ylabel('BER');
    title(sprintf('M6-%s: BER vs SNR（不同雨率）', cfg.rain_band_labels{bi}));
end
sgtitle('雨衰对 NR-NTN 下行链路影响（ITU-R P.618 / TR 38.811 §6.6.5）');
drawnow; saveas(fh, fullfile(results_dir, 'M6_rain_BER.png'));
end

function plot_compensation(r, cfg, results_dir)
fh = figure('Name','M7 补偿前后对比','NumberTitle','off','Position',[120 120 1200 480]);
subplot(1,2,1);
semilogy(r.snr_dB, max(r.ber_nocomp,1e-5), 'r-s', 'LineWidth', 1.8, 'MarkerSize', 6); hold on;
semilogy(r.snr_dB, max(r.ber_comp,  1e-5), 'b-o', 'LineWidth', 1.8, 'MarkerSize', 6);
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
function plot_curve_set(x, y_mat, labels, y_min)
markers = {'-o','-s','-d','-^','-v','-x','-*'};
colors  = lines(max(size(y_mat,1), 5));
for i = 1:size(y_mat,1)
    semilogy(x, max(y_mat(i,:), y_min), markers{mod(i-1,length(markers))+1}, ...
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
rx_syms = rx_syms(:);
k       = log2(M);
N_sym   = length(rx_syms);
llr     = zeros(N_sym*k, 1);
sigma2  = max(noise_var, 1e-10);
switch M
    case 2
        llr = 2 * real(rx_syms) / sigma2;
    case 4
        llr(1:2:end) = 2*sqrt(2)*real(rx_syms)/sigma2;
        llr(2:2:end) = 2*sqrt(2)*imag(rx_syms)/sigma2;
    otherwise
        c = qammod_const(M);
        c = c / sqrt(mean(abs(c).^2));
        for si = 1:N_sym
            r = rx_syms(si);
            for b = 1:k
                all_idx = 0:M-1;
                bit_b   = floor(all_idx / 2^(k-b));
                mask0 = mod(bit_b,2)==0; mask1 = ~mask0;
                d0 = min(abs(r - c(mask0)).^2);
                d1 = min(abs(r - c(mask1)).^2);
                llr((si-1)*k + b) = (d1 - d0) / sigma2;
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
