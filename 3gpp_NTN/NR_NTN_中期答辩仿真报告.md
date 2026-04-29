# NR 非地面网络（NTN）下行链路仿真报告

> 本科毕业设计中期答辩 · 基于 3GPP TR 38.811
>
> 作者：周思灵 · 学号：U202242040 · 日期：2026-04-29

---

## 1 课题背景

非地面网络（Non-Terrestrial Network, NTN）是 5G NR 在卫星与空中平台场景下的扩展。3GPP 在 **TR 38.811** 给出了 NTN 信道建模方法，在 **TR 38.821** 给出了 NTN 下 NR 物理层与协议层的解决方案。LEO（低轨）卫星（典型轨道高度 500~1500 km）是中国"星网/商用低轨星座"重点部署的场景，具有下列典型挑战：

| 挑战项 | 典型量级 | 对物理层影响 |
|--------|---------|-------------|
| 高多普勒频移 | LEO-600 @2 GHz 最大 ±40 kHz | 子载波间干扰，信道估计失准 |
| 长传播时延 | 单向 2~10 ms | 超出 NR 标准 TA 范围，需 TA 预补偿 |
| 大尺度路径损耗 | FSPL @2GHz, 600km ≈ 154 dB | 链路预算紧张 |
| 莱斯 vs 瑞利 | K 因子随仰角 -5~20 dB | 低仰角接近瑞利，性能差 |
| 雨衰（高频段） | Ka 频段 25 mm/h 可达 10+ dB | 高频段部署受限 |

本报告在已搭建的 NR-NTN 链路级仿真器上，围绕 **6 个核心问题**做仿真分析，为后续毕业设计完整论文提供数据支撑。

---

## 2 系统模型与仿真工具

### 2.1 链路框图

```
  info bits -> LDPC 编码 -> QAM 调制 -> 资源映射 -> OFDM 调制
                                                         |
                                               NTN 信道（莱斯 + 多普勒
                                                       + 阴影 + 路损
                                                       + 雨衰 + 时延）
                                                         |
             +--- GNSS 频偏补偿 <--- AWGN <--- 残余时延 <--+
             |
         OFDM 解调 -> 信道估计 -> ZF 均衡 -> 软解映射 -> LDPC 解码
```

### 2.2 NR 参数（3GPP TS 38.211）

| 参数 | 值 | 说明 |
|------|-----|------|
| Numerology μ | 1 | 30 kHz SCS（NTN 典型配置，容多普勒） |
| RB 数 | 25 | 约 10 MHz 信道带宽 |
| FFT 点数 | 512 | 采样率 15.36 MHz |
| CP 长度 | 36 采样点 / 2.34 μs | Normal CP |
| DMRS | 符号 3/8/12, 梳状-2 | TS 38.211 Type1 简化 |
| 调制 | QPSK | 答辩默认值，可改 |
| 编码 | NR LDPC BG2, 码率 1/3, SMS 迭代 10 | TS 38.212 标准基图 |

### 2.3 NTN 信道模型（TR 38.811）

- **莱斯衰落**：K(dB) ≈ -12 + 0.3·仰角（TR 38.811 §6.7.2 S 频段城市近似）
- **路径损耗**：FSPL + 大气吸收（≤6 GHz 取 0.1 dB/sinθ）+ 雨衰（ITU-R P.618）
- **阴影衰落**：对数正态，默认 σ=2 dB（答辩对比图已关闭以便可复现）
- **多普勒**：LEO-600 @2 GHz 最大 ±40 kHz，可由外部施加 `extra_fd_hz`
- **时延**：斜距 / 光速；搭配 **TA 预补偿**（TR 38.821 §6.3）

### 2.4 代码架构

```
3gpp_NTN/
├── ntn_sim_main.m                   # 主入口（本报告所有图由此生成）
├── ntn_channel_model.m              # NTN 信道模型（支持雨衰）
├── ntn_rain_attenuation.m           # ITU-R P.618 雨衰计算
├── ntn_ta_precomp.m                 # TA 预补偿（TR 38.821 §6.3）
├── nr_ldpc_encode.m / decode.m      # NR LDPC 编解码（TS 38.212）
├── nr_ldpc_bg_tables.m              # BG1/BG2 基图移位系数表
├── nr_chest_dmrs_standalone.m       # DMRS 信道估计（LS/MMSE/2D_WIENER）
└── results/                         # 仿真结果（*.png, *.mat）
```

### 2.5 运行方法

```matlab
%% 在 MATLAB 命令行：
>> cd D:\study\lab\5G-system-and-link-level-simulator\LTE-上行-下行链路\3gpp_NTN

%% （A）一键跑全部 7 个模块
>> cfg.modules = {'all'};
>> ntn_sim_main

%% （B）只跑指定模块，例如只看雨衰
>> cfg.modules = {'rain','compensation'};
>> ntn_sim_main

%% （C）调整 SNR 范围和帧数（参数可在 ntn_sim_main.m §1 直接改）
>> cfg.snr_dB_vec = 0:2:16;
>> cfg.n_frames   = 120;
>> cfg.modules    = {'kfactor'};
>> ntn_sim_main
```

---

## 3 仿真结果与分析

### 3.1 M1 · 理想 AWGN BER（校准基线）

**目的**：关闭所有信道损伤，仿真 AWGN 下的 QPSK 未编码 BER，与理论 Q 函数对比，**校准链路仿真器的正确性**。

**仿真设置**：

| 项 | 值 |
|----|----|
| 信道 | 纯 AWGN（关 fading/shadow/doppler） |
| 信道估计 | 跳过（直接 H=1） |
| 编码 | 关闭（与未编码 BER 理论对比） |
| SNR | -5 : 3 : 20 dB |

**理论公式**（QPSK, AWGN, Gray）：

$$
\text{BER} = Q\!\left(\sqrt{\frac{E_s/N_0}{1}}\right) = Q(\sqrt{\text{SNR}})
$$

**结果**：`results/M1_ideal_AWGN_BER.png`

**分析**：仿真曲线与理论曲线应在 0.5~1 dB 以内重合（OFDM DMRS 占位引入轻微功率折扣）。此图证明链路仿真器的发射—AWGN—接收流程正确，为后续带衰落曲线的判读提供了基线。

---

### 3.2 M2 · Rician K 因子对 BER / BLER 的影响

**目的**：体现 LEO NTN 小尺度衰落特性。K 因子反映 LOS 分量对散射分量的功率比，随仰角上升而变大（TR 38.811 §6.7.2）。

**仿真设置**：

| 项 | 值 |
|----|----|
| K (dB) | {-∞ (Rayleigh), 0, 7, 15} |
| 多普勒 | 0（排除多普勒干扰） |
| LDPC | 开，1/3 码率 |

**结果**：`results/M2_Kfactor_BER_BLER.png`

**分析要点**：

- 瑞利（无 LOS）曲线下降最慢，BER≈10⁻² 需要 SNR ~15 dB
- K=15 dB 曲线接近 AWGN，LOS 分量主导
- K=0 dB 与 K=7 dB 差距约 3 dB，符合 Rician CDF 的理论预期
- **意义**：仰角 > 30° 时 K ≈ -3 dB 以上，NTN 下行质量显著优于瑞利多径环境

---

### 3.3 M3 · LEO 多普勒频移对 BER 的影响

**目的**：展示未补偿多普勒时 OFDM 子载波正交性被破坏的程度。

**背景**：LEO-600 卫星线速度 ~7.56 km/s。@2 GHz 用户仰角 45° 时多普勒 ≈ 40 kHz (极端值)；过顶附近多普勒变化率最大达 300 Hz/s。

**仿真设置**：

| 项 | 值 |
|----|----|
| 额外频偏 fd | {0, 500, 2 k, 7 k} Hz |
| K | 7 dB（典型 LOS 场景） |

**结果**：`results/M3_doppler_BER.png`

**分析要点**：

- 子载波间隔 30 kHz，fd=500 Hz ≈ 1.7% SCS，几乎无影响
- fd=2 kHz ≈ 6.7% SCS，信道估计已能跟踪，BER 损失 ~1 dB
- fd=7 kHz ≈ 23% SCS，ICI 严重，BER 在高 SNR 出现 floor
- **结论**：NTN 下行必须在 UE 侧做多普勒预补偿（GNSS 反演），见 M7

---

### 3.4 M4 · 残余时延影响与 TA 预补偿效果

**目的**：LEO 下行到 UE 单向时延 2~10 ms，远超 NR 原生 TA 能力，必须做 **UE 侧 TA 预补偿**。此实验展示残余时延（TA 补偿不彻底）对 BER 的危害。

**仿真设置**：

| 项 | 值 |
|----|----|
| 残余时延 | {0, 1, 5, 10, 30} μs |
| TA 补偿曲线 | 启用 ntn_ta_precomp，理论残余 ≈ 0 |
| K | 7 dB |

CP 长度 2.34 μs（30 kHz SCS）是时延门限：

- 残余 ≤ CP：仅相位旋转，可被信道估计吸收
- 残余 > CP：产生 ISI，BER floor
- 残余 ≈ N_fft·Ts (33 μs)：FFT 窗完全错位，BER ≈ 0.5

**结果**：`results/M4_delay_TA.png`

**分析要点**：

- 残余 1 μs（< CP）：几乎与 0 μs 重合，健壮
- 残余 5 μs（约 2× CP）：高 SNR BER floor 出现
- 残余 30 μs：BER 趋近 0.5，链路完全失效
- 启用 TA 预补偿（黑虚线 ★）：完美贴合 0 μs 曲线，验证 TR 38.821 §6.3 方案有效

---

### 3.5 M5 · 仰角与路径损耗

**目的**：理解不同仰角下各损耗分量的贡献，并说明为何 NTN 链路预算以"最恶劣仰角"设计。

**仿真设置**：

| 项 | 值 |
|----|----|
| 仰角 | 10° : 5° : 90° |
| 频段 | S (2 GHz) / Ka (20 GHz) |
| 雨率 | 25 mm/h（中雨） |

**结果**：`results/M5_pathloss_elevation.png`

**分析要点（左图）**：

- FSPL 在 10° 仰角较 90° 增加 ~10 dB（斜距增长因子 1/sinθ）
- 大气吸收 S 频段极小（<1 dB）
- S 频段雨衰 @25 mm/h 仅 ~0.3 dB
- Ka 频段雨衰可达 10+ dB，**雨衰 ≫ 大气吸收**，是 Ka NTN 的主要挑战

**分析要点（右图）**：固定 SNR=10 dB（接收端），BER 随仰角下降（因 K 因子上升，LOS 分量增强），仰角 < 20° 时 BER > 10⁻²。

---

### 3.6 M6 · 雨衰对 BER 的影响（S vs Ka）

**目的**：ITU-R P.618 模型量化不同雨率对 BER 的退化，直观对比 S 与 Ka 频段在雨天的性能差异。

**模型**（TR 38.811 §6.6.5 引用 ITU-R P.618-13）：

$$
A_\text{rain} = k\cdot R^{\alpha}\cdot L_\text{eff}\ \text{[dB]}
$$

其中 k, α 按 P.838-3 与频率关系插值，L_eff 考虑雨层几何路径长度。

**仿真设置**：

| 项 | 值 |
|----|----|
| 频段 | S (2 GHz), Ka (20 GHz) |
| 雨率 | {0, 5, 25, 50} mm/h |
| 仰角 | 45° |

**结果**：`results/M6_rain_BER.png`

**分析要点**：

- S 频段：50 mm/h 暴雨雨衰仅 ~0.3 dB，几乎不影响 BER
- Ka 频段：50 mm/h 雨衰可达 **~20 dB**，BER 曲线右移明显
- **结论**：NR-NTN 高频段部署必须部署站点分集 / 自适应链路预算 / UPC（上行功率控制），低频段（S/L）天生抗雨

---

### 3.7 M7 · 补偿算法对比（GNSS 频偏补偿 & 信道估计）

**目的**：展示在 LEO 大多普勒场景下两类补偿算法的价值：

1. **GNSS 频偏预补偿**：UE 用 GNSS 定位 + 星历反演多普勒，在接收端做前馈频偏纠正
2. **信道估计算法**：从最简 LS_INTERP 到时频联合 2D_WIENER，补偿多普勒导致的信道时变

**仿真设置**：

| 项 | 值 |
|----|----|
| 多普勒 | fd = 7 kHz（高难度场景） |
| GNSS 残差 | 10 Hz（LEO-GNSS 星历典型精度） |
| 估计算法 | {LS_INTERP, MMSE, 2D_WIENER} |

**结果**：`results/M7_compensation.png`

**分析要点（左图，GNSS 补偿前后）**：

- 无补偿：7 kHz 残留多普勒使 BER ~10⁻¹，高 SNR 显示 floor
- GNSS 补偿后（残余 10 Hz）：BER 迅速回到接近 0 Hz 水平，增益 >>5 dB
- **GNSS 反演 + 前馈补偿是 LEO NTN 最基础的必备机制**

**分析要点（右图，估计算法对比）**：

- LS_INTERP 在高 SNR 仍因时域跨符号插值误差出现 floor
- MMSE 用导频局部平滑，性能改善
- 2D_WIENER 联合时频二维滤波，性能最好，**推荐用于 LEO 高多普勒场景**

---

## 4 结论与下阶段工作

### 已完成（v1/v2 + 本次中期）

- [x] 基于 TR 38.811 的 NTN 信道模型（Rician + 多普勒 + 阴影 + 路损 + **雨衰**）
- [x] NR LDPC 编解码（TS 38.212 标准 BG2）
- [x] DMRS 信道估计（LS / MMSE / 2D_WIENER）
- [x] TA 预补偿（TR 38.821 §6.3）
- [x] GNSS 频偏补偿前后对比
- [x] 6 类核心分析图全部完成

### 下阶段（中期 → 终期）

1. **双极化 MIMO**：引入 2×2 DP-MIMO，对比 SISO 性能
2. **随机接入 PRACH**：NTN 大时延下的 PRACH 增强（TR 38.821 §6.6）
3. **HARQ**：NTN 长 RTT 下 HARQ 失效问题（TR 38.821 §6.5）对吞吐量影响
4. **切换**：卫星过顶切换导致的短时链路中断对 BLER 的影响
5. **系统级评估**：多用户/多波束干扰场景（若时间允许）

---

## 5 参考标准与文献

1. 3GPP TR 38.811 V15.4.0, *Study on New Radio (NR) to support non-terrestrial networks*
2. 3GPP TR 38.821 V16.2.0, *Solutions for NR to support non-terrestrial networks (NTN)*
3. 3GPP TS 38.211 / 212 / 213, *NR Physical Channels, Coding, Procedures*
4. ITU-R P.618-13, *Propagation data and prediction methods required for the design of Earth-space telecommunication systems*
5. ITU-R P.838-3, *Specific attenuation model for rain for use in prediction methods*

---

## 附录 A · 参数速查表

```matlab
% === ntn_sim_main.m §1 参数区 ===
cfg.modules            = {'all'};            % 模块开关
cfg.mu                 = 1;                   % 30 kHz SCS
cfg.n_rb               = 25;                  % 10 MHz 带宽
cfg.mod_order          = 4;                   % QPSK
cfg.coding             = true;                % NR LDPC
cfg.ldpc_rate          = 1/3;
cfg.ldpc_iter          = 10;

cfg.altitude_m         = 600e3;               % LEO-600
cfg.fc_hz              = 2e9;                 % S 频段
cfg.elevation_deg      = 45;
cfg.enable_shadow      = false;

cfg.snr_dB_vec         = -5:3:20;             % 9 个 SNR 点
cfg.n_frames           = 60;                  % 快速模式

cfg.K_factor_list_db   = [-Inf, 0, 7, 15];
cfg.doppler_list_hz    = [0, 500, 2000, 7000];
cfg.delay_list_us      = [0, 1, 5, 10, 30];
cfg.rain_rate_list     = [0, 5, 25, 50];
cfg.rain_bands_hz      = [2e9, 20e9];
cfg.chest_methods      = {'LS_INTERP','MMSE','2D_WIENER'};
cfg.pathloss_elev_deg  = 10:5:90;
```

## 附录 B · 交付物清单

| 文件 | 类型 | 说明 |
|------|------|------|
| `ntn_sim_main.m` | 主函数 | 用户入口，7 个仿真模块 |
| `ntn_rain_attenuation.m` | 模块 | ITU-R P.618 雨衰模型 |
| `ntn_channel_model.m` | 修改 | 扩展雨衰支持 |
| `results/M1_ideal_AWGN_BER.png` | 图 | 理想 AWGN 基线 |
| `results/M2_Kfactor_BER_BLER.png` | 图 | K 因子影响 |
| `results/M3_doppler_BER.png` | 图 | 多普勒影响 |
| `results/M4_delay_TA.png` | 图 | 时延与 TA 补偿 |
| `results/M5_pathloss_elevation.png` | 图 | 仰角与路径损耗 |
| `results/M6_rain_BER.png` | 图 | 雨衰影响 |
| `results/M7_compensation.png` | 图 | 补偿算法对比 |
| `results/ntn_sim_main_results.mat` | 数据 | 所有仿真原始数据 |
