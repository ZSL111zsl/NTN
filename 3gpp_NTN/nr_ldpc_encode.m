%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   filename:    nr_ldpc_encode.m
%   description: NR LDPC 编码器，基于 3GPP TS 38.212 §5.3.2。
%                H 矩阵使用 nr_ldpc_bg_tables 提供的标准基图移位系数构造。
%
%   编码流程：
%       1. 根据 K 和 R 选择 BG 和 Zc
%       2. 从 nr_ldpc_bg_tables 获取标准基图，按 Zc 展开 H
%       3. 回代法（back-substitution）求校验比特
%       4. 速率匹配（环形缓冲，打孔前2*Zc系统位）
%
%   input:
%       info_bits   [K×1]   信息比特（0/1）
%       R           scalar  目标码率
%       mode        string  'BG1'|'BG2'|'auto'
%
%   output:
%       coded_bits  [E×1]   编码比特（速率匹配后）
%       enc_params  struct  编码参数
%
%   update note:
%       2026-04-15  created by wangzl
%       2026-04-26  updated by wangzl  (替换为标准BG移位系数表)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [coded_bits, enc_params] = nr_ldpc_encode(info_bits, R, mode)

if nargin < 3, mode = 'auto'; end
info_bits = info_bits(:);
K = length(info_bits);

% ---- BG 选择（TS 38.212 §5.2.2）----
if strcmpi(mode, 'auto')
    if K > 3840 || R > 2/3
        bg = 'BG1';
    else
        bg = 'BG2';
    end
else
    bg = upper(mode);
end

% ---- BG 参数 ----
if strcmpi(bg, 'BG1')
    Kb = 22; mb = 46; Nb = 68;
else
    Kb = 10; mb = 42; Nb = 52;
end

% ---- 选择 Zc（TS 38.212 Table 5.3.2-1）----
Zc_set = sort(unique([2,4,8,16,32,64,128,256,384, ...
                      3,6,12,24,48,96,192,384, ...
                      5,10,20,40,80,160,320, ...
                      7,14,28,56,112,224, ...
                      9,18,36,72,144,288, ...
                      11,22,44,88,176,352, ...
                      13,26,52,104,208]));
Zc = Zc_set(find(Kb * Zc_set >= K, 1, 'first'));
if isempty(Zc)
    Zc = max(Zc_set);
    warning('K=%d 超过 BG%s 支持范围，使用最大 Zc=%d', K, bg(3), Zc);
end

% ---- 尺寸 ----
K_ldpc = Kb * Zc;
N_ldpc = Nb * Zc;
M_ldpc = mb * Zc;

% 填充
if K < K_ldpc
    info_padded = [info_bits; zeros(K_ldpc - K, 1)];
else
    info_padded = info_bits(1:K_ldpc);
end

% ---- 构造标准 H 矩阵 ----
H = nr_ldpc_expand_H(bg, Zc, mb, Kb, Nb);

% ---- 回代法编码 ----
H_s = H(:, 1:K_ldpc);
H_p = H(:, K_ldpc+1:end);
syndrome = mod(H_s * info_padded, 2);
parity   = nr_ldpc_gf2_solve(H_p, syndrome);
codeword = [info_padded; parity];

% ---- 速率匹配 ----
E = round(K / R);
E = max(E, 1);
E = min(E, N_ldpc - 2*Zc);

RM_buffer = codeword(2*Zc+1:end);
if E <= length(RM_buffer)
    coded_bits = RM_buffer(1:E);
else
    reps = ceil(E / length(RM_buffer));
    coded_bits = repmat(RM_buffer, reps, 1);
    coded_bits = coded_bits(1:E);
end

% ---- 输出参数 ----
enc_params.bg     = bg;
enc_params.Zc     = Zc;
enc_params.K      = K;
enc_params.K_ldpc = K_ldpc;
enc_params.N_ldpc = N_ldpc;
enc_params.E      = E;
enc_params.R      = K / E;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   展开 H 矩阵：从基图移位系数表构造完整稀疏 H
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function H = nr_ldpc_expand_H(bg, Zc, mb, Kb, Nb)

[V, ~] = nr_ldpc_bg_tables(bg);

% 按实际 Zc 换算移位值：V_actual = mod(V_base, Zc)
% 分配稀疏矩阵
H = sparse(mb * Zc, Nb * Zc);

I_Zc = speye(Zc);
for row = 1:mb
    for col = 1:Nb
        v = double(V(row, col));
        if v < 0, continue; end         % -1 = 零块，跳过
        shift = mod(v, Zc);
        % 循环置换矩阵 I(shift)：列循环右移 shift 位
        sub = circshift(I_Zc, shift, 2);
        r_idx = (row-1)*Zc+1 : row*Zc;
        c_idx = (col-1)*Zc+1 : col*Zc;
        H(r_idx, c_idx) = sub;
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   GF(2) 线性方程组求解：H_p * p = s (mod 2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = nr_ldpc_gf2_solve(H_p, s)

[m, n] = size(H_p);
A = full(mod([H_p, s], 2));
pivot = zeros(1, n);

col = 1;
for row = 1:m
    if col > n, break; end
    found = false;
    for r = row:m
        if A(r, col) == 1
            if r ~= row
                A([row,r], :) = A([r,row], :);
            end
            found = true;
            break;
        end
    end
    if ~found
        col = col + 1;
        continue;
    end
    pivot(col) = row;
    for r = 1:m
        if r ~= row && A(r, col) == 1
            A(r, :) = mod(A(r, :) + A(row, :), 2);
        end
    end
    col = col + 1;
end

p = zeros(n, 1);
for c = 1:n
    if pivot(c) > 0
        p(c) = A(pivot(c), end);
    end
end

end
