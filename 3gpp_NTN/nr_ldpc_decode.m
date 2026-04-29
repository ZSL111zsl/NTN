%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   filename:    nr_ldpc_decode.m
%   description: NR LDPC 解码器（Log-BP / Min-Sum），H 矩阵使用与编码侧完全
%                相同的标准基图（nr_ldpc_bg_tables + nr_ldpc_expand_H）。
%
%   算法：
%       Sum-Product (SP)、Min-Sum (MS)、Scaled Min-Sum (SMS, alpha=0.75)
%
%   input:
%       llr_in       [E×1 或 N_ldpc×1]  接收 LLR（正值倾向0，负值倾向1）
%       enc_params   struct              来自 nr_ldpc_encode 的编码参数
%       max_iter     scalar              最大 BP 迭代次数（默认20）
%       decoder_type string              'SP'|'MS'|'SMS'（默认'SMS'）
%
%   output:
%       decoded_bits [K×1]   解码信息比特
%       crc_pass     logical  是否校验通过（H*c=0）
%       n_iter       scalar   实际迭代次数
%
%   update note:
%       2026-04-15  created by wangzl
%       2026-04-26  updated by wangzl  (使用标准BG移位系数表，去除随机H矩阵)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [decoded_bits, crc_pass, n_iter] = nr_ldpc_decode(llr_in, enc_params, max_iter, decoder_type)

if nargin < 3 || isempty(max_iter),     max_iter     = 20;    end
if nargin < 4 || isempty(decoder_type), decoder_type = 'SMS'; end

llr_in = llr_in(:);

bg     = enc_params.bg;
Zc     = enc_params.Zc;
K      = enc_params.K;
K_ldpc = enc_params.K_ldpc;
N_ldpc = enc_params.N_ldpc;
E      = enc_params.E;

alpha_ms = 0.75;

% ---- BG 参数 ----
if strcmpi(bg, 'BG1')
    Kb = 22; mb = 46; Nb = 68;
else
    Kb = 10; mb = 42; Nb = 52;
end

% ---- 速率解匹配 ----
if length(llr_in) == N_ldpc
    llr_full = llr_in;
else
    llr_full = zeros(N_ldpc, 1);
    llr_full(1:2*Zc) = 0;   % 打孔位无先验
    buffer_len = N_ldpc - 2*Zc;
    if E <= buffer_len
        llr_full(2*Zc+1 : 2*Zc+E) = llr_in;
    else
        temp = zeros(buffer_len, 1);
        for i = 1:E
            idx = mod(i-1, buffer_len) + 1;
            temp(idx) = temp(idx) + llr_in(i);
        end
        llr_full(2*Zc+1:end) = temp;
    end
end

% ---- 构造标准 H 矩阵 ----
H = nr_ldpc_expand_H(bg, Zc, mb, Kb, Nb);

% ---- BP 初始化 ----
[M, N] = size(H);
[row_idx, col_idx] = find(H);

L_ch    = llr_full;
L_cn2vn = zeros(nnz(H), 1);

vn_neighbors = cell(N, 1);
for v = 1:N
    vn_neighbors{v} = find(H(:, v));
end
cn_neighbors = cell(M, 1);
for c = 1:M
    cn_neighbors{c} = find(H(c, :));
end

edge_idx = zeros(M, N);
edge_idx(sub2ind([M,N], row_idx, col_idx)) = 1:length(row_idx);

crc_pass = false;
n_iter   = 0;

% ---- BP 迭代 ----
for iter = 1:max_iter
    n_iter = iter;

    % VN→CN 消息
    L_vn2cn = zeros(length(row_idx), 1);
    for v = 1:N
        nbrs = vn_neighbors{v};
        if isempty(nbrs), continue; end
        e_all = edge_idx(nbrs, v);
        sum_all = L_ch(v) + sum(L_cn2vn(e_all));
        for nb_c = nbrs'
            e = edge_idx(nb_c, v);
            L_vn2cn(e) = sum_all - L_cn2vn(e);
        end
    end

    % CN→VN 消息
    for c = 1:M
        nbrs = cn_neighbors{c};
        if isempty(nbrs), continue; end
        msgs = L_vn2cn(edge_idx(c, nbrs));

        switch upper(decoder_type)
            case 'SP'
                prod_tanh = prod(tanh(msgs / 2));
                prod_tanh = max(-1+1e-10, min(1-1e-10, prod_tanh));
                for i = 1:length(nbrs)
                    ti = tanh(msgs(i)/2);
                    if abs(ti) < 1-1e-10
                        pt_excl = prod_tanh / ti;
                    else
                        pt_excl = prod(tanh([msgs(1:i-1); msgs(i+1:end)]/2));
                    end
                    pt_excl = max(-1+1e-10, min(1-1e-10, pt_excl));
                    L_cn2vn(edge_idx(c, nbrs(i))) = 2 * atanh(pt_excl);
                end

            case {'MS','SMS'}
                alpha = 1.0;
                if strcmpi(decoder_type,'SMS'), alpha = alpha_ms; end
                abs_msgs  = abs(msgs);
                sign_msgs = sign(msgs);
                sign_msgs(sign_msgs == 0) = 1;
                prod_sign = prod(sign_msgs);
                [sorted_abs, sort_idx] = sort(abs_msgs);
                min1 = sorted_abs(1);
                min2 = sorted_abs(2);
                for i = 1:length(nbrs)
                    sgn_excl = prod_sign * sign_msgs(i);
                    min_excl = min1 + (sort_idx(1)==i) * (min2 - min1);
                    L_cn2vn(edge_idx(c, nbrs(i))) = alpha * sgn_excl * min_excl;
                end
        end
    end

    % 后验 LLR 与硬判决
    L_post = L_ch;
    for v = 1:N
        nbrs = vn_neighbors{v};
        L_post(v) = L_ch(v) + sum(L_cn2vn(edge_idx(nbrs, v)));
    end
    c_hat = double(L_post < 0);

    % 校验收敛
    if all(mod(H * c_hat, 2) == 0)
        crc_pass = true;
        break;
    end
end

decoded_bits = c_hat(1:K);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   展开 H 矩阵（与 nr_ldpc_encode 共用逻辑，独立副本避免跨文件调用）
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function H = nr_ldpc_expand_H(bg, Zc, mb, Kb, Nb)

[V, ~] = nr_ldpc_bg_tables(bg);
H = sparse(mb * Zc, Nb * Zc);
I_Zc = speye(Zc);

for row = 1:mb
    for col = 1:Nb
        v = double(V(row, col));
        if v < 0, continue; end
        shift = mod(v, Zc);
        sub = circshift(I_Zc, shift, 2);
        r_idx = (row-1)*Zc+1 : row*Zc;
        c_idx = (col-1)*Zc+1 : col*Zc;
        H(r_idx, c_idx) = sub;
    end
end

end
