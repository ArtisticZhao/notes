function [Hest] = channel_estimate(rxgrid, refgrid, dmrs_k_idx, dmrs_slot_idx, expand, interp_method, use_cir, nfft, min_cp)
%CHANNEL_ESTIMATE 信道估计，给我有信号的部分,1 slot
%   Args:
%     - rxgrid: 一个slot内的接收时频网格
%     - refgrid: 参考网格
%     - dmrs_k_idx: dmrs频率上的索引位置
%     - dmrs_slot_idx: 在当前slot内DMRS的位置
%     - expand: 每个符号copy那个DMRS的估计结果；如导频在[2 8 12], Nsym=14;
%              expand = [2 -1 2 2 2 2 8 -1 8 8 12 -1 12 12], 其中-1表示不copy
%     - interp_method: 插值方法，可选'linear', 'spline', 'polar_linear', 'polar_spline'
%     - use_cir: 是否使用CIR噪声消除技术
%   Returns:
%     - Hest: 估计完的全slot 的信道估计， 均衡用 equ = rxgrid .* Hest;

%% 入参预处理
if nargin < 7
    use_cir = false; % 如果不提供use_cir 参数， 则不启用
end
dmrs_k_idx = dmrs_k_idx(:);
% ------------------------- 信道估计， 先求出包含RS的符号的H

if ~use_cir % 普通估计
    Hest = zeros(size(rxgrid));
    % 有导频的地方，采用最小二乘估计
    Hest(dmrs_k_idx, dmrs_slot_idx) = rxgrid(dmrs_k_idx, dmrs_slot_idx) ./ refgrid(dmrs_k_idx, dmrs_slot_idx);
    % DMRS所在符号内没有导频的地方插值
    idx = 1:size(rxgrid, 1);
    for i=1:length(dmrs_slot_idx)
        c_H = Hest(dmrs_k_idx, dmrs_slot_idx(i));
        Hest(:, dmrs_slot_idx(i)) = interp1_cart_polar(dmrs_k_idx, c_H, idx, interp_method);
    end
else % LS估计后， 采用CIR技术消除额外噪声
    ERB = 2; % 在执行CIR防止边缘起飞，多2个RB一边
    ERE = ERB * 12;
    % 0. 计算冲击相应窗
    % Calculate minimum cyclic prefix length in terms of a DFT of size K+eK
    cp = floor(min_cp / nfft * (size(rxgrid, 1)+ERE*2));
    
    % Create time-domain windowing function for CIR denoising
    w = raised_cosine_window(cp*2,cp);
    w = [w; zeros([size(rxgrid, 1)+ERE*2-length(w) 1])];
    w = circshift(w,-cp-floor(cp/2));


    HestExtend = zeros(ERE + size(rxgrid, 1) + ERE, size(rxgrid, 2));
    dmrs_k_idx_offset = dmrs_k_idx + ERE;

    % 有导频的地方，采用最小二乘估计, 放在扩展后的位置
    HestExtend(dmrs_k_idx_offset, dmrs_slot_idx) = rxgrid(dmrs_k_idx, dmrs_slot_idx) ./ refgrid(dmrs_k_idx, dmrs_slot_idx);
    % 虚拟的导频，每6个SC一个
    Extend_DMRS_idx_left = (1:6:ERE).';
    Extend_DMRS_idx_right = (ERE + size(rxgrid, 1) + 1 : 6 : size(rxgrid, 1) + ERE).';
    for i=1:length(dmrs_slot_idx) % 按符号计算CIR降噪的信道H
        % ------ 1.首先将H的结果向左右各扩展2个RB
        % 1.1 left 虚拟RS
        left_H_idx = dmrs_k_idx_offset(1:6);
        left_H = HestExtend(left_H_idx, dmrs_slot_idx(i));
        VP_left = polar_interp_virtual_pilots(left_H, left_H_idx, Extend_DMRS_idx_left);
        HestExtend(Extend_DMRS_idx_left, dmrs_slot_idx(i)) = VP_left;
        % 1.2 right 虚拟RS
        right_H_idx = dmrs_k_idx_offset(end-5:end);
        right_H = HestExtend(right_H_idx, dmrs_slot_idx(i));
        VP_right = polar_interp_virtual_pilots(right_H, right_H_idx, Extend_DMRS_idx_right);
        HestExtend(Extend_DMRS_idx_right, dmrs_slot_idx(i)) = VP_right;

        % ------ 2. 插值！
        all_rs_idx = [Extend_DMRS_idx_left; dmrs_k_idx_offset ;Extend_DMRS_idx_right];
        H_sym = HestExtend(all_rs_idx,  dmrs_slot_idx(i));
        all_idx = (1:size(HestExtend, 1)).';
        H = interp1_cart_polar(all_rs_idx, H_sym, all_idx, interp_method);

        % ------ 3. 信道CIR计算
        h = ifft(H);
        % Apply time domain windowing function to denoise CIR
        h = h .* w;
        % back to H
        H = fft(h);

        HestExtend(:,  dmrs_slot_idx(i)) = H;
    end
    Hest = HestExtend(ERE+1:end-ERE, :);

end


% ------------------------- expand to slot
Nsym = size(rxgrid, 2); % 认为输入是一个slot的数据
for i=1:Nsym
    if expand(i) == -1
        continue
    end
    Hest(:, i) = Hest(:, expand(i));
end

end

% interp1, 支持按IQ插值或者按极坐标插值
function [Yq] = interp1_cart_polar(X, V, Xq, method)
    % method 可选 'linear', 'spline', 'polar_linear', 'polar_spline'
    is_polar = contains(method, 'polar');
    if  contains(method, 'linear') 
        interp_method = 'linear';
    else
        interp_method = 'spline';
    end
    if is_polar
        V_mag = abs(V);
        V_ang = unwrap(angle(V));
        Y_mag = interp1(X, V_mag, Xq, interp_method, 'extrap');
        Y_ang = interp1(X, V_ang, Xq, interp_method, 'extrap');
        Yq = Y_mag .* exp(1j * Y_ang);
    else
        Yq = interp1(X, V, Xq, interp_method, 'extrap');
    end
end


function [Vps] = polar_interp_virtual_pilots(H, H_idx, Extend_idx)
% 按模和辐角插值， 并且向带外插值，采用一次拟合的方法。
    % 向带外扩展的时候，采用多项式拟合，才能利用上全部的点
    left_H_mag = abs(H);
    left_H_ang = unwrap(angle(H));
    % 1.1 mag
    p = polyfit(H_idx, left_H_mag, 1); % 一次函数拟合
    virtual_H_mag = Extend_idx * p(1) + p(2);
    % 1.2 angle
    p = polyfit(H_idx, left_H_ang, 1); % 一次函数拟合
    virtual_H_ang = Extend_idx * p(1) + p(2);
    
    Vps =  virtual_H_mag .* exp(1j*virtual_H_ang);
end

% Raised cosine window creation; creates a window function of length n+w
% with raised cosine transitions on the first and last 'w' samples.
function p = raised_cosine_window(n,w)
    
    p = 0.5*(1-sin(pi*(w+1-2*(1:w).')/(2*w)));
    p = [p; ones([n-w 1]); flipud(p)];
    
end