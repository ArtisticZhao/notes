function [tm] = ofdmmod(grid, cp_len, nfft, void_DC)
%ofdmmod OFDM 调制
% Args:
%       - grid：  时频网格，Nsc*Nsym，Grid 自动拆成左右两个部分，放在频带正中间
%       - cp_len: cp长度，如果有多种cp长度配置（LTE，NR）输入最少循环节
%       - nfft
%       - void_DC: fft映射时避开直流位置（wlan） default = false；
if nargin<4
    void_DC = false;
end
if void_DC
    void_DC = 1;
else
    void_DC = 0;
end

nsc = size(grid, 1);
n_left = ceil(nsc / 2);
n_right = nsc - n_left;
tm = [];
cps = repmat(cp_len(:), floor(size(grid, 2)/length(cp_len))+1, 1); %复制足够长的CP
for i=1: size(grid, 2)
    cp = cps(i);
    fsym = grid (:, i);
    pre_fft = zeros (nfft, 1);
    left_sc = fsym(1:n_left);
    right_sc = fsym(n_left+1 : end);
    pre_fft(1+void_DC:n_right+void_DC) = right_sc;
    pre_fft(end-n_left+1:end) = left_sc;
    tm_0 = ifft (pre_fft) * sqrt(nfft);
    tm = [tm; tm_0(end-cp+1:end); tm_0];
end