function [grid] = ofdmdemod(tm, cp_len, nfft, nsc, void_DC, cp_fraction)
%ofdmmod OFDM 调制
%   Args:
%       - grid: 时频网格，Nsc*Nsym，Grid 自动拆成左右两个部分，放在频带正中间
%       - cp_len: cp长度， 按照输入符号的个数，要为每个符号提供一个cp
%       - nfft
%       - nsc: 使用的有效子载波个数
%       - cp_fraction: 防止混叠技术，默认0.55
if nargin<6
    cp_fraction = 0.55;
end
if nargin<5
    void_DC = false;
end

if void_DC
    void_DC = 1;
else
    void_DC = 0;
end

n_left = ceil(nsc / 2);
n_right = nsc - n_left;

idx_left = nfft-n_left+1:nfft;
idx_right = 1+void_DC:n_right+void_DC;

nsym = length(cp_len);    

start = 1;
grid = zeros(nsc, nsym);
for i=1:nsym
    cp_offset = floor(cp_len(i) * cp_fraction);
    this_sym_tm = tm(start: start + cp_len(i) + nfft - 1);

    phase_correction = conj(exp(-1j*2*pi*(cp_len(i)-cp_offset)/nfft * (0:nfft-1).'));
    freq_sym = fft(this_sym_tm(1+cp_offset:cp_offset+nfft)) .* phase_correction ./ sqrt(nfft);
    
    left_sc = freq_sym(idx_left);
    right_sc = freq_sym(idx_right);
    grid(:, i) = [left_sc; right_sc];

    start = start + cp_len(i) + nfft;
end