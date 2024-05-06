% @file generate_demo_OFDM2.m
% @brief 生成一个OFDM信号，用于测试各种信道和解调性能；
% @details In this version, we introduced a preamble in full subcarriers
%          dedicate for channel estimate. In the next simulation assume
%          that the channel is time constantly.
% @author zhao.yuhao

%% OFDM
scs = 15e3;
nfft = 1024;
fs = scs * nfft;
cp_len = 256; % length of cyclic prefix
nsym = 14; % 生成OFDM符号数
M_data = 4; % modulation order of data
M_preamble = 2; % modulation order of preamble
nsc_all = 600; % all activate subcarrier|
nsc_data = nsc_all - 0; % data sc;
nsym_preamble = 2; % 两个preamble給发射分集的almouti码用，在接收分集时只用第一个符号的preamble
nsym_data = nsym - nsym_preamble;
%
idx_data_freq = (1:nsc_all).';
idx_data_sym = (nsym_preamble+1:nsym).';
idx_preamble_sym = 1:nsym_preamble;
% ------ generate data in frequency domain
grids = zeros(nsc_all, nsym);

% data all
x_data = randi([0 M_data-1], nsc_data*nsym_data, 1);
y = qammod(x_data, M_data, 'UnitAveragePower', true);
x_data_bin = double(dec2bin(x_data, M_data)) - double('0');
grids(idx_data_freq, idx_data_sym) = reshape(y, nsc_data, nsym_data);
% preamble
x_preamble = randi([0 M_preamble-1], nsc_all*nsym_preamble, 1);

preamble = qammod(x_preamble, M_preamble, 'UnitAveragePower', true);
preamble = reshape (preamble, nsc_all, nsym_preamble);
grids (:, idx_preamble_sym) = preamble;
%
figure;
plot_scatterIQ(grids); title("Tx all constellation");
% ------ OFDM mod
tm = azcomm.ofdmmod(grids, cp_len, nfft);