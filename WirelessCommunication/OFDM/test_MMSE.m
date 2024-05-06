
close all

generate_demo_OFDM2;
cSNR = 20;

% channel .....
rx = awgn(tm, cSNR, 'measured');

% demod
freq = azcomm.ofdmdemod(rx, repmat(cp_len, nsym, 1), nfft, nsc_all);

% channel est and equ
%% 1. LS
H = channel_est(freq, grids, (1:nsc_all), idx_preamble_sym);

%% 2. MMSE matlab
% [H, n] = nrChannelEstimate(freq, (1:nsc_all*2), freq(1:nsc_all*2));
equ = nrEqualizeMMSE(freq(:), H(:), 0);
equ = reshape(equ,[],14);

figure;
plot_scatterIQ(equ)

function [H] = channel_est(rx_grid, raw_grid, rs_k, rs_sym)
    N_sym = size(rx_grid,2);
    H = rx_grid(rs_k, rs_sym) ./ raw_grid(rs_k, rs_sym);
    % Expand H to whole grid, using symbol direction average the H
    H = repmat(mean(H, 2), 1, N_sym); 
end