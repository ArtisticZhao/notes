% @file main_MRC_rayleigh
% @brief MRC receiver under rayleigh channel.
% @author yuhao.zhao
% @details New structure of MIMO simulation. To get the experiment the MIMO
%          channel effect and MIMO techniques
close all

%% Configure of simulation
cfg = struct;
cfg.n_test = 5; % simulate 5 times in each SNR.
cfg.snr = 0:2:20; % dB 0:2:20

n_rx_ant = [1 2 4];  % define rx diversity
%% Generate Tx signal
% This scripts purpose a OFDM signal with single preamble to estimate
% channel, and assume the whole simulation under the time-const channel
% model.
generate_demo_OFDM2;

ber = cell(length(n_rx_ant), 1);
for n = 1 : length(n_rx_ant)
    n_rx = n_rx_ant(n);

    %% define process function
    % -------------- 1. multipath propogate
    multi_path = @(tx) repmat(tx, 1, n_rx);
    % -------------- 2. rayleigh channel
    pathDelays = [0 10 20 30 ]*1e-9; % in seconds|
    avgPathGains0 = [0 -3 -10 -8 ];  % dB
    fD = 0; % Max Doppler shift in Hz
    rayleigh_channels = cell(n_rx, 1);
    for i=1:n_rx
        rayleigh_channels{i} = comm.RayleighChannel('SampleRate', fs, ...
                                            'PathDelays',pathDelays, ...
                            'AveragePathGains',avgPathGains0, ...
                            'MaximumDopplerShift',fD, ...
                            'Visualization', 'Off'); %
    end
    channel_func = @(tx) channel(tx, rayleigh_channels);
    post_demod = @(rx) rx; % bypass
    % -------------- 3. demodulation
    demod_func = @(rx) MRC_demodulate(rx, cp_len, nsym, nfft, nsc_all, preamble);
    % -------------- 4. measure BER
    measure_func = @(freq_equ) measure_BER(freq_equ, idx_data_freq, idx_data_sym, M_data, x_data_bin);
    %% Simulation!
    % snr_playground(cfg, tx, pre_process, channel, demodulate, measurement)
    ber{n} = snr_playground(cfg, tm, multi_path, channel_func, post_demod, demod_func, measure_func);
end

%% plot 
figure;
lb = [];
for i=1:length(n_rx_ant)
    ber_n = ber{i};
    ber_p = mean (ber_n, 2);

    semilogy(cfg.snr, ber_p,'-*'); hold on;
    lb = [lb sprintf("N_{rx} = %d", n_rx_ant(i))];
end
legend(lb)
title("Rayleigh信道的OFDM-MRC接收机的接收分集性能")
grid on;
xlabel('SNR (dB)', 'FontSize',12);
ylabel('BER', 'FontSize', 12);
%EOF


%% function defination
function [rx] = channel(tx, channels)
    n_path = size(tx, 2);
    rx = zeros(size(tx), 'like', tx);
    for i=1:n_path
        ch = channels{i};
        rx(:, 1) = ch(tx(:, 1));
    end
end

function [freq_equ] = MRC_demodulate(rx, cp_len, nsym, nfft, nsc_all, preamble)
    nrx = size(rx,2);
    if nrx == 1 % SISO
        % ------ OFDM demod
        freq = azcomm.ofdmdemod (rx, repmat(cp_len, nsym, 1), nfft, nsc_all);
        % ------ channel est
        rx_preamble = freq(:, 1);
        Hest = rx_preamble ./ preamble(:, 1);
        Hest = repmat (Hest, 1, nsym);
        freq = freq./ Hest;
        freq_equ = freq;
    else % MRC
        all_freq = zeros(nsc_all, nsym, nrx);
        for r=1:nrx
            %------ OFDM demod
            freq = azcomm.ofdmdemod(rx(:, r), repmat(cp_len, nsym, 1), nfft, nsc_all);
            % ------ channel est
            rx_preamble = freq(:, 1);
            Hest = rx_preamble ./ preamble(:, 1);
            Hest = repmat(Hest, 1, nsym);
            freq = conj(Hest).* freq; %！！这里利用估计的衰落信道的估计值，热行MRC
            all_freq(:,:,r) = freq;
        end
        % MRC with equ first, than you must equ AGAIN!!
        MRC_freq = sum(all_freq, 3);
        rx_preamble = MRC_freq(:, 1);
        Hest = rx_preamble ./ preamble(:, 1);
        Hest = repmat (Hest, 1, nsym);
        freq_equ = MRC_freq./ Hest;
    end
end