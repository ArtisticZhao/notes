
close all

generate_demo_OFDM2;

%% rayleigh channel model
pathDelays = [0 10 20 30 ]*1e-9; % in seconds
avgPathGains0 = [0 -3 -10 -8 ]; % dB

fD = 0; % Max Doppler shift in Hz
rayleigh = comm.RayleighChannel('SampleRate',fs, ...
        'PathDelays' ,pathDelays,...
        'AveragePathGains',avgPathGains0,...
        'MaximumDopplerShift',fD, ...
        'Visualization', 'Off');

measure_func = @(freq_equ) measure_BER(freq_equ, idx_data_freq, idx_data_sym, M_data, x_data_bin);

n_test = 20; % simulate 5 times in each SNR.
snr = 0:0.5:20; % dB 0:2:20
snr = 0:2:20;
snr = 20;
ber_LS = zeros(length(snr), n_test);
ber_MMSE = zeros(length(snr), n_test);
ber_LS_MMSE = zeros(length(snr), n_test);
ber_LS_Book = zeros(length(snr), n_test);
for i=1:length(snr)
    cSNR = snr(i);
    fprintf('Current simulate at SNR=%.2f\n', cSNR);
    for t=1:n_test
        % channel .....
        rx = rayleigh(tm);
        rx = awgn(rx, cSNR, 'measured');

        % demod
        freq = azcomm.ofdmdemod(rx, repmat(cp_len, nsym, 1), nfft, nsc_all);

        % channel est and equ
        %% 1. LS
        H_LS = channel_est_LS(freq, grids, (1:nsc_all), idx_preamble_sym);
        equ_LS = method_ls(freq, H_LS);
        ber_LS(i, t) = measure_func(equ_LS);

        %% 2. MMSE matlab
        [H_Mat, n] = nrChannelEstimate(freq, (1:nsc_all*2), grids(1:nsc_all*2));
        H_Mat = mean_expand_H(H_Mat(:, idx_preamble_sym), size(freq, 2));
        equ_Mat = nrEqualizeMMSE(freq(:), H_Mat(:), n);
        equ_Mat = reshape(equ_Mat,[],14);
        ber_MMSE(i, t) = measure_func(equ_Mat);

        %% 3. LS + matlab equalizer
        % H = channel_est_LS(freq, grids, (1:nsc_all), idx_preamble_sym);
        equ_LS_mat = nrEqualizeMMSE(freq(:), H_LS(:), 0);
        equ_LS_mat = reshape(equ_LS_mat,[],14);
        ber_LS_MMSE(i, t) = measure_func(equ_LS_mat);

        %% 4. book MMSE channel estimate
        % N_sym = size(freq,2);
        % H_est1 = MMSE_CE(freq(:,1).',grids(:,1).',1:600,600, 1, cSNR);
        % H_est2 = MMSE_CE(freq(:,2).',grids(:,2).',1:600,600, 1, cSNR);
        % H_book = [H_est1.' H_est2.'];
        % H_book = mean(H_book, 2);
        % H_book = repmat(H_book, 1, N_sym);
        H_book = channel_est_MMSE(freq, grids, (1:nsc_all), idx_preamble_sym, cSNR);
        equ_book = method_ls(freq, H_book);
        ber_LS_Book(i, t) = measure_func(equ_book);
    end
end

figure;
semilogy(snr, mean(ber_LS, 2),'-*'); hold on;
semilogy(snr, mean(ber_MMSE, 2),'-*'); hold on;
semilogy(snr, mean(ber_LS_MMSE, 2),'o'); hold on;
semilogy(snr, mean(ber_LS_Book, 2),'--*'); hold on;
grid on;
legend('LS', 'MMSE matlab', 'LS + matlab MMSE', 'Book')
xlabel('SNR (dB)', 'FontSize',12);
ylabel('BER', 'FontSize', 12);

% Hest = repmat(Hest, 1, size(freq,2));

% equ = nrEqualizeMMSE(freq, Hest, 0);

% equ = Hest' * ((Hest*Hest') \ freq(:, 1));

% ber{n} = snr_playground(cfg, tm, bypass, bypass, bypass, demod_func, measure_func);
%% compare EVM level
evmer = comm.EVM;
evm_ls = evmer(equ_LS, grids);

release(evmer);
evm_mat = evmer(equ_Mat, grids);

release(evmer);
evm_book = evmer(equ_book, grids);

figure;
subplot(131)
plot_scatterIQ(equ_LS)
title(sprintf("LS EVM: %.2f", mean(evm_ls)))

subplot(132)
plot_scatterIQ(equ_Mat)
title(sprintf("mat EVM: %.2f", mean(evm_mat)))

subplot(133)
plot_scatterIQ(equ_book)
title(sprintf("book EVM: %.2f", mean(evm_book)))


function [H] = channel_est_LS(rx_grid, raw_grid, rs_k, rs_sym)
    N_sym = size(rx_grid,2);
    H = rx_grid(rs_k, rs_sym) ./ raw_grid(rs_k, rs_sym);
    % Expand H to whole grid, using symbol direction average the H
    H = mean_expand_H(H, N_sym);
end

function [H] = channel_est_MMSE(rx_grid, raw_grid, rs_k, rs_sym, cSNR)
    N_sym = size(rx_grid,2);
    Nsc = size(rx_grid,1);
    nsym_preamble = length(rs_sym);
    H_rs = zeros(Nsc, nsym_preamble);

    for i=1:nsym_preamble
        l_rs = rs_sym(i);
        H_est1 = MMSE_CE(rx_grid(:,l_rs).',raw_grid(:,l_rs).', rs_k, Nsc, 1, cSNR);
        H_rs(:, i) = H_est1;
    end

    H = mean_expand_H(H_rs, N_sym);
end

function [H] = mean_expand_H(H, N_sym)
    H = mean(H, 2);
    H = movmean(H, 3);
    H = repmat(H, 1, N_sym); 
end

function [equ] = method_ls(rx_grid, H_grid)
    % for LS equalization
    equ = rx_grid ./ H_grid;
end
