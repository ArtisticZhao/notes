% @file main_Alamouti_rayleigh.m
% @brief 仿真发射分集，采用Alamouti编码方式
% @author yuhao.zhao 
close all

%% Configure of simulation
cfg = struct;
cfg.n_test = 5; % simulate 5 times in each SNR.
cfg.snr = 0:2:20; % dB

n_tx_ant = [1 2];

%% generation
% For alamouti coding, only use the frequency part, and make alamouti 
generate_demo_OFDM2;

ber = cell(length(n_tx_ant), 1);

% -------------- 2. rayleigh channel in const
pathDelays = [0 10 20 30 ]*1e-9; % in seconds
avgPathGains0 = [0 -3 -10 -8 ]; % dB

fD = 0; % Max Doppler shift in Hz
rayleigh_channels = cell(max(n_tx_ant), 1);
for i=1:max(n_tx_ant)
    rayleigh_channels{i} = comm.RayleighChannel('SampleRate',fs, ...
        'PathDelays' ,pathDelays,...
        'AveragePathGains',avgPathGains0,...
        'MaximumDopplerShift',fD, ...
        'Visualization', 'Off'); %
end

for n = 1 : length(n_tx_ant)
    n_tx = n_tx_ant(n);

    % ------ alamouti coding
    grid_alamouti = Alamouti_coding(grids, n_tx);
    % ------ OFDM mod
    tm = zeros (nsym* (nfft+cp_len), n_tx);
    for i=1:n_tx
        tm(:, i) = azcomm.ofdmmod (grid_alamouti(:, :, i), cp_len, nfft);
    end
    
    %% define process function
    % -------------- 1. multipath propogate
    % multi_path = @(tx) tx ./ sqrt(n_tx); % 终T发射总功率与单T拉齐
    multi_path = @(tx) tx;    % 多T发射总功率与单T拉齐
    
    rayleigh_channels2 = cell(n_tx, 1);
    for i=1:n_tx
        rayleigh_channels2{i} = rayleigh_channels{i};
    end
    channel_func = @(tx) channel(tx, rayleigh_channels2);

    post_demod = @(rx) rx;
    % -------------- 3. alamouti decode
    demod_func = @(rx) demodulate(rx, n_tx, cp_len, nsym, nfft, nsc_all, preamble, idx_preamble_sym, idx_data_sym);
    % -------------- 4. measure BER
    measure_func = @(freq_equ) measure_BER(freq_equ, idx_data_freq, idx_data_sym, M_data, x_data_bin);
    
    %% Simulation!!
    ber{n} = snr_playground (cfg, tm, multi_path, channel_func, post_demod, demod_func, measure_func);
end

%% plot
figure;
lb = [];
for i=1: length (n_tx_ant)
    ber_n = ber{i};
    ber_p = mean(ber_n, 2);
    semilogy(cfg.snr, ber_p,'-*'); hold on;
    lb = [lb sprintf("N_{tx} = %d", n_tx_ant(i))];
end
legend (lb)
title("Rayleigh信道的OFDM-A1amouti 的发射分集性能")
grid on ;
xlabel('SNR (dB)', 'FontSize',12);
ylabel('BER', 'FontSize', 12);
%EOF

%% functions
function [grid_alamouti] = Alamouti_coding(grid, n_tx)
    if n_tx == 1
        grid_alamouti = grid;
        return;
    end
    nsym = size(grid, 2);
    k = size(grid, 1);
    if mod(nsym, n_tx) ~=0
        error ("alamouti code must has multiple n_tx");
    end
    grid_alamouti = zeros(k, nsym, n_tx);
    for i=1:nsym/n_tx
        seg_g = grid(:, (i-1) * n_tx + (1:n_tx));
        if n_tx == 2
            %     t0 t1
            % AO: s0 s1*
            % A1: s1 -s0*
            grid_alamouti(:, (i-1) * n_tx + (1:n_tx), 1) = [seg_g(:, 1) conj(seg_g(:, 2))];
            grid_alamouti(:, (i-1) * n_tx + (1:n_tx), 2) = [seg_g(:, 2) -conj(seg_g(:, 1))];
        else
            error("Ntx=%d, alamouti code not supported!", n_tx);
        end
    end
end

function [H] = Alamouti_channel_estimate (grid, preamble, n_tx)
    if n_tx == 2
        y0 = grid(:, 1);
        y1 = grid(:, 2);
        c0 = preamble (:, 1);
        c1 = preamble (:, 2);

        H0 = conj(c0) .* y0 + c1 .* y1;
        H1 = conj(c1) .* y0 - c0 .* y1;

        factor = c0 .* conj(c0) + c1 .* conj(c1);
        H = [H0 ./ factor, H1 ./ factor];
    else
        error("Ntx=%d, alamouti channel estimate not supported!", n_tx);
    end
end

function [decode] = Alamouti_decode(grid, H, n_tx)
    k = size(grid, 1);
    nsym = size(grid, 2);
    decode = zeros (k, nsym);
    for i=1:nsym/n_tx
        rx_seg = grid(:, (i-1) * n_tx +(1:n_tx));
        if n_tx == 2
            r0 = rx_seg(:, 1);
            r1 = rx_seg(:, 2);
            H0 = H(:, 1);
            H1 = H(:, 2);

            X0 = conj(H0) .* r0 - H1 .* conj(r1);
            X1 = conj(H1) .* r0 + H0 .* conj(r1);
            factor = H0 .* conj(H0) + H1 .* conj(H1);
            decode (:, (i-1) * n_tx + (1:n_tx)) = [X0 ./ factor, X1 ./ factor];
        else
            error("Ntx-%d, alamouti channel estimate not supported!", n_tx);
        end
    end
end

function [freq_equ] = demodulate(rx, n_tx_ant, cp_len, nsym, nfft, nsc_all, preamble, idx_preamble_sym, idx_data_sym)
    if n_tx_ant == 1 % SISO % ------ OFDM demod
        freq = azcomm.ofdmdemod(rx, repmat (cp_len, nsym, 1), nfft, nsc_all);
        % ------ channel est
        rx_preamble = freq(:, 1);
        Hest = rx_preamble ./ preamble(:, 1);
        Hest = repmat (Hest, 1, nsym);
        freq = freq./ Hest;
        freq_equ = freq;
    else % Alamouti demod
        % --------------- 1. channel estimate|
        freq = azcomm.ofdmdemod(rx, repmat (cp_len, nsym, 1), nfft, nsc_all);
        rx_dmrs = freq(:,idx_preamble_sym);
        H = Alamouti_channel_estimate(rx_dmrs, preamble, n_tx_ant);
        % --------------- 2. demod
        rx_data = freq(:, idx_data_sym);
        data = Alamouti_decode (rx_data, H, n_tx_ant);
        freq_equ = [preamble data];
    end
end

function [rx] = channel(tx, channels)
    n_path = size(tx, 2);

    rx = zeros(size(tx), 'like', tx);
    for i=1:n_path
        ch = channels{i};
        rx(:, i) = ch(tx(:, i));
        % rx(:, i) = tx(:, i);
    end
    rx = sum(rx,2); % 接收分集 合并接收
    rx = rx ./sqrt(2);
    pow = mag2db(mean (abs (rx).^2));
    fprintf('pow=%.2f dB\n', pow);
end