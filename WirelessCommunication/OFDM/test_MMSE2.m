close all

generate_demo_OFDM2;
cSNR = 20;

%% channel .....
if 0
pathDelays = [0 10 20 30 ]*1e-9; % in seconds
avgPathGains0 = [0 -3 -10 -8 ]; % dB

fD = 0; % Max Doppler shift in Hz
rayleigh = comm.RayleighChannel('SampleRate',fs, ...
        'PathDelays' ,pathDelays,...
        'AveragePathGains',avgPathGains0,...
        'MaximumDopplerShift',fD, ...
        'Visualization', 'Impulse and frequency responses');
output1 = rayleigh(tm);
release(rayleigh)
output2 = rayleigh(tm);

% **if we call the channel multiple times , return different responses.**
isequal(output1, output2)
end

%% find a way to have same response in multiple time
if 0
pathDelays = [0 10 20 30 ]*1e-9; % in seconds
avgPathGains0 = [0 -3 -10 -8 ]; % dB

fD = 0; % Max Doppler shift in Hz
rayleigh = comm.RayleighChannel('SampleRate',fs, ...
        'PathDelays' ,pathDelays,...
        'AveragePathGains',avgPathGains0,...
        'MaximumDopplerShift',fD, ...
        'RandomStream','mt19937ar with seed', ...
        'Seed',22, ...  
        'Visualization', 'Off');
output1 = rayleigh(tm);
release(rayleigh);
rayleigh.RandomStream = 'Global stream';
rng(22);

output2 = rayleigh(tm);
isequal(output1, output2)
end
%% channel H
pathDelays = [0 10 20 30 ]*1e-9; % in seconds
avgPathGains0 = [0 -3 -10 -8 ]; % dB

fD = 0; % Max Doppler shift in Hz
rayleigh = comm.RayleighChannel('SampleRate',fs, ...
        'PathDelays' ,pathDelays,...
        'AveragePathGains',avgPathGains0,...
        'MaximumDopplerShift',fD, ...
        'RandomStream','mt19937ar with seed', ...
        'Seed',22, ...  
        'Visualization', 'Off');
rx = rayleigh(tm);  % get signal through channel
release(rayleigh);
rayleigh.RandomStream = 'Global stream';
rng(22);

impulse = [1; zeros(nfft-1, 1)];
h = rayleigh(impulse);  % get the impulse response on the same channel
H = fft(h);

fax = (-nfft/2:nfft/2-1) * (fs / nfft);

figure;
plot(fax, mag2db(abs(fftshift(H)))); hold on;
title('frequency response H')



% demod
freq = azcomm.ofdmdemod(rx, repmat(cp_len, nsym, 1), nfft, nsc_all);

H_hat_LS = channel_est_LS(freq, grids, (1:nsc_all), idx_preamble_sym);

% plot estimate by LS
plot(fax((nfft-nsc_all)/2+1:(nfft-nsc_all)/2+nsc_all), mag2db(abs(H_hat_LS)), 'rx');

function [H] = channel_est_LS(rx_grid, raw_grid, rs_k, rs_sym)
    N_sym = size(rx_grid,2);
    H = rx_grid(rs_k, rs_sym) ./ raw_grid(rs_k, rs_sym);
    % Expand H to whole grid, using symbol direction average the H
    H = repmat(mean(H, 2), 1, N_sym); 
end