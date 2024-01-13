% @file main_om.m
% @brief 利用O&M SLN算法进行开环符号同步，采集数据验证
% @author zhao.yuhao
close all
is_ideal = true;

if is_ideal
    %% ideal data generation
    % -- generate TX -- RX
    % Parameters
    M = 256; % BPSK modulation, so M = 2
    sps = 4; % Samples per symbol (sps parameter)
    dataLength = 1000; % Length of the data
    osf = 64; % over sample factor
    delay_offset = round(osf * 0.38*3);
    fprintf("Fraction Timeoffset = %f Ts\n", delay_offset/osf);
    % Generate random binary data
    data = randi([0 M-1], dataLength, 1);
    
    modData = qammod(data, M, 'UnitAveragePower', true);
    figure;
    scatterIQ(modData)
    title(sprintf("constellation of %d-order modulation", M))
    
    % -- pulse shaping
    % Root-raised Cosine Filter parameters
    rolloff = 0.35; % Rolloff factor of the filter
    span = 10; % Filter span in symbols
    rcosFilter_os = rcosdesign(rolloff, span, sps*osf, 'sqrt');
    rcosFilter = rcosdesign(rolloff, span, sps, 'sqrt');
    % fvtool(rcosFilter);
    
    txSignal = upfirdn(modData, rcosFilter_os/sum(rcosFilter_os)*sps*osf, sps*osf);
    txSignal = txSignal(groupDelay(rcosFilter_os)+1:groupDelay(rcosFilter_os)+sps*osf*dataLength);
    
    % -- RX + channel
    % downsample to sps
    rxSignal = txSignal(1+delay_offset:osf:end);
    rxSignal = upfirdn(rxSignal, rcosFilter / sum(rcosFilter), 1);
    rxSignal = rxSignal(groupDelay(rcosFilter)+1: groupDelay(rcosFilter)+sps*dataLength);
    demod_Nsym = dataLength;
    % -- phase error
    phase_err = 30;
    rxSignal = rxSignal .* exp(1j * deg2rad(phase_err));
else % load from samples
    %% parameter % load
    beta = 0.3;
    max_demod_Nsym = 800;
    % load data
%     [~, rx_buffer] = compare2('./data/qpsk_3M072_ref_r0.3'); sps = 40;
%     [~, rx_buffer] = compare2('./data/qpsk_30M72_ref_r0.3'); sps = 4;
    % [~, rx_buffer] = compare2('./data/qam256_30M72_ref_r0.3'); sps = 4;
    % [~, rx_buffer] = compare2('./data/qam256_30M72_ref_r0.75'); sps = 4; beta = 0.75;
    [~, rx_buffer] = compare2('./data/qam256_3M072_ref_r0.3'); sps = 40;
    % [~, rx_buffer] = compare2('./data/qam1024_30M72_ref_r0.3'); sps = 4;
%     [~, rx_buffer] = compare2('./data/qam1024_3M072_ref_r0.75'); sps = 40; beta = 0.75;
%     [~, rx_buffer] = compare2('./data/qam1024_7M68_ref_r0.75'); sps = 16; beta = 0.75;
    data_Nsym = floor(length(rx_buffer)/sps); % calculate symbols contain in data.
    
    demod_Nsym = min(max_demod_Nsym, data_Nsym-1);
    % clip data from buffer
    rx = rx_buffer(1:(demod_Nsym+1) * sps);
    rx = normalize(rx);
    % rx rcosfilter
    rcosFilter = rcosdesign(beta, 10, sps, 'sqrt');
    rxSignal = upfirdn(rx, rcosFilter / sum(rcosFilter), 1);
    rxSignal = rxSignal(groupDelay(rcosFilter)+1: groupDelay(rcosFilter)+sps*(demod_Nsym+1));
    
    rxSignal = rxSignal(1+sps:end); % clip the first symbol to avoid effect by filter.
end
%% O&M algorithm
epsilon = OM_timing_error_estimate(rxSignal, sps, true);

%% resample
decision_idx = 1:sps:length(rxSignal);
rxModSync = spline(1:length(rxSignal), rxSignal, decision_idx +epsilon).';
% evmer.reset()
% evmSync = evmer(rxModSync, modData);
% fprintf("EVM with time sync %.2f %%\n", evmSync);

figure;
% plot(real(txSignal(1:osf:end)))
% hold on
stem(1:sps:sps*demod_Nsym, real(rxModSync))
xlim([0 sps*10])

figure;
scatterIQ(rxModSync, 'bo');
title("symbol sync")
hold on
% scatterIQ(modData, 'r+');

%% viterbi phase estimate
z = rxModSync;
vM = 4;
l=16;
Im = 0;
Re = 0;
for i=1:demod_Nsym
    tmp = abs(z(i)).^l * exp(1j*vM*angle(z(i)));
    Im = Im + imag(tmp);
    Re = Re + real(tmp);
end
theta = 1/vM * atan2(Im,Re) + pi/4; % viterbi get pi/4 rotation
fprintf("phase error = %d deg\n", rad2deg(theta));
figure;
scatterIQ(rxModSync.*exp(-1j*theta), 'bo');
title("fix phase error")

function [d] = groupDelay(coef)
    d = (length(coef)-1)/2;
end