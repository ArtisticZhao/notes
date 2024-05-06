close all

%% generate signal

M = 4; % QAM order
N = 10000;
sps = 4;
rrc_beta = 0.3;

bb = qammod(randi([0 M-1], N, 1), M, 'gray','InputType','integer', 'UnitAveragePower',true);

% rcoef = rcosdesign(rrc_beta, 10, sps);

phznoise = comm.PhaseNoise;

bb_noise = phznoise(bb);

figure;
plot_scatterIQ(bb_noise)