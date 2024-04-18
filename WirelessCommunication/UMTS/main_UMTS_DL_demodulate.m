close all

sps = 16;

%% Generation of UMTS DL
preconfigParams = umtsDownlinkReferenceChannels('TM1_4DPCH'); % Get H-Set parameters
preconfigParams.TotFrames = 2;
preconfigParams.FilterTpe = 'RRC';
preconfigParams.OversamplingRatio = sps;
frcWaveform = umtsDownlinkWaveformGenerator (preconfigParams); % Generate H-Set waveform

%% RX: match filter
rrc = rcosdesign(0.22, 10, sps, "sqrt");
rc = rcosdesign(0.22, 10, sps, "normal");
[frcWaveform, axx] = filter_sync(rrc, frcWaveform, 1);

%%  downsample to 3.84M * 4
phi = 1;
frcWaveform = frcWaveform(1+phi:sps/4:end); sps = sps/4;
eps = OM_timing_error_estimate(frcWaveform, sps, true);
frcWaveform = timing_error_correct(frcWaveform, eps, sps, 38400*2); sps=1;
figure; plot_scatterIQ(frcWaveform)

%% PSC search

psc_full = ch.PSCH();
if sps > 1
    psc_full = filter_sync(rc, psc_full, sps);
end
[r, l] = xcorr(frcWaveform(1:38400), psc_full);
[~, ml] = max(abs(r));
dly = l(ml);
fprintf('PSC find at %d\n', dly);
figure;
plot (l, abs(r))
title('PSC sync')

frcWaveform = frcWaveform(1+dly:sps:dly+38400*sps); % downsample

%% SSC search
ssc_rho = zeros(64, 1);
for i=1:64
    ssc_rho(i) = abs(sum(frcWaveform.*conj(ch.SSCH(i-1))));
end
[~,PrimaryCodeGroup] = max(ssc_rho);
PrimaryCodeGroup = PrimaryCodeGroup - 1;
fprintf('PrimaryCodeGroup: %d\n', PrimaryCodeGroup);
figure;
stem((0:63).', ssc_rho)

%% scrambling code determine
scrambling_rho = zeros(8, 1);
pc =[];
for i=1:8 % 主扰码组已经通过SSC确定了，之后在该组内确定具体的主扰码号 0-7
    pcpich_1frame = ch.pcpich(PrimaryCodeGroup, i-1);
    pc = [pc pcpich_1frame];
    scrambling_rho(i) = abs (sum(frcWaveform(1:38400).*conj(pcpich_1frame)));
end
figure;
stem((0:7).', scrambling_rho)
[~, l] = max(scrambling_rho);
PrimaryScramblingCodeIdx = l-1;

% xcorr way
figure;
tiledlayout(2, 4); 
axx=[];
for i=1:8
    [r, l] = xcorr(frcWaveform, pc(:, i));
    ax = nexttile;
    plot(l, abs(r));
    axx=[axx ax];
end
linkaxes (axx, 'xy')

%% demodulate
% 解扰
SC = code.scrambling_code(16 * (8*PrimaryCodeGroup+PrimaryScramblingCodeIdx), 38400);% 主犹码
wave = frcWaveform(1:38400).*conj(SC);
% 解扩
ovsf_pcpich = code.OVSF(256, 0);
c_pcpich = dechannelization(wave, ovsf_pcpich);
c_pcpich_m = mean(c_pcpich);
h = c_pcpich_m ./ (1+1j);
c_pcpich = c_pcpich ./ h;
figure;
plot_scatterIQ(c_pcpich)
evm = comm.EVM;
e = evm(c_pcpich, repmat(1+1j, 150, 1));
fprintf('EVM of P-CPICH: %.2f %%\n', e);


