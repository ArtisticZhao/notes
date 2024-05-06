close all
% @brief: TM4 only have PSCH, SSCH, PCCPCH.
%         here we will discuss the demodulation of PCCPCH
%         Additionally, I also implement a simple frequency offset
%         estimation via PSCH/SSCH signals.

%% Generation of UMTS DL
preconfigParams = umtsDownlinkReferenceChannels('TM4'); % Get H-Set parameters
preconfigParams.TotFrames = 1;
preconfigParams.FilterType = 'Off';
preconfigParams.OversamplingRatio = 1;
frcWaveform = umtsDownlinkWaveformGenerator(preconfigParams); % Generate H-Set waveform

PrimaryCodeGroup = fix(preconfigParams.PrimaryScramblingCode/8);
PrimaryScramblingCodeIdx = mod(preconfigParams.PrimaryScramblingCode,PrimaryCodeGroup*8);

%% frequency offset estimate
% assume the frcWaveform have perfect sychnorized and corrected.
df = 200; % Hz
fs = 3.84e6; % 1sps

frcWaveform = nco(frcWaveform, df, fs);

% frequency offset estimate by PSC and SSC

psc_full = ch.PSCH();
ssc_full = ch.SSCH(PrimaryCodeGroup);
ideal_psc_ssc = psc_full + ssc_full;

sync = frcWaveform .* conj(ideal_psc_ssc);
sync = reshape(sync, 2560, []);

sync_const = sum(sync);

% calc freq offset
A_sync = unwrap(angle(sync_const.'));

figure;
subplot(121)
plot_scatterIQ(sync_const)
subplot(122)
plot((0:14), A_sync)
xlabel('#slot')

d_angle = diff(A_sync);
f_offset = d_angle / 2 / pi / 2560 * fs;
f_offset = mean(f_offset);
fprintf('Frequency offset: %.2f Hz\n', f_offset);
% correction
frcWaveform = nco(frcWaveform, -f_offset, fs);

%% scrambling code for P-CCPCH using PSC

% Primary scrambling code
SC = code.scrambling_code(16 * (8*PrimaryCodeGroup+PrimaryScramblingCodeIdx), 38400);

wave = frcWaveform(1:38400).*conj(SC);

%% channelisation for P-CCPCH use Cch,256,1 ref. 25.213 ss5.2.1
pccpch_ovsf = code.OVSF(256, 1);

mask_pccpch = [zeros(256, 1); ones(2560-256, 1)];
mask_pccpch = repmat(mask_pccpch, 15, 1);
% mask for removing the PSCH/SSCH
wave = wave.*mask_pccpch;
c_pccpch = dechannelization(wave, pccpch_ovsf);

figure;
plot_scatterIQ(c_pccpch)
