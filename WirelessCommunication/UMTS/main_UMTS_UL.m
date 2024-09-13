close all

%% generate waveform
frc = 'FRC1';     % FRC number
preconfigParams = umtsUplinkReferenceChannels(frc);        % Get FRC parameters
preconfigParams.OversamplingRatio = 1;
preconfigParams.FilterType = 'Off';
preconfigParams.DPDCH.Enable = 'Off';
preconfigParams.HSUPA.Enable = 'Off';

frcWaveform = umtsUplinkWaveformGenerator(preconfigParams);% Generate FRC waveform

%% channel
rx = frcWaveform .* exp(1j * pi/6);

%% phase error estimate by pilots in DPCCH



ref_pilots_slot0 = ch.UL_DPCCH(0, 0);
scrambling_code = code.scrambling_code(preconfigParams.ScramblingCode, 38400);

rx_1 = frcWaveform .* conj(scrambling_code);
figure;
scatterIQ(frcWaveform)

test_1 = m_DPCCH  .* conj(scrambling_code);