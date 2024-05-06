close all


%% Generation of UMTS DL
preconfigParams = umtsDownlinkReferenceChannels('TM1_4DPCH'); % Get H-Set parameters
preconfigParams.TotFrames = 2;
preconfigParams.FilterTpe = 'Off';
preconfigParams.OversamplingRatio = 1;
frcWaveform = umtsDownlinkWaveformGenerator (preconfigParams); % Generate H-Set waveform

%% sync
psc_full = ch.PSCH();
ssc_full = ch.SSCH(PrimaryCodeGroup);
ideal_psc_ssc = psc_full + ssc_full;
[r, l] = xcorr(frcWaveform, ideal_psc_ssc);
figure;
plot(l, abs(r));

frcWaveform = frcWaveform(1+16:16+38400);

SC = code.scrambling_code(16 * (8*0+0), 38400);% 主犹码
w  = frcWaveform.*conj(SC);

available_SF = [512, 256, 128, 64, 32, 16, 8, 4];

cdp = zeros(512, length(available_SF));

for s=1:length(available_SF)
    sf = available_SF(s);
    for i=1:sf
        k = i-1;
        [L, U] = CDP_mapping(sf, k);
        cd = code.OVSF(sf, k);
        r = dechannelization(w, cd);
        mean_p = mean(sqrt(real(r).^2 + imag(r).^2));
        cdp((L:U)+1, s) = mean_p;
    end
end


figure;
tiledlayout(2, 4);
axx = [];
for s=1:length(available_SF)
    ax = nexttile;
    axx = [axx ax];
   
    stem(cdp(:, s))
    title(sprintf('SF=%d', available_SF(s)))
end
linkaxes(axx, 'xy')




function [L, U] = CDP_mapping(SF, k)

cSymbolRate = 3.84e6 / SF;

DF = cSymbolRate / 7.5e3;

L = k*DF;
U = L + DF - 1;

end