function [y] = nco(x, Freqoffset, fs, phase)
%NCO 
%   fs: default = 1
%   phase: default = 0.0, rad
if nargin < 4
    phase = 0;
end
if nargin < 3
    fs = 1;
end
% time vector
N = length(x);
t = (0:N-1).' / fs;              % time in seconds, for every sample

frequencyShift = exp(1j * 2 * pi * Freqoffset * t) .* exp(1j * phase);


y = x .* frequencyShift;

end

