function [x] = timing_error_correct(y, epsilon, sps, len)
%timing error Correct and dounsample
% Args:
%   - y: input time sequence
%   - epsilon: timing error in samples
%   - sps: downsample factor, default 1
%   - len: length of output
if nargin < 4
    len = floor (length(y)/sps);
end
if nargin < 3
    sps = 1;
    len = length(y);
end

idx = (1:sps:len*sps);
x = spline((1:length(y)), y, idx+epsilon);
x = x(:);
end