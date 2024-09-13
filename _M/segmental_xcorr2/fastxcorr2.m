close all
rng(0);  % Set random seed for reproducibility
N = 1228800;  % Length of the input sequence
M = 8192+20;  % Length of the reference sequence with offset
offset = 400;  % Offset for generating the reference sequence

% Generate random input and reference sequences
in = randi(9, N, 1);  % Input sequence: random integers from 1 to 9
ref = in(1+offset:M+offset);  % Reference sequence: part of input with offset

% Example for testing with shorter sequences
% in = [3 1 4 1 5 9 2 6 5 3 5 8 9 7 9 3 2 3 8 9 6].';  % Example input
% ref = [1 4 1].';  % Example reference

%% MATLAB's xcorr function to compute cross-correlation
tic
r_matlab = xcorr(in, ref);  % Compute cross-correlation using built-in xcorr
toc  % Display time elapsed

%% Custom implementation of cross-correlation
% 0. parameters calculation
tic
n_in = length(in);  % Length of input
n_ref = length(ref);  % Length of reference
n_fft = findTransformLength(2*n_ref);  % Find transform length (probably a custom function)

% Segment size for FFT processing
n_seg_size = n_fft + 1 - n_ref;  
n_fft_shift_back = n_ref-1;  % **Shift for FFT in the back
n_seg = ceil(n_in/n_seg_size);  % Number of segments
n_padding = n_seg*n_seg_size-n_in;  % Padding needed
n_overlap = n_fft - n_seg_size;  % Overlap between segments
% Display calculated values
fprintf("n=%d, m=%d, l=%d, nseg=%d\n", n_seg_size, n_ref, n_fft, n_seg);


%% 1. Reshape the input sequence and apply zero padding
in_ = reshape([in; zeros(n_padding, 1)], n_seg_size, n_seg);  % Input reshaped into segments with padding

ref_f_conj = conj(fft(ref, n_fft));  % Take FFT of the reference and conjugate it
offset0 = n_in-n_ref;  % Offset calculated based on input and reference lengths

% FFT of the input sequence
t_in_f = fft(in_, n_fft);  % FFT of each segment of input
res_ = ifft(t_in_f .* ref_f_conj);  % Element-wise multiplication in the frequency domain and inverse FFT

% Adjust the result to handle FFT shift
res_ = [res_(end-n_fft_shift_back+1:end, :); res_(1:end-n_fft_shift_back, :)];  

%% 2. Accumulate results for each segment to construct whole result
% ---- Copy part of result first
R = [zeros(offset0, 1); ...
     res_(:, 1); ...
     reshape(res_(1+n_overlap:n_fft, 2:end), [], 1)];
% ---- Overlap-add method for handling overlapping segments
for n=1:n_seg
    % Add the overlapping part from the previous segment
    if n ~= 1
        R((1:n_overlap)+offset0+n_seg_size*(n-1)) = R((1:n_overlap)+offset0+n_seg_size*(n-1)) + res_(1:n_overlap, n);
    end
end

toc  % End timer for the custom implementation


%% compare result
figure;
ax1 = subplot(211);  % Create the first subplot
plot(abs(r_matlab));  % Plot the absolute value of r_matlab (result of MATLAB's xcorr)
ax2 = subplot(212);  % Create the second subplot
plot(abs(R));  % Plot the absolute value of the custom implementation's result R
linkaxes([ax1 ax2], 'xy');  % Link the axes of both subplots for panning and zooming

%% EOF
function m = findTransformLength(m)
    m = 2*m;  % Multiply the input by 2
    while true
        r = m;
        for p = [2 3 5 7]  % Loop over prime factors 2, 3, 5, 7
            while (r > 1) && (mod(r, p) == 0)  % While r is divisible by p
                r = r / p;  % Remove the factor from r
            end
        end
        if r == 1  % If r is fully factorized
            break;  % Exit the loop
        end
        m = m + 1;  % Increment m and try again
    end
end