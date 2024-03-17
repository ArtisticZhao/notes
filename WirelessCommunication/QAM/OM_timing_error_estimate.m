function epsilon = OM_timing_error_estimate(rxSignal, sps, debug)
% OM_timing_error_estimate OERDER & MEYR 非线性平方率 开环采样误差估计
% @Refer Oerder, M.; Meyr, H. Digital Filter and Square Timing Recovery. IEEE Transactions on Communications 1988, 36 (5), 605–612. https://doi.org/10.1109/26.1476.
% @Args:
%   - rxSignal 接收信号，需先做接收匹配滤波
%   - sps samples per symbol
%   - debug default false;
% @returns:
%   - epsilon: 采样偏差， 相对**采样时钟**
if nargin < 3
    debug = false;
end
demod_Nsym = floor(length(rxSignal) / sps);
x = abs(rxSignal).^2;
% x = (rxSignal).^2;
X = 0;
for i=1:length(rxSignal)
    X = X + x(i) * exp(-1j*2*pi*(i-1) / sps);
end
epsilon = -angle(X)/2/pi*sps;
% 算法本质是平方后， 在1/Tsym处存在一个符号定时时钟的分量，
% 看该分量的相位值即可得到，采样时钟偏差 epsilon
if debug
    Fx = fft(x);
    F_axis = (0:1/demod_Nsym:sps-1/demod_Nsym).';
    figure;
    sgtitle("spectrum of x^2")
    ax1=subplot(211);
    plot(F_axis, abs(Fx));
    grid on
    ax2=subplot(212);
    plot(F_axis, -angle(Fx)/2/pi*sps);
    ylabel("\tau Ts");
    grid on
    linkaxes([ax1 ax2], 'x')
    
    epsilon2 = -angle(Fx(demod_Nsym+1))/2/pi*sps;
    fprintf("Epsilon = %f, Epsilon from FFT = %f\n", epsilon, epsilon2);
end

end