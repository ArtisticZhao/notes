% 参数定义
tau_max = 5e-6;  % 最大延时5微秒
fD_max = 100;    % 最大多普勒频移100Hz
sigma2 = 1;      % 假设单位功率
tau_corr = 1e-6; % 延迟相关性时间
nu_corr = 10;    % 多普勒相关性时间

% 离散化
tau = linspace(0, tau_max, 100);
nu = linspace(-fD_max, fD_max, 200);
[Tau, Nu] = meshgrid(tau, nu);

% 计算RHH
RHH_2 = sigma2 * exp(-2*pi*1i*Nu.*Tau) .* exp(-abs(Tau)/tau_corr - abs(Nu)/nu_corr);

% 可视化
figure;
mesh(Tau, Nu, abs(RHH_2));
xlabel('Delay (s)');
ylabel('Doppler Shift (Hz)');
zlabel('|R_{HH}|');
title('Magnitude of R_{HH} in a Uniform Scattering Environment');

figure;
mesh(abs(RHH));