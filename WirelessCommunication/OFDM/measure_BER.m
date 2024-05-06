function [BER] = measure_BER (freq_equ, idx_data_freq, idx_data_sym, M_data, x_data_bin)
%MEASURE_BER 测筆结果BER
% Args:
%   - freq_equ：完整的时频grid，包含RS信号
%   -idx_data_freq: data 频域位置
%   - idx_data_sym:  data 时域位置
%   - M_data： data调制阶数
%   - x_data_bin：发射的真实比特，比较用
% Returns:
%   -BER:bit员码率
    all_data = freq_equ(idx_data_freq, idx_data_sym);
    all_data = all_data(:);
    x_hat_data = qamdemod (all_data, M_data, 'UnitAveragePower', true);
    x_hat_data_bit = double(dec2bin(x_hat_data, M_data)) - double('0');
    err = sum(sum(abs(x_hat_data_bit-x_data_bin)));
    BER = err/numel(x_hat_data_bit)/M_data;
end