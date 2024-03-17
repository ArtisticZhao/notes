function [measure_result] = snr_playground(cfg, tx, pre_process, channel, post_demod, demodulate, measurement)
%SNR_PLAYGROUND 算法性能测试
% Args:
%       - cfg 测试配置结构体
%       - cfg.n_test：蒙特卡洛循环次数
%       - cfg.snr: snr array
%       - tx：发射机发射信号
%       - pre_process：预处理函数，如额外的mimo或者 ostbc编码，输出內空域数据
%       - channe1：信道模拟函数，衰落信道模拟，不附加额外的SNR，信道比在此函数处理后添加
%       - post_demod：按收预处理函数，用于合并接收
%       - demodulate：解调函数，解调，解mimo等
%       - measurement:EVM 或者 BER
% Returns:
%       -measure_result: EVM 或者 BER
snr = cfg.snr;
n_test = cfg.n_test;

measure_result = zeros(length(snr), n_test);

for i=1:length(snr)
    cSNR = snr(i);
    fprintf('Current simulate at SNR=%.2f\n', cSNR);
    for t=1:n_test
        pre_rx = pre_process(tx);
        rx = channel(pre_rx);
        rx= awgn(rx, cSNR, 'measured');
        rx = post_demod(rx) ;
        freq = demodulate (rx);
        ber = measurement (freq);

        measure_result(i, t) = ber;
    end
end