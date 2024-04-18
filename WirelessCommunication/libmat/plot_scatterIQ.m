function [] = plot_scatterIQ(data, plot_cfg)
    if nargin ~= 2
        plot_cfg = 'b.';
        cfg = false;
    else
        cfg = true;
    end
    if ~cfg
        scatter(real(data), imag(data), '.');
    else
        plot (real (data), imag(data), plot_cfg);
    end
    grid on;
    % auto-scale.
    max_x = max(abs (real (data))) * 1.2;
    max_y = max(abs(imag(data))) * 1.2;
    xlim ([-max_x, max_x]);
    ylim ([-max_y, max_y]);
end