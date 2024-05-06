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
    max_x = max(max(abs(real(data))));
    max_y = max(max(abs(imag(data))));

    mm = max(max_x, max_y) * 1.2;
    xlim ([-mm, mm]);
    ylim ([-mm, mm]);
end