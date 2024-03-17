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
end