function [y, axx] = filter_sync(coef, x, interp, debug)
%FILTER_SYNC Summary of this function goes here
%   Detailed explanation goes here
    if nargin < 4
        debug = false;
    end
    if nargin < 3
        interp = 1;
    end
    
    assert(size(x, 1) == 1 || size(x, 2) == 1, "input x must 1D!");
    x = x(:);
    
    if interp > 1
        % upsample with 0
        x_ = [x.'; zeros(interp-1, length(x))];
        x_ = x_(:);
    else
        x_ = x;
    end
    
    dly = groupDelay(coef);
    % filter
    y_f = filter(coef, abs(sum(coef)) / interp, [x_; zeros(dly, 1, "like", x_)]);
    y = y_f(1+dly: dly+length(x)*interp);

    axx = [];
    if debug
        figure;
        sgtitle('debug of sync filter')
        ax=subplot(211);
        plot(real(y)); hold on;
        plot(real(x))
        legend('filtered', 'interpreted');
        ax1=subplot(212);
        plot(imag(y)); hold on; 
        plot(imag(x_))
        legend('filtered', 'interpreted');
        axx = [ax ax1];
        linkaxes(axx, 'xy');
    end
end

function [d] = groupDelay(coef)
    d = (length(coef)-1)/2;
end
