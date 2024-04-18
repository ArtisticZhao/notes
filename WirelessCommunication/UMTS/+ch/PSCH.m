function [psc] = PSCH(slot)
%PSCH Summary of this function goes here
%   Detailed explanation goes here
    if nargin<1
        slot = (0:14);
    end
    if isscalar(slot)
        assert(0<=slot && slot<=14, "slot num in [0,14]")
    else
        if ~any(slot == (0:14))
            warning("slot for vector must is (0:14), aka full frame")
            slot = (0:14);
        end
    end
    is_full = ~isscalar(slot);


    psc = code.primary_sync_code();
    if is_full
        % repmat for 15 slots.
        psc = repmat([psc; zeros(2560-256, 1)], 15, 1);
    end
end

