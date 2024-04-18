function [x] = dechannelization(y, ovsf)
    sf = length(ovsf);
    nsym = length(y) / sf;
    x = zeros(nsym, 1, 'like', y);
    for i=1:nsym
    idx = (1:sf) + (i-1)*sf;
    x(i) = y(idx).' * ovsf;
    end
    x = x ./ sf;
end