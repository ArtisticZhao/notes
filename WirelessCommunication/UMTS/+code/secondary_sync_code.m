function [c_ssc] = secondary_sync_code(k)
%SSC UMTS secondary synchronization codes, SSC
% k for 1-16
assert(k>=1 && k<=16, "SSC idx in [1, 16]");
a = [1, 1, 1, 1, 1, 1, -1, -1, 1, -1, 1, -1, 1, -1, -1, 1];
b = [a(1:8) -a(9:16)];
z = [b, b, b, -b, b, b, -b, -b, b, -b, b, -b, -b, -b, -b, -b];

H = hadamard(8);

m = 16*(k-1);
hm = H(m+1, :);

c_ssc = (1+1j) * (hm .* z);
c_ssc = c_ssc(:);
end


function [H] = hadamard(k)
    if k==0
        H = 1;
    else
        Hk = hadamard(k-1);
        H = [Hk Hk;
             Hk -Hk];
    end
end