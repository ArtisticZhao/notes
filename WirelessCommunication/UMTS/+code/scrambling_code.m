function [scr_code] = scrambling_code (N, len)
%SCRAMBLING_CODE 
L = len;

x = code.LFSR([1; zeros(17, 1)], [0 7]);
y = code.LFSR(ones (18, 1), [0 5 7 10]);

for i=1:N
    x = step(x);
end

X_0 = zeros(L, 1, 'int8');
X_m = zeros(L, 1, 'int8');
Y_0 = zeros(L, 1, 'int8');
Y_m = zeros(L, 1, 'int8');

for i=1:L
    X_0(i) = x.LSB();
    X_m(i) = x.read([4 6 15]);
    x = step(x);
    Y_0(i) = y.LSB();
    Y_m(i) = y.read([5 6 8 9 10 11 12 13 14 15]);
    y = step(y);
end
I = -(double(xor(X_0, Y_0))*2 - 1);
Q = -(double(xor(X_m, Y_m))*2 - 1);

scr_code = I + 1j * Q;


end