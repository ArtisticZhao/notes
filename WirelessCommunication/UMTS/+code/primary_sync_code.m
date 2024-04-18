function [c_psc]=primary_sync_code()
%PSC UMTS primary synchronization codes, PSC
a = [1, 1, 1, 1, 1, 1, -1, -1, 1, -1, 1, -1, 1, -1, -1, 1];
c_psc = (1+1j) * [a, a, a, -a, -a, a, -a, -a, a, a, a, -a, a, -a, a, a];
c_psc = c_psc(:);
end