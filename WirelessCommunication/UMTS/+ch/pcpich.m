function [pcpich] = pcpich(n_main_group, k)
%PCPICH 
% n_main_group: primary scrambling code group, 0-63
% k: primary scrambling code idx, 0-7

assert(0<=n_main_group && n_main_group<=63, "PCPICH scrambling group with primary code! 0-63");
assert(0<=k && k<=7, "PCPICH scrambling with primary code!0-7");

ovsf_pcpich = comm.OVSFCode('SpreadingFactor', 256, 'Index', 0, 'SamplesPerFrame',256);
spread_code_pcpich = ovsf_pcpich();

pcpich = repmat(spread_code_pcpich, 10, 1); % pcpich 1slot
pcpich = repmat(pcpich, 15, 1); % pcpich 1frame

pcpich = pcpich + 1j*pcpich;
pcpich = pcpich ./ 2;

n_psc = 16 * (8*n_main_group+k);
psc = code.scrambling_code(n_psc, 38400);
pcpich = pcpich .* psc;
end