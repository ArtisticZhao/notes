function [code] = OVSF(SF, m)
%OVSF 此处显示有关此函数的摘要% 此处显示详细说明
ov = comm.OVSFCode ('SpreadingFactor', SF, 'Index', m, 'SamplesPerFrame',SF);
code = ov();
end