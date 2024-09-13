function [DPCCH] = UL_DPCCH(slotformat, i_slot, TPCData, TFCI, FBIData)
%UL_DPCCH 
%   Detailed explanation goes here

param = DPCCH_fields(slotformat);
pilots = code.DPCCH_pilots(param.N_pilot, i_slot);
TFCI = zeros(param.N_TFCI, 1);
FBI = zeros(param.N_TFBI, 1);
TPC = zeros(param.N_TPC, 1);
all_data = [pilots; TFCI; FBI; TPC];

all_data = all_data * 2 - 1;
ovsf = code.OVSF(256, 0);

dpcch = all_data.' .* ovsf;
DPCCH = dpcch(:);

end


function [pdcch_param] = DPCCH_fields(slotformat)
% refer to TS25.211 ss5.2.1.1 Table 2
if(slotformat~=0)
    error("Only support format 0")
end
    pdcch_param = struct;
    pdcch_param.N_pilot = 6;
    pdcch_param.N_TPC = 2;
    pdcch_param.N_TFCI = 2;
    pdcch_param.N_TFBI = 0;
end