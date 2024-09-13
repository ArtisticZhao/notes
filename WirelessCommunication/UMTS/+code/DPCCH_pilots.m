function [pilots_slot] = DPCCH_pilots(Npilot, i_slot)
%DPCCH_PILOTS get DPCCH pilots
%   @ref TS 25.211 ss5.2.1.1 table 3 and 4

switch (Npilot)
    case 6
        pilots = [
            1 1 1 1 1 0;
            1 0 0 0 1 0;
            1 0 1 1 1 0;
            1 0 1 1 1 0;
            1 1 0 1 1 0;
            1 1 1 1 1 0;
            1 1 1 1 1 0;
            1 0 1 1 0 0;
            1 0 0 1 1 0;
            1 0 1 1 0 1;
            1 0 0 0 0 1;
            1 0 1 1 1 1
        ];
    otherwise
        error('unimplemented data!')
end
pilots_slot = pilots(i_slot+1, :);
pilots_slot = pilots_slot(:);
end

