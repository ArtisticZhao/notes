function [code] = OVSF(SF, m, debug)
%OVSF
if nargin < 3
    debug = false;
end

code = my_ovsf(SF, m);
if debug
    ov = comm.OVSFCode ('SpreadingFactor', SF, 'Index', m, 'SamplesPerFrame',SF);
    code_2 = ov();
    
    err = sum(abs(code - code_2));
    assert(err==0, "mismatch OVSF!")
end
end


function [code] = my_ovsf(SF, k)
    assert(log2(SF) == round(log2(SF)), "SF must be power of 2")
    if (SF==1)
        code = 1;
        return
    end
    father_code = my_ovsf(SF/2, floor(k/2));
    if(mod(k,2)==0)
        code = [father_code; father_code];
    else
        code = [father_code; -father_code];
    end

end

