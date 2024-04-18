classdef LFSR
    properties 
        buffer;
        update;
    end
    methods
        function obj = LFSR(init, update)
            obj.buffer = int8(init);
            obj.update = update + 1;
        end
        function obj = step(obj)
            s = obj.buffer(2:end);
            u = mod (sum(obj.buffer(obj.update)), 2);
            obj.buffer = [s; u];
        end
        function x0 = LSB(obj)
            x0 = obj.buffer(1);
        end
        function xx = read(obj, mask)
            xx = mod(sum(obj.buffer(mask+1)), 2);
        end
    end
end