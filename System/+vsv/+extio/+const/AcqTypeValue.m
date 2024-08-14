classdef AcqTypeValue
%ACQTYPEVALUE Summary of this class goes here
%   Detailed explanation goes here
    
    enumeration
        OnDemand("on-demand");
        Sampled("sampled");
    end
    
    properties
        String (1,1) string;
    end
    
    methods
        function this = AcqTypeValue(val)
            this.String = val;
        end
    end
    
end

