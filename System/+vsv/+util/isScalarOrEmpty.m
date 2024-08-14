function isSoE = isScalarOrEmpty( obj )
%ISSCALAROREMPTY Summary of this function goes here
%   Detailed explanation goes here

    isSoE = isempty(obj) || isscalar(obj);
end

