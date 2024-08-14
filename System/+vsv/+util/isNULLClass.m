function isnull = isNULLClass(obj)
% Checks whether given handle object is null or not
% 
%   this tests whether object is an empty arrya of class(obj). this is not
%   using empty function because some collections are implementing isempty
%   function. An object is null if either its an empty object
%
%
%  NOTE: This function is undocumented and may change in future releases
%
% Parameters:
%   @param val - the value to check (generic)
%
% Return values:
%   @return is0 - @type logical, which is true if val is NULL false otherwise
%
% Version 1.0 | 2019-09-23 
% $Author: Dr. Daniel Rohrbach 
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
%

    isnull = false;
    if isa(obj, 'handle')
        isnull = ~isvalid( obj );
        if isempty(isnull)
            isnull = true;
        end
    end
    
end

