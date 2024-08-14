function is0 = isNULL(val)
% Checks whether given value is null or not
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

    is0 = isnumeric( val ) && isempty( val);
   
end

