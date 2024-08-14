function mustBeNonInf(val)
% MUSTBENONINF throws an error 
%
% Version 1.0 | 2020-10-05 
% $Author: Dr. Daniel Rohrbach
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
% 

    if any( isinf(val) )
        error( 'vsv:mustBeNonInf', 'value cannot be Inf' );
    end

end

