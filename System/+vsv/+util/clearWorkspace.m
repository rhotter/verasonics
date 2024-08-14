function clearWorkspace( clearMex )
% Clears all base workspace variables without clearing functions
%
% Clear all is not recommended because it also deletes all function handles
% and just-in-time compiled code that is used by Matlab to optimize
% function calls. This can cause significant performance issues. Therefore
% Verasonics recommends using clearWorkspace function instead
%
% @Usage
%   vsv.util.clearWorkspace();
%
% Version 1.0 | 2020-04-02
% $Author: Dr. Daniel Rohrbach
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
%

    if nargin < 1 || isempty(clearMex)
        clearMex = true;
    end
    if clearMex
        evalin( 'base', 'clear(''mex'');');
    end
    evalin( 'base', 'clear(''variables'');');


end
