function isAny = mustBeA( variable, type)
% Will throw an error if variable is not of any of the given types
%
% Note this function is to support older Matlab versions that do not have
% this function. This will be available in R2020b
% 
% Version 1.0 | 2020-10-05 
% $Author: Dr. Daniel Rohrbach
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
% 

    arguments
        variable
        type (1,:) string;
    end
    
    numTypes = numel(type);
    isAny = false;
    
    for i = 1:numTypes
        isAny = isAny | isa(variable, type(i));
    end
    if ~isAny
        error( 'mustBeA:noMatchingType', 'Given variable does not match any of the types' );
    end
    
end

