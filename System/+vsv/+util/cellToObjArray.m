function objArr = cellToObjArray(cellArray)
%cellToObjArray converts the given cell array to an object array
%  
% All elements of cellArray must be of the same data type or implement
% matlab.mixin.Heterogeneous 
%
% @param cellArray - @type cell array with all objects being of the same
%                    type
%
% Version 1.0 | 2020-10-05 
% $Author: Dr. Daniel Rohrbach
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
% 

    arr = [ cellArray{:}];
    objArr = reshape( arr, size(cellArray) );

end

