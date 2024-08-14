function mustBeCellString(cell)
%MUSTBECELLSTRING will throw an error if cell is not a cell string
%   A cell string is a cell array that only contains char arrays
%
% Version 1.0 | 2020-05-15 
% $Author: Dr. Daniel Rohrbach
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
% 

    if ~iscellstr(cell) %#ok 
        error('vsv:mustBeCellString', 'Given cell array is not a cellstring')
    end
    
end
