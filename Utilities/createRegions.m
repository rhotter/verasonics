function [Region,numRegs] = createRegions(PData)
%
% Copyright 2001-2017 Verasonics, Inc.  All world-wide rights and remedies under all intellectual property laws and industrial property laws are reserved.  Verasonics Registered U.S. Patent and Trademark Office.
%
% [Region,numRegs] = createRegions(PData)
% For the Vantage 3.0 SW release, the new computeRegions function replaces
% createRegions.  For backwards compatibility with older scripts written
% for earlier Vantage releases, this function simply calls the new
% computeRegions.  The 'numRegs' field is no longer used, so a dummy value
% of 1 is returned.

[Region] = computeRegions(PData);
numRegs = 1;

end
