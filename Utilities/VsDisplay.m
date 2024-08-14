function VsDisplay(dwindex,DisplayData)
%VsDisplay - updates DisplayWindow(dwindex) with DisplayData
% Copyright 2001-2017 Verasonics, Inc.  All world-wide rights and remedies under all intellectual property laws and industrial property laws are reserved.  Verasonics Registered U.S. Patent and Trademark Office.
%
%   This function is called by runAcq to update the DisplayWindow defined
%   in the Resource structure.  The handle to the DisplayWindow image is
%   obtained from Resource.DisplayWindow(dwindex).imageHandle.
%
%   If the user desires to manage a custom display, they can replace this
%   function (keeping the same name) with their own. To provide runAcq with
%   the correct size information, the following Resource.DisplayWindow
%   attributes should be set in the SetUp script:
%       mode: '2d' or '3d'  (if not set, defaults to '2d')
%       Orientation: 'xz','yz','xy'  (defaults to 'xz')
%       Position  (the width and height values are required to set size.)
%       numFrames  (if a cineloop buffer is desired)
%       ReferencePt (required to set location with respect to PData)
%       pdelta
%   The dwindex value passed to this function by runAcq will be determined
%   by the imageDisplay Process structure.

Resource = evalin('base','Resource');
hFig = handle(Resource.DisplayWindow(dwindex).imageHandle);
hFig.CData=DisplayData;

end

