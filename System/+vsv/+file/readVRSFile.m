function [data, vrsFileHandle] = readVRSFile(file)
%READVRSFILE reads the vrs file and returns formatted data and file handle
%
% Usage:
%   file =  'TestFile030.vrs'
%   [data, vrsFileHandle] = vsv.file.readVRSFile(file);
%
%   % plot the image data
%   figure; imagesc( log(  abs( hilbert( data) ) )  ); colormap gray;
%
%   % print information about the file
%   vrsFileHandle
%   vrsFileHandle.header
%
% @param file (optional) a char file name of a vrs file to open. if empty
% or not given an open dialog will be opened
% @return data the formatted data as captured during acquisition
% @return vrsFileHandle the vrs file handle object to retreive more
% information
%
% $Author: Dr. Daniel Rohrbach,
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.

    if nargin < 1
        file = [];
    end

    vrsFileHandle = vsv.file.VRSBinaryFile(file);
    data = vrsFileHandle.formatData();

end

