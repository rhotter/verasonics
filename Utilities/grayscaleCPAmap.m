function [cmap] = grayscaleCPAmap()
%
% Copyright 2001-2017 Verasonics, Inc.  All world-wide rights and remedies under all intellectual property laws and industrial property laws are reserved.  Verasonics Registered U.S. Patent and Trademark Office.
%
% Function to create colormap for CPA.

M = hot(256);
cmap = zeros(256,3);

% load linear greyscale in 1:128
cmap(1:128,1) = (0:(1/127):1)';
cmap(1:128,2) = (0:(1/127):1)';
cmap(1:128,3) = (0:(1/127):1)';

% load color map in 129:256
cmap(129:256,:) = M(1:2:256,:);

return

