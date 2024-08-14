function Media = createDopplerMedia()
% Create Media structure for simulating Doppler string phantom
%
% Parameters:
%
% Return value:
%   @return Media - @type Media struct to be use by runAcq()
%
% Copyright 2021 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.


DEPTH =  40;  %set depth of string
REFLECTIVITY = 0.001;  %set reflectivity coefficient of string

Media.Model = 'PointTargets1';

% Generate line of scatterers
Xpoints = (-45:5:45);
for a = 1:length(Xpoints)
    Media.MP(a,:) = [Xpoints(a), 0, DEPTH, REFLECTIVITY];
end

Media.numPoints = size(Media.MP,1);
Media.function = 'vsv.gpu.media.moveDopplerPoints';  %use moveDopplerPoints to give points velocity