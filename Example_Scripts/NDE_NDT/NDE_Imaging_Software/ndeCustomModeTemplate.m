function [plT, plR, W] = ndeCustomModeTemplate(x, z, elementIndex, imagingParameters)

% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use.
%
% File name: ndeCustomModeTemplate.m - Template for custom mode
% look-up-table calculation function for use with the NDE Imaging
% Software
%
% Description: An example processing function for use with the "custom
% look up table" class of custom mode in the NDE Imaging
% Software. Outputs transmit and receive path lengths in uS for simple
% direct mode TFM.
%
% Last update:
% 05/12/2021

% Inputs x and z are vectors of all pixel locations values in mm,
% elementIndex is element number index.  imagingParameters is a data
% structure containing parameters related user specified arrays, geometry,
% inspection configuration.  Outputs, plT and plR are vectors of the
% transmit and acoustic path lengths from the element to the pixel
% locations. W is a vector of amplitude weighting coefficients for each
% pixel location.
                                

ElementPosXMm = imagingParameters.ElementPosXMm(elementIndex);%Element x position for index "elementIndex"
ElementPosZMm = imagingParameters.ElementPosZMm(elementIndex);%Element z position for index "elementIndex"

velocity = imagingParameters.LongitudinalVelMps; %Material longitudinal velocity

%computes path length in uS
plT = sqrt((x - ElementPosXMm).*(x - ElementPosXMm) + (z-ElementPosZMm).*(z-ElementPosZMm))./velocity; %travel from pixel to element in uS
plR = plT;

%Weighting value
W=ones(size(plT));

end

