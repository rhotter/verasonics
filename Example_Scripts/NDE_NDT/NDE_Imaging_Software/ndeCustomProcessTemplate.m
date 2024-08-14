function [imgData] = ndeCustomProcessTemplate( RData, imagingParameters)

% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use.
%
% File name: ndeCustomProcessTemplate.m - Template for custom
% processing functions for use with the NDE Imaging Software
%
% Description: A template processing function for use with the "custom
% look up table" class of custom mode in the NDE Imaging
% software. Outputs random image data in correct format.
%
% Last update:
% 05/12/2021

imgData = rand(size(imagingParameters.XCoordsMm));
 
end
