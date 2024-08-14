function [RDataOut] = ndeCustomRFTemplate( RData, imagingParameters )

%imagingParameters

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
% look up table" class of custom mode in the NDE Imaging software. Outputs
% RF data in correct format. 
%
% Last update:
% 05/12/2021

RDatain = double(RData);

RDataOut = int16(RDatain);



end
