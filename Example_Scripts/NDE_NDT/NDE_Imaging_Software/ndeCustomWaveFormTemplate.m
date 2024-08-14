function [States] = ndeCustomWaveFormTemplate( imagingParameters)

% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use.
%
% File name: ndeCustomWaveFormTemplate.m - Template for custom
% waveform functions for use with the NDE Imaging Software
%
% Description: A template processing function for use with custom waveform
% acquisition setting in the NDE Imaging software. Outputs example
% 5MHz pulse "States" data in correct format.
%


States = [0 -1 0 1 0 -1 0 1 30; 2 9 12 17 8 17 12 9 0].';


end
