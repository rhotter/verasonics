classdef SampleMode
%SAMPLEMODE Summary of this class goes here
%   Detailed explanation goes here
%
%
% Version 1.0 | 2020-12-29 
% $Author: Dr. Daniel Rohrbach 
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
% 
    
    enumeration
        
        % Nyquist sampling (200% bandwidth) of demodFrequency. ‘NS200BWI‘
        NS200BW;
        
        % 2-1 interleaved sampling (requires 2 acquisitions). 
        NS200BWI;
        
        % 100% bandwidth sampling of demodFrequency.
        BS100BW;
        
        % 67% bandwidth sampling of demodFrequency.
        BS67BW;
        
        % 50% bandwidth sampling of demodFrequency.
        BS50BW;
        
        % sample rate set by decimSampleRate
        custom;
        
    end
    
   
end

