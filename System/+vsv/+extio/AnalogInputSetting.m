classdef (Sealed) AnalogInputSetting < vsv.extio.AbstractInputSetting
%ANALOGINPUTSETTING defines the settings for a set of analog inputs
%
% The number of analog input channels is specified by 'NumInput'. 'Mode'
% specfies connection mode, which may be "differential" or "single-ended"
% and will default to "single-ended" is unspecified.  Range is peak voltage
% range of input (V). 'AcqType' specifies the type of acquisition; "OnDemand"
% for scalar reads, "sampled" for sampled, internally buffered vector
% reads. AcqType will default to "OnDemand" if unspecified.
%
% On demand acquisitions are software controlled acquisitions of the scalar
% voltage of each channel.  Sampled acquisitions are buffered internally on
% the IO device at a specified rate and may be hardware triggered.
%
% For sampled acquisitions, the follow additional settings apply. 'Samples'
% defines number of samples to be acquired. 'SampleRate' is the acquisition
% rate in Hz. SampledAcqTrigger specifies how the acquisition is triggered.
% If "softwareImmediate" the acquisitiion will be software trigger by VSX
% prior to the Vantage sequence.  "softwareEvent" allows the acquisition to
% be software triggered within the Vantage event sequence. "hardwareSingle"
% will start the acquisition on detection of a hardware trigger pulse,
% subsequent triggers will be ignored. "hardwareRetriggerable" will allow
% the acquisition to be restarted after completion by subsequent hardware
% triggers.
%
% @Example:
%   There are three different ways of creating an analog input settings object     
%
%   % create an analog setting that defines a single single-ended, on demand-input with 
%   %range of 10V.
%   AISetting = vsv.extio.AnalogInputSetting(); 
%
%   create an analog input setting with three single single-ended, on demand-inputs 
%   AISetting = vsv.extio.AnalogInputSetting(3); 
%
%   create analog input setting for sampled inputs with hardware trigger
%    AISetting = vsv.extio.AnalogInputSetting('NumInput', 3, 'AcqType', "sampled", ...
%                                              'SampledAcqTrigger', "hardwareSingle", ...
%                                              'Samples', 100, 'SampleRate', 100); 
%
% Version 1.0 | 2020-06-16
% $Author: Dr. Jack Potter, Dr. Daniel Rohrbach
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
%
    
    properties
        
        % number of analog input channels @type int32 must be > 0
        NumInput (1,1) int32 {mustBeGreaterThan(NumInput, 0)} = 1;
        
        % mode which can be differential or single-ended
        Mode (1,1) string {mustBeMember(Mode, ["differential", "single-ended" ]) } ...
            = "single-ended";
        
        % Range setting
        Range (1,1) double {mustBeGreaterThan(Range, 0), mustBeLessThanOrEqual(Range, 10)} = 10;
        
        % acquisition type can be on-demand or sampled
        AcqType (1,1) vsv.extio.const.AcqTypeValue = "OnDemand"
        
        % sampled acquisition trigger can be "software" or "hardware"
        SampledAcqTrigger (1,1) vsv.extio.const.AcqTriggerValue = "softwareImmediate"
        
        % Number of samples in AI Acquisition
        Samples (1,1) int32 {mustBeGreaterThan(Samples, 1)} = 100;
        
        % Sample rate fn AI Acquisition
        SampleRate (1,1) int32 {mustBeGreaterThan(SampleRate, 0)} = 100;
        
    end
    
    %% setter and getter
    methods
        
        function this = AnalogInputSetting( varargin)
            % Create an instance of analog input setting
            %
            % @Usage
            %   % create an analog setting that defines a single single-ended, on demand-input with 
            %   %range of 10V.
            %   AISetting = vsv.extio.AnalogInputSetting(); 
            %
            %   create an analog input setting with three single single-ended, on demand-inputs 
            %   AISetting = vsv.extio.AnalogInputSetting(3); 
            %
            %   create analog input setting for sampled inputs with hardware trigger
            %    AISetting = vsv.extio.AnalogInputSetting('NumInput', 3, 'AcqType', "sampled", ...
            %                                              'SampledAcqTrigger', "hardwareSingle", ...
            %                                              'Samples', 100, 'SampleRate', 100); 

            if nargin == 1 && isnumeric( varargin{1} )
                varargin = { 'NumInput', varargin{1} };
            end
            
            this@vsv.extio.AbstractInputSetting( varargin{:} );
        end
        
    end
  
    methods(Access=protected, Sealed)
        
        function varargin = preInit(~, varargin)
        end
    end
end

