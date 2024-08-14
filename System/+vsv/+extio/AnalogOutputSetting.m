classdef (Sealed) AnalogOutputSetting < vsv.extio.OutputSetting
%ANALOGOUTPUTSETTING defines the settings for a extended IO analog outputs
%   
% 'NumOutputs' specifies the number of analog output channels.  'Range'
% specifies the maximum peak output voltage. 'Samplerate' is the rate in Hz
% that output samples will be written. 'Trigger' specifies how the
% acquisition is triggered. If "softwareImmediate" the output will be
% software triggered by VSX prior to the Vantage sequence.  "softwareEvent"
% allows the output to be software triggered within the Vantage event
% sequence. "hardwareSingle" will start the output on detection of a
% hardware trigger pulse, subsequent triggers will be ignored.
% "hardwareRetriggerable" will allow the output to be restarted after
% completion by subsequent hardware triggers ("hardwareRetriggerable"
% cannot be used for "continuous" outputs).
%
% The output waveforms are specified by the 'Waveform' field (the maximum
% value should not exceed the 'Range' value).  This should be a numeric
% array with first dimension equal to the number of outputs and the second
% being the number of output samples.  The 'Repeat' setting may be "single"
% which will output the specified waveform once or "continuous" which will
% repeat the waveform continuously.
%
% @Example:
%
%   % creates two analog outputs with single repeat and hardware trigger
%   as = vsv.extio.AnalogOutputSetting( 'NumOutputs', 2, 'Range', 10, ...
%                                  'SampleRate', 100, 'Waveform', rand(2,
%                                  100), 'Repeat', "single", 'Trigger',
%                                  "hardwareSingle");
%    
% Version 1.0 | 2020-06-16 
% $Author: Dr. Jack Potter, Dr. Daniel Rohrbach 
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
%    
        
    properties
        
        % @type double 
        % Maximum output voltage
        Range (1, 1) double {mustBeGreaterThan(Range, 0)} = 10;
        
        % @type double 
        % Output waveforms, first dimension should equal NumOutputs
        Waveform (:, :) double {mustBeReal} = zeros(1,100);
    
    end
    
    methods(Static)
        
        function this = loadobj(this)
            if isstruct(this)
                try
                    this = vsv.extio.AnalogOutputSetting( this );                
                catch
                    warning( 'AnalogOutputSetting:loadobj:invalidData', ...
                        'Unable to read AnalogOutput from file, invalid data');
                end
                    
            end
        end
        
    end
    
    methods
        
        function set.Waveform(this, value)
        % setter for Waveform

            if this.isValidSize(value)
                this.Waveform = value;
            else
                this.errInvalidSize('Waveform');
            end
        end
        
    end
    
    methods(Access=protected)
        
        function handleNumOutputsChanged(this)
        % create a new wave form if number of outputs changed
        % 
        
            this.Waveform = repmat( defaultWaveform(), this.NumOutputs, 1 );

        end
    end
    
end

function defwf = defaultWaveform()
% returns the default value for CountPerPulse 
%
% @return defwf - @type double the default value zeros(1,100);

    defwf = zeros(1, 100);

end

