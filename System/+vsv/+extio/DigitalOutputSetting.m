classdef (Sealed) DigitalOutputSetting < vsv.extio.OutputSetting
%DIGITALOUTPUTSETTING defines the settings for a extended IO digital outputs
%   
% 'NumOutputs' specifies the number of digital output channels.
% 'Samplerate' is the rate in Hz that output samples will be written.
% 'Trigger' specifies how the acquisition is triggered. If
% "softwareImmediate" the output will be software triggered by VSX prior to
% the Vantage sequence.  "softwareEvent" allows the output to be software
% triggered within the Vantage event sequence. "hardwareSingle" will start
% the output on detection of a hardware trigger pulse, subsequent triggers
% will be ignored. "hardwareRetriggerable" will allow the output to be
% restarted after completion by subsequent hardware triggers
% ("hardwareRetriggerable" cannot be used for "continuous" outputs).
%
% The output digital waveforms are specified by the 'Waveform' field.  This
% should be a binary array with first dimension equal to the number of
% outputs and the second being the number of output samples.  The 'Repeat'
% setting may be "single" which will output the specified waveform once or
% "continuous" will repeat the waveform continuously.

% @Example:
%   % creates two digital outputs. ds = vsv.extio.DigitalOutputSetting(
%                                  'NumOutputs', 2, ...
%                                  'SampleRate', 100, ...
%                                  'Waveform', ones(2, 100), ...
%                                  'Repeat',  "single", 'Trigger', "hardwareSingle");

%    
% Version 1.0 | 2020-06-16 
% $Author: Dr. Jack Potter, Dr. Daniel Rohrbach 
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
%    
    
    properties
        
        % @type int32 
        % Output waveforms, first dimension should equal NumOutputs
        Waveform (:, :) int32 {mustBeReal, mustBeGreaterThanOrEqual(Waveform, 0)} = zeros(1,100);
          
    end
    
    methods(Static)
        
        function this = loadobj(this)
            if isstruct(this)
                try
                    this = vsv.extio.DigitalOutputSetting( this );                
                catch
                    warning( 'DigitalOutputSetting:loadobj:invalidData', ...
                        'Unable to read DigitalOutput from file, invalid data');
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

            this.Waveform = repmat( defaultWaveform(), this.NumOutputs, 1 );

        end
    end
    
end

function defwf = defaultWaveform()
%
% @return defwf - @type double the default value zeros(1,100);

    defwf = zeros(1,100, 'int32');

end

