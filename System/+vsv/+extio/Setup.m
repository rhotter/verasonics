classdef Setup < handle
%Setup ExtIO setup structure defining settings for inputs and outputs 
%   
%
%
% Version 1.0 | 2020-11-04 
% $Author: Dr. Jack Potter, Dr. Daniel Rohrbach 
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
%
   
    properties
        
        % encoder settings as defined in the script
        Encoder       vsv.extio.EncoderSetting       {mustBeScalarOrEmpty};
        
        % analog input settings as defined in the script
        AnalogInput   vsv.extio.AnalogInputSetting   {mustBeScalarOrEmpty};
         
        % digital input settings as defined in the script
        DigitalInput  vsv.extio.DigitalInputSetting  {mustBeScalarOrEmpty};
         
        % analog output settings as defined in the script
        AnalogOutput  vsv.extio.AnalogOutputSetting  {mustBeScalarOrEmpty};
         
        % digital output settings as defined in the script
        DigitalOutput vsv.extio.DigitalOutputSetting {mustBeScalarOrEmpty};
        
    end

    methods
        
        function createEncoder(this, varargin)
        % creats a new encoder settings 
        %
        % @Example:
        %   setup = vsv.extio.Setup;
        %   There are three different ways of creating an encoder settings object     
        %
        %   % create an encoder setting that defines a single encoder with 
        %   % X1 as Decoding and Zindex = 0
        %   setup.newEncoder();
        %
        %   % creates an encoder setting for defining 3 encoders  
        %   setup.newEncoder(3);
        %   
        %   % creates an encoder setting for defining 4 encoders by using
        %   % parameter value pairs
        %   setup.newEncoder( 'NumEncoder', 4, ...
        %                     'Decoding', ["X1" "X2", "X1", "X4"], ...
        %                     'Zindex', [100 10 1 0], 'ZIndexPhase', ["AHighBLow" "AHighBHigh" "AHighBHigh" "AHighBHigh"]); 
        %
        %   % the user can modify parameters after creations 
        %   setup.Encoder.Decoding(3) = 'X1';  
        %   setup.Encoder.Zindex = ones(1, es.NumEncoder) .* 12;
        %   
            this.Encoder = vsv.extio.EncoderSetting(varargin{:});
        end
                
        function createAnalogInput(this, varargin)
            this.AnalogInput = vsv.extio.AnalogInputSetting(varargin{:});
        end
                
        function createDigitalInput(this, varargin)
            this.DigitalInput = vsv.extio.DigitalInputSetting(varargin{:});
        end

        function createAnalogOutput(this, varargin)
            this.AnalogOutput = vsv.extio.AnalogOutputSetting(varargin{:});
        end
        
        function createDigitalOutput(this, varargin)
            this.DigitalOutput = vsv.extio.DigitalOutputSetting(varargin{:});
        end
        
        
    end
    
end



function mustBeScalarOrEmpty(val)
    vsv.util.mustBeScalarOrEmpty(val);
end