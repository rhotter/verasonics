classdef DigitalInputSetting< vsv.extio.AbstractInputSetting
%DIGITALINPUTSETTING defines the settings for a set of digital inputs
%The number of encoders is defined by NumInput.  No other settings
%required.
%
% @Examples:
%
%   % create an digital input setting that defines a single input
%   es = vsv.extio.DigitalInputSetting();
%
%   % creates an digital input setting for defining 3 inputs
%   es = vsv.extio.DigitalInputSetting(3);
%
%   % creates an digital input setting for defining 4 inputs
%   es = vsv.extio.DigitalInputSetting('NumInput', 4);
%
%
% Version 1.0 | 2020-06-16
% $Author: Dr. Jack Potter, Dr. Daniel Rohrbach
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
%
    
    properties
        
        NumInput (1,1) int32 {mustBeGreaterThan(NumInput, 0), mustBeLessThanOrEqual(NumInput, 8)} = 1; %Number of digital input lines
        
    end
    
    
    methods(Static)
        
        function defnumin = defaultNumInput()
            % returns the default value for NumInput
            % @return defcp - @type int32 the default NumInput = 1
            defnumin = 1;
        end
        
    end
    
    methods
        
        function this = DigitalInputSetting( varargin)
            % Create an instance of digital input setting
            %
            % @Usage
            %   % create an digital input setting that defines a single input
            %   es = vsv.extio.DigitalInputSetting();
            %
            %   % creates an digital input setting for defining 3 inputs
            %   es = vsv.extio.DigitalInputSetting(3);
            %
            %   % creates an digital input setting for defining 4 inputs
            %   es = vsv.extio.EncoderSetting( 'NumInput', 4);
            
            
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

