classdef (Abstract) OutputSetting < vsv.extio.mixin.ParamImportable
% OUTPUTSETTING defines the settings for a extended IO digital outputs
%   
%
% @Example:
%   % creates two digital outputs
%   es = vsv.extio.OutputSetting( 'NumOutputs', 2, ...
%                                  'SampleRate', 100, 'Waveform', ones(2, 100), 'Repeat', "single"); 

%    
% Version 1.0 | 2020-06-16 
% $Author: Dr. Jack Potter, Dr. Daniel Rohrbach 
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
%    
    properties(SetAccess=private)
        
        % number of encoders must be a positive scalar integer value 
        % note if a double value is set this will be converted to int32
        % This value cannot be changed after the user created the settings
        % object. NumEncoder must be specified during Setting Creation
        NumOutputs (1,1) int32 {mustBeGreaterThan(NumOutputs, 0), mustBeLessThanOrEqual(NumOutputs, 8)} = 1;
        
    end
    
    properties
        
        % Sample rate of output
        SampleRate (1, 1) int32 {mustBeGreaterThan(SampleRate, 0)} = 100;
        
        
        % @type string 
        %ouput start trigger can be "software" or "hardware"
        Trigger (1,1) vsv.extio.const.AcqTriggerValue = "softwareImmediate";
        
        % @type string 
        % output is "single", "continuous"
        Repeat (1,1) vsv.extio.const.RepeatValue = "continuous";
    end
    
   
    %% setter and getter
    methods(Sealed)
        
        function this = OutputSetting( varargin)
        % Create an instance of analogue output setting
        %
        % @Usage
        %   % creates two digital outputs
        %   es = vsv.extio.DigitalOutputSetting( 'NumOutputs', 2, ...
        %                                  'SampleRate', 100, 'Waveform', ones(2, 100), 'Repeat', "single"); 
        
            
            if nargin == 1 && isnumeric( varargin{1} )
                varargin = { 'NumOutputs', varargin{1} };
            end
            
            this@vsv.extio.mixin.ParamImportable( varargin{:} );
        end
    end
    
    %% getter and setter
    methods
        
        function set.NumOutputs(this, value)
        % setter for number of encoder values
        %
        
            oldV = this.NumOutputs;
            this.NumOutputs = value;
            try
                % make sure the values are updated correctly
                this.handleNumOutputsChanged();
            catch err
                % make sure if there is an error that we restore the value.
                % this should never happen if the other parameter are
                % always correct
                this.NumOutputs = oldV;
                rethrow(err);
            end
            
        end
         
    end
    
    
    %% helper functions 
    methods(Sealed)
        
        function isv = isValidSize( this, array )
        % Determines whether the given element array has valid size
        %  
        % an array has valid size if its first dimension fits the number
        % of elements defined by this settings
        %
        % @param array - @type numeric, the array to check
        % @return isv - @type logical true if size(array, 1) ==
        %               this.NumOutputs
        
            isv = size(array, 1) == this.NumOutputs;
            
        end
        
    end
    
    methods(Access=protected, Abstract)
        handleNumOutputsChanged(this);
    end
    
    methods(Access=protected)
                
        function errInvalidSize(~, varName)
            error( 'OutputSetting:invalidLength', ...
                [ 'Error number of elements in ' varName ' does not match the number of outputs of this DigitalOutputSetting' ])
        end
        
    end
    
    methods(Access=protected, Sealed)
        
        function varargin = preInit(this, varargin)
            
            if nargin == 2 && isstruct(varargin{1})
                str = varargin{1};
                if isfield(str, 'NumOutputs' )
                    this.NumOutputs = str.NumOutputs;
                    str = rmfield(str, 'NumOutputs' );
                    varargin = {str}; 
                end
            else
                ind = find( strcmp('NumOutputs', varargin) );
                if ~isempty(ind)
                    if ind+1 <= numel(varargin)
                        this.NumOutputs = varargin{ind+1};
                    else
                        error('OutputSetting:invalidvalueparam', ...
                            'Input for output settings must be value param pairs. Cannot set num outputs');
                    end
                    varargin( ind:ind+1 ) = []; 
                end
            end
            
        end
        
    end
end


