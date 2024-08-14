classdef (Sealed) EncoderSetting < vsv.extio.AbstractInputSetting
%ENCODERSETTING defines the settings for a set of encoders 
%   
%   The number of encoders is defined by NumEncoder. For each encoder the
%   user can specify several parameters such as Decoding, Zindex and
%   other parameters.
%
%   The parameters that are specific for each encoder must be a vector that
%   has the same number of elements as EncoderSetting.NumEncoder. The
%   EncoderSetting class guarantees that encoder specific parameters 
%   (e.g., Decoding and Zindex) are vectors with a number of elements
%   as this.NumEncoder. The class must be initialized with the number of 
%   encoders (if the user want to specify  more than one encoder). After
%   the encoder setting was generated the NumEncoder value cannot be modfied.
%
%   The parameter Decoding specfies the pulse decoding method which may be
%   "X1", "X2" or "X4". The default value is "X1" if unspecified.
%   Zindex is the value in encoder ticks to which to reset the measurement 
%   when signal Z is high.  If omitted or zero, the z index is disabled.
%
%   ZIndexPhase specifies the states at which signal A and signal B must be
%   while signal Z is high to reset the measurement. This may be
%   "AHighBHigh", "AHighBLow", "ALowBHigh" or "ALowBLow".  The default is
%   "AHighBHigh" if unspecified.
%   
% @Example:
%   There are three different ways of creating an encoder settings object     
%
%   % create an encoder setting that defines a single encoder with 
%   % X1 as Decoding and Zindex = 0
%   es = vsv.extio.EncoderSetting(); 
%
%   % creates an encoder setting for defining 3 encoders  
%   es = vsv.extio.EncoderSetting(3); 
%
%   % creates an encoder setting for defining 4 encoders
%   es = vsv.extio.EncoderSetting( 'NumEncoder', 4, 'Decoding', ["X1" "X2", "X1", "X4"], ...
%                                  'Zindex', [100 10 1 0], 'ZIndexPhase', ["AHighBLow" "AHighBHigh" "AHighBHigh" "AHighBHigh"]); 
%
%   % the user can modify parameters after creations 
%   es.Decoding(3) = 'X1';  
%   es.Zindex = ones(1, es.NumEncoder) .* 12;
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
        NumEncoder (1,1) int32 {mustBeGreaterThan(NumEncoder, 0), mustBeLessThanOrEqual(NumEncoder, 4)} = 1;
        
    end
    
    properties
        
        % @type vsv.extio.const.CountPerPulseValue which can be X1, X2, X4
        % the length of countPerPulse must match 
        Decoding (1, :) vsv.extio.const.CountPerPulseValue ....
                          = vsv.extio.EncoderSetting.defaultDecoding();
        
        % @type double must be of same length as NumEncoder 
        Zindex (1, :) double {mustBeGreaterThanOrEqual(Zindex, 0)} ...
                          = vsv.extio.EncoderSetting.defaultZindex();
                      
        % @type vsv.extio.const.ZIndexPhase which can be AHighBHigh,
        % AHighBLow, ALowBHigh, ALowBLow
        % must be of same length as NumEncoder
        ZIndexPhase (1, :) vsv.extio.const.ZIndexPhase ....
                          = vsv.extio.EncoderSetting.defaultZIndexPhase();
                      
        % buffer settings are part of all ExtIO settings 
        % @type vsv.extio.BufferSetting
        PopSetting vsv.extio.PopSetting ...
                        {vsv.util.mustBeScalarOrEmpty(PopSetting)} ...
                            = vsv.extio.PopSetting.empty();             
        
    end
    
    

    methods(Static)
        
        function defdc = defaultDecoding()
        % returns the default value for CountPerPulse 
        % @return defcp - @type string the default CountPerPulse value = X1
            defdc = vsv.extio.const.CountPerPulseValue.X1;
        end
        
        
        function defz = defaultZindex()
        % returns the default value for Zindex 
        % @return defz - @type double the default zindex == 0
            defz = 0;
        end
        
        function defzip = defaultZIndexPhase()
        % returns the default value for ZIndexPhase 
        % @return defzip - @type string the default ZIndexPhase value = AHighBHigh
            defzip = vsv.extio.const.ZIndexPhase.AHighBHigh;
        end
        
    end
    
    %% setter and getter
    methods
        
        function this = EncoderSetting( varargin)
        % Create an instance of encoder setting
        %
        % @Usage
        %   % create an encoder setting that defines a single encoder with 
        %   % X1 as CountPerPulse and Zindex = 0
        %   es = vsv.extio.EncoderSetting(); 
        %
        %   % creates an encoder setting for defining 3 encoders  
        %   es = vsv.extio.EncoderSetting(3); 
        %
        %   % creates an encoder setting for defining 4 encoders
        %   es = vsv.extio.EncoderSetting( 'NumEncoder', 4, 'Decoding', ["X1" "X2", "X1", "X4"], ...
        %                                  'Zindex', [1 0 1 0], 'ZIndexPhase', ["AHighBLow" "AHighBHigh" "AHighBHigh" "AHighBHigh"]); 
        
        
            
            if nargin == 1 && isnumeric( varargin{1} )
                varargin = { 'NumEncoder', varargin{1} };
            end
            
            this@vsv.extio.AbstractInputSetting( varargin{:} );
            
            this.PopSetting = vsv.extio.PopSetting();

        end
        
        
        function set.NumEncoder(this, value)
        % setter for number of encoder values
        %
        
            oldV = this.NumEncoder;
            this.NumEncoder = value;
            try
                % make sure the values are updated correctly
                this.handleEncoderNumChanged();
            catch err
                % make sure if there is an error that we restore the value.
                % this should never happen if the other parameter are
                % always correct
                this.NumEncoder = oldV;
                rethrow(err);
            end
            
        end
        
        function set.Decoding(this, value)
        % setter for count per pulse
        %
        
            if this.isValidNumElements(value)
                this.Decoding = value;
            else
                this.errInvlaidLength('Decoding');
            end
        end
        
        function set.Zindex(this, value)
        % setter for count per pulse
        %
        
            if this.isValidNumElements(value)
                this.Zindex = value;
            else
                this.errInvlaidLength('Zindex');
            end
        end
        
        function set.ZIndexPhase(this, value)
            % setter for Z index phase
            %
            if this.isValidNumElements(value)
                this.ZIndexPhase = value;
            else
                this.errInvlaidLength('ZIndexPhase');
            end
            
        end
        
    end
    
    
    methods
        
        function isv = isValidNumElements( this, vector )
        % Determines whether the given element vector has valid size
        %  
        % a vector has valid size if its number of elements fits the number
        % of elements defined by this settings
        %
        % @param vector - @type numeric, the vector to check
        % @return isv - @type logical true if numel(vector) ==
        %               this.NumEncoder
        
            isv = numel(vector) == this.NumEncoder;
            
        end
        
        
    end
    
    methods(Access=protected, Sealed)
        
        function varargin = preInit(this, varargin)
            
            if nargin == 2 && isstruct(varargin{1})
                str = varargin{1};
                if isfield(str, 'NumEncoder' )
                    this.NumEncoder = str.NumEncoder;
                    str = rmfield(str, 'NumEncoder' );
                    varargin = {str}; 
                end
            else
                ind = find( strcmp('NumEncoder', varargin) );
                if ~isempty(ind)
                    if ind+1 <= numel(varargin)
                        this.NumEncoder = varargin{ind+1};
                    else
                        error('EncoderSetting:invalidvalueparam', ...
                            'Input for encoder settings must be value param pairs. Cannot set num outputs');
                    end
                    varargin( ind:ind+1 ) = []; 
                end
            end
            
        end
        
    end
    
    
    methods(Access=private)
        
        
        
        function handleEncoderNumChanged(this)
        
            import vsv.extio.EncoderSetting;
            
            defDecoding  = EncoderSetting.defaultDecoding();
            defZindex = EncoderSetting.defaultZindex();
            defZindexPhase =  EncoderSetting.defaultZIndexPhase();

            this.Decoding = this.updateArrayParameter( this.Decoding, defDecoding);
            this.Zindex        = this.updateArrayParameter( this.Zindex, defZindex);
            this.ZIndexPhase   = this.updateArrayParameter( this.ZIndexPhase, defZindexPhase);

        end
        
        function narray = updateArrayParameter(this, arrayParam, defParam)
        % if num encoder changes this parameter needs to be updated
        %
        
            nEncoder = this.NumEncoder;
            nCount   = numel(arrayParam);
            
            if nEncoder ~= nCount
                
                if nCount > nEncoder
                    
                    % old array fits into narray 
                    narray(1, 1:nEncoder) = arrayParam(1:nEncoder); 
                    
                else % nCount < nEncoder cannot be equal
                    
                    narray(1, 1:nCount) = arrayParam(1:nCount); 
                    narray(1, nCount+1:nEncoder) = defParam;
                end
            else
                narray = arrayParam;
            end
        
        end
        
        
        function errInvlaidLength(~, varName)
            error( 'EncoderSetting:invalidLength', ...
                [ 'Error number of elements in ' varName ' does not match the number of elements of this EncoderSettings' ])
        end
        
    end
    
end

