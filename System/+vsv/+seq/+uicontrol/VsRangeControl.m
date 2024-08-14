classdef VsRangeControl <   vsv.seq.uicontrol.UserUiControl  ...
                    & vsv.seq.uicontrol.mixin.HasLabel ...
                    & vsv.ctrl.mixin.CreatePropertyControlSetting ...
                    & vsv.events.EnableEventsMixin ...
                    & vsv.seq.storage.HasStorageParameter ...
%VSRANGECONTROl describes a range control that can be defined by the user 
%   
%   A range control defines a value that can be set within a certain range.
%
% Attributes are RangeMin, RangeMax. Constrains on these attributes are 
% RangeMin < RangeMax.  Use setRange function to set min and max values at
% the same time. 
%
% The slider state is represented by the CurrentValue which is a numeric
% value, which is > RangeMin and < RangeMax. Setting a CurrentValue outside
% of the range will cause an error.
%
% The class also provides some functions to validate whether a value is
% within the acceptable range:
%
%   % create a range control
%   range = vsv.seq.uicontrol.VsRangeControl();
%
%   % display range min, max, and  value
%   range.RangeMin
%   range.RangeMax
%   % will display [1 100]
%   range.getRange()
%   % will display true
%   range.isValueInLimits( 2)
%   % will display false
%   range.isValueInLimits( 0)
% 
% Version 1.0 | 2020-01-24 
% $Author: Dr. Daniel Rohrbach
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
%

    properties(SetAccess=private)
        
        % @type numeric range minimum value, setting the current value 
        RangeMin    = 1;
        
        % @type numeric three dimensional array
        RangeMax    = 100;
        
        % @type char a format string
        ValueFormat = '%3.0f';
        
    end
    
    properties
       
        % @type function handle, a user defined function that will be
        % executed before setting the current value property. e.g. @round
        ValueFormatFnc function_handle = function_handle.empty();
    end
    
    
    %% event indicating value changed
    events
        CurrentValueChanged;
    end
    
    %% Current value property that can change after creation
    properties(SetAccess=private)
        
        % @type this is the current value
        CurrentValue = 1;
        
    end
    
    %% public methods
    methods
        
        function this = VsRangeControl( varargin )
        % constructor of the range control    
        %
        % @param rangeValue - @type numeric the range initial value
        % @param UI - @type the UI control structure
        
            
            this@vsv.seq.uicontrol.UserUiControl(varargin{:});
            this.setStyle( 'RangeControl' );
        end
        
        % @override vsv.seq.storage.HasStorageParameter
        function ID = getStorageID(this)
            ID = this.Label;
        end
        
                
        function strUI = convertToStructUIControl(this)
        % converts this ui control to the old fashion Struct representation
        % of this UI control
        %
        % @return strUI - the struct represention of this UI control
        
            strUI.Control  = {this.LocationCode, ...
                              'Style',          this.Style, ...
                              'Label',          this.Label, ...
                              'RangeMinMaxVal', [ this.RangeMin this.RangeMax this.CurrentValue], ...
                              'ValueFormat',    this.ValueFormat};
            strUI.Callback = this.Callback;
            if ~vsv.util.isNULL( this.Statement )
                strUI.Statement = this.Statement;
            end
        end
        
        function isin = isValueInLimits(this, val)
        % test whether given value is within the ranges
        %   
        %  A value is in range if its larger or equal then RangeMin and
        %  smaller or equal then RangeMax
        % @param val - @type numeric the value to compare
        
            isin = val >= this.RangeMin && val <= this.RangeMax;
            
        end
        
        function success = setCurrentValueCallback( this, val)
        % will set the current value of the range control and exec. Callback
        %
        % This will set the current value of the range control using the
        % given value val. the function returns true if setting was
        % successful
        % 
        % The function also informs listener, listeners will called before
        % executing the callback
        % 
        % @param val - @type numeric, the current value
        % @return val - @type numeric, @type logical, true if successful
        %               false otherwise
        
            ov = this.CurrentValue;
            success = this.setCurrentValue(val);
            % we cannot put this object as the Source!!!!! the first
            % argument must be [], there are scripts that expect the first
            % argument to be a graphics handle which will cause an error if
            % the graphcis handle is actually a range control
            this.executeCallbackControlSave( [], 'CurrentValue', ov);
            
        end
        
        function success = setCurrentValue( this, val)
        % will set the current value of the range control and NOT exec. Callback
        %
        % This will set the current value of the range control using the
        % given value val. the function returns true if setting was
        % successful
        % 
        % The function also informs listener
        % 
        % @param val - @type numeric, the current value
        % @return val - @type numeric, @type logical, true if successful
        %               false otherwise
        
            success = false; 
            val = this.applyValueFormatFnc(val);
            
            if this.isValueInLimits(val)
                this.CurrentValue = val;
                this.notifyEnabled( 'CurrentValueChanged' );
                success = true;
            end
            
        end
        
    end
   
       
    %% setters and getters 
    methods
        
        function set.ValueFormatFnc(this, fnc)
            this.ValueFormatFnc = fnc;
            
            if ~isempty(this.ValueFormatFnc)
                % make sure its getting applied
                this.setCurrentValue(this.CurrentValue);%#ok
            end
        end
        
        function range = getRange(this)
        % returns the range, i.e., [ this.RangeMin, this.RangeMax ]
        % @return range - @type numeric, 2 element
        %
            range = [ this.RangeMin, this.RangeMax ];
        end
        
        function setRange(this, range)
        % set the range value for this range control
        %
        % @param range - @type numeric the range, should be a 2 element
        %               array
        
            if isnumeric(range) && numel(range) == 2
                if range(1) <= range(2)

                    this.RangeMin = range(1);
                    this.RangeMax = range(2);

                    if this.CurrentValue > this.RangeMax 
                        this.CurrentValue = this.RangeMax;
                    end
                    
                    if this.CurrentValue < this.RangeMin
                        this.CurrentValue = this.RangeMin;
                    end
                    
                else
                    error('VsSlider:setRange:invalidRange', ...
                        'Min range is larger then Max range');
                end
            else
                error('VsSlider:setRange:invalidArgument', ...
                        'range must be a numeric 2 element array');
            end
        end
        
        function setMinMaxVal(this, mmv)
        % set the min and maximum value mmv
        %
        %@param mmv - @type numeric the min max value, should be a 3
        %             element array
        
            if length(mmv) >= 3
                
                rangeMin = this.RangeMin;
                rangeMax = this.RangeMax;
                
                this.setRange( mmv(1:2) );
                success = this.setCurrentValue( mmv(3));
                if ~success
                    % RangeMin old state
                    this.RangeMin = rangeMin;
                    this.RangeMax = rangeMax;
                    error( 'VsSlider:setSliderMinMaxVal:invalidArgument' ,...
                        'Value is outside of limits');
                end
               
            else
                error( 'VsSlider:setSliderMinMaxVal:invalidArgument' ,...
                    'sliderMinMaxVal must be a three dimensional numeric array');
            end
        end
        
        function setRangeMin(this, rMin)
        % setter for the range min value
        % 
        %   this will throw an error if rMin is larger or equal range Max 
        %
        % @param rMin - @numeric must be the minimum and smaller then range
        %               max
        
            if rMin <= this.RangeMax
                this.RangeMin = rMin;
                if this.CurrentValue < tthis.RangeMin
                    this.setCurrentValue( this.RangeMin );
                end
            else
                error('VsSlider:setRangeMin:invalidRange', ...
                    'Min range is larger then Max range');
            end
            
        end
        
        function setRangeMax(this, rMax)
        % setter for the range max value
        % 
        %   this will throw an error if rMax is smaller or equal RangeMin 
        %
        % @param rMin - @numeric must be the minimum and smaller then range
        %               max
        
            if rMax >= this.RangeMin
                this.RangeMax = rMax;
                if this.CurrentValue > this.RangeMax
                    this.setCurrentValue( this.RangeMax );
                end
            else
                error('VsSlider:setRangeMin:invalidRange', ...
                    'Min range is larger then Max range');
            end
            
        end
        
        
        function setValueFormat(this, vf)
        % value display format
        %
        %@param vf - @type char the value format
            if isstring(vf)
                vf = convertStringsToChars(vf);
            end
            
            if ischar(vf)
                this.ValueFormat = vf;
            else
                error('VsSlider:setValueFormat:invalidArgument' ,...
                    'ValueFormat must be a char');
            end
                
        end
        
        function isvalid = isValidStyle(~, style)
            isvalid = ~isempty(style) && ( strcmp( style, 'VsSlider' ) || strcmp( style, 'RangeControl' ) );  
        end
        
        
        function settings = getPropertyControlSettings(this)
            settings = this.createPropertyControlSetting('CurrentValue');
        end
        
    end
    
    methods(Access=protected)
        
        function val =  applyValueFormatFnc(this, val)
            
            if ~isempty( this.ValueFormatFnc )
                
                try
                    val = this.ValueFormatFnc(val);
                catch err
                    
                    excp = MException( 'VsRangeControl:applyValueFormatFnc:errorFnc', ...
                                       ' The user defined value format function "ValueFormatFnc" of the VsRangeControl errored! ' );
                    excp.addCause(err);
                    excp.throw();
                end
                
            end
            
        end
        
        % @overriden vsv.seq.storage.HasStorageParameter
        function [paramName, values] = getStorageParameter(this)
        
            paramName = {'CurrentValue'};
            values    = {this.CurrentValue};
        end
        
        % @overriden vsv.seq.storage.HasStorageParameter
        function Class = getSupportedClass(~)
            Class = 'vsv.seq.uicontrol.VsRangeControl'; 
        end
        
        % @overriden vsv.seq.storage.HasStorageParameter
        function success = importSingleStorageParameter(this, parameterName, value) 
            
            switch parameterName
                case 'CurrentValue'
                    % only import if value are different, this save some
                    % hick ups
                    if this.CurrentValue ~= value
                        success = this.setCurrentValueCallback(value);
                    else
                        success = true;
                    end
                otherwise
                    error('Unable to set the parameter value');
            end
            
        end
        
        
        function success = applyParamValue(this, param, value)
            
            success = true;
            switch param
                case 'Range'
                    this.setRange( value );
                case 'RangeMin'
                    this.setRangeMin( value );
                case 'RangeMax'
                    this.setRangeMax( value );
                case 'MinMaxVal'
                    this.setMinMaxVal(value);
                case 'ValueFormat'
                    this.setValueFormat(value);
                case 'Label'
                    this.setLabel(value);
                case 'ValueFormatFnc'
                    this.ValueFormatFnc = value;
                otherwise 
                    success = this.applyParamValue@vsv.seq.uicontrol.UserUiControl( param, value );
            end
        end
    end
    
    methods(Access ={?vsv.ctrl.mixin.CreatePropertyControlSetting, ...
                               ?vsv.ctrl.AbstractPropertyController})
        
        function setting = createPropertyControlSetting(this, property)
        % create a vsv.ctrl.PropertyControllerSettings for the given
        % property name
        %
        %   the default setting to describe access rights for a given
        %   property is direct property access with post set properties. In
        %   case the inheriting class uses different getters and setters
        %   this can be used to create the settings object that best
        %   describes the access rights for the given property
        %
        % @param property - @type char the name of the property
        % @return setting - @type vsv.ctrl.PropertyControllerSettings
            
            
            persistent settingTemplate;
            
            if strcmp(property, 'CurrentValue') 
                
                if isempty(settingTemplate)
                % the settings object
                    settingTemplate = vsv.ctrl.PropertyControllerSettings(...
                                    'PropertyName',        property , ...
                                    'PropertyType',        'numeric', ...
                                    'PropertySetter',      'setCurrentValue', ...
                                    'PropertyLimits',      [ this.RangeMin this.RangeMax ], ...
                                    'PropertyChangeEvent', 'CurrentValueChanged');
                end
                setting = settingTemplate.newCopy();
                setting.setPropertyLimits( [ this.RangeMin this.RangeMax ] );
                setting.setPropertyName( property );
                
                if ~isempty( this.ValueFormat )
                    format = vsv.ctrl.format.EditNumberFormat( this.ValueFormat );
                    setting.setPropFormatter( format );
                end
                
            else
                error('VsRangeControl:createPropertyControlSetting:unknownProperty', ...
                    'Given Property is unknown');
            end
        end
        
    end

    
end