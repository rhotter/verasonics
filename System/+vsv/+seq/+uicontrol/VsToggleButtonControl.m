classdef VsToggleButtonControl <   vsv.seq.uicontrol.VsButtonControl ....
                                 & vsv.ctrl.mixin.CreatePropertyControlSetting ...
                                 & vsv.events.EnableEventsMixin ...
                                 & vsv.seq.storage.HasStorageParameter
%VsToggleButtonControl uicontrol defining a toggle button with a state
%   
% A Toggle button control inherits from a VsButton Control and extends it
% with a State field
%
% -----------
% Usage:
%
%       tb = vsv.seq.uicontrol.VsToggleButtonControl( ...
%                                     'LocationCode',   'UserC3', ...
%                                     'Style',          'VsToggleButton', ... 
%                                     'Label',          'HWFilter');
%       tb.State
%
% Version 1.0 | 2020-04-01 
% $Author: Dr. Daniel Rohrbach
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
%    

    properties( SetAccess=private )
        
        % the state of the toggle button @type logical
        State = false;
        
        % @type char, a lable that can be used for changing if the state
        % changes to true
        OffLabel = '';
        
    end
    
    %% event indicating value changed
    events
        StateChanged;
    end
    
    %% public methods
    methods
        
        function this = VsToggleButtonControl( varargin )
        % constructor of the VSButtonControl    
        %
        % @param UI - @type the UI control structure 
        
          
            this@vsv.seq.uicontrol.VsButtonControl( varargin{:} );
            this.setStyle( 'VsToggleButton' );
            
        end
        
        % @override vsv.seq.storage.HasStorageParameter
        function ID = getStorageID(this)
            ID = this.Label;
        end
        
        function setStateCallback(this, state)
        % set the sate of the toggle button control, will execute callback
        %
        %   the function supports numeric values in which case the state
        %   gets converted to logical. State must be a scalar value (i.e.,
        %   of length 1) and it must be not nan
        %
        %
        % @error VsToggleButtonControl:setState:invalidType - if state is
        %        invalid
        %
        % @param state - @type logical or numeric, if numeric the state is
        %               converted to logical
        
            ov = this.State;
            this.setState(state);
            this.executeCallback( this, 'State', ov);
            
        end
        
        function setOffLabel(this, label)
        % setter for off label 
        % @param label - @type char the off label string
            if ischar( label )
                this.OffLabel = label;
            else
                error('VsToggleButtonControl:setOffLabel:invalidType', ...
                    'off label must be a char');
            end
            
        end
        
        function setState(this, state)
        % set the sate of the toggle button control will NOT execute CallBk
        %
        %   the function supports numeric values in which case the state
        %   gets converted to logical. State must be a scalar value (i.e.,
        %   of length 1) and it must be not nan
        %
        % @error VsToggleButtonControl:setState:invalidType - if state is
        %        invalid
        %
        % @param state - @type logical or numeric, if numeric the state is
        %               converted to logical
            if isnumeric(state) && ~isnan(state)
                state = logical(state);
            end
            if islogical(state) && numel(state) == 1
                this.State = state;
                this.notifyEnabled( 'StateChanged' );
            else
                error('VsToggleButtonControl:setState:invalidType', ...
                    'Type of given state is invalid, must be a scalar logical or numeric and not nan');
            end
        end
        
        function isvalid = isValidStyle(~, style)
        % Test whether style is a valid state
        %
        % @overwritten
        % @return isvalid - @type logical, true if style == 'VsToggleButton'
            isvalid = ~isempty(style) && strcmp( style, 'VsToggleButton' ) ;  
        end
        
        function isvalid = isValidControl(this, UI)
        % check whether the UI control is valid
        %
        % A ui control should be valid if it is a isUserUiControl
        %
        % @overwritten  
        % @param UI - the ui control struct
        %
            isvalid = this.isValidControl@vsv.seq.uicontrol.UserUiControl( UI );
            if isvalid
                Style = vsv.seq.uicontrol.UserUiControl.getUIStyleProtected(UI);
                isvalid = ~isempty(Style) && strcmp( Style, 'VsToggleButton' );
            end
            
        end
        
        function settings = getPropertyControlSettings(this)
            settings = this.createPropertyControlSetting('State');
        end
        
    end
    
    methods(Access ={?vsv.ctrl.mixin.CreatePropertyControlSetting, ...
                               ?vsv.ctrl.AbstractPropertyController})
        
        function setting = createPropertyControlSetting(~, property)
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
            
            if strcmp(property, 'State') 
                
                if isempty(settingTemplate)
                % the settings object
                    settingTemplate = vsv.ctrl.PropertyControllerSettings(       ...
                                    'PropertyName',        property ,            ...
                                    'PropertyType',        'logical',            ...
                                    'PropertySetter',      'setState', ...
                                    'PropertyChangeEvent', 'StateChanged');
                end
                setting = settingTemplate.newCopy();                
            else
                error('VsToggleButtonControl:createPropertyControlSetting:unknownProperty', ...
                    'Given Property is unknown');
            end
        end
        
    end

    
    methods(Access=protected)
        
        % @overriden vsv.seq.storage.HasStorageParameter
        function [paramName, values] = getStorageParameter(this)
        
            paramName = {'State'};
            values    = {this.State};
            
        end
        
        % @overriden vsv.seq.storage.HasStorageParameter
        function success = importSingleStorageParameter(this, parameterName, value) 
            success = this.applyParamValue(parameterName, value);
        end
        
        function success = applyParamValue(this, param, value)
        % @overwritten
        %
            success = true;
            switch param
                case 'State'
                    this.setState(value);
                case 'Label'
                    if iscell(value) 
                        this.setLabel(value{1});
                        if length( value ) > 1
                            this.setOffLabel( value{2} );
                        end
                    else
                        this.setLabel(value)
                    end
                otherwise 
                    success = this.applyParamValue@vsv.seq.uicontrol.VsButtonControl( param, value );
            end
        end
        
    end
    
    
end

