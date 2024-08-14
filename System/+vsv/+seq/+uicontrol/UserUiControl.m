classdef (Abstract) UserUiControl <  vsv.seq.uicontrol.mixin.HasStatement ...
                                   & vsv.seq.uicontrol.mixin.HasCallback
% UserUiControl: A user control defines a control object in a vantage
% application  
%  
% This is the basic class for all UserUIControls such as
% VsSliderControls, VsButtons, VsButtonGroupControl, or
% VsToggleButtonControls
%
% 
% A user control consists of a Style a LocationCode and a callback. 
%
% To facilitate creating user controls on the Verasonics GUI, the Vantage 
% software provides some built-in UI functions that can be
% positioned at certain predefined locations. 
%
% The currently available controls that can be selected are:
% VsSliderControl, VsPushButton, VsToggleButton, and VsButtonGroup. 
%
% Each UIControl defines a LocationCode which is a string in the form of 
% 'UserA1', 'UserA2', 'UserB1', ...'UserB8', 'UserC1', ..., 'UserC8'.
% Depending on the application this code defines the location where to put
% the custom defined UIControl. 
%
% We recommend to use these UIControls to define user defined Controls
% independent of the GUI. This way the script can also be used by different
% GUI applications. 
%
% The user can also provide a callback that will be executed if the value
% of the UIControl is changing. The callback should be a function handle: 
% use @()myFunction(); to create a function handle. the function should
% accept 3 input parameter for all value-based ui controls (slider, toggle
% button, button group) and should be provided as @(src, evt,
% UIValue)myFunction(src, evt,UIValue). If fewer arguments are given the
% function will be called with only the arguments that fit. Lets say if the
% user only defines one input argument the function will only be called
% with the src. 
% The src is the source graphics object that triggered the event. We
% provide access to that handle object (hObject) to support older scripts.
% However, we do not recommend to modify GUI components directly. As
% mathworks is developing new GUI technologies we may not be able to
% support the matlab uicontrols in future releases. The evt object provides
% some means to access the calling application which also gives access to
% the other controls. 
% 
%   
%
% Version 1.0 | 2020-02-04 
% $Author: Dr. Daniel Rohrbach
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
%
        
    properties(SetAccess=private)
        
        % the ui control style @type char, this will be filled in
        % automatically by the appropriate class. A VsSliderControl will
        % fill this in as 'VsSlider'
        Style = '';
        
        % the ui control location code @type char some UIs define a
        % location code where the UI can be placed.
        LocationCode = '';
        
        % @type vsv.apps.AbstractVSXApp the parent application
        ParentApplication;
        
    end
    
    
    methods(Abstract)
        
        UI = convertToStructUIControl(this);
        % converts this ui control to the old fashion Struct representation
        % of this UI control
        %
        % @return strUI - the struct representation of this UI control
        
        isvalid = isValidStyle(this, style);
        % test whether the given style is valid
        %
        % @param style - @type char the style to check for 
        % @return isvalid - @type logical true if style is valid, false
        % otherwise
        
    end
    
    
    %% static helper functions
    methods(Static=true)
    
        function isUserUI = isUserUiControl(UI)
        % Checks whether the given UI is a user defined UI control
        %
        % @param UI - struct with a control field and a callback field
            
            isUserUI = isfield(UI,'Control') && ~isempty(UI.Control) && ...
                strncmp(UI.Control{1},'User',4);
        end
        
        function Style = getUIStyle(UI)
        % parses the ui Style if present
        % @return Style - the style set in the ui control, if no style is
        %                 found the function returns NULL i.e. []
            if vsv.seq.uicontrol.UserUiControl.isUserUiControl(UI)
               
                Style = vsv.seq.uicontrol.UserUiControl.getUIStyleProtected(UI);
            else
                Style = '';
            end
            
        end
        
      
    end
    
    methods(Static,Access=protected)
        
        function Style = getUIStyleProtected(UI)
        % parses the ui Style if present
        % @return Style - the style set in the ui control, if no style is
        %                 found the function returns NULL i.e. []
            
            control = UI.Control;
            index = find( strcmp('Style',control) );

            if ~isempty(index)
                Style = control{index+1};
            else
                Style = [];
            end
            
        end
    end
    
    
    methods
             
        function this = UserUiControl(UI, varargin)
        % Constructor to create the user UI control
        %   
            
            if nargin >= 1
               this.initUserUiControl( UI, varargin{:} );
            else
               this.initUserUiControl( );
            end
        end
        
        function isvalid = isValidControl(this, UI)
        % check whether the UI control is valid
        %
        % A ui control should be valid if it is a isUserUiControl
        %
        % @param UI - the ui control struct
        %
            isvalid = this.isUserUiControl( UI );
            
        end
        
        
        
    end
    
    %% setter and getter
    methods
        
        function setParentApplication(this, pa)
        % setter for parent application
            if isa(pa, 'vsv.apps.AbstractVSXApp') || vsv.util.isNULL(pa)
                this.ParentApplication = pa;
            else
                error( 'UserUIControl:invalidParentApplication', ...
                    'Parent application must be a vsv.apps.AbstractVSXApp' );
            end
        end
        
        function pa = getParentApplication(this)
            pa = this.ParentApplication;
        end
        
        function setLocationCode(this, code)
        % setter for the location code
        %
        % @param code - @type char the code spcifying the location
        
            if vsv.util.isStringScalar(code)
                code = convertStringsToChars(code);
            end
            
            if ischar(code)
                this.LocationCode = code;
            else
                error('UserUiControl:LocationCode:invalidArgument', ...
                    'Code must be a char or string');
            end
        end
        
        function setCallback(this, fnc)
        % make sure that the callback be converted to a functionDef 
        %
        % vsv.seq.function.ExFunctionDef will be cleared after they run
        % out of scope. So if this userUiControl is running out of scope
        % the function def will run out of scope and its function space
        % will be cleared. This is important for persistent variables that
        % should be cleared for each run of a script
        % 
        % @param fnc - @type function_handle or [] or vsv.seq.function.AbstractFunctionDef  
            
            if isa( fnc, 'vsv.seq.function.AbstractFunctionDef' )
                this.setCallback@vsv.seq.uicontrol.mixin.HasCallback(fnc);
            elseif vsv.util.isNULL(fnc)
                this.setCallback@vsv.seq.uicontrol.mixin.HasCallback(fnc);
            else
                fname = [ this.LocationCode 'Callback' ];
                
                exFunc = vsv.seq.function.ExFunctionDef( fname, fnc );
                this.setCallback@vsv.seq.uicontrol.mixin.HasCallback(exFunc);
                
            end
        end
        
    end
    
    methods(Access=protected)
        
        function initUserUiControl(this, UI, varargin)
        % Constructor to create the user UI control
        %   
            
            if nargin == 2 && ~isempty(UI)
                this.parseUserUi(UI);
            elseif nargin > 2
                varargin = [ {UI} varargin ];
                this.applyControls( varargin{:} );
            end
        end
        
        
        function success = applyParamValue(this, param, value)
        % set the following param value pairs
        %
        
            success = true;
            switch param
                case 'Style'
                    if this.isValidStyle( value )
                        this.setStyle(value);
                    else
                        success = false;
                    end
                case 'LocationCode'
                    this.setLocationCode(value);
                case 'Callback'
                    this.setCallback(value);
                case 'Statement'
                    this.setStatement(value);
                otherwise 
                    success = false;
            end
        end
    
        
        function parseUserUi(this, UI)
        % will parse the user
        %
            % we know it starts with a user tag and has a control field
            if this.isValidControl(UI)
                
                this.setLocationCode( UI.Control{1} );
                if isfield(UI, 'Callback') 
                    this.setCallback(UI.Callback);
                end
                if isfield(UI, 'Statement')
                    this.setStatement( UI.Statement);
                end
                
                this.applyControls( UI.Control{2:end});
                
            else
                error( 'UserUiControl:parseUserUi:invalidControl', ....
                        'given UI control is not a valid VS ui control');
            end
            
        end
        
        function applyControls(this, varargin)
        % used by the parseUserUi to apply the controls
        %
            params = varargin(1:2:end);
            values = varargin(2:2:end);
            
            nParams = numel(params);
            for i = 1:nParams
                success = this.applyParamValue(params{i}, values{i}); 
                if ~success
                    if ischar( params{i} ) || isstring( params{i} )
                        error( 'UserUiControl:applyControls:invalidParam', ...
                            ['Cannot import given parameter value pairs ' ...
                            'because of unknown parameter or invalid ' ...
                            'parameter value: ' params{i} ]);
                    else
                        error( 'UserUiControl:applyControls:invalidParam', ...
                            [' Parameter name at index ' num2str(i) ' is not a char or string' ] );
                    end
                end
            end
        end
       
        
        function setStyle(this, style)
        % setter for the style property
        
            if isstring(style)
                style = convertStringsToChars(style);
            end
            
            if ischar(style)
                this.Style = style;
            else
                error('UserUiControl:setStyle:invalidArgument', ...
                    'Label must be a char or string');
            end
        end
    end
    
end

