classdef VsButtonControl  <   vsv.seq.uicontrol.UserUiControl  ...
                    & vsv.seq.uicontrol.mixin.HasLabel
%VSBUTTONCONTROL describes a use ui button control
%
% A button control does not extend any property of the
% vsv.seq.uicontrol.UserUiControl. It predefines the style to be
% VsPushButton 
%
% -----------
% Usage:
%
%       b = vsv.seq.uicontrol.VsButtonControl( ...
%                                     'LocationCode',   'UserC3', ...
%                                     'Label',          'HWFilter');
%       b.setLabel( 'Test' );
%       b.Label
%
%
% Version 1.0 | 2020-04-01 
% $Author: Dr. Daniel Rohrbach
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
%
    
    
    methods
        
        function this = VsButtonControl( varargin )
        % constructor of the VSButtonControl    
        %
        % @param UI - @type the UI control structure 
        
                        
            this@vsv.seq.uicontrol.UserUiControl( varargin{:} );
            this.setStyle( 'VsPushButton' );
            
        end
        
        function strUI = convertToStructUIControl(this)
        % converts this ui control to the old fashion Struct representation
        % of this UI control
        %
        % @return strUI - the struct represention of this UI control
        
            strUI.Control  = {this.LocationCode, ...
                              'Style',          this.Style, ...
                              'Label',          this.Label };
            strUI.Callback = this.Callback;
            
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
                isvalid = this.isValidStyle(Style);
            end
            
        end
        
        function isvalid = isValidStyle(~, style)
            isvalid = ~isempty(style) && strcmp( style, 'VsPushButton' );
        end
        
    end
    
    methods(Access=protected)
        
        function success = applyParamValue(this, param, value)
            
            success = true;
            switch param
                case 'Label'
                    this.setLabel(value);
                otherwise 
                    success = this.applyParamValue@vsv.seq.uicontrol.UserUiControl( param, value );
            end
        end
        
    end
end

