classdef VsSliderControl < vsv.seq.uicontrol.VsRangeControl
%VsSliderControl defines a slider control handle that can be added to a GUI 
%   
%   A VsSliderControl inherits from the VsRangeControl, which defines a
%   range control object in a GUI environment. 
%
%   to get more information type 
%   >> help vsv.seq.uicontrol.VsRangeControl 
%
%   To get general information about the UIControl type 
%   >> help vsv.seq.uicontrol.UserUIControl
%
%   The UserUIControl is the basic class for all ui controls and provides
%   information about common attributes such as Callback, and LocationCodes
%
%   Usage:
%   ------
%      import vsv.seq.uicontrol.VsSliderControl;
%
%      vsSlider1 = VsSliderControl( 'LocationCode',    'UserB7', ...
%                                   'Label',           'Sens. Cutoff',   ...
%                                   'SliderMinMaxVal', [ 0,1.0,Recon(1).senscutoff], ...
%                                   'SliderStep',      [0.025,0.1], ...
%                                   'ValueFormat',     '%1.3f');
%   
%
%      % use setter functions to modify some of the fields
%      vsSlider1.setCallback( @(s,e)disp('Test') );
%      vsSlider1.setLabel( 'Sensitivity Cutoff' ); 
%      vsSlider1.setRange( [0 5] );  
%       
%
% Version 1.0 | 2020-01-24 
% $Author: Dr. Daniel Rohrbach
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
% 
    
    properties(SetAccess=protected)
        
        % @type numeric two dimensional array provides step for the slider
        SliderStep = [0.01,0.1];
        
    end
        
    %% setters and getters 
    methods
        
        function this = VsSliderControl( varargin )
        % constructor of the VsSliderControl    
        %
        % @param rangeValue - @type numeric the range initial value
        % @param UI - @type the UI control structure 
        
            this@vsv.seq.uicontrol.VsRangeControl(varargin{:});
            this.setStyle( 'VsSlider' );
            
        end
        
        function strUI = convertToStructUIControl(this)
        % converts this ui control to the old fashion Struct representation
        % of this UI control
        %
        % @return strUI - the struct represention of this UI control
        
            strUI.Control  = {this.LocationCode, ...
                              'Style',          this.Style, ...
                              'Label',          this.Label, ...
                              'SliderMinMaxVal', [ this.RangeMin this.RangeMax this.CurrentValue],...
                              'SliderStep',     this.SliderStep, ...
                              'ValueFormat',    this.ValueFormat};
            strUI.Callback = this.Callback;
            if ~vsv.util.isNULL( this.Statement )
                strUI.Statement = this.Statement;
            end
            
        end
        
        function setSliderStep(this, ss)
        % setter for the slider step
        %
        
            if isnumeric(ss) && length(ss) == 2
                this.SliderStep = ss;
            else
                error( 'VsSlider:setSliderStep:invalidArgument' ,...
                    'SliderStep must be a 2 dimensional numeric array');
            end
        end
        
        function isvalid = isValidControl(this, UI)
        % check whether the UI control is valid
        %
        % A ui control should be valid if it is a isUserUiControl
        %
        % @overwritten  
        % @param UI - the ui control struct
        %
         
            
            isvalid = this.isValidControl@vsv.seq.uicontrol.VsRangeControl( UI );
            if isvalid
                Style = vsv.seq.uicontrol.UserUiControl.getUIStyleProtected(UI);
                isvalid = ~isempty(Style) && strcmp( Style, 'VsSlider' );
            end
            
        end
       
    end
    
    methods(Access=protected)
        
        function success = applyParamValue(this, param, value)
        %
        %
        
            success = true;
            switch param
                case 'SliderStep'
                    this.setSliderStep(value);
                case 'SliderMinMaxVal'
                    this.setMinMaxVal(value);    
                otherwise 
                    success = this.applyParamValue@vsv.seq.uicontrol.VsRangeControl(param, value);
            end
        end
    end
    
end

