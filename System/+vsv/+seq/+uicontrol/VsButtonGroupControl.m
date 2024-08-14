classdef VsButtonGroupControl < vsv.seq.uicontrol.VsCaseControl 
%VSBUTTONGROUPCONTROL a group of selectable options using radio buttons 
%
%   A VsButtonGroupControl inherits from the VsCaseControl, which defines a
%   control that allows to select cases in a GUI environment. 
%
%   to get more information type 
%   >> help vsv.seq.uicontrol.VsCaseControl 
%
%   To get general information about the UIControl type 
%   >> help vsv.seq.uicontrol.UserUIControl
%
%   % ----------
%   Usage:
%       bg =  vsv.seq.uicontrol.VsButtonGroupControl( ...
%                                 'LocationCode',    'UserA1', ...
%                                 'Title',           'Mode', ...
%                                 'PossibleCases',   {'Velocity','Power'});
%
%       bg.setSelectedCase( 'Power' );
%       bg.SelectedIndex
%       bg.setSelectedIndex( 2 );
%       bg.SelectedCase
%
%
% Version 1.0 | 2020-01-28 
% $Author: Dr. Daniel Rohrbach
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
%
      
      
    methods
        
       
        function this = VsButtonGroupControl(varargin)
        % Constructor for the VsButtonGroupControl class
        %
        %   the button group control must be initialized with the possible
        %   cases that can be selected
        %
        % Errors
        %   @error VSButtonGroupContro:setPossibleCases:invalidArgument if
        %          posc is not a cellstring
        %
        % @param posc - possible cases @type cellstring
            
        
            if nargin == 1
                
                % we have an UI struct
                UI = varargin{1};
                if this.isValidControl(UI)
                    % replace Labels as possible cases to user super init
                    
                    varargin = [ { 'LocationCode', UI.Control{1} } ...
                                 UI.Control(2:end) ];
                                 
                    if isfield(UI, 'Callback')
                        varargin = [ varargin ...
                                     {'Callback', UI.Callback }];
                    end
                        
                else
                    error('VsButtonGroupControl:invalidUI', ...
                        [ 'If one argument is given, given argument must be a valid UIControl of type' ... 
                          'VsButtonGroup' ] );
                end                
            end
            
            if numel(varargin) >= 2
                
                params = varargin( 1:2:end );
                
                varargin(1:2:end) = regexprep( params, '^Labels$', 'PossibleCases' );    
                ind = strcmp( varargin,  'PossibleCases' );
                
                if ~any(ind)
                    ind = find( strcmp( varargin,  'NumButtons' ) );
                    if numel(ind) == 1
                        
                        nB = varargin{ind+1};
                        possibleCases = cell( 1, nB);
                        for i = 1:nB
                            possibleCases{i} = num2str(i);
                        end
                        
                        varargin = [ {'PossibleCases', possibleCases}, varargin ];
                    end
                end
            end
            this.initVsCaseControl(varargin{:});
            this.setStyle( 'VsButtonGroup' );
        end
        
        function isvalid = isValidStyle(~, style)
            isvalid = ~isempty(style) && strcmp( style, 'VsButtonGroup' );
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
                isvalid = ~isempty(Style) && strcmp( Style, 'VsButtonGroup' );
            end
            
        end
        
    end
    
    
    methods(Access=protected)
        
        
        function success = applyParamValue(this, param, value)
        %
        %
        
            success = true;
            switch param
                case 'NumButtons'
                    % just do not do anything its for backward comp
                case 'Labels'
                    this.setPossibleCases(value);
                otherwise 
                    success = this.applyParamValue@vsv.seq.uicontrol.VsCaseControl(param, value);
            end
            
        end
        
    end
    
    
end

