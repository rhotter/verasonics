classdef VsDropDownControl < vsv.seq.uicontrol.VsCaseControl 
%VSDROPDOWNCONTROL defines a drop down menue that can be added to a gui
%   
% The drop down control inherits from a VSsCaseControl and can be used in
% similar cases as the button group control.
%
% ----------
%   Usage:
%       dd =  vsv.seq.uicontrol.VsDropDownControl( ...
%                                 'LocationCode',    'UserA1', ...
%                                 'Title',           'Mode', ...
%                                 'PossibleCases',   {'Velocity','Power'});
%
%       dd.setSelectedCase( 'Power' );
%       dd.SelectedIndex
%       dd.setSelectedIndex( 2 );
%       dd.SelectedCase
%
% Version 1.0 | 2020-04-01 
% $Author: Dr. Daniel Rohrbach
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
%    
    
    methods
        
        function this = VsDropDownControl(varargin)
            this@vsv.seq.uicontrol.VsCaseControl( varargin{:} );
            this.setStyle( 'VsDropDown' );
        end
        
        function isvalid = isValidStyle(~, style)
            isvalid = ~isempty(style) && strcmp( style, 'VsDropDown' );
        end
        
    end
    
end

