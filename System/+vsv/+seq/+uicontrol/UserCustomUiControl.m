classdef (Sealed) UserCustomUiControl <  vsv.seq.uicontrol.mixin.HasStatement ...
                                   & vsv.seq.uicontrol.mixin.HasCallback
%USERCUSTOMUICONTROL for backward compatibiliy
%
% this is 
%
% Version 1.0 | 2020-02-07 
% $Author: Dr. Daniel Rohrbach
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
%

    properties
        
        Control;
        
        ParentAplication;
        
    end
    
    methods(Static=true)
        
        function ID = getCustomID(startID)
        % returns a unique Statement ID within the program live cycle. 
        %
        % calling clear all will delete the ID and restart counting
        %
        % @return ID - @type char a unique ID within a live cycle 
        %
        
            persistent usedID 
            
            if isempty(usedID)
                if nargin < 1
                    usedID = 0;
                else
                    usedID = startID;
                end
                    
            end
            ID = [ 'CustomControl' num2str( usedID ) ];
            
        end
        
    end
    
    methods
        
        function this = UserCustomUiControl(UI, varargin)
            if nargin >= 1
               this.initUserUiControl( UI, varargin{:} );
            else
               this.initUserUiControl( );
            end
        end
        
        function app = getParentApplication(this)
            app = this.ParentAplication;
        end
        
        function UI = convertToStructUIControl(this)
        % converts this ui control to the old fashion Struct representation
        % of this UI control
        %
        % @return strUI - the struct represention of this UI control
            
            
            UI = struct();
            if this.hasControl()
                UI.Control = this.Control;
            end
            
            if this.hasCallback()
                UI.Callback = this.Callback;
            end
            
            if this.hasStatement()
                UI.Statement = this.Statement;
            end
            
        end
        
        
       
        function hasCntr = hasControl(this)
            hasCntr = ~vsv.util.isNULL(this.Control);
        end
        
        function setControl(this, cntr)
        % setter for control
        %
        
            if vsv.util.isNULL(cntr) || iscell(cntr)
                this.Control = cntr;
            else
                error('UserCustomUiControl:setControl:invalidControl', ...
                'Given controller is invalid');
            end
        end
        
    end
    
    methods(Access=protected)
        
        function parseUserUi(this, UI)
        % will parse the user
        %
            % we know it starts with a user tag and has a control field
            
                
            
            if isfield(UI, 'Callback') 
                this.setCallback(UI.Callback);
            end
            
            if isfield(UI, 'Statement')
                this.setStatement( UI.Statement);
            end
            
            if isfield(UI, 'Control')
                this.setControl( UI.Control);
            end
            
        end
        
        function initUserUiControl(this, UI, varargin)
        % Constructor to create the user UI control
        %   
            
            if nargin == 2 && ~isempty(UI)
                this.parseUserUi(UI);
            elseif nargin > 2
                varargin = [ {UI} varargin ];
                success = this.applyControls( varargin{:} );
                if ~success
                    error('UserCustomUiControl:initUserUiControl:invalidInput', ...
                        'Only Control, Callback and Statement are allowed for UserCustomUiControl')
                end
            end
        end
        
        function success = applyControls(this, param, value)
        % set the following param value pairs
        %
        
            success = true;
            switch param
                case 'Control'
                    this.setCallback
                case 'Callback'
                    this.setCallback(value);
                case 'Statement'
                    this.setStatement(value);
                otherwise 
                    success = false;
            end
        end
    
        
    end
end

