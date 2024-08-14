classdef VsStatementOnly < vsv.seq.uicontrol.UserUiControl
%VSSTATEMENTONLY Summary of this class goes here
%   Detailed explanation goes here
    
    properties
        
    end
    
    methods(Static=true)
        
        function ID = getStatementID(startID)
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
            ID = [ 'Statement' num2str( ID ) ];
            
        end
        
    end
    
    methods
        
        
        function UI = convertToStructUIControl(this)
        % converts this ui control to the old fashion Struct representation
        % of this UI control
        %
        % @return strUI - the struct represention of this UI control
            
           UI.Statement = this.Statement;
        
        end
        
        function isvalid = isValidStyle(~, ~)
        % test whether the given style is valid
        %
        % @param style - @type char the style to check for 
        % @return isvalid - @type logical true if style is valid, false
        % otherwise
            isvalid = true;
        end
        
    end
    
end

