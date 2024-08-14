classdef (Abstract) ParamImportable < handle
%STRUCTIMPORTABLE Provides function to import struct
%   Detailed explanation goes here
%
% Version 1.0 | 2020-11-04 
% $Author: Dr. Daniel Rohrbach 
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
% 
    
    methods
    
        function this = ParamImportable(varargin)
        % Allows input of a struct with the fields matching the class
        % fields, or a param value pairs if desired. Param value pairs is
        % recommended 
        %
        
            varargin = this.preInit(varargin{:});
            if nargin == 1 && isstruct(varargin{1})
                this.importStruct(varargin{1});
            else
                this.setParamValue(varargin{:});
            end
        end
        
        function importStruct(this, str)
        % allows to import a struct whoes fields matches the props of the
        % class
        %
        % This is not recommended but allows the user to use structs if
        % that is preferred
        %
        % @param str - @type struct a struct to import the parameter, the
        %              struct fields must match exactly the fields of this
        %              class properties
        % 
        
            fields = fieldnames(str);
            vals   = struct2cell(str);
            
            var = cell( 1, numel(fields)*2 );
            var( 1:2:end ) = fields;
            var( 2:2:end ) = vals;
            
            this.setParamValue( var{:} );
            
        end
        
        
        function setParamValue(this, varargin)
            import vsv.util.mixin.ParamValuePair;
            [params, vals, msg] = ParamValuePair.prepareParamVals(varargin{:});
            % if msg is empty means no errors
            if isempty(msg)
                
                % apply the parameter value pairs 
                msg = ParamValuePair.setParamValuePair(this, params, vals);
                if ~isempty(msg)
                    error('ParamImportable:setParamValue:cannotImport', ...
                    [ 'Cannot import parameter', newline ...
                      msg ]);
                end
                
            elseif nargin == 1
                error('ParamImportable:setParamValue:invalidNumberInput', ...
                    [ 'Number of input parameter must be even or 0', newline ...
                      msg ]);
            end
            
        end
        
                
    end
    
    methods(Access=protected, Abstract)
        
        varargin = preInit(this, varargin);
        
    end
    
    
end

