classdef (Abstract) AbstractFunctionDef < handle
%ABSTRACTFUNCTIONDEF Raw model for function definitions
%   
% 
% Version 1.0 | 2020-01-03 
% $Author: Dr. Daniel Rohrbach 
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
%     


    properties(SetAccess=protected)
        
        % @type function_handle the name of the function
        FunctionName;
        
        % @type char the name of the function
        FunctionHandle;
    end
    
    methods
        
        isvalid = isValidFunction(this)
        % returns true if this function is valid and can be executed
        %
        % @return isvalid - @type logical true if this can be executed
        %                   false otherwise
        
    end
    
    methods
        
        function varargout = executeFunction(this, varargin )
        % executes the function handle with the given arguments
        % @param varargin - gerneric input arguments accepted by the
        %                   underlying function handle
        % @param varargout - gerneric output defined by the underlying
        %                    function
        
            if this.isValidFunction()
                varargout = cell( 1, nargout); 
                [varargout{:}] = this.FunctionHandle( varargin{:} );
            else
                error('ExFunctionDef:functionInvalid', ...
                     'This function is not valid to be executed, isValidFunction returned false');
            end
            
        end
        
    end
    
    methods(Access=protected)
        
        function clearFunction(this)
            
            if ~isempty(this.FunctionHandle) 
                fh = functions( this.FunctionHandle );
                if isfield( fh, 'within_file_path')
                    clear(fh.within_file_path);
                end
                if isfield( fh, 'file')
                    % this causes an annoying warning in R2018a in the unit
                    % tests otherwise 
                    %@ToDo remove with R2018a not supported
                    info  = warning( 'query', 'MATLAB:lang:cannotClearExecutingFunction');
                    warning('off', 'MATLAB:lang:cannotClearExecutingFunction');
                    clear(fh.file);
                    warning( info.state, 'MATLAB:lang:cannotClearExecutingFunction');
                end
                if ~isempty(this.FunctionName ) && ischar(this.FunctionName )
                    clear( this.FunctionName );
                end
% % keep this in case this is necessary, it will drop performance but we may  
% % need to do this to ensure all links to the function are deleted
%                 if isfield( fh, 'workspace') && ~isempty(fh.workspace)
%                     % clear all objects of the workspace 
%                     if iscell( fh.workspace )                        
%                         workspace = fh.workspace;
%                         nw = numel(workspace);
%                         for i = 1:nw
%                             wi = workspace{i};
%                             if isstruct(wi)
%                                 fi = fields(wi);
%                                 for k = 1:length(fi)
%                                      obj = wi.( fi{k} );
%                                      clname = class(obj);
%                                      
%                                      elements = regexp(clname, '\.', 'split');
%                                      clear( elements{end});
%                                 end
%                             end
%                         end
%                     end
%                 end
            end
            
        end
        
        function setFunctionName(this, functionName)
        % setter for function name
        %
        % protected because function name should not change after created
        
            if ischar(functionName)
                functionName = regexprep(functionName, '\.m$', '');
            end
            
            if ischar(functionName) ...
                    && ~isempty(functionName) ...
                    && isempty( regexp(functionName, '^[0-9_]+|[^a-zA-Z0-9_]', 'once') )
                this.FunctionName = functionName;
            else
                 error('ExFunctionDef:invalidArgument', ...
                     'functionName must be a char and must not contain special characters');
            end
            
        end
        
        function setFunctionHandle(this, functionHandle)
        % setter for function handle
        %
        % protected because function name should not change after created
        
            if isa(functionHandle, 'function_handle')
                this.FunctionHandle = functionHandle;
            else
                 error('ExFunctionDef:invalidArgument', ...
                     'functionHandle must be a function_handle');
            end
            
        end
    end
    
    
end

