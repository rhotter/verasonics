classdef (Sealed) ExFunctionDef < vsv.seq.function.AbstractFunctionDef
% EXFUNCTIONDEF definition of an external function using a function name
% and handle
%   
%  The ExFunctionDef can be understood as a tuple of function name and
%  function handle ( i.e., the class assigns a function name to the
%  function handle ). Otherwise there is no easy way to extract the
%  function name from a function handle if the function handle is an
%  anonymous function. An anonymous function handle assigns a one line
%  function to a variable without a function name.
%
%  It is also possible to use direct function handles. The difference
%  between a function handle and an function handle to an anonymous
%  function is the following
%
%  Anonymous function handle:
%  fnc = @(X,Y,Z)aFunction(X,Y,Z)
%       The anonymous function handle (indicated by @) creates a handle to
%       a new function that has no name but, in the example above, takes two
%       input parameter X,Y, and Z and passes them to a function CALLED
%       INSIDE the anonymous function (aFunction, in that example). The
%       advantage of anonymous functions is that they allow to create new
%       functions that redefine input parameter without the need of
%       creating a new m-file. For example aFunction in the example above
%       takes 3 input parameter. But lets say we want to pass the function
%       handle to a place where it will be called only with two input
%       parameter. this is the case for many callback functions in
%       Matlab graphical user interfaces. In this case you can do the
%       following:
%
%  MyY = 'My Test Valye'; 
%  fnc = @(X,Z)aFunction(X,MyY,Z);
%       In this case MyY will be stored with the function handle and will
%       be passed to it every time the function is called. Now fnc defines
%       a variable to an anonymous function handle that only takes 2 input
%       parameter (X and Z). You can call it with
%
%  fnc(val1, val2);
%
%       This call will be evaluated to  aFunction(val1,MyY, val2);
% 
%  A function direct function handle save a handle to a defined function
%  directly, which doesn't allow to define new input parameter but has the
%  advantage that the user doesn't need to define it themselves:
%
%  fnc = @aFunction 
%  Now fnc points directly to aFunction and will act exactly like that when
%  called fnc(X,Y,Z)
%
%  In previous Vantage software versions we allowed the user to specify
%  char arrays as function names that then will be executed using feval.
%  However, we believe that the function handle approach provides more
%  advantages. For example using a char can be ambiguous when the
%  function is evaluated and could change inbetween calls. Furthermore, the
%  char array approach does not allow to generate function handles to
%  nested functions defined in scripts. 
%
%
%
% @Dependency
%  The class inherits from vsv.seq.function.AbstractFunctionDef, which
%  defines the two main properties of the class:
%
%         % @type char the name of the function
%         FunctionName;
%
%         % @type function_handle the name of the function
%         FunctionHandle;
%
% @Usage:
%
%      % using anonymous functions and constructor input parameter with
%      % function name and anonymous function handle 
%      exf = vsv.seq.function.ExFunctionDef('myProcFunction', ...
%                                   @(RData)myProcFunction(RData));
%
%      % using anonymous functions and a struct as input parameter with
%      % function name and anonymous function handle
%      str.FunctionName   = 'myProcFunction';
%      str.FunctionHandle = @(RData)myProcFunction(RData)
%      exf = vsv.seq.function.ExFunctionDef( str );
%  
%      % same with direct function handles 
%      exf = vsv.seq.function.ExFunctionDef('myProcFunction', ...
%                                              @myProcFunction);
%
%      % using a struct as input parameter with
%      % function name and function handle
%      str.FunctionName   = 'myProcFunction';
%      str.FunctionHandle = @myProcFunction
%      exf = vsv.seq.function.ExFunctionDef( str );
%
% @Note
%  The class is designed to be an immutable class which means, once an
%  instance is created the function handle and name cannot be changed
%
% Version 1.2 | 2020-01-03 
%   - 2020-22-05 update documentation and moved to public folder
% $Author: Dr. Daniel Rohrbach 
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
%     
    properties(SetAccess=protected)
      
        
        % @type logical true if the function should clear the function
        % handle this will ensure that persistent variables are cleared
        % default is true
        IsClearDeletedFunction = true;
        
    end
    
    methods(Access=protected)
        
        function initExFunctionDef(this, functionName, functionHandle )
        % initialize the function definition allowed input are a single
        % struct or the function name and function handle
        %
            
            if nargin == 2 && isstruct(functionName)
                this.initExFunctionDef( functionName.FunctionName, functionName.FunctionHandle );
            else
                this.setFunctionHandle( functionHandle );
                this.setFunctionName( functionName );
            end
            
        end
 
    end
    
    methods
        
        function this = ExFunctionDef(varargin)
        % Creates a External function definition
        %   
        % Usage:
        %      exf = vsv.seq.function.ExFunctionDef('myProcFunction', ...
        %                                   @(RData)myProcFunction(RData));
        %
        %
        %      str.FunctionName   = 'myProcFunction';
        %      str.FunctionHandle = @(RData)myProcFunction(RData)
        %      exf = vsv.seq.function.ExFunctionDef( str );
        %  
        %
        % @param functionName - @type char the name of the function
        % @param functionHandle - @type function_handle the function handle
        % 
        
            this.initExFunctionDef( varargin{:} );
        end
        
        function isvalid = isValidFunction(~)
        % this is always true because the class can only be initialized
        % with a valid function handle
        
            isvalid = true;
        end
        
        function enableIsClearDeletedFunction(this)
        % this will ENABLE to clear function handle when ExFunctionDef is
        % deleted
        %
        %   this is the default behavior
        %
            this.IsClearDeletedFunction = true;
        end
        
        function disableIsClearDeletedFunction(this)
        % this will DISABLE to clear function handle when ExFunctionDef is
        % deleted
        %
            this.IsClearDeletedFunction = false;
        end

        
        function delete(this)
        % Deconstructor 
        %
        % This will clear the calling function handle or its parent file.
        % Some external functions make use of persistent variables which
        % need to be cleared before use.
        %
        
            if this.IsClearDeletedFunction
                % if the constructor throws and error the following can
                % throw an error
                %disp([ 'Clear Ex Function Def: ' this.FunctionName ]);
                this.clearFunction();
            end
        end
        
    end
    
    
end

