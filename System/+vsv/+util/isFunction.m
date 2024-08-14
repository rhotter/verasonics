function [ res, ID ] = isFunction(fun)
%ISFUNCTION tests whether the given fun is a function or something else
%   fun can be a char that is funciton in the seach path, or a function
%   handle
% 
%   Also returns an identifyer ID indicating what fun may is. 
%       
%     ID = -1 ; % script
%     ID = -2 ; % probably another type of file, or it does not exist
%     ID = -3 ; % probably a handle, but not to a function
%     ID = -4 ; % probably a variable or an array
%     ID = 0 ; % unknown cause for error
% 
% Version 1.0 | 2020-05-18 
% $Author: Dr. Daniel Rohrbach
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
% 

    % @ToDo check whether there is a better way to do this. The problem is
    % that a try catch for testing something is not a good solution, will
    % cost performance for sure. 
    
    try
        nargin(fun); % nargin errors when FUNNAME is not a function
        ID = 1;
    catch err
        % catch the error of nargin
        switch (err.identifier)        
            case { 'MATLAB:nargin:isScript', 'MATLAB:persistentNotInFunction' }
                ID = -1 ; % script
            case 'MATLAB:narginout:notValidMfile'
                ID = -2 ; % probably another type of file, or it does not exist
            case 'MATLAB:narginout:functionDoesnotExist'
                ID = -3 ; % probably a handle, but not to a function
            case 'MATLAB:narginout:BadInput'
                ID = -4 ; % probably a variable or an array
            otherwise
                ID = 0 ; % unknown cause for error
        end
    end
    
    res = ID > 0;
    
end

