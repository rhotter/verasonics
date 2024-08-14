function [resultTable, results, benchmarks, workspace] = vsProfile( cmd, varargin )
%VSPROVFILE profile execution times of scripts
%   
% Examples
%
%    % profile a script and review results in a table 
%    vsProfile start
%    VSX
%    SetUpL11_5Flash
%    vsProfile stop
%    vsProfile review          
% 
%    % all commands can be executed as string inputs. Commands will always
%    % be executed on last profiled script. To get the results, result
%    % table, and benchmarks from the last profile you can use
%    [resultTable, results, benchmarks] = vsProfile( 'review');
% 
%    % scripts can be profiled in a single line call
%    [resultTable, results, benchmarks] = vsProfile( 'script', 'SetUpL11_5gHAcquireRF' );
%    
%    % the profiler also comes with a simple visualization of the profiler
%    % results 
%    vsProfile view
% 
%
% Version 1.0 | 2021-12-24
% $Author: Bryan Cunitz, Dr. Daniel Rohrbach 
% Copyright 2001-2021 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.
%

    timing = vsv.profile.TimingTool();
    
    switch(cmd)
        case "review"
           [resultTable, results, benchmarks] = timing.(cmd)(varargin{:});
           warningNoReturn(resultTable)

        case "script"
           [resultTable, results, benchmarks, workspace] = timing.profileScript(varargin{:});
           warningNoReturn(resultTable)

        otherwise
           timing.(cmd);
           if nargout > 0
                error('vsProfile:nargout', [ 'given command ' cmd  ' Does not support return arguments']);
           end
    end
    
    
end

function warningNoReturn(resultTable)
    if isempty(resultTable)
        warning('vsProfile:emptyProfile','Profile review is empty.  This is most likely due to vsProfile not being initiated with "vsProfile start"');
    end
end

