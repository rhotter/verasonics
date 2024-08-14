function VVT(varargin)
%-------------------------------------------------------------------------%
%  Runs VVT across the secondary systems, either individually or sequentially %
%
%  usage:
%    vsv.multi.VVT  by itself it defaults to sequetially running VVT on all
%                   the salve systems and creating a single log file under
%                   Hardwaretest
%
%    vsv.multi.VVT(secondaryNum) runs VVT remotely on the secondary system and creates
%                         a single log file under the primary Hardwaretest
%                         directory as well as the secondary's hardwaretest
%                         directroy
%  Example:
%     vsv.multi.VVT(2)    runs VVT on secondary #2
%-------------------------------------------------------------------------%

% Copyright 2001-2021 Verasonics, Inc.  All world-wide rights and remedies
% under all intellectual property laws and industrial property laws are
% reserved.  Verasonics Registered U.S. Patent and Trademark Office.

basedir = [getenv('VERASONICS_VPF_ROOT') '/'];
testLogsDir = [basedir '/HardwareTest/TestResultLogs/'];
numSecondaries = vsv.multi.getNumSecondaries();
switch numSecondaries
    case 4
        IPlist = {'10.10.0.2', '10.10.1.2', '10.10.2.2', '10.10.3.2'};  % 5 system configuration
    case 3
        IPlist = {'10.10.1.2', '10.10.2.2', '10.10.3.2'}; % 4 system configuration
    case 2
        IPlist = {'10.10.1.2', '10.10.2.2'}; % 3 system configuration
    case 1
        IPlist = {'10.10.1.2'}; % 2 system configuration
    otherwise
        error('The system network IP addresses are in a non-standard configuration.  Please contact Verasonics support.')
end

disp('######################################')
disp('###   Verasonics MultiSystem       ###')
disp('###     Verification Test          ###')
disp('######################################')


if length(varargin)>1
    error('incorrect number of inputs, just a number 0 to 3 for the secondary or leave empty')
elseif isempty(varargin)
    disp('== Performing VVT on all systems in parallel ==');
    for secondarynum = 1:numSecondaries
        fprintf('--- Perform VVT on Secondary %i ---\n',secondarynum);
        cmd = sprintf('ssh %s "killall -9 MATLAB"', IPlist{secondarynum}); %stop all running matlab sessions
        [~,~]=unix(cmd);
        cmd = sprintf(['ssh %s -X ''matlab -sd %s -nosplash -r \'...
            '"run activate.m;'....
            ' result=VVT(false,false,true,true);'...
            ' exit\"'' & '],...
            IPlist{secondarynum}, basedir); %start VVT on all systems
        [~,~]=unix(cmd);
    end
    if numSecondaries~=4
        fprintf('--- Perform VVT on Primary ---\n');
        [~]=VVT(false,false,true,true); %just do the digital test
    end
    %----------- check if remote systems are still running ---------------%
    for secondaryNum = 1:numSecondaries
        fprintf('waiting for system %i to finish\n', secondaryNum)
        for timeoutCount = 1:30 %timeout of 30 seconds for remote VVT
            cmd = sprintf('ssh %s "ps -e |grep MATLAB"', IPlist{secondaryNum}); % probe remote systems and look to see if Matlab is still running
            [~,result]=unix(cmd);
            if  isempty(result)
                break
            end
            pause(1);
        end
        if timeoutCount == 30 %30 second timeout trying each system
            warning('VVT timed out on system')
        end
    end
    
else  % perform VVT remotely on just 1 secondary
    fprintf('--- S Secondary %i ---\n',varargin{1});
    cmd = sprintf(['ssh 10.10.%i.2 -X ''matlab -sd %s -nosplash -r \'...
        '"run activate.m;'...
        ' result=VVT(false,false,true,true);'...
        ' exit\"'' & '],...
        varargin{1}, basedir, varargin{1},testLogsDir); %start VVT on all systems
    
    [~,~]=unix(cmd);
end

fprintf('Results of vsv.multi.VVT can be found in: %s\n', testLogsDir)
