function VVTmultiSys(varargin)
%-------------------------------------------------------------------------%
%  Runs VVT across the secondary systems, either individually or sequentially %
%
%  usage:
%    VVTmultiSys    by itself it defaults to sequetially running VVT on all
%                   the salve systems and creating a single log file under
%                   Hardwaretest
%
%   VVTmultiSys(secondaryNum) runs VVT remotely on the secondary system and creates
%                         a single log file under the primary Hardwaretest
%                         directory as well as the secondary's hardwaretest
%                         directroy
%  Example:
%     VVTmultiSys(2)    runs VVT on secondary #2
%-------------------------------------------------------------------------%
warning('VVTmultiSys:deprecated', 'VVTmultiSys.m has been deprecated and should be replaced with vsv.multi.VVT');

c = fix(clock);
setenv('LD_LIBRARY_PATH', '/usr/lib/x86_64-linux-gnu/libexpat.so.1') %this supresses a libcrypto warning
basedir = [getenv('VERASONICS_VPF_ROOT') '/'];
testLogsDir = [basedir '/HardwareTest/TestResultLogs/'];
datetag = [ num2str(c(1)) num2str(c(2),'%02i') num2str(c(3),'%02i') '_' num2str(c(4),'%02i') num2str(c(5),'%02i') num2str(c(6),'%02i') ];
numSecondaries = 3; %TODO: autodetect
secondaryList = 4-[numSecondaries:-1:1];
disp('######################################')
disp('###   Verasonics MultiSystem       ###')
disp('###     Verification Test          ###')
disp('######################################')


if length(varargin)>1
    error('incorrect number of inputs, just a number 0 to 3 for the secondary or leave empty')
elseif isempty(varargin)
    disp('--- Performing VVT on all systems in parallel ---');
    logFile = [testLogsDir 'VVTmultiSys_' datetag '_log.txt'];
    for secondarynum = secondaryList
        fprintf('--- Perform VVT on Secondary %i ---\n',secondarynum);
        cmd = sprintf('ssh 10.10.%i.2 "killall -9 MATLAB"',secondarynum); %stop all running matlab sessions
        [~,result]=unix(cmd);
        cmd = sprintf('ssh 10.10.%i.2 -X ''matlab -sd %s -nosplash -r \"run activate.m;result=VVT(false,false,true);exit\"'' & ',secondarynum,basedir); %start VVT on all systems
        [~,result]=unix(cmd);
    end
    if secondaryList(1)~=0
        fprintf('--- Perform VVT on Primary ---\n');
        result=VVT(false,false,true); %just do the digital test
    end
    %----------- check if remote systems are still running ---------------%
    for a = secondaryList
        waiting = 1;
        timeoutCount = 1;
        fprintf('waiting for system %i to finish\n',a)
        while waiting
            cmd = sprintf('ssh 10.10.%i.2 "ps -e |grep MATLAB"',a); % probe remote systems and look to see if Matlab is still running
            [~,result]=unix(cmd);
            if  ~isempty(result)
                waiting = 0;
            end
            if timeoutCount == 20 %try 20 probes and if matlab is still running, give up and check the next system
                warning('VVT timed out on system')
                waiting = 0;
            end
            timeoutCount = timeoutCount + 1;
        end
    end
    
else  % perform VVT remotely on just 1 secondary
    logFile = [testLogsDir 'VVT_Secondary' num2str(varargin{1}) '_' datetag '_log.txt'];
    fid=fopen(logFile,'a+');
    fprintf('--- S Secondary %i ---\n',varargin{1});
    cmd = sprintf('ssh 10.10.%i.2 -X ''matlab -sd %s -nosplash -r \"run activate.m;result=VVT(false,false,true);exit\"''  ',varargin{1},basedir); %just do the digital test
    [~,result]=unix(cmd);
    fprintf(fid,'============  VVT on Secondary %i ===========\n',varargin{1});
    fwrite(fid,result);
    fclose(fid);
end

fprintf('Results of VVTmultiSys can be found: %s%s\n', testLogsDir, logFile)
