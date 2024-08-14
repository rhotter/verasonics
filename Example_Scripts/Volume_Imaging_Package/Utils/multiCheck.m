function secondaryReturn = multiCheck(varargin)
%function multiCheck(mode)
%
%  mode 0 = simple check (just before a script runs, just check that all VDAS boards are online)
%  mode 1 = thorough check (when matlab session starts)
%
% MultiSystem Pre-check
%
% Performs all neccesary HW checks before starting a multisystem script to
% make sure that all the HW is configured and working properly
%
% returns the secondary# configuration 3 or 4 or 0 if error
warning('multiCheck:deprecated', 'multiCheck.m has been deprecated and should be replaced with vsv.multi.configCheck');

% x. Primary check:
%   . linux kernel 4.12.8?
%   . what is the version of the mellanox drivers
%   . # of mellanox cards connected?
%   . what is their firmware?
%   --- primary specific ---
%   . Are cards UP?
%   . check tmux version
%   -----------------------
%   . card temperature?
%   . is there an HCA connected to the Primary?
%   . check license
%   . check matlab version
%   . --> determine configuration
%
% x. Secondary Checks:
%   . can we ping the secondary?
%   . can we ssh to the secondary?
%   . linux kernel 4.12.8?
%   . check mellanox cards: FW? driver versions
%   . can we communicate with the HCA and Vantage?
%   . check SW&HW versions
%   . check matlab version
%   . check license
% x. If everything is up to date and working secondaryReturn =1
%   .
if length(varargin) ==1
    mode = varargin{1};
else
    mode = 0;
end
if mode
    disp('---> Full MultiCheck <---')
else
    disp('MultiSystem Pre-Check')
end

secondaryReturn = 1;  %start with the return flag true
VERBOSE = 1;
PASS = 'verasonics';  %need to have root password for some checks (remove in future)
SWDIRECTORY = getenv('VERASONICS_VPF_ROOT');
%----
MINLINUXKERNEL = '4.15.0';         %
MAXLINUXKERNEL = '4.19.1';         % as of 5/29/2020 this is the newest linux kernel supported by WinDriver
MINMATLABVERSION = '9.7.0';        % Matlab 2019b
MINOFEDVERSION =  '5.0-2.1.8.0';   % as of 5/29/2020 this is the latest available Mellanox driver
MINMLXFWVERSION = '16.27.2008';    % as of 5/29/2020 this is the latest available FW
MLXMAXTEMP = 55;                   % max temperature (degC)of mellanox card (from ConnectX®-5 VPI Adapter Cards User Manual)
SWVERSION = [4 4 0];
DRIVERVERSION =  '12.9.0';         % as of 5/29/2020 this is the latest WinDriver Version supported
TMUXVERSION = 'next-3.1';
%----
ipList = {'10.10.0.2','10.10.1.2','10.10.2.2','10.10.3.2'};
warning off backtrace;

%=== Primary Check ===%
fprintf('==================================================================\n');
[~,result]=unix('ibv_devinfo -l');
nSecondarys = sscanf(result, '%d HCAs');  %get the number of secondarys
fprintf('Checking Primary System and %i RDMA connections\n',nSecondarys)
if mode
    %--- Check Linux kernel ---%
    [~,result] = unix('uname -r');
    if (verStr2verDec(result)>= verStr2verDec(MINLINUXKERNEL))&&(verStr2verDec(result)<= verStr2verDec(MAXLINUXKERNEL))
        if VERBOSE, fprintf('Linux Kernel %s is within compatible range %s to %s\n', strtrim(result),MINLINUXKERNEL,MAXLINUXKERNEL); end
    else
        warning('Incorrect Linux Kernel.  You are using: %s but it need to be between %s and %s.  Check GRUB configuration.',result, MINLINUXKERNEL, MAXLINUXKERNEL)
    end
    %--- Check Mellanox Drivers (mode 1) ---%
    [~,result] = unix('ofed_info -n');
    if verStr2verDec(result)>= verStr2verDec(MINOFEDVERSION)
        if VERBOSE, fprintf('Mellanox OFED driver %s is greater than or equal to minimum version %s\n', strtrim(result),MINOFEDVERSION); end
    else
        warning('OFED version mismatch.  You have version: %s but version %s should be installed.',result, MINOFEDVERSION)
    end
    %-.-.-.-.- Primary Specific Checks -.-.-.-.-.-%
    %--- Tmux version check ---%
    [~, result]=unix('tmux -V');
    tmuxVersion = strsplit(result);tmuxVersion=tmuxVersion{2};
    if strcmp(tmuxVersion,TMUXVERSION)
        fprintf('Tmux Version is up to date: %s\n',tmuxVersion);
    else
        warning('Incorrect Tmux Version: %s, but it should be %s\n',tmuxVersion, TMUXVERSION)
    end
    
    %--- LosecondaryReturn for Mellanox Cards with ibvdev2netdev (mode 1) ---%
    [~,result] = unix('ibdev2netdev');
    nPorts = length(strfind(result,'mlx'));
    portsUp = length(strfind(result,'Up'));
    if nPorts == portsUp
        if VERBOSE, fprintf('%d Mellanox Ports detected and all are active.\n%s',nPorts,result); end
    else
        warning(['%d Mellanox Ports detected but only %d are active\n' result],nPorts,portsUp)
    end
    
    
    %--- Check FW on each one (mode 1) ---%
    [~,result] = unix('ibv_devinfo | grep fw_ver:');
    if (length(strfind(result,MINMLXFWVERSION)))==nPorts
        if VERBOSE, fprintf('All Mellanox cards firmware is up-to-date: %s\n', MINMLXFWVERSION); end
    else
        warning('Mellanox card firmware is not all up to date should be %s:\n %s', MINMLXFWVERSION, result)
    end
    
    %--- Check card temperature (mode 0) ---%
    [~,result] = unix('ls /dev/mst | grep -v "\.1"'); % get card list
    cardList = strsplit(result);cardList=cardList(1:end-1);% edit card list
    for a = 1:length(cardList)
        cmd = sprintf('sudo mget_temp -d /dev/mst/%s',cardList{a});
        [~, result] = unix(cmd);
        cardTemp(a) = sscanf(result,'%d');
        if cardTemp(a)<MLXMAXTEMP
            if VERBOSE,        fprintf('%s temperature: %d degC\n',cardList{a}, cardTemp(a)); end
        else
            warning('%s temperature: %d degC exceeds recommended limit of %d degC\n',cardList{a}, cardTemp(a), MLXMAXTEMP);
        end
    end
    
    %--- Matlab version check ---%
    matlabVersion = strsplit(version);matlabVersion = matlabVersion{1};
    if (verStr2verDec(matlabVersion)>= verStr2verDec(MINMATLABVERSION))
        fprintf('Matlab Version %s is supported (%s or newer).\n',matlabVersion,MINMATLABVERSION);
        primaryMatlabVersion = matlabVersion;
    else
        warning('Incorrect Matlab Version: %s, but it should be %s or newer.\n', matlabVersion, MINMATLABVERSION)
        secondaryReturn = 0 ;
    end
    
else
    %    hardwareReset;
end

%-------------------------------------%
SysConfig = hwConfigCheck(1,1);
%--- HW config check ---%
if SysConfig.SWconfigFault
    warning('SW Config Fault Detected')
elseif SysConfig.HWconfigFault
    warning(['HW Config Fault Detected.  please run VVT on on this system with the command ''VVTmultiSys(' num2str(a) ')'''])
elseif SysConfig.FPGAconfigFault
    warning('FPGA Config Fault Detected')
end
%--- License Check ---%
if (~SysConfig.multisystemRDMA)
    error('Multisystem license is invalid or not found')
end

if SysConfig.VDAS; fprintf('VDAS Connected with %d Boards\n',sum(SysConfig.AcqSlots));
    if contains(SysConfig.DriverVersion,DRIVERVERSION)
        fprintf('Correct WinDriver Version: %s\n', strtrim(SysConfig.DriverVersion));
    else, warning('Incorrect WinDriver version %s  Should be %s\n',SysConfig.DriverVersion, DRIVERVERSION );
    end
else
    warning('VDAS not enabled.  Must run in Simulate Mode');
end
optical8725 = pexEeproms('8725','check','./Hal/lib/pex/PEX8725_optical.bin');
optical8748 = pexEeproms('8748','check','./Hal/lib/pex/PEX8748_optical.bin');
if strcmpi(optical8725,'success')&&strcmpi(optical8748,'success')
    disp('System Configured for Optical')
else
    warning('System Configured for Copper - Please reconfigure for Optical')
end


%% Secondary Check
%NOTE: 'unset LD_LIBRARY_PATH &&' is needed before every ssh to suppress a
%warning: ssh: /usr/local/MATLAB/R2017b/bin/glnxa64/libcrypto.so.1.0.0: no version information available (required by ssh)
for a = (5-nSecondarys):length(ipList)
    fprintf('========================================================================================\n');
    %--- Ping address ---
    [~,result]=unix(['ping -c 1 ' ipList{a} ' >&1 >/dev/null; echo $?']);
    notReachable=sscanf(result,'%d');
    if ~notReachable %if reachable
        fprintf('Connecting to secondary:%d, IP: %s\n',a -1, ipList{a})
        %==== Start Checking Secondary n ===%
        if mode %run in mode 1, full check
            %--- Check Linux kernel ---%
            [~,result] = unix(['unset LD_LIBRARY_PATH && ssh ' ipList{a} ' uname -r;']);
            if (verStr2verDec(result)>= verStr2verDec(MINLINUXKERNEL))&&(verStr2verDec(result)<= verStr2verDec(MAXLINUXKERNEL))
                if VERBOSE, fprintf('Linux Kernel %s is within compatible range %s to %s\n', strtrim(result),MINLINUXKERNEL,MAXLINUXKERNEL); end
            else
                warning('Incorrect Linux Kernel.  You are using: %s but it need to be between %s and %s.  Check GRUB configuration.',result, MINLINUXKERNEL, MAXLINUXKERNEL)
            end
            %---  Check NFS ---%
            [~,result] = unix(['unset LD_LIBRARY_PATH && ssh ' ipList{a} ' df -h |grep  /home/verasonics/cloud']);
            if ~isempty(result)
                if VERBOSE, fprintf('NFS is setup and working.\n'); end
            else
                warning('NFS is not functioning properly or not setup. Please restart the system.  If the problem persists please contact Verasonics Support.')
            end
            %---  Check Mellanox Drivers ---%
            fprintf('----------------------------------------------------------------------------------------\n');
            [~,result] = unix(['unset LD_LIBRARY_PATH && ssh ' ipList{a} ' ofed_info -n']);
            if verStr2verDec(result)>= verStr2verDec(MINOFEDVERSION)
                if VERBOSE, fprintf('Mellanox OFED driver %s is greater than or equal to minimum version %s\n', strtrim(result),MINOFEDVERSION); end
            else
                warning('OFED version mismatch.  You have version: %s but version %s should be installed.',result, MINOFEDVERSION)
            end
            
            %--- Check Mellanox FW version ---%
            [~,result] = unix(['unset LD_LIBRARY_PATH && ssh ' ipList{a} ' ibv_devinfo | grep fw_ver:']);
            result = strtrim(result(10:end));
            if verStr2verDec(result)>= verStr2verDec(MINMLXFWVERSION)
                if VERBOSE, fprintf('Mellanox card firmware is up to date. Version: %s\n', MINMLXFWVERSION); end
            else
                warning('Mellanox card firmware is not all up to date:\n %s', result)
            end
            
            %--- Check card Temperature ---%
            cmd = sprintf('unset LD_LIBRARY_PATH && ssh %s "echo %s | sudo mget_temp -d /dev/mst/mt4119_pciconf0"', ipList{a}, PASS);
            [~, result] = unix(cmd);
            cardTemp = sscanf(result,'%d');
            if ~isempty(cardTemp)
                if cardTemp<MLXMAXTEMP
                    if VERBOSE,        fprintf('%s temperature: %d degC\n',cardList{1}, cardTemp); end
                else
                    warning('%s temperature: %d degC exceeds recommended limit of %d degC\n',cardList{1}, cardTemp, MLXMAXTEMP);
                end
            else
                warning('could not run mget_temp on remote system')
            end
            
            %--- Matlab version check ---%
            fprintf('----------------------------------------------------------------------------------------\n');
            disp('Starting Matlab remotely... please wait')
            cmd = sprintf(['unset LD_LIBRARY_PATH && ssh %s "cd %s; matlab -r \\"'...
                'run activate.m;'...
                'disp(version);'...
                'SysConfig = hwConfigCheck(1,1);'...
                'disp(SysConfig);'...
                'com.verasonics.hal.hardware.Hardware.closeHardware;',...
                'exit \\""'], ipList{a},SWDIRECTORY);
            [~,result] = unix(cmd);
            temp = strsplit(result(strfind(result,SWDIRECTORY)+length(SWDIRECTORY+2):end));
            matlabVersion = temp{2};
            if (verStr2verDec(matlabVersion)>= verStr2verDec(MINMATLABVERSION))
                fprintf('Matlab Version %s is %s or newer.\n',matlabVersion,MINMATLABVERSION);
            elseif (primaryMatlabVersion~=matlabVersion)
                warning('Secondary system is not running the same version of Matlab as the Primary system.  This could cause unpredictable errors.');
            else
                warning('Incorrect Matlab Version: %s, but it should be %s or newer.\n',matlabVersion, MINMATLABVERSION)
                secondaryReturn = 0 ;
            end
            
            %--- HW config check ---%
            SWconfigFault=     eval(result(strfind(result,'SWconfigFault:')+15));
            HWconfigFault=     eval(result(strfind(result,'HWconfigFault:')+15));
            FPGAconfigFault=   eval(result(strfind(result,'FPGAconfigFault:')+17));
            SWversion =        eval(result(strfind(result,'SWversion:')+[11:20]));
            DriverVersion = strtrim(result(strfind(result,'DriverVersion:')+[16:21]));
            VDAS =          str2num(result(strfind(result,'VDAS:')+[6]));
            AcqSlots =         eval(result(strfind(result,'AcqSlots:')+[10:26]));
            %-------------------------------------%
            if SWconfigFault
                warning('SW Config Fault Detected')
            elseif HWconfigFault
                warning(['HW Config Fault Detected.  Please run VVT on on this system with the command ''VVTmultiSys(' num2str(a) ')'''])
            elseif FPGAconfigFault
                warning('FPGA Config Fault Detected.  Please run ''reprogramHardware'' from the matlab command line on each system')
            end
            %-------------------------------------%
            if SWversion(1)*4^3+SWversion(2)*4^2+SWversion(3) == SWVERSION(1)*4^3+SWVERSION(2)*4^2+SWVERSION(3)
                fprintf('Correct SW Version: %d.%d.%d\n', SWversion);
            else
                warning('Incorrect SW version %d.%d.%d  Should be %d.%d.%d\n',SWversion(1), SWversion(2), SWversion(3), SWVERSION(1), SWVERSION(2), SWVERSION(3) );
                secondaryReturn = 0;
            end
            if VDAS
                fprintf('VDAS Connected with %d Boards\n',sum(AcqSlots));
                if contains(DriverVersion,DRIVERVERSION)
                    fprintf('Correct WinDriver Version: %s\n', strtrim(DriverVersion));
                else
                    warning('Incorrect WinDriver version %s  Should be %s\n',DriverVersion, DRIVERVERSION );
                    secondaryReturn=0;
                end
            end
            disp('Checking PCIe:Copper/Optical settings.  Starting Matlab remotely... please wait')
            cmd = sprintf(['unset LD_LIBRARY_PATH && ssh %s "cd %s; matlab -r \\"'...
                'activate;'...
                'secondaryReturn=pexEeproms(''8725'', ''check'', ''./Hal/lib/pex/PEX8725_optical.bin'');'...%                'pexEeproms(''8748'',''check'',''./Hal/lib/pex/PEX8748_optical.bin'');'...%'disp([''PCIeConfig==Optical:'' secondaryReturn ]);',...
                'fprintf(''8725:'');',...
                'disp(secondaryReturn);',...
                'secondaryReturn=pexEeproms(''8748'', ''check'', ''./Hal/lib/pex/PEX8748_optical.bin'');'...%                'pexEeproms(''8748'',''check'',''./Hal/lib/pex/PEX8748_optical.bin'');'...%'disp([''PCIeConfig==Optical:'' secondaryReturn ]);',...
                'fprintf(''8748:'');',...
                'disp(secondaryReturn);',...
                'exit \\""'], ipList{a},SWDIRECTORY);
            [~,result] = unix(cmd);
            if (contains(result,'8725:Success'))&&(contains(result,'8748:Success'))
                disp('System Configured for Optical')
            else
                warning('System Configured for Copper - Please reconfigure for Optical')
            end
        else
            % quick check only tests
            [~,result] = unix(['unset LD_LIBRARY_PATH && ssh ' ipList{a} ' ' SWDIRECTORY '/System/hwstatus.linux64 GetBoardStatus BKP']);
            if contains(result,'Error')
                warning(result)
                secondaryReturn = 0;
            else
                secondaryReturn = secondaryReturn + 1;
                disp('Vantage System Online');
            end
        end
        
    else
        secondaryReturn = 0;
        warning('Secondary:%d, IP: %s NOT REACHABLE\n',a -1, ipList{a});
    end
end

function verDec = verStr2verDec(versionString)
%convert Major.Minor.Revision to number for comparison
[verArray]=sscanf(versionString,'%i.');
N=min(length(verArray),3); %only do Major, Minor, rev (3)
verDec=0;
for a = N:-1:1
    verDec = verArray(a)*4^a + verDec;
end
return
