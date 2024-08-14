function result = configCheck(varargin)
% function result = vsv.multi.configCheck(testLevel)
%
% testLevel:
%     0 or empty is a quick test to make sure that all hardware is online (default)
%     1 is a full configuration test for all systems and will
%     produce an output log file for debugging
%
% returns:
%   vector with length of number of connected systems, 0 if offline, 1 if
%   system is ready

% Copyright 2001-2021 Verasonics, Inc.  All world-wide rights and remedies
% under all intellectual property laws and industrial property laws are
% reserved.  Verasonics Registered U.S. Patent and Trademark Office.

%-------------------------------------------------------------------------%
config.MINLINUXKERNEL = '4.15.0';         % as of 3/9/2021 this is the oldest linux kernel supported by WinDriver
config.MAXLINUXKERNEL = '4.19.1';         % as of 3/9/2021 this is the newest linux kernel supported by WinDriver
config.OFEDVERSION =  '5.0-2.1.8.0';      % as of 3/9/2021 this is the latest available Mellanox driver
config.MINTMUXVERSION = 3.1;              % as of 3/9/2021 this is the oldest version of tmux that has been tested with the system
config.MAXTMUXVERSION = 3.3;              % as of 3/9/2021 this is the oldest version of tmux that has been tested with the system
config.MINMLXFWVERSION = '16.27.2008';    % as of 3/9/2021 this is the latest available FW16.27.2008
config.MLXMAXTEMP = 55;                   % max temperature (degC) of mellanox card (from ConnectX®-5 VPI Adapter Cards User Manual)
config.WINDRIVERVERSION =  1500;          % as of 3/9/2021 this is the latest WinDriver Version supported

%-------------------------------------------------------------------------%

if nargin==0
    testLevel=0;
else
    testLevel=varargin{1};
end

[nPorts,portsUp] = checkPorts('localhost');
if nPorts~=portsUp
    warning('vsv:multi:configCheck','!!! The number of ports connected (%i) does not match the number of ports configured (%i) !!!\n Verify that all secondary systems are properly connected.', portsUp, nPorts);
end
numSecondaries = vsv.multi.getNumSecondaries();
fprintf('Testing system with primary + %i Secondaries\n',numSecondaries);

switch numSecondaries
    case 4
        ipList = {'localhost', '10.10.0.2', '10.10.1.2', '10.10.2.2', '10.10.3.2'};  % 5 system configuration
    case 3
        ipList = {'localhost', '10.10.1.2', '10.10.2.2', '10.10.3.2'}; % 4 system configuration
    case 2
        ipList = {'localhost', '10.10.1.2', '10.10.2.2'}; % 3 system configuration
    case 1
        ipList = {'localhost', '10.10.1.2'}; % 2 system configuration
    otherwise
        error('The system network IP addresses are in a non-standard configuration.  Please contact Verasonics support.')
end
sshTestResult = ones(1,numSecondaries+1);
vantageStatuResult = ones(1,numSecondaries+1);
configurationResult = ones(1,numSecondaries+1);
matlabTestResult = ones(1,numSecondaries+1);

%% ALL TESTS %%
% 1. Check Communications between systems
disp('=== Check SSH connection to systems ===')
for a = 1:length(ipList)
    sshTestResult(a) = sshTest(ipList{a});
end

if ~all(sshTestResult)
    error('Cannot ssh to one or more systems.  Check that all systems have booted fully.')
end

% 2. Quick check to make sure systems are online
disp('=== Perform Quick Vantage Test ===')
if (numSecondaries == 4)  %5 system configuration
    HWTestList = 2:length(ipList);
else
    HWTestList = 1:length(ipList);
end

% check nfs on secondaries
for a = HWTestList(2:end)
    nfsStatuResult(a) = checkNfsStatus(ipList{a},nPorts);
end

for a = HWTestList
    vantageStatuResult(a) = checkHardwareStatus(ipList{a},nPorts);
end

if testLevel>0
    % Open file for logging
    fname = sprintf('multisystemConfigCheck%s.log',datestr(now,30));
    fid = fopen([getenv('VERASONICS_VPF_ROOT') '/' fname],'w');
    cleanup = onCleanup(@()cleanupFunc(fid));
    fprintf('opening %s for logging output\n', fname);
    
    % 4. Perform Configuration Tests
    for a = 1:length(ipList)
        assertReport(true,sprintf('=== Testing Configuration of %s ===',ipList{a}),'',fid);
        configurationResult(a)=checkConfiguration(ipList{a},config,fid);
    end
    
    % 5. Check hwconfig() through Matlab
    for a = HWTestList
        assertReport(true,sprintf('=== Testing Hardware of %s ===',ipList{a}),'',fid);
        matlabTestResult(a)=matlabTests(ipList{a},fid);
    end
    
end

result = sshTestResult & vantageStatuResult & configurationResult & matlabTestResult;
if ~all(result)
    warning('vsv:multi:configCheck','\nmultisystem config/hardware check failed on %s', ipList{find(~result)})
    if testLevel>0
        warning('vsv:multi:configCheck','Check log file ''%s'' for more details and send to support@verasonics if you need additional help debugging your system', fname)
    end
end

end

%% Tests
% Test functionality to ssh to all secondary systems
function sshTestResult = sshTest(ipAddress)
[status,~] = unix(sprintf('ssh -q -o BatchMode=yes -o StrictHostKeyChecking=no -o ConnectTimeOut=5 %s ''exit 0'' && echo $?', ipAddress));
if status==0
    fprintf('connection to %s successful\n', ipAddress)
    sshTestResult= true;
else
    warning('ssh connection to %s failed', ipAddress);
    sshTestResult = false;
end
end

% Test NFS Status (bash)
function nfsStatusResult = checkNfsStatus(ipAddress,nPorts)
cmd = sprintf('df -h |grep cloud');
resultStr = multiCommand(ipAddress, cmd, 15); % 15 second timeout
if isempty(resultStr)
    warning('NFS not online for secondary system.  Restablishing connection.')
    cmd = sprintf('"echo verasonics | sudo -S mount -t nfs %s\b1:/home/verasonics/cloud ~/cloud"',ipAddress);
    resultStr = multiCommand(ipAddress, cmd, 15); % 15 second timeout
    %try again
    cmd = sprintf('df -h |grep cloud');
    resultStr = multiCommand(ipAddress, cmd, 15); % 15 second timeout
    if isempty(resultStr)
            nfsStatusResult = false;
    end
else
    nfsStatusResult = true;
    fprintf('%s - NFS System Online\n', ipAddress);
end
end


% Test Hardware Status (bash)
function vantageStatusResult = checkHardwareStatus(ipAddress,nPorts)
cmd = sprintf('%s/System/hwstatus.linux64 GetBoardStatus BKP', getenv('VERASONICS_VPF_ROOT'));
resultStr = multiCommand(ipAddress, cmd, 15); % 15 second timeout
if isempty(resultStr)||contains(resultStr,'Error')
    warning('vsv:multi:configCheck',resultStr)
    vantageStatusResult = false;
    if (strcmp(ipAddress,'localhost')&&(nPorts~=5)) %for 5 system config, localhost is not connected to a vantage
        warning('vsv:multi:configCheck:checkHardwareStatus','%s - Vantage System Off-line.  If this is a 5 system configuration, ignore this warning.\n', ipAddress);
    end
else
    vantageStatusResult = true;
    fprintf('%s - Vantage System Online\n', ipAddress);
end
end

% Test System Configuration
function configurationResult = checkConfiguration(ipAddress,config,fid)
if strcmp(ipAddress,'localhost') % Primary only test
    % Check Tmux version
    tmuxVersion = getTmuxVersion();
    tmuxVersionOK = (tmuxVersion>=config.MINTMUXVERSION)&&(tmuxVersion<=config.MAXTMUXVERSION);
    assertReport(tmuxVersionOK,...
        sprintf('Tmux version within range.  It is: %.1f', tmuxVersion),...
        sprintf('Tmux version outside range.  It is: %.1f but should be between %.1f and  %.1f',tmuxVersion, config.MINTMUXVERSION, config.MAXTMUXVERSION),fid);
else
    tmuxVersionOK=true;
end



% Check card temperature
mlxCardTemp = getMlxCardTemp(ipAddress);
mlxCardTempOK = all(mlxCardTemp<config.MLXMAXTEMP);
for a = 1:length(mlxCardTemp)
    assertReport(mlxCardTemp(a)<config.MLXMAXTEMP, ...
        sprintf('Mellanox Card Temperature OK: %.1f ', mlxCardTemp(a)), ...
        sprintf('Mellanox Card Temperature too high: %.1f ', mlxCardTemp(a)),fid);
end

% Check ports
[nPorts,portsUp] = checkPorts(ipAddress);
portsOK = nPorts==portsUp;
assertReport(portsOK,...
    sprintf('All ports functional'),...
    sprintf('Not all ports are functional.  %i are available but only %i are up ', nPorts, portsUp), fid);

% Check mellanox fw version
mlxFwVersionStr = getmlxFwVersion(ipAddress);
mlxFwVersionNum = reshape(sscanf(mlxFwVersionStr,'fw_ver: %i.%i.%i\n'), 3, nPorts)';

for b = 1:nPorts
    mlxFwVersionOK(b) = verArr2verDec(mlxFwVersionNum(b,:))>=verStr2verDec(config.MINMLXFWVERSION);
    assertReport(mlxFwVersionOK(b),...
        sprintf('Mellanox firmware version up to date. It is: %i.%i.%i',mlxFwVersionNum(b,:)'),...
        sprintf('Mellanox firmware version incorrect. It is: %i.%i.%i it should be greater than %s ',mlxFwVersionNum(b,:), config.MINMLXFWVERSION),fid);
end

% Check Ubuntu kernel version
kernelVersionStr = getKernelVersion(ipAddress);
kernelVersionNum = verStr2verDec(kernelVersionStr);
kernelVersionOK = (kernelVersionNum>verStr2verDec(config.MINLINUXKERNEL)) && (kernelVersionNum<verStr2verDec(config.MAXLINUXKERNEL));
assertReport(kernelVersionOK,...
    sprintf('Ubuntu kernel within range.  It is: %s',kernelVersionStr),...
    sprintf('Ubuntu kernel outside range.  It is: %s but should be less than %s and greater than %s',kernelVersionStr, config.MAXLINUXKERNEL, config.MINLINUXKERNEL),fid);

% Check OFED driver version
ofedDriverVersion=getOFEDDriverVersion(ipAddress);
ofedDriverVersionOK = strcmp(ofedDriverVersion,config.OFEDVERSION);
assertReport(ofedDriverVersionOK,...
    sprintf('OFED Driver version up to date.  It is: %s ', ofedDriverVersion),...
    sprintf('OFED Driver version incorrect.  It is: %s, it should be greater than %s', ofedDriverVersion, config.OFEDVERSION),fid);

% Check windriver version
winDriverVersion = getWinDriverVersion(ipAddress);
winDriverVersionOK = winDriverVersion==config.WINDRIVERVERSION;
assertReport(winDriverVersionOK,...
    sprintf('WinDriver version correct.  It is: %i', winDriverVersion),...
    sprintf('WinDriver version incorrect.  It is: %i, it should be %i', winDriverVersion, config.WINDRIVERVERSION),fid);


configurationResult = all([tmuxVersionOK mlxCardTempOK portsOK mlxFwVersionOK kernelVersionOK ofedDriverVersionOK winDriverVersionOK]);
end

% Test hwConfigCheck()
function matlabResults = matlabTests(ipAddress,fid)
persistent primaryMatlabVersion
fprintf('----------------------------------------------------------------------------------------\n');
fprintf('Starting Matlab hwConfigCheck() tests on system %s ... please wait this can take up to 40 seconds\n', ipAddress)
cmd = sprintf(['"cd %s; matlab -nodesktop -r \\"'...
    'activate;',...
    'vsv.multi.remoteMatlabTest();',...
    'exit \\""'], getenv('VERASONICS_VPF_ROOT'));

resultStr = multiCommand(ipAddress, cmd, 40); %40 second timeout
disp(resultStr);%BWC
remoteCmdOk = ~isempty(resultStr);
assertReport(remoteCmdOk,'','Remote commands to Matlab Failed',fid);

if remoteCmdOk
    %--- get Matlab Version ---%
    matlabVersionStr = readString(resultStr,'MatlabVersion');
    matlabVersionNum =  verStr2verDec(matlabVersionStr);
    
    %--- get configurations ---%
    SWconfigFault=readScalar(resultStr,'SWconfigFault');
    HWconfigFault=readScalar(resultStr,'HWconfigFault');
    FPGAconfigFault=readScalar(resultStr,'FPGAconfigFault');
    %SWversion=readVector(resultStr,'SWversion');
    %WinDriverVersion=readString(resultStr,'DriverVersion');
    VDAS=readScalar(resultStr,'VDAS');
    AcqSlots =readVector(resultStr,'AcqSlots');
    
    %-- Get license status
    multiSystemRDMA=readScalar(resultStr,'multisystemRDMA');
    hwLicense=readScalar(resultStr,'hwLicense');
    reconLicense=readScalar(resultStr,'reconLicense');
    
    %-- Get PCIE configuration
    IOPanelType = readString(resultStr,'IOPanelType');
    IOPanelTypeOK = contains(IOPanelType,'CopperOptical');
    hcaOptical = readScalar(resultStr,'hcaOptical');
    ioPanelOptical = readScalar(resultStr,'ioPanelOptical');

    %=== Actual Tests ===%
    %-- Check that the same version of Matlab is being run on all systems
    if strcmp(ipAddress,'localhost')
        primaryMatlabVersion = matlabVersionNum;
        matlabVersionOK=true;
    else
        matlabVersionOK = primaryMatlabVersion==matlabVersionNum;
        assertReport(matlabVersionOK,...
            'Secondary system is running the same version of Matlab as the Primary system.',...
            'Secondary system is not running the same version of Matlab as the Primary system.  This could cause unpredictable errors.',fid)
    end
    
    %-- Check for any Faults
    assertReport(~SWconfigFault, '','--- SW Config Fault Detected ---',fid);
    if ~isnan(HWconfigFault)
        assertReport(~HWconfigFault, '','--- HW Config Fault Detected ---',fid);
    else
        HWconfigFault = 0;
    end
    if ~isnan(FPGAconfigFault)
        assertReport(~FPGAconfigFault, '','--- FPGA Config Fault Detected ---',fid);
        assertReport(~FPGAconfigFault, '','An out-of-date or unrecognized version of FPGA code has been found in the HW system.',fid);
    else
        FPGAconfigFault = 0;
    end
    
    %-- Check that all necessary licenses are installed/detected
    assertReport(multiSystemRDMA, '','No valid multisystem license detected', fid);
    assertReport(hwLicense, '','No valid hardware license detected', fid);
    assertReport(reconLicense, '','No valid reconstruction license detected', fid);
    
    %-- Check that VDAS is functional
    assertReport(VDAS, ...
        sprintf('VDAS Connected with %d Boards',sum(AcqSlots)), ...
        'VDAS not detected', fid);
    assertReport(true, resultStr,'', fid);
    
    %-- Check that the PCIE is configured for optical cable
    fprintf(fid,'\n');
    assertReport(IOPanelTypeOK,...
        'Host Controller Card is CopperOptical',...
        sprintf('Host Controller Card is NOT CopperOptical but %s', IOPanelType'),fid);
    assertReport(hcaOptical,...
        'Host Controller Card Configured for Optical',...
        'Host Controller Card NOT Configured for Optical',fid);
    assertReport(ioPanelOptical,...
        'ioPanel Configured for Optical',...
        'ioPanel NOT Configured for Optical',fid);

    matlabResults = all([matlabVersionOK ~SWconfigFault ~HWconfigFault ~FPGAconfigFault multiSystemRDMA hwLicense reconLicense VDAS IOPanelTypeOK hcaOptical ioPanelOptical]);
else
    matlabResults = false;
end

end

% Check Symbolic Link for license
function symbolicLink = checkSymbolicLink()
symbolicLink = false;
[a, b]=system('stat -c ''%F'' licenseMgr.p');
if (a==0)
    symbolicLink = contains(b,'symbolic');
end
end
%% Get Values Functions

% Obtain the Linux Kernel version and return it as a numeric for comparison
function kernelVersionStr = getKernelVersion(ipAddress)
kernelVersionStr = multiCommand(ipAddress, 'uname -r',5); % 5 second timeout
end

% obtain the OFED driver version and return as a string
function ofedVer = getOFEDDriverVersion(ipAddress)
ofedVer = multiCommand(ipAddress, 'ofed_info -n',5); % 5 second timeout
end

% Obtain the Tmux version and return a string (only for primary)
function tmuxVersion = getTmuxVersion()
[~, result]=unix('tmux -V');
tmuxVersion = strsplit(result);
tmuxVersion = sscanf(tmuxVersion{2},'%f');
end

% Get Mellanox card firmware version
function mlxFwVersion = getmlxFwVersion(ipAddress)
mlxFwVersion = multiCommand(ipAddress, 'ibv_devinfo | grep fw_ver:',5); % 5 second timeout
end

% Check Mellanox card temperature
function mlxCardTemp = getMlxCardTemp(ipAddress)
result = multiCommand(ipAddress, 'ls /dev/mst | grep -v "\.1"', 5); % 5 second timeout
cardList = strsplit(result);
mlxCardTemp = zeros(1,length(cardList));
for a = 1:length(cardList)
    result = multiCommand(ipAddress, sprintf('sudo mget_temp -d /dev/mst/%s',cardList{a}),5); % 5 second timeout
    mlxCardTemp(a) = str2double(result);
end
end

% Check WinDriver Version
function winDriverVersion = getWinDriverVersion(ipAddress)
result = multiCommand(ipAddress, 'lsmod |grep windrvr',5); % 5 second timeout
winDriverVersion = sscanf(result,'windrvr%i');
end

%% --- General Functions ------------------------------------------------%%


% Send command to ip address. Return the same result whether it be used on a local system or remote
function returnStr = multiCommand(ipAddress, commandStr, timeoutVal)
[status,result] = unix(sprintf('timeout %i ssh %s %s;', timeoutVal, ipAddress, commandStr));
if (status~=0)
    if isempty(result)
        warning('vsv:multi:configCheck:timeout','timeout trying to run: %s\n', commandStr)
    else
        warning('vsv:multi:configCheck:err','error trying to run: %s\n %s', commandStr, result)
    end
    returnStr = '';
else
    returnStr = strip(result);
end
end

% Convert Major.Minor.Revision to number for comparison
function verDec = verStr2verDec(versionString)
[verArray]=sscanf(versionString,'%i.');
N=min(length(verArray),3); %only do Major, Minor, rev (3)
verDec=0;
for a = N:-1:1
    verDec = verArray(a)*4^(N-a) + verDec;
end
end

% Convert Major.Minor.Revision to number for comparison
function verDec = verArr2verDec(versionArray)
N=min(length(versionArray),3); %only do Major, Minor, rev (3)
verDec=0;
for a = N:-1:1
    verDec = versionArray(a)*4^a + verDec;
end
end

% Check the available ports on the Mellanox cards
function [nPorts,portsUp] = checkPorts(ipAddress)
result = multiCommand(ipAddress, 'ibdev2netdev',1);
nPorts = length(strfind(result,'mlx'));
portsUp = length(strfind(result,'Up'));
end

% Use regexp to find scalars in the returned string
function value = readScalar(inputStr, param)
value=str2double(regexp(inputStr, sprintf('(?:%s:\\s*)(\\d*)+',param),'once','tokens'));
end

% Use regexp to find vectors in the returned string
function value = readVector(inputStr, param)
value=str2num(cell2mat(regexp(inputStr, sprintf('(?<=%s:\\s*\\[)((\\d.\\s*)+\\d)',param),'once','tokens')));
end

% Use regexp to find string in the returned string
function value = readString(inputStr, param)
value=cell2mat(regexp(inputStr, sprintf('(?:%s: )(.*)',param),'once','tokens','dotexceptnewline'));
end

% Custom assert function to report only warning and not stop script
function assertReport(test, okMsg, errMsg, fid)
% test = boolean True/False
% okMsg = string of succesful test
% errMsg = string of the error result

if test==false
    warning('vsv:multi:configCheck',errMsg);
    fprintf(fid, [errMsg newline]);
else
    disp(okMsg);
    fprintf(fid, [okMsg newline]);
end
end

% On function close close any open files
function cleanupFunc(fid)
fclose(fid);
end
