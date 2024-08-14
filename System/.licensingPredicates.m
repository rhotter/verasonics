% The following is a 'dummy' function provided only so that
% saving this file does not generate a warning from MATLAB
% mode in Emacs.
function licensingPredicates()
    assert(false, 'licensingPredicates: Not expecting to be called.')
end

%
% License code debugging is enabled by setting an environment variable,
% VS_LICENSING_DEBUG.
%
% This predicate returns true if license code debugging is enabled,
% and false if it is not.
%
function result = isLicensingDebugEnabled()
  result = not(isempty(getenv('VS_LICENSE_DEBUG')));
end

%
% Is macAddr a MAC address associated with this machine?
%
function result = lpMac(macAddr)
    persistent macMap;

    if isempty(macMap)
        macMap = containers.Map();
    end

    if isKey(macMap, macAddr)
        result = macMap(macAddr);
    else
        result = lpMacActual(macAddr);
        macMap(macAddr) = result;
    end

    if isLicensingDebugEnabled
        fprintf('lpMac: result = %d\n', result);
    end
end

function result = lpMacActual(macAddr)
    assert(length(macAddr) == 12, 'lpMac: Wrong length for a MAC address')
    result = false;
    interfaces = java.net.NetworkInterface.getNetworkInterfaces;
    while (interfaces.hasMoreElements)
        systemMacAddr = interfaces.nextElement.getHardwareAddress;
        systemMacAddr = sprintf('%.2X', typecast(systemMacAddr, 'uint8'));
        if strcmpi(macAddr, systemMacAddr)
            result = true;
        end
    end
end

%
% Is chassisId the chassis ID of the attached Vantage system?
%
% Note: There is no means to actually check this at this time,
%       so this predicate always succeeds.
%
function result = lpChassisId(chassisId)
    result = true;
end

%
% Does the host controller have the specified hwConnected status?
%
function result = lpHwConnected(hwConnected)
    persistent hwConnectedMap;

    if isempty(hwConnectedMap)
        hwConnectedMap = containers.Map();
    end

    if isKey(hwConnectedMap, hwConnected)
        result = hwConnectedMap(hwConnected);
    else
        result = lpHwConnectedActual(hwConnected);
        hwConnectedMap(hwConnected) = result;
    end

    if isLicensingDebugEnabled
        fprintf('lpHwConnected: result = %d\n', result);
    end
end

function result = lpHwConnectedActual(hwConnected)
    [a, ~] = getHardwareInfo('SHI', 'AllAttributes');
    result = false;
    if (strcmp('Success', a) && strcmp('true', hwConnected))
        result = true;
    end
    if (~strcmp('Success', a) && strcmp('false', hwConnected))
        result = true;
    end
end

%
% Is hwId the hardware ID of the attached Vantage system?
%
function result = lpHwId(hwId)
    persistent hwIdMap;

    if isempty(hwIdMap)
        hwIdMap = containers.Map();
    end

    if isKey(hwIdMap, hwId)
        result = hwIdMap(hwId);
    else
        result = lpHwIdActual(hwId);
        hwIdMap(hwId) = result;
    end

    if isLicensingDebugEnabled
        fprintf('lpHwId: result = %d\n', result);
    end
end

function result = lpHwIdActual(hwId)
    persistent myHwId;
    if isempty(myHwId)
      [a, b] = getHardwareInfo('SHI', 'AllAttributes');
      if strcmp('Success', a)
          myHwId = deblank(b(7, 1:end));
      else
          myHwId = '';
      end
    end
    result = strcmp(myHwId, hwId);
end

%
% Is swVer the software version of the currently running software
% release?
%
% Note: This function assumes that the isCapabilityLicensed()
%       function in licenseCommon.c has set the value of the
%       variable vs_swVer__ in the MATLAB workspace.
%
function result = lpSwVer(swVer)
    result = strcmp(evalin('base', 'vs_swVer__'), swVer);

    if isLicensingDebugEnabled
        fprintf('lpSwVer: result = %d\n', result);
    end
end

%
% Is the current date within the specified date range?
%
% Note: The endDate is conceptually extended by a 30-day
% grace period.
%
function result = lpDateIn(startDate, endDate)
    persistent warningIssued;

    [result, daysGraceRemaining] = lpDateInAux(startDate, endDate);

    if isLicensingDebugEnabled
        fprintf('lpDateIn: result = %d\n', result);
    end

    if (isempty(warningIssued))
        if (0 < daysGraceRemaining && daysGraceRemaining <= 30)
            warningIssued = true;

            warnBacktrace = warning('query', 'backtrace');
            warning('off', 'backtrace');
            warning(['The license associated with this computer is ' ...
                     'within the defined grace period.']);
            warning(['There are %d calendar days remaining before license ' ...
                     'expiration.'], daysGraceRemaining);
            warning(['Please contact licensing-request@verasonics.com ' ...
                     'for assistance in renewing software licenses.']);
            warning(warnBacktrace.state, 'backtrace');
        end
    end

    if isLicensingDebugEnabled
        fprintf('lpSwVer: result = %d\n', result);
    end
end

%
% How many days remain in the current grace period?
%
function daysGraceRemaining = lpDaysGraceRemaining(startDate, endDate)
    [~, daysGraceRemaining] = lpDateInAux(startDate, endDate);
end

%
% Helper function for lpDateIn. Computes whether the current date is
% within startDate, endDate, or, if not, if within the 30 day grace
% period window that follows endDate.
%
function [result, daysGraceRemaining] = lpDateInAux(startDate, endDate)

    % Get the current date
    c = clock();
    cYear = c(1);
    cMonth = c(2);
    cDay = c(3);

    result = true;
    daysGraceRemaining = 0;

    if isequal(startDate, '')
        sYear = cYear;
        sMonth = cMonth;
        sDay = cDay;
    else
        sYear = str2double(startDate(7:10));
        sMonth = str2double(startDate(1:2));
        sDay = str2double(startDate(4:5));
    end

    if isequal(endDate, '')
        eYear = 9999;
        eMonth = 12;
        eDay = 31;
    else
        eYear = str2double(endDate(7:10));
        eMonth = str2double(endDate(1:2));
        eDay = str2double(endDate(4:5));
    end

    cDatenum = datenum(sprintf('%02d/%02d/%04d', cMonth, cDay, cYear), 'mm/dd/yyyy');
    eDatenum = datenum(sprintf('%02d/%02d/%04d', eMonth, eDay, eYear), 'mm/dd/yyyy');
    sDatenum = datenum(sprintf('%02d/%02d/%04d', sMonth, sDay, sYear), 'mm/dd/yyyy');

    if (cDatenum < sDatenum || cDatenum > eDatenum)
        result = false;
    end

    if (not(result))
        if (cDatenum > eDatenum)
            daysPastEndDate = round(cDatenum - eDatenum);
            gracePeriodInDays = 30;
            daysGraceRemaining = gracePeriodInDays - daysPastEndDate;
            if (daysGraceRemaining < 0)
                daysGraceRemaining = 0;
            else
                result = true;
            end
        end
    end
end

% Return a string that describes the currently-licensed capabilities.
function result = showLicensed()
    result = ['On this system, at this time, the status of capabilities ' ...
              'that may be licensed is as follows: \n\n'];

    fmt = '%s 1. The hardware capability %s licensed.\n';
    if licenseMgr('isCapabilityLicensed', 'hardware')
        result = sprintf(fmt, result, 'is');
    else
        result = sprintf(fmt, result, 'is not');
    end

    fmt = '%s 2. The simulation capability %s licensed.\n';
    if licenseMgr('isCapabilityLicensed', 'simulation')
        result = sprintf(fmt, result, 'is');
    else
        result = sprintf(fmt, result, 'is not');
    end

    fmt = '%s 3. The arbwaveTk capability %s licensed.\n';
    if licenseMgr('isCapabilityLicensed', 'arbwaveTk')
        result = sprintf(fmt, result, 'is');
    else
        result = sprintf(fmt, result, 'is not');
    end

    fmt = '%s 4. The arbwaveTx capability %s licensed.\n';
    if licenseMgr('isCapabilityLicensed', 'arbwaveTx')
        result = sprintf(fmt, result, 'is');
    else
        result = sprintf(fmt, result, 'is not');
    end

    fmt = '%s 5. The extendedIO capability %s licensed.\n';
    if licenseMgr('isCapabilityLicensed', 'extendedIO')
        result = sprintf(fmt, result, 'is');
    else
        result = sprintf(fmt, result, 'is not');
    end

    fmt = '%s 6. The extendedTx capability %s licensed.\n';
    if licenseMgr('isCapabilityLicensed', 'extendedTx')
        result = sprintf(fmt, result, 'is');
    else
        result = sprintf(fmt, result, 'is not');
    end

    fmt = '%s 7. The recon capability %s licensed.\n';
    if licenseMgr('isCapabilityLicensed', 'recon')
        result = sprintf(fmt, result, 'is');
    else
        result = sprintf(fmt, result, 'is not');
    end

    fmt = '%s 8. The trigger capability %s licensed.\n';
    if licenseMgr('isCapabilityLicensed', 'trigger')
        result = sprintf(fmt, result, 'is');
    else
        result = sprintf(fmt, result, 'is not');
    end

    fmt = '%s 9. The FUS 2D capability %s licensed.\n';
    if licenseMgr('isCapabilityLicensed', 'hifuPlex')
        result = sprintf(fmt, result, 'is');
    else
        result = sprintf(fmt, result, 'is not');
    end

    fmt = '%s10. The FUS Elite 3000 capability %s licensed.\n';
    if licenseMgr('isCapabilityLicensed', 'hifuPlexPlus')
        result = sprintf(fmt, result, 'is');
    else
        result = sprintf(fmt, result, 'is not');
    end

    fmt = '%s11. The thermalStrainImaging capability %s licensed.\n';
    if licenseMgr('isCapabilityLicensed', 'thermalStrainImaging')
        result = sprintf(fmt, result, 'is');
    else
        result = sprintf(fmt, result, 'is not');
    end

    fmt = '%s12. The zMeas capability %s licensed.\n';
    if licenseMgr('isCapabilityLicensed', 'zMeas')
        result = sprintf(fmt, result, 'is');
    else
        result = sprintf(fmt, result, 'is not');
    end

    fmt = '%s13. The multisystemRDMA capability %s licensed.\n';
    if licenseMgr('isCapabilityLicensed', 'multisystemRDMA')
        result = sprintf(fmt, result, 'is');
    else
        result = sprintf(fmt, result, 'is not');
    end

    fmt = '%s14. The ultrasoundCurriculum capability %s licensed.\n';
    if licenseMgr('isCapabilityLicensed', 'ultrasoundCurriculum')
        result = sprintf(fmt, result, 'is');
    else
        result = sprintf(fmt, result, 'is not');
    end

    fmt = '%s15. The ndeResearch capability %s licensed.\n';
    if licenseMgr('isCapabilityLicensed', 'sonivueNDE')
        result = sprintf(fmt, result, 'is');
    else
        result = sprintf(fmt, result, 'is not');
    end

    fmt = '%s16. The gpuToolkit capability %s licensed.\n';
    if licenseMgr('isCapabilityLicensed', 'gpuToolkit')
        result = sprintf(fmt, result, 'is');
    else
        result = sprintf(fmt, result, 'is not');
    end

    fmt = '%s17. The gpuDirect capability %s licensed.\n';
    if licenseMgr('isCapabilityLicensed', 'gpuDirect')
        result = sprintf(fmt, result, 'is');
    else
        result = sprintf(fmt, result, 'is not');
    end
end

%
% Create a 'systemInfo.txt' file that the customer can
% transmit to Verasonics for help with support issues.
%
function [] = systemInfo()

    [vs_ret__, vs_hostName__] = system('hostname');
    if vs_ret__ ~= 0
       if ispc
          vs_hostName__ = getenv('COMPUTERNAME');
       else
          vs_hostName__ = getenv('HOSTNAME');
       end
    end
    vs_hostName__ = strtrim(lower(vs_hostName__));

    vs_fid__ = fopen('systemInfo.txt', 'w');
    if vs_fid__ == -1
      error('Unable to open text file \"systemInfo.txt\" for writing.');
    end

    function printSeparator
        fprintf(vs_fid__, '-----------------------------------------------------------------------------------------------------\n');
    end

    fprintf(vs_fid__, 'Verasonics Product Software System Information Report\n');
    fprintf(vs_fid__, '=====================================================================================================\n');
    fprintf(vs_fid__, '\n');
    fprintf(vs_fid__, 'This text file contains information describing the computer\n');
    fprintf(vs_fid__, 'system on which it was generated. This information may be used\n');
    fprintf(vs_fid__, 'for providing technical support and/or for licensing the use of\n');
    fprintf(vs_fid__, 'Verasonics, Inc., product software.\n');
    fprintf(vs_fid__, '\n');
    fprintf(vs_fid__, 'This file is created whenever the systemInfo command is executed\n');
    fprintf(vs_fid__, 'from within the MATLAB console window.\n');
    fprintf(vs_fid__, '\n');
    fprintf(vs_fid__, 'If the Verasonics product software indicates that it is unlicensed,\n');
    fprintf(vs_fid__, 'or that it is executing with a defined grace period, you are asked\n');
    fprintf(vs_fid__, 'to execute the systemInfo command from within a MATLAB console window\n');
    fprintf(vs_fid__, 'and to then send the contents of this text file (systemInfo.txt) to\n');
    fprintf(vs_fid__, 'the Verasonics, Inc., technical support group by electronic mail,\n');
    fprintf(vs_fid__, 'as described below.\n');
    fprintf(vs_fid__, '\n');

    fprintf(vs_fid__, 'Contacting Verasonics by Electronic Mail\n');
    printSeparator
    fprintf(vs_fid__, '\n');
    fprintf(vs_fid__, 'Please attach this file to an electronic mail message. Requests for technical\n');
    fprintf(vs_fid__, 'support unrelated to licensing should be addressed to support@verasonics.com.\n');
    fprintf(vs_fid__, 'Requests for assistance with product licensing issues should be addressed to\n');
    fprintf(vs_fid__, 'licensing-request@verasonics.com.\n');
    fprintf(vs_fid__, '\n');
    fprintf(vs_fid__, 'Please include your full name, the name of your institution or organization,\n');
    fprintf(vs_fid__, 'your telephone number(s), your preferred email address, and any additional\n');
    fprintf(vs_fid__, 'contact information you wish to include in the body of your email message.\n');
    fprintf(vs_fid__, '\n');
    fprintf(vs_fid__, 'You will be contacted promptly thereafter for assistance with your technical\n');
    fprintf(vs_fid__, 'support or product software licensing request.\n');
    fprintf(vs_fid__, '\n');

    printSeparator

    fprintf(vs_fid__, 'Local Time: %s\n', datestr(now));
    fprintf(vs_fid__, 'Hardware ID: ');

    [vs_status__, vs_hwSerNo__] = getHardwareInfo('SHI', 'AllAttributes');
    if strcmp('Success', vs_status__)
        vs_hwSerNo__ = deblank(vs_hwSerNo__(7, 1:end));
    else
        vs_hwSerNo__ = '';
    end

    if isempty(vs_hwSerNo__)
      fprintf(vs_fid__, 'No attached hardware detected.\n');
    else
      fprintf(vs_fid__, '%s\n', vs_hwSerNo__);
    end

    vs_interfaces__ = java.net.NetworkInterface.getNetworkInterfaces;
    vs_systemMACAddr__ = {};
    while (vs_interfaces__.hasMoreElements)
      macAddr = vs_interfaces__.nextElement.getHardwareAddress;
      macAddr = sprintf('%.2X', typecast(macAddr, 'uint8'));
      if numel(macAddr) == 12
        vs_systemMACAddr__{end+1} = macAddr;
      end
    end

    vs_systemMACAddr__ = unique(vs_systemMACAddr__);
    for i=1:length(vs_systemMACAddr__)
        fprintf(vs_fid__, 'MAC Address: %s\n', vs_systemMACAddr__{i});
    end

    fprintf(vs_fid__, 'Hostname: %s\n', vs_hostName__);
    fprintf(vs_fid__, 'SW Version: %s\n', evalin('base', 'vs_fullSwVer__'));
    vs_loginName__ = char(java.lang.System.getProperty('user.name'));
    fprintf(vs_fid__, 'Login Name: %s\n', vs_loginName__);
    fprintf(vs_fid__, 'LP Memoization: Enabled\n');
    vs_osVersion__ = split(evalc('ver'), newline);
    fprintf(vs_fid__, '%s\n', vs_osVersion__{[2 4:end-1]});

    printSeparator
    fprintf(vs_fid__, '\n');

    fprintf(vs_fid__, licenseMgr('explain'));

    printSeparator

    fprintf(vs_fid__, licenseMgr('showLicensed'));

    printSeparator

    vs_map__ = java.lang.System.getenv();
    vs_map__ = containers.Map(upper(cell(vs_map__.keySet.toArray)), cell(vs_map__.values.toArray));
    vs_keys__ = vs_map__.keys;

    for vs_i__=1:numel(vs_keys__)
        fprintf(vs_fid__, '%s=%s\n', vs_keys__{vs_i__}, vs_map__(vs_keys__{vs_i__}));
    end

    printSeparator

    vs_startups__ = which('startup', '-all');
    fprintf(vs_fid__, 'MATLAB startup files:\n');
    fprintf(vs_fid__, ' %s\n', vs_startups__{:});
    fprintf(vs_fid__, 'MATLAB path:\n');
    vs_paths = strsplit(path(), pathsep);
    fprintf(vs_fid__, ' %s\n', vs_paths{:});

    fclose(vs_fid__);
end
