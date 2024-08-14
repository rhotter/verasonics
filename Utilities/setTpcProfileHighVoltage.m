function [result, hvSet] = setTpcProfileHighVoltage(hvReq, tpcNum, isVsxInit)
% setTpcProfileHighVoltage  Change transmit voltage
%   Checks that the requested voltage value is within the allowed range.
%   The allowed range is determined by the Vantage min TPC voltage and
%   the high voltage limit values specified in the TPC structure. The
%   Control structure is updated to command runAcq to set the TPC voltage.
%   The external power supply is also commanded to change voltage when
%   in use.
%
%   Argument definitions:
%       hvReq:  Matlab double, desired transmit voltage in Volts
%       tpcNum: TPC profile being set; Matlab double with integer value
%               in range 1:5
%       hvSet:  Set equal to hvReq for result = 'Success'; otherwise set
%               to empty
%       result: Matlab string set to 'Success' or an appropriate error message
%       isVsxInit: logical, set true during VSX initialization and false 
%               during run-time.  This argument is optional, and defaults
%               to false if not provided.

% Revision history
% Aug 7 2020 VTS-1796 optional input argument isVsxInit

% VTS-1796 set default if isVsxInit not provided
if nargin < 3
    isVsxInit = false; % default is to assume run time
end

minTpcVoltage = evalin('base', 'minTpcVoltage');
try
    TPC = evalin('base', 'TPC');
catch
    error('setTpcProfileHighVoltage: TPC structure required but not found.');
end
tpcNum = round(tpcNum);
if tpcNum < 1 || tpcNum > 5 || TPC(tpcNum).inUse == 0
    result = 'Out of range TPC profile index';
    hvSet = [];
    return
end
% restrict hvReq to stay within TPC limits (VTS-2509)
TPChvMax = min(TPC(tpcNum).maxHighVoltage, TPC(tpcNum).highVoltageLimit);
hvReq = min(hvReq, TPChvMax);
if hvReq < minTpcVoltage
    result = 'Out of range hv Value';
    hvSet = TPC(tpcNum).hv;
    return
end

VDAS = evalin('base', 'VDAS');
if VDAS
    extps = TPC(tpcNum).inUse == 2;
    if tpcNum == 5 && extps
        % Call external power supply control function.
        PSerrorflag = extPwrCtrl('SETV', hvReq);
        if PSerrorflag
            result = 'ERROR! Communication error with external power supply.';
            hvSet = TPC(tpcNum).hv;
            return
        end
    end
end

% VTS-1796 Check VSX initialization status.  If VSX is in initialization
% don't create set&Run command for TPC hv; Sequence Load will automatically
% pick up the value that has been written to TPC(tpcNum).hv
if  ~isVsxInit
    % get Control structure from base workspace, or set to empty if not
    % found and then create the structure
    Control = vsv.seq.getBaseComp('Control');
    
    if isempty(Control)
        Control = struct( 'Command', [], 'Parameters', [] );
    end
    
    if isempty(Control(1).Command)
        n = 1;
    else
        n = length(Control)+1;
    end
    Control(n).Command = 'set&Run';
    Control(n).Parameters = {'TPC', tpcNum, 'hv', hvReq};
    assignin('base', 'Control', Control);
end

TPC(tpcNum).hv = hvReq;
assignin('base', 'TPC', TPC);
result = 'Success';
hvSet = hvReq;

