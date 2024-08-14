function rc = setProfile5VoltageLimits(P5min, P5max)
% setProfile5VoltageLimits.m  Utility function for user control of TPC
% Profile 5 HV monitor limits
%
%   Created Jan 12, 2020 for 4.2.0 SW release (VTS-1583)
%   This function allows the user to enable and set the thresholds for TPC
%   hardware fault monitoring of the the TPC profile 5 transmit voltage.
%   This monitoring feature is automatically disabled when VSX loads and
%   begins running a user script; the user can enable monitoring whenever
%   desired by calling this function with P5min and P5max set to the
%   desired fault thresholds in Volts, below and above the nominal transmit
%   HV level that has been set.  (Note these limits must be updated whenever
%   the HV level is changed.)  Monitoring can also be disabled whenever
%   desired, by simply setting P5min equal to or greater than P5max.  The
%   return argument "rc" will be a string variable set to 'Success' if the
%   command was completed successfully, or to an error message string
%   reporting the problem if there was an error.

rc = 'Success'; % initialize return status argument
% fetch TPC structure and update fields
TPC = evalin('base', 'TPC');
TPC(5).HVmonitorLimits = [P5min, P5max]; % updated TPC structure to new values
% clip user-supplied values to <= 100 Volts
TPC(5).HVmonitorLimits = min(TPC(5).HVmonitorLimits, 100);
% clip user-supplied values to => 0 Volts
TPC(5).HVmonitorLimits = max(TPC(5).HVmonitorLimits, 0);

% return updated TPC structure to base workspace
assignin('base', 'TPC', TPC);
% set action flag for VSX to process at next return to matlab
evalin('base', 'action = ''setTPC5hvMonitorLimits'';');

end

