function TpcProfileMissed(profileNumber)
%TpcProfileMissed - Issues warning when a TPC level transition is missed.
%   If a transmit occurs while the TPC is ramping from one voltage level to
%   another, the TPC sets an error flag that is detected by runAcq and
%   runAcq then calls this function to warn the user of the condition.  For errors
%   due to changes in level caused by HV slider movement, we want to supress this
%   warning, and so we check the time elapsed from the time in the variable,
%   'tStartHvSldr', which was set when the slider was moved.  If the time is
%   less than 0.5 seconds, no warning is issued.
%
% Copyright 2001-2017 Verasonics, Inc.  All world-wide rights and remedies under all intellectual property laws and industrial property laws are reserved.  Verasonics Registered U.S. Patent and Trademark Office.

tElapsed = 1.0;  % default elapsed time to something greater than 0.5
if evalin('base','exist(''tStartHvSldr'',''var'');')
    tStart = evalin('base','tStartHvSldr');
    tElapsed = toc(tStart);
end


% Return if elapsed time less than 0.5
if tElapsed < 0.5, return, end

% Output warning
fprintf('Warning: A transmit occurred while a transition to high voltage profile %d was in progress.\n',profileNumber);


end
