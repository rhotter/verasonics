function [] = TXEventCheck(isVsxInit)
% Copyright 2001-2021 Verasonics, Inc. Verasonics Registered U.S. Patent and Trademark Office.
%
% Analyze TX events from a setup script .mat file after VSX initialization,
% and after any 'update&Run' command that affects the HW event sequence.
% Check for gate driver supply operating limits, Transmit device
% temperature rise limits, and per-channel output current limits.  Estimate
% total system power dissipation due to HIFU transmit activity (but this
% does not generate a limit, since temperature rise is monitored actively
% during operation).  Also check for push capacitor peak current limit, and
% TXbus distribution maximum current limit per acquisition board.
%
% The input argument "isVsxInit" is a logical, set by VSX to true during
% initialization and false while VSX is in its run-time while loop.  This
% argument can be used within TXEventCheck for situations where the
% behavior needs to be different between initialization and runtime.  For
% backward compatibility, it is optional: if it is not present it will be
% created and set to false.

% TXEventCheck finds all active transmit events in the event sequence as
% defined by the user's setup script, estimates the operating state
% currents, dissipation levels, etc. and compares them to the limits that
% have been defined to ensure safe, reliable operation of the system.
% TXEventCheck is automatically invoked by VSX while in the process of
% loading a user setup script and also while executing an 'update&Run'
% command that could modify the transmit operating state. If any parameter
% is found to exceed the limit, TXEventCheck reports that specific failure
% and automatically reduces the HV maximum voltage to stay within the
% limits.  For those cases where reducing HV will not resolve the problem,
% TXEventCheck will cause VSX to exit with an error message, with the
% intent that the user would then modify the setup script to correct the
% problem and then run VSX again, iterating in this fashion until no
% overlimit conditions are found.  The limit evaluations are applied to all
% transmit events for peak output current, transformer FLux saturation
% limit, and other constraints that applly to all transmit activity.
% Transmit events using TPC profile 5 (or an imaging profile with a small
% number of channels) at high total burst energy levels (referred to as
% "HPEB" events, for High Power Extended Burst) will have additional
% thermal limits evaluated due to the high power disspiation levels that
% are possible with the HIFU or Extended Transmit options.  For the
% dissipation limits non-HPEB events are simply ignored, but their duration
% is taken into account when calculating the effective PRI for the HPEB
% events.

% If any overlimit condition is detected that can be corrected by reducing
% the maximum transmit voltage level, TXEVentCHeck will automatically
% reduce the associated TPC Profile max voltage limit, and the current HV
% setting if necessary, to enforce the new reduced maximum.  In this case,
% the script will be allowed to continue running but at the reduced voltage
% level.

% TXEventCheck produces an array TXEvent(:) with an entry for each HPEB
% transmit event, and saves this array to the matlab workspace.  Each entry
% is a structure containing all of the estimated operating state parameters
% for that event, to facilitate analysis and debugging activities.  A
% separate structure array TXthermal(:) is created, with the thermal
% characteristics and estimated temperature rise for individual transmit
% devices.

% TXEventCheck synthesizes a "TXEvent" sequence that will be accurate for
% all combinations of loops, calls, etc. from the user's Event sequence.
% From this sequence plus the characteristics of each trasmit event, total
% duty cycles, power levels, etc. can be estimated. The approach is to
% simply follow the sequence control commands as written, for a single pass
% through the event sequence and thus create a new TXEvent every time a
% HPEB event is encountered.

% As TXEventCheck steps through the event sequence, it accumulates noops
% and non-HPEB events to determine the actual PRI of the previous HPEB
% Event. Therefore the PRI is not actually known until the next TXEvent has
% been reached.  To fill in the PRI for the last TXEvent in a repeating
% sequence, the jump back to start repeating must be followed until the
% next HPEB event is reached.  At that jump, a "finalHPEB" flag is set so
% the routine will know it is time to quit when the first HPEB event is
% reached after the jump.

% IDENTIFYING THE END OF THE SEQUENCE:  A key assumption in this routine is
% that there is only one way a sequence can repeat- by a jump command back
% to a previous event, and never while in a loop or subroutine call.  In
% addition, TXEventCheck currently assumes the ONLY case in which there
% would be a jump to an earlier event is the start of the repeating
% sequence.  We may need to add additional logic to handle cases where this
% assumption is false (e.g. a jump to an earlier event that has never been
% executed, since it is earlier than the startEvent point).

% The PRI for each transmit event is determined by evaluating the duration
% of all events plus any noop, timeToNextAcq, or timeToNextEB instructions.
% The estimated PRI will not include other sources of delay such as waiting
% for a trigger input, or pausing to wait for DMA's to complete, or pausing
% to wait for processing in a synchronous script.  The effect of any such
% additional delay will be to reduce the power and dissipation levels from
% those estimated by TXEventCheck and thus the limit calculations will
% reflect a 'worst case' operating state.

% TXEventCheck identifies TPC profile transition commands in the event
% sequence, so it can apply the correct profile to each transmit event. But
% it does not evaluate the transition time allowed between the set TPC
% profile command and the following TX event, to determine if sufficient
% transition time has been provided for the TPC to complete the transition.
% Any missed TPC profile transition time while the system is running will
% result in an error or warning message to the matlab command line,
% however.

% The limits calculated by TXEventCheck are highly dependent on the actual
% transducer load impedance presented to the system, as set by the user in
% the parameter Trans.impedance.  If this value is not correct, the limits
% evaluated here will be meaningless and may allow the system to
% operate in a self-destructive state.  If the transducer element impedance
% varies significantly from element to element, it may be necessary to
% modify TXEventCheck to specify impedance on a per-element basis.  Note
% also that the Trans.impedance parameter is intended to be the complex
% impedance presented to the system as a function of frequency.  The limit
% algorithms assume the transducer tuning network (which is usually a
% series inductor) will result in a load impedance much higher at the
% harmonics of the transmit frequency than at the fundamental, and thus
% power consumption at the harmonics is assumed to be insignificant.

% Note that all of the limit parameters and associated values as defined
% within TXEventCheck are based on the capabilites of the VDAS system by
% itself, up to the transducer connector.  In many cases, the transmit
% power and/or transmit current that can be safely handled by the
% transducer itself will be far more restrictive than the system's transmit
% output capability. TXEventCheck is only intended to protect the system
% from self-destructive operating states; it is up to the user to impose
% additional limits as needed to protect the transducer and/or to stay
% within regulatory limits for acoustic output, transducer temperature,
% etc.

% Revision History
% Aug 7 2020 VTS-1796 optional input argument isVsxInit
% Nov 4 2019 VTS-1419 add checks for SeqControl structure not defined, or
%   optional SeqControl fields not defined
% Sep 1 2019 VTS-1416 add support for Active Clamp configuration
% Updated April 13, 2018 for 3.4.2 issue VTS-476, to check for and enforce
%   the optional transducer maximum average power limit given by
%   Trans.maxAvgPower (output power going to the transducer, averaged over
%   entire event sequence).

%% get copies of required variables from base workspace

% VTS-1796 check for optionial input argument
if nargin == 0
    % default to runtime status if isVsxInit not provided
    isVsxInit = false;
end

% VTS-1419: SeqControl may not exist
if evalin('base', 'exist(''SeqControl'', ''var'')')
    SeqControl = evalin('base', 'SeqControl'); % get SeqControl structure for access below
    % create dummies for optional fields if they don't exist
    if ~isfield(SeqControl, 'condition')
        SeqControl(1).condition = [];
    end
    if ~isfield(SeqControl, 'argument')
        SeqControl(1).argument = [];
    end
else
    % SeqControl not found; create a dummy in case event sequence tries to
    % use it
    SeqControl.command = [];
    SeqControl.condition = [];
    SeqControl.argument = [];
end
TPC = evalin('base', 'TPC'); % get the TPC structure from base workspace
% TPC(5).inUse = 0:  Access to TPC Profile 5 is disabled
% TPC(5).inUse = 1:  Extended Transmit Option is installed and in use
% TPC(5).inUse = 2:  HIFU configuration is installed and in use
Trans = evalin('base', 'Trans'); % and Trans

% import sequence structures that will be referenced:
TX = evalin('base', 'TX');
TW = evalin('base', 'TW');
if evalin('base', 'exist(''Receive'', ''var'')')
    Receive = evalin('base', 'Receive');
end
Event = evalin('base', 'Event');

numBoards = evalin('base', 'numBoards'); % number of acquisition boards present in HW system

Resource = evalin('base', 'Resource'); % get the Resource structure from base workspace


%% initialize variables and constants for this run of TXEventCHeck

% enable or disable text warnings to command window based on
% Resource.HIFU.verbose parameter value.
if (isfield(Resource,'HIFU')&&isfield(Resource.HIFU,'verbose')&&~isempty(Resource.HIFU.verbose))
    verbose = Resource.HIFU.verbose;
else
     % default to state being used by VSX if not specified by user
    verbose = evalin('base', 'Resource.Parameters.verbose');
end

% Find which profiles are in use.
profile = 1; % default to profile 1 if there are no profile change
% for j=1:size(TPC,2)
%     TPC(j).inUse = 0;
% end
for j=1:size(SeqControl,2)
    if strcmp(SeqControl(j).command,'setTPCProfile')
        % At least one setTPCProfile command exists, so we will start with
        % the profile set to zero to check for a set profile prior to first
        % transmit event.
        profile = 0;
        % Never demote a profile that was previously in use, but add any
        % now being used that were not previously
        TPC(SeqControl(j).argument).inUse = max(TPC(SeqControl(j).argument).inUse, 1);
    end
end
% Profile 1 is active by default, if none other active.
if ~any([TPC.inUse])
    TPC(1).inUse = 1;
end
if profile == 0 && TPC(5).inUse == 0
    % multiple profiles are being used, but profile 5 is not.  So in this
    % case set an initial default profile of 1, so the requirement to
    % initialize the profile selection before the first transmit event will
    % not be enforced.  Failure to initialize a TPC profile selection will
    % only trigger an error if profile 5 is being used in the script.
    profile = 1;
elseif TPC(5).inUse
%     TPC(5).inUse = Resource.SysConfig.p5req;
end

% reset highVoltageLimit in each TPC structure, since this run of
% TXEventCheck may result in a higher limit than was set previously
% Create record of previous highVoltageLimit values, so we will be able to
% determine if any have changed
PreviousHvLimit = zeros(1, 5);
for i = 1:5
    if TPC(i).inUse
        % note unused profiles may have empty values and we will leave them
        % that way.
        if isempty(TPC(i).maxHighVoltage)
            % user did not initialize, so do it for them
            TPC(i).maxHighVoltage = Trans.maxHighVoltage;
            TPC(i).highVoltageLimit = TPC(i).maxHighVoltage;
        end
        PreviousHvLimit(i) = TPC(i).highVoltageLimit; % capture existing value before resetting it
        TPC(i).highVoltageLimit = TPC(i).maxHighVoltage;
    end
end


% check for the "voltageTrackP5" variable and set default if needed.
if TPC(5).inUse && (isfield(Resource,'HIFU') && isfield(Resource.HIFU,'voltageTrackP5') ...
        && ~isempty(Resource.HIFU.voltageTrackP5))
    voltageTrackP5 = Resource.HIFU.voltageTrackP5;
else
     % default to disabled state if not specified by user, or if P5 not
     % being used
    voltageTrackP5 = 0;
end

if TPC(5).inUse && (~isfield(TPC, 'currentLimit') || isempty(TPC(5).currentLimit))
    TPC(5).currentLimit = 1; % initialize if necessary, to 1 Amp level
end

% get the TXindex value from the Resource structure
TXindex = Resource.SysConfig.TXindex;
% A check for unrecognized TXindex values will be made later, when the
% limit values that depend on configuration are defined.

P5HVperBdindex = Resource.SysConfig.P5HVperBdindex;
% indicates whether backplane per-board limits on Profile 5 HV current are
% programmable, or fixed at 40 Amps

startEvent = Resource.Parameters.startEvent; % first event in the sequence (usually set to 1)

numEvents = size(Event,2); % total number of events in the script

% characteristics of the system used in evaluating limits:
clampV = 0.8; % Voltage drop across the transmit clamp at receive path input, when transmit current is flowing

PRIovh = 5; % typical overhead time needed for system to finish one tx-rx event and start another, in usec

P5rqdHVmin = 5; % threshold level in Volts
% If the required maxHV setting to allow use of an imaging
% profile is below this level, the sysExtenBL flag will be
% set.

HPEBprofiles = zeros(1, 5); % to track which profiles are driving HPEB events

%% Read through Event List, decode SeqControl instructions, find Profiles & find PRI's


% --- initialize some added parameters in each event
for eventnum = 1:numEvents
    Event(eventnum).repeat = 0;    % zero indicates event has not been used
    Event(eventnum).HPEB = 0;    % flag for HIFU burst (profile 5 with extended burst length)
    Event(eventnum).callID = []; % empty means no call command in this event
end

% % clear sysExtendBL in all TX structures
% for txnum = 1:length(TX)
%     TX(txnum).sysExtendBL = 0;
% end

nextprofiledelayed = profile; % used for 'next' condition in non-tr event
loopcountset = -ones(8,1); % minus one state signifies it has not been initialized by a loopCnt command (required for looping)
activeCall = 0; % indicates we have not made a subroutine call
numcalls = 0; % will increment for each call command found
LoopCtrLst = zeros(1,8);  % identifies which loop counters are actually being used

%% While-loop to interpret/ decode Event Sequence:

% The while loop below will step through the event sequence in the same way
% the HW system will, to derive the values (such as effective PRI for each
% HPEB event) needed to evaluate transmit limits.  The following steps are
% performed within this while loop:
    % 1: Decode sequence control instructions, and assign active TPC profile
    % number to each Event
    % 2: Assign total PRI values for each HPEB event based on non-HPEB
    % events, noops found in other events, etc.
    % 3: Create a TXEvent structure for every HPEB event encountered in the
    % event sequence

eventnum = startEvent; % first event to be evaluated.

TEnum = 0;
finalHPEB = 0; % flag to determine when we have returned for a second pass
activePRI = 0;
activettna = 0;
fluxLimStr1 = []; % warning string for flux limit HV reduction

while eventnum > 0   % eventnum will be set to zero when we reach the end of the sequence,
    % or the point where it jumps back and starts repeating

    nextprofile = profile;      % default to existing profile if there are no profile changes
    TXnum = Event(eventnum).tx; % pointer to the TX structure number for this event (or zero if no TX in the event)
    if TXnum > 0
        TWnum = TX(TXnum).waveform; % pointer to the TW structure number for this event's TX
    else
        TWnum = 0; % or set to zero if there is no TX in this event
    end

    % --- record profile for this event
    Event(eventnum).profile = profile;

    % --- detect T or R in current event
    if TXnum ~= 0 || Event(eventnum).rcv ~= 0
        tr = 1;
        nextprofile = nextprofiledelayed;   % since this is a tr event, any delayed next profile will take effect
                                            % at the end of this event
    else
        tr = 0; % no transmit or receive activity in this event
        PRI = 0; % default value for time taken in event
    end
    Event(eventnum).tr = tr; % set to one if event has TX or RCV

    % --- set repeat variable to record that this event has been
    % executed (was previously initialized to zero for all events)
    Event(eventnum).repeat = Event(eventnum).repeat + 1;

    % check flux limit and imaging burst total charge limit for the
    % TXstructure and TPC profile being used in each TX event
    if TXnum > 0 && TX(TXnum).VDASBistDriveEnable == 0
        % don't make these checks if BIST driver is being used
        if profile == 0 % First check if TPC profile has been initialized
            fprintf(2,['TXEventCheck Error: This script is using multiple TPC profiles, but a profile was not selected prior to the first transmit in Event ',...
                 num2str(eventnum), ',\n']);
            fprintf('Add a setTPCProfile command at the start of the event sequence, and then run the SetUp script and VSX again.\n\n');
            error('TXEventCheck: Exiting due to error condition identified above.')
        end
        if TW(TWnum).fluxHVlimit < TPC(profile).highVoltageLimit
            TPC(profile).highVoltageLimit = floor(TW(TWnum).fluxHVlimit);
            if verbose
                fluxLimStr1 = ['TXEventCheck Warning: reducing TPC(', num2str(profile), ').highVoltageLimit to ', ...
                    num2str(TPC(profile).highVoltageLimit), ' Volts due to\n'];
                fluxLimStr2 = ['transformer flux limit of TW(', num2str(TWnum), ') in Event(', num2str(eventnum), ').\n'];
                fluxLimStr3 = ' \n';
            end
        end
    end

    % --- set default (increment) value for next event to execute, if no jumps etc. encountered
    if eventnum < numEvents
        Event(eventnum).nextevent = eventnum + 1;
    else
        Event(eventnum).nextevent = 0; % zero indicates we have reached the end of the sequence
    end

    % --- set default values of timing command arguments
    % These three parameters will be updated with actual value if the
    % associated seq. control command is present in this event.
    noop = 0;
    ttna = 0;
    ttnEB = 0;

    % --- manage multiple sequence control commands
    % this seqcmd flag lets us detect if this event contains more than one event sequencing control command
    %  (loopCnt, loopTst, stop, jump, call, return), and declare an error if so
    seqcmd = 0;

    % --- 1: START Decode sequence control instructions ------------------------------ 1 -----------------------------
    for j=1:size(Event(eventnum).seqControl,2) % start decoding SeqControl commands in this event
        if Event(eventnum).seqControl(j)>0 % skip the decoding if SeqControl value is zero

            % We have a non-zero seqControl index so read in the field
            % values from the indexed SeqControl structure
            command = SeqControl(Event(eventnum).seqControl(j)).command;
            condition = SeqControl(Event(eventnum).seqControl(j)).condition;
            argument = SeqControl(Event(eventnum).seqControl(j)).argument;

            switch command
                case 'setTPCProfile'        % ----------------------
                    if strcmp(condition, 'immediate')
                        if tr == 1
                            fprintf(2, 'TXEventCheck: setTPCProfile immediate cannot be used in an event with active transmit or receive. \n') %#ok<*PRTCAL>
                            fprintf('correct the error and then run the SetUp script and VSX again.\n\n');
                            error('TXEventCheck: Exiting due to error condition identified above.')
                            return
                        end
                        % At this point, we have a valid setTPCProfile 'immediate' command, so put new profile # in nextprofile
                        nextprofile = argument;
                        nextprofiledelayed = argument;
                    else % either 'next' or nothing which defaults to next
                        nextprofiledelayed = argument;
                        if tr == 1
                            nextprofile = argument; % if tr = 1, takes effect at end of this event
                        end
                    end

                case 'timeToNextAcq'     % ----------------------
                    if tr == 1
                        ttna = argument;
                    else
                        fprintf(2, 'TXEventCheck: timeToNextAcq command cannot be used in event with no TX or Receive. \n')
                        fprintf(['correct the error in event ', num2str(eventnum), '  and then run the SetUp script and VSX again.\n\n']);
                        error('TXEventCheck: Exiting due to error condition identified above.')
                        return
                    end

                case 'timeToNextEB'      % ----------------------
                    if TXnum > 0 && profile == 5
                        ttnEB = argument;
                        TX(TXnum).sysExtendBL = 1; % Force to true regardless of waveform, so TTNEB will function properly
                    else
                        fprintf(2, 'TXEventCheck: timeToNextEB command can be used only in extended burst profile 5 TX events. \n')
                        fprintf(['correct the error in event ', num2str(eventnum), '  and then run the SetUp script and VSX again.\n\n']);
                        error('TXEventCheck: Exiting due to error condition identified above.')
                        return
                    end

                case 'loopCnt'           % ----------------------
                    seqcmd = seqcmd + 1; % increment the count of event sequencing commands

                    % --- Determine which counter is being used
                        if isempty(condition) % if present, condition value will specify loop counter ID number
                            CounterID = 1;  % default to counter # 1 if ID not specified
                            LoopCtrLst(CounterID) = 1; % note that this counter is being used
                        else
                             CounterID = str2num(condition(8));
                             LoopCtrLst(CounterID) = 1;
                            % condition value should be 'counterN' where N is an
                            % integer from 1:8.  So just convert 8th character in
                            % string to a number to get counter ID.
                        end
                    loopcountset(CounterID) = argument; % number of times we are to execute the loop: ?needs to equal the Event(n).repeat?
                    loopcounter(CounterID) = argument;  % if argument = 0, the associated loopTst will not jump back to repeat the loop


                case 'loopTst'           % ----------------------
                    seqcmd = seqcmd + 1; % increment the count of event sequencing commands

                    % --- set CounterID and do error check
                        if isempty(condition) % if present, condition value will specify loop counter ID number
                            CounterID = 1;  % default to counter # 1 if ID not specified
                        else
                             CounterID = str2num(condition(8));                 %#ok<*ST2NM>
                            % condition value should be 'counterN' where N is an
                            % integer from 1:8.  So just convert 8th character in
                            % string to a number to get counter ID.
                        end

                        if loopcountset(CounterID) == -1 % error if no loopcount value
                            fprintf(2, ['TXEventCheck error: Event ', num2str(eventnum), ' has loopTst command but loop count has not been set.\n\n']);
                            error('TXEventCheck: Exiting due to error condition identified above.')
                            return
                        end

                    % --- mimic the loop
                    if loopcounter(CounterID) > 0
                        % We are executing the loop, and need to keep track of repeats, (see Event(n).repeat)
                        % for correct number of passes through event sequence.
                        % the next event will be the top of the loop.
                        loopcounter(CounterID) = loopcounter(CounterID) -1; % decrement the counter
                        Event(eventnum).nextevent = argument; % jump to top of loop as next event

                        % now we need to test if this is an "endless" loop,
                        % where the loopCnt command is being used within
                        % the loop to keep it from ever ending.  Detect
                        % this state by checking for repeat count of the
                        % loopTst event greater than the loopCbnt initial
                        % value
                        if Event(eventnum).repeat > loopcountset(CounterID) +3
                            % add three to insure we've already gone around
                            % the loop enough times to have filled in
                            % complete PRI for any HPEB events within the
                            % loop.  Once we reach that point just quit, by
                            % setting nextevent to zero.
                            Event(eventnum).nextevent = 0;
                        end

                    elseif loopcounter(CounterID) == 0 % DONE with loop
                        % either we just finished filling in the loop, or the
                        % loop count is zero so no filling is needed.  In
                        % either case, we are done handling this loopTst command.
                        % now the sequencer will just go to event following the
                        % loopTst event, or quit if we are at the end of the
                        % event list
                        if eventnum < numEvents % go to next event if there is one
                            Event(eventnum).nextevent = eventnum + 1;
                        else
                            Event(eventnum).nextevent = 0; % zero indicates we are done with the sequence
                        end
                    end

                case 'call'             % ----------------------
                    seqcmd = seqcmd + 1; % increment the count of event sequencing commands
                    if isempty(Event(eventnum).callID) % if empty, we haven't processed this call command before
                        numcalls = numcalls + 1; % increment the running total of call commands
                        Event(eventnum).callID = numcalls; % callID value for this event
                        returnEvent(numcalls) = eventnum + 1; % the return destination from this call is the following event
                        if verbose > 2
                            disp( ['Call to Event # ' num2str(argument) ' ...'])
                            disp( ['... Return to Event # ' num2str(eventnum + 1)])
                        end
                    end
                    activeCall = Event(eventnum).callID;
                    Event(eventnum).nextevent = argument; % argument of the call command is the subroutine starting event

                case 'rtn'               % ----------------------
                    seqcmd = seqcmd + 1; % increment the count of event sequencing commands
                    if activeCall == 0 % return without a previous call is an error
                        fprintf(2, ['TXEventCheck error: Event ', num2str(eventnum), ' has a return command but there was no previous call.\n\n']);
                        error('TXEventCheck: Exiting due to error condition identified above.')
                        return
                    end
                    Event(eventnum).nextevent  = returnEvent(activeCall); % the return tells us what the next event will be
                    activeCall = 0; % set back to zero since we are now returning from the call

                case 'stop'              % ----------------------
                    seqcmd = seqcmd + 1; % increment the count of event sequencing commands
                    Event(eventnum).nextevent = 0; % zero indicates we are done with the sequence

                case 'jump'              % ----------------------
                    seqcmd = seqcmd + 1; % increment the count of event sequencing commands
                    Event(eventnum).nextevent = argument;
                    if argument < eventnum
                        % this is a jump backwards, so we must be at the
                        % end of the sequence and starting to repeat the
                        % endless run-time loop
                        finalHPEB = 1;
                        % We will follow the jump and keep processing
                        % events to fill in the PRI for last HPEB event in
                        % sequence; the finalHPEB flag indicates that it is
                        % time to quit when the next HPEB is reached after
                        % the jump back
                        if TEnum == 0 || Event(eventnum).repeat > 1
                            % there are no HPEB events in the sequence so
                            % just quit instead of jumping back; also quit
                            % if we've already been here more than once, to
                            % avoid an endless loop
                            Event(eventnum).nextevent = 0;
                        end
                    end

                case 'noop'             % ----------------------
                    % if this noop is not in a tr event, then set tr to 2 so it
                    % will be counted as if it was part of a non-HPEB tr event
                    noop = noop + argument/5; % noop is in 0.2 usec units

                case {'pause', 'triggerIn', 'multiSysSync'}     % external trigger input commands
                    if verbose > 1
                        fprintf( ['TXEventCheck: A "pause", "triggerIn", or "multiSysSync" command for external trigger input was detected in Event number ' num2str(eventnum) '. \n']);
                        fprintf('   TXEventCheck will not account for any delay that may be introduced while waiting for the trigger input. \n');
                    end

                case {'markTransferProcessed', 'returnToMatlab', 'triggerOut' 'waitForTransferComplete' 'encoderWait'}
                    % These commands are ignored by TXEventCheck since they have
                    % no impact on the event sequence timing.

                case {'transferToHost', 'sync', 'setRcvProfile'}
                    % These commands are ignored by TXEventCheck based on
                    % the assumption that in normal operation they will have
                    % no impact on the HW event sequence timing.

                    % Note however that these commands will add delay to
                    % the event sequence if the HW sequencer has to pause
                    % to wait for a DMA to complete or for the SW sequence
                    % to catch up, in synchronous operation.

                case 'cBranch'
                    % This command is ignored by TXEventCheck, meaning it
                    % will never take the conditional branch and instead
                    % will just proceed to the following event.  Thus the
                    % Event sequence code that would be executed when the
                    % branch is taken will not be evaluated by
                    % TXEventCheck.  This needs to be corrected in a future
                    % release, by providing some means for TXEventCheck to
                    % evaluate the sequence for both alternatives at the
                    % cBranch point.

                case 'rdmaWrite'
                    % placeholder for rdmawrite check

                case 'rdmaSync'
                    % placeholder for rdmawrite check

                otherwise
                    % an unrecognized command has been found, so report an
                    % error and quit

                    fprintf(2,['TXEventCheck Error: Unrecognized or unsupported SeqControl command "', command, '"\n']);
                    fprintf(2,['was detected in Event number ', num2str(eventnum),  '.\n\n']);
                    error('TXEventCheck: Exiting due to error condition identified above.')
                    return

            end
        end
    end % finished with seqControl array for this event
    % --- END 1 ---------- End of for-loop to decode SeqControl commands for this event --------------- 1 -------------------

        % Test for error if multiple sequencing commands in one event
        if seqcmd > 1
            fprintf(2,['TXEventCheck Error: Multiple event sequencing commands detected in Event ', num2str(eventnum),  ' .\n']);
            fprintf(2, '(loopCnt, loopTst, stop, jump, call, or rtn commands.)\n');
            fprintf('Correct the error and then run the SetUp script and VSX again.\n\n');
            error('TXEventCheck: Exiting due to error condition identified above.')
            return
        end

    % --- 2: Assign atomic PRI values ('atomic' meaning each event taken by itself for HPEB events). -------- 2 --------------
    %   ... find the PRI if this is a tr event
    if tr > 0 % find the PRI since this is a TR event
        PRI = 0; % start out at zero for PRI (either TX or Receive may be absent)
        if TXnum > 0  % this event does have a TX?
            if TX(TXnum).VDASBistDriveEnable == 0
                % We have an active TX not using the BIST driver
                % Check for presence of profile 5 vs imaging profiles
                % 1:4
                if (profile == 5 && TX(TXnum).imgProfileMaxHv < TPC(profile).maxHighVoltage) ...
                        || (nnz(TX(TXnum).Apod) > 0 && nnz(TX(TXnum).Apod) < 17)
                    % This TX event is using TPC profile 5, with max
                    % voltage level that exceeds capacity of imaging supply
                    % so this is a "HPEB" event to which the thermal and
                    % HIFU current limits must be applied (otherwise, we
                    % view this as an 'imaging' transmit done with TPC
                    % profile 5 but with energy level low enough that it
                    % does not contribute significantly to device
                    % dissipation).
                    % (VTS-2509) Alternatively, if less than 17 channels
                    % are active we need to evaluate thermal limits
                    % regardless of which profile is active since for small
                    % active channel counts the imaging HV supply can
                    % provide enough power to damage a channel
                    Event(eventnum).HPEB = 1; % flag for Profile 5 Extended Burst
%                     TX(TXnum).sysExtendBL = 1; % Don't force to 1 if user
%                     didn't set it, but already is using profile 5 (to
%                     avoid potentially changing how any TTNEB commands
%                     will behave)
                end
                % (VTS-2509) regardless of whether HPEB flag got set above,
                % we now need to check separately for use of an imaging
                % provile and limit HV if beyond the imaging supply
                % capacity
                if profile < 5
                    % This TX event is using one of the imaging profiles 1
                    % through 4; check if the TX exceeds capacity of
                    % imaging supply
                    if TX(TXnum).imgProfileMaxHv < TPC(profile).highVoltageLimit
                        if TX(TXnum).imgProfileMaxHv > P5rqdHVmin
                            % This burst can be executed using the imaging
                            % TPC, if we reduce the highVoltageLimit
                            TPC(profile).highVoltageLimit = floor(TX(TXnum).imgProfileMaxHv);
                            if verbose
                                % warn the user if verbose > 0
                                fluxLimStr1 = ['TXEventCheck Warning: reducing TPC(', num2str(profile), ').highVoltageLimit to ', ...
                                    num2str(TPC(profile).highVoltageLimit), ' Volts due to \n'];
                                fluxLimStr2 = ['imaging transmit power supply capacity limit for TW(', num2str(TWnum), ') in Event(', num2str(eventnum), ').\n'];
                                fluxLimStr3 = 'TPC Profile 5 must be used for this Transmit Event if a higher Voltage is desired.\n';
                            end
                        else
                            % this burst requires imaging HV to be less
                            % than P5rqdHVmin; report an error condtion
                            % that TPC profile 5 is required
                            fprintf(2,['TXEventCheck Error: In Event ', num2str(eventnum), ', TX(', num2str(TXnum),'), using TPC profile ', num2str(profile),...
                                ' the transmit waveform and active\n']);
                            fprintf(2, 'aperture exceed the capacity of the imaging transmit profiles; TPC profile 5 must be used.\n');
                            fprintf('Correct the error and then run the SetUp script and VSX again.\n\n');
                            TX_Limits.Event = Event;
                            assignin('base','TX_Limits',TX_Limits);
                            TX(TXnum).sysExtendBL = 1;
                            assignin('base','TX',TX);
                            error('TXEventCheck: Exiting due to error condition identified above.')
                            return
                        end
                    end
                end
            end
            PRI = TW(TWnum).Bdur  +  max(TX(TXnum).Delay)/Trans.frequency;
            % Burst duration in usec based on TX.Bdur plus max TX.Delay
            % based on Trans.frequency
        end

        if Event(eventnum).rcv > 0
            PRI = max(PRI, Receive(Event(eventnum).rcv).endDepth*2/Trans.frequency); % or duration of receive interval if longer
        end

        % since this is a TR event, find the actual PRI from the previous
        % one
        lastPRI = max(activePRI, activettna);

        % set new active PRI and ttna values for the current TR event, so we can find the biggest when
        % we reach the next PRI
        activePRI = PRI + PRIovh + noop; % add per-event overhead and any noop delays found in this event
        activettna = ttna; % save any ttna value in the event structure to refer to later

    else
        % This is not a TR event so add any noop from this event to the
        % activePRI value
        activePRI = activePRI + noop;
        lastPRI = 0;
    end
    % -------- END 2 -------------- finished finding PRI for this tr event ------------------------------ 2 --------------


            if (verbose > 2)
                fprintf(['Event#= ' num2str(eventnum)...
                    '\t Prof= ' num2str(profile)...
                    '\t tr= '  num2str(tr)...
                    '\t PRI= '  num2str(max(PRI, ttna))...
                    '\t startEvent= ' num2str(startEvent)...
                    '\t\t ' Event(eventnum).info ...
                    '\n' ] )
            end



    % --- 3: Create TXEvent if this is a HPEB Event ----------- 4 --------------
    if Event(eventnum).HPEB == 1
        if TEnum == 0
            % there has not been a previous HPEB event, so we can't do
            % anything with current "active" values from previous tr events
            % or noops
        else
            % find the actual cumulative PRI for the previous TXEvent,
            % since we are now starting a new one
            TXEvent(TEnum).cumPRI = TXEvent(TEnum).cumPRI + lastPRI; % add PRI from previous event
            TXEvent(TEnum).cumPRI = max(TXEvent(TEnum).cumPRI, activettnEB); % use ttnEB if greater
        end
        if finalHPEB
            % we are at the beginning of the second pass through the entire
            % sequence; PRI for the last TXEvent has just been updated so
            % now we can quit
            Event(eventnum).nextevent = 0; % zero indicates we are done parsing through the sequence
        else
            % not finished yet, so create the initial values for the new TXEvent
            TEnum = TEnum + 1;
            activettnEB = ttnEB;
            TXEvent(TEnum).Event = eventnum;
            TXEvent(TEnum).tx = TXnum;
            TXEvent(TEnum).tw = TWnum;
            TXEvent(TEnum).profile = profile;
            TXEvent(TEnum).cumPRI = 0;
            Bdur(TEnum) = TW(TWnum).Bdur; % burst duration in usec
            Numpulses(TEnum, :) = TX(TXnum).Numpulses;
            CumOnTime(TEnum, :) = TX(TXnum).CumOnTime;
            Zload(TEnum) = TW(TWnum).Zload;
            Zsource(TEnum) = TW(TWnum).Zsource;
            HPEBprofiles(profile) = 1;

            % Now create the per-channel ChIpk1V array, representing a sine
            % wave at estimatedAverageFrequency with amplitude based on a
            % square wave with relative pulse width set by CumOnTime/Bdur
            TXEvent(TEnum).ChIpk1V = 4*sin(CumOnTime(TEnum, :)*pi/(2*Bdur(TEnum)))/(pi*abs(Zload(TEnum) + Zsource(TEnum)));

            % fill in placeholders for other values to be added later
            TXthermal(TEnum).ChPrGdP = 0;
            TXthermal(TEnum).BdGdI = 0;
            TXthermal(TEnum).ChPrGdTempRise = 0;
            TXthermal(TEnum).BdGdVdroop = 0;
            TXthermal(TEnum).BdGdPSTempRise = 0;
            TXthermal(TEnum).ChFETTempRise = 0;
            TXthermal(TEnum).ChXfmrTempRise = 0;
            TXthermal(TEnum).BdTxBusTempRise = 0;
        end
    elseif TEnum>0
        % This is not a HPEB event, but there has been a previoius one so
        % add PRI from the current event.
        TXEvent(TEnum).cumPRI = TXEvent(TEnum).cumPRI + lastPRI;
    end

    % --- 3: update the profile number to be used in the next event  ----------- 3 --------------
    profile = nextprofile;

    eventnum = Event(eventnum).nextevent; % go to the next event
    % (or drop out of the while loop if nextevent = 0)

end
% finished with the pass through event list, to decode SeqControl, find TR
% event PRI's, and build the TXEvent list of HPEB events

assignin('base', 'TX', TX); % write TX back to base workspace with updated sysExtendBL values

TX_Limits.Event = Event; % save stuff added to Event structure

if finalHPEB == 0 && TEnum > 0
    % If finalHPEB is zero, that means we never made a jump backwards &
    % thus this was a one-shot event sequence.  So the final TXEvent, if there was
    % one, still needs to get updated to reflect a ttnEB if that was longer
    % than the cumulative PRI's: (on a repeating sequence we would follow
    % the jump, and this update would happen when we encountered the first
    % HPEB event and then jumped out of the while loop)
    if TXEvent(TEnum).cumPRI == 0
        % This was a one-shot sequence that had no other TR events after
        % the last TXEvent, so it currently has a cumPRI of zero.  Add the
        % current active values, representing it's own PRI.
        TXEvent(TEnum).cumPRI = max(activePRI, activettna);
    end
    TXEvent(TEnum).cumPRI = max(TXEvent(TEnum).cumPRI, activettnEB); % use ttnEB if greater
end

% Check for highVoltageLimit values less than the current hv setting for
% each imaging profile
hv2GUIprofile = evalin('base', 'hv2GUIprofile');
for i = 1:5
    if TPC(i).inUse
        % note unused profiles may have empty values and we will leave them
        % that way.
        if TPC(i).hv > TPC(i).highVoltageLimit
            % current setting exceeds the new limit, so reduce it.
            newhv = TPC(i).highVoltageLimit-.01;
            if i == 1
                % TPC 1 HV slider needs to be updated to show new voltage setting
                set(findobj('Tag','hv1Sldr'),'Value',newhv);
                set(findobj('Tag','hv1Value'),'String',num2str(newhv,'%.1f'));
            elseif i == hv2GUIprofile
                % 2nd HV slider needs to be updated to show new voltage setting
                set(findobj('Tag','hv2Sldr'),'Value',newhv);
                set(findobj('Tag','hv2Value'),'String',num2str(newhv,'%.1f'));
            end
            TPC(i).hv = newhv;
        end
    end
end

% print warning strings if they were created
if ~isempty(fluxLimStr1)
    fprintf(2, fluxLimStr1);
    fprintf(2, fluxLimStr2);
    fprintf(2, fluxLimStr3);
end

% ================================================================================================================
% Print out a message to say we've been here
if TEnum>0
    if verbose > 2
        disp(['TXEventCheck event parsing complete; generated ', num2str(TEnum), ' TXEvent structures.']);
    end
    numTXEvents = TEnum;
else
    assignin('base','TPC', TPC);
    % also clear the TXEvent structure TXEventCheck normally creates
    evalin('base','TX_Limits = [];');
    if verbose > 2
        disp('TXEventCheck event parsing complete; no profile 5 extended burst events are present.');
    end
    return
end

limProfiles = find(HPEBprofiles);
%% Fill in cumulative totals from each TXEvent Structure and check limits

% HVLim must be set on a per-event basis for the TPC Profile being used by
% that event (VTS-2509; previously it was assumed that limits needed to be
% evaluated only for events using TPC profile 5, at power levels beyond the
% imaging supply capacity)

% Gate Driver power supply output current is directly proportional to
% number of channels active and Burst frequency in MHz. According to Larry,
% all 64 channels active at 4 MHz would produce a total per-board load
% current of 4 Amps.
gdIperCHF = 1/64; % gate driver supply current in Amps, per channel at 1 MHz

gdVGD = 12; % Gate drive power is 12 V. times supply current

HPEBchIpkLim = 2; % maximum allowed peak per-channel output current in Amps for HPEB events
% note FET switching loss model now increases dramatically at currents
% above 1.5 Amps.  This allows the maximum current limit to be raised from
% 1.5 Amps back to 2 Amps.
% This limit is 2 Amps for non-extended burst transmits.


% Here we calculate both cumulative and per-TXEvent values that will be
% utilized repeatedly in the following limit evaluation loops through the
% sequence.  Initialize variables to zero, then accumulate over all TXEvents
totalPRI = 0;
idleT = zeros(1, size(TXEvent, 2)); % idle time per event (cumPRI - Bdur)

TX_Limits.TXEvent = TXEvent;

% create array indices for summing pairs of channels sharing a gate driver;
% within each AFE group of 8 channels, channels N and N+4 are using the
% same gate driver package.
incr4 = [1, 2, 3, 4];
GDindxA = incr4;
for i = 1:(numBoards*8 - 1)
    GDindxA = [GDindxA, (incr4+8*i)];
end
GDindxB = GDindxA + 4;
TX_Limits.GDindxB = GDindxB;
assignin('base', 'TX_Limits', TX_Limits);


bdCumGdPsLd = 0;
for TEnum = 1:numTXEvents
    HVLim = TPC(TXEvent(TEnum).profile).highVoltageLimit - clampV;
    % HVLim is max transmit voltage for this profile after subtracting the clamp diode drop
    totalPRI = totalPRI + TXEvent(TEnum).cumPRI; % accumulate PRI over all TXEvents
    idleT(TEnum) = TXEvent(TEnum).cumPRI - Bdur(TEnum);

    % find per-channel gate driver supply current based on average
    % cycle rate over the burst interval, and also the total gate driver
    % supply current per board by summing over all channels on each board
    TXthermal(TEnum).ChGdIperCh = gdIperCHF * Numpulses(TEnum, :)/(2*Bdur(TEnum));
    ChPrGdI = TXthermal(TEnum).ChGdIperCh(GDindxA) + TXthermal(TEnum).ChGdIperCh(GDindxB);
    for bdnum = 1:numBoards
        chindex = (64*(bdnum-1) + (1:64));
        BdGdI(bdnum) = sum(TXthermal(TEnum).ChGdIperCh(chindex));
    end
    bdCumGdPsLd = bdCumGdPsLd + Bdur(TEnum)*BdGdI; % total gate driver supply load in Amp-usec
    TXthermal(TEnum).BdGdI = BdGdI;
    TXthermal(TEnum).ChPrGdP = ChPrGdI*gdVGD;


    % find peak per channel Xmt current and enforce the limit
    maxChIpk1V = max(TXEvent(TEnum).ChIpk1V); % max current over all channels
    if maxChIpk1V*HVLim > HPEBchIpkLim
        % Over the peak current limit, so reduce HVLim as needed
        newHVLim = HPEBchIpkLim/maxChIpk1V;
        newHVLim = max(floor(newHVLim + clampV), 2) - clampV; % round hvLimit down to integer value

        % Now check newHVLim against Larry's peak current plus peak voltage
        % limit, but only for standard freq. system
        if TXindex == 1
            Vscale = (1/95 + maxChIpk1V/1.5);
            if (newHVLim * Vscale) > 1.5
                newHVLim = 1.5/Vscale;
                newHVLim = max(floor(newHVLim + clampV), 2) - clampV; % round hvLimit down to integer value
            end
        elseif TXindex == 5
            % active clamp standard frequency board
            Vscale = (1/95 + maxChIpk1V/1.5);
            if (newHVLim * Vscale) > 1.5
                newHVLim = 1.5/Vscale;
                newHVLim = max(floor(newHVLim + clampV), 2) - clampV; % round hvLimit down to integer value
            end
        end

        % we reduced HVLim; only warn the user if verbose is enabled
        if verbose > 0 % display warning only if enabled
            disp(' ')
            disp(['Peak transmit output current of ',num2str(maxChIpk1V*HVLim,'%.2f'), ' Amps exceeds ',num2str(HPEBchIpkLim,'%.2f'), ' Amp limit.']);
            disp(['Reducing hvLimit to ', num2str(newHVLim+clampV, '%.0f'), ' Volts to conform to the limit.']);
            disp(' ')
        end
        HVLim = newHVLim; % Update HVLim with the new limiting value
        TPC(TXEvent(TEnum).profile).highVoltageLimit = HVLim + clampV;
    end

end



%% Loops to evaluate the thermal model-based limits

% ***************** GATE DRIVER DEVICE DISSIPATION LIMIT *******************
Lim.ChPrGdTempRise = 35; % 125 degree absolute die limit, with 90 degree PCB temperature
Rth = 23; % Thermal resistance of the gate driver package, in degrees C. per Watt
% Note that this thermal resistance applies to the entire package, serving
% two channels.  The gate driver dissipation values for the pair of
% channels served by one package have been summed, so this limit is
% evaluated for each gate driver package not each channel.
tau = 100e3; % 100 msec thermal time constant scaled to usec for gate driver package
prevValue = 0; % stores the value found in previous event; initialized to zero
maxValue = 0; % highest value encountered so far
oldmaxValue = 0; % stores the maxValue from previous pass throgh while loop, to check if we have reached steady state
rising = 1; % controls when to exit the while loop
while rising
    for TEnum = 1:numTXEvents
        TXthermal(TEnum).ChPrGdTempRise  = TXthermal(TEnum).ChPrGdP*Rth  -  (TXthermal(TEnum).ChPrGdP*Rth - prevValue)*exp(-Bdur(TEnum)/tau);
        maxValue = max(maxValue, max(TXthermal(TEnum).ChPrGdTempRise));
        % this gives us temperature at the end of the current event's transmit
        % burst, based on prevTemp from the previous TXEvent
        prevValue = TXthermal(TEnum).ChPrGdTempRise*exp(-idleT(TEnum)/tau);
    end
    if maxValue <= oldmaxValue || ((maxValue - oldmaxValue)/maxValue < .01) || maxValue > Lim.ChPrGdTempRise
        % exit from while loop immediately if we have already passed the
        % limit, otherwise wait until maxValue settles within 1 percent.
        rising = 0;
    end
    oldmaxValue = maxValue;
end

TX_Limits.TXthermal = TXthermal;

if maxValue > Lim.ChPrGdTempRise
    disp(' ');
    fprintf(2,'ERROR:\n');
    fprintf(2,['Gate driver device estimated temperature rise of ',num2str(maxValue,'%.1f'),...
        ' degrees C. exceeds ', num2str(Lim.ChPrGdTempRise,'%.1f'), ' degrees maximum limit.\n\n']);
    assignin('base','TX_Limits',TX_Limits); % save the TXEvent structure completed thus far
    error('TXEventCheck: Exiting due to error condition identified above.')
    return
end

% ***************** GATE DRIVER POWER SUPPLY DROOP LIMIT *******************
gdPsImax = 4; % maximum gate driver supply current in Amps
Cgdps = 1100; % power supply storage capacitor in uF
maxVdroopLimit = 0.25; % maximum allowed droop in volts at power supply storage capacitor

% First check for total gate driver supply load (in Amp-usec) exceeding
% its total capacity over all events for this event sequence, so droop will
% increase with each pass through event sequence
if max(bdCumGdPsLd)>gdPsImax*totalPRI
    disp(' ');
    fprintf(2,'ERROR:\n');
    fprintf(2,['Maximum per-board total Gate driver supply load of ', num2str(max(bdCumGdPsLd),'%.1f'),...
        ' Amp-usec exceeds power supply capacity of ', num2str(gdPsImax*totalPRI,'%.1f'), ' Amp-usec.\n\n']);
    assignin('base','TX_Limits',TX_Limits); % save the TXEvent structure completed thus far
    error('TXEventCheck: Exiting due to error condition identified above.')
    return
end

% now that we know droop will reach a stable value, loop through event
% sequence to find event with highest droop and check against limit.
prevVdroop = 0;
maxVdroop = 0;
oldmaxVdroop = 0;
rising = 1; % controls when to exit the while loop
while rising
    for TEnum = 1:numTXEvents
        TXthermal(TEnum).BdGdVdroop = prevVdroop + Bdur(TEnum) * max(0, TXthermal(TEnum).BdGdI - gdPsImax) / Cgdps;
        maxVdroop = max(maxVdroop, max(TXthermal(TEnum).BdGdVdroop));
        prevVdroop = max(0, TXthermal(TEnum).BdGdVdroop - gdPsImax*idleT(TEnum)/Cgdps);
    end
    if maxVdroop<=oldmaxVdroop || maxVdroop > maxVdroopLimit
        % exit from while loop immediately if we have already passed the
        % limit, otherwise wait until maxVdroop settles.
        rising = 0;
    end
    oldmaxVdroop = maxVdroop;
end

TX_Limits.TXthermal = TXthermal;
if maxVdroop > maxVdroopLimit
    disp(' ');
    fprintf(2,'ERROR:\n');
    fprintf(2,['Gate driver power supply estimated preregulator droop of ',num2str(maxVdroop,'%.2f'),...
        ' Volts exceeds ', num2str(maxVdroopLimit,'%.2f'), ' Volt maximum limit.\n\n']);
    assignin('base','TX_Limits',TX_Limits); % save the TXEvent structure completed thus far
    error('TXEventCheck: Exiting due to error condition identified above.')
    return
end


% ************* GATE DRIVER POWER SUPPLY DISSIPATION LIMIT ***************
Lim.BdGdPSTempRise = 35; % 125 degree absolute die limit, with 90 degree PCB temperature
Rth = 28; % Thermal resistance of the power supply device package, in degrees C. per Watt
Rdiss = 0.2; % effective series resistance in power supply that leads to the dissipation
tau = 100e3; % 100 msec thermal time constant scaled to usec for gate driver package
prevValue = 0; % stores the value found in previous event; initialized to zero
maxValue = 0; % highest value encountered so far
oldmaxValue = 0; % stores the maxValue from previous pass throgh while loop, to check if we have reached steady state
rising = 1; % controls when to exit the while loop
while rising
    for TEnum = 1:numTXEvents
        % find Isquared R dissipation level for this burst
        bdGdPSdiss = Rdiss * TXthermal(TEnum).BdGdI.*TXthermal(TEnum).BdGdI;
        TXthermal(TEnum).BdGdPSTempRise  = bdGdPSdiss*Rth  -  (bdGdPSdiss*Rth - prevValue)*exp(-Bdur(TEnum)/tau);
        maxValue = max(maxValue, max(TXthermal(TEnum).BdGdPSTempRise));
        % this gives us temperature at the end of the current event's transmit
        % burst, based on prevTemp from the previous TXEvent
        prevValue = TXthermal(TEnum).BdGdPSTempRise*exp(-idleT(TEnum)/tau);
    end
    if maxValue <= oldmaxValue || (maxValue - oldmaxValue)/maxValue < .01 || maxValue > Lim.BdGdPSTempRise
        rising = 0;
    end
    oldmaxValue = maxValue;
end

TX_Limits.TXthermal = TXthermal;
if maxValue > Lim.BdGdPSTempRise
    disp(' ');
    fprintf(2,'ERROR:\n');
    fprintf(2,['Gate driver Power Supply estimated temperature rise of ',num2str(maxValue,'%.1f'),...
        ' degrees C. exceeds ', num2str(Lim.BdGdPSTempRise,'%.1f'), ' degrees maximum limit.\n\n']);
    assignin('base','TX_Limits',TX_Limits); % save the TXEvent structure completed thus far
    error('TXEventCheck: Exiting due to error condition identified above.')
    return
end

% ***************** TRANSFORMER DISSIPATION LIMIT *******************
% set limits and transformer parameters based on system transmit configuration

% Coefficients for the FET switching are also set here, for each TXindex
% value.  FET switching loss is based on a model that switching loss
% represents four terms all of which scale with the "cycle rate" in MHz
% where a "cycle" is two pulses or four FET switching transitions.  Thus
% cycle rate in MHz is equal to Numpulses/(2 * Bdur) with Bdur in usec.
% The four terms are associated with coefficients from the Cswl array as
% listed below:

    % Cswl(1): constant term independent of HV, due to dissipation
    % driven by the gate driver power supply.

    % Cswl(2): losses directly proportional to HV

    % Cswl(3): losses proportional to HV squared

    % Cswl(4): losses proportional to transmit power, or HV * Iout

    % Cswl(5): additional losses proportional to transmit power, or HV *
    % Iout, when peak Iout exceeds 1.5 Amps  (i.e. Cswl(4) gets larger when
    % Iout > 1.5 A., since in this situation with a reactive load
    % cross-conduction can start to occur whith the H-bridge, dramatically
    % increasing dissipation in the FETs.)

    % Cswl(6): This coefficient is not used to estimate FET
    % switching loss; the FET switching loss model is based on
    % estimated FET dissipation for a 'worst-case' individual FET
    % package and thus the sum of this estimate over all channels
    % would overstate the expected total system switching loss.
    % Instead, Cswl(6) is used to estimate typical switching loss
    % per channel driven by the HV supply, based on this
    % coefficient times the cycle rate (not proportional to HV?).  At
    % presetn Cswl(6) is set to zero, and total system dissipation
    % estimate for switching loss is based on total gate driver power
    % consumption, with no regard to how it is partitioned between FETs and
    % gate drivers.


switch TXindex
    case 1
        % system has the standard frequency transformer DA2319:
        Lim.ChXfmrTempRise = 30; % 120 degree absolute limit, with 90 degree PCB temperature
        Rth = 130; % Thermal resistance of the Transformer, in degrees C. per Watt
        tau = 50e3; % 50 msec thermal time constant scaled to usec for transformer
        windingR = 1.45; % total transformer winding resistance in Ohms from the DA2319 data sheet
        % coefficients for Transmit FET switching loss
        % (see definitions in the FET section)
        Cswl = [0.0695, 0, 0, 0, 0.03, 0]; % These values from Larry's analysis Aug. 3,
        % and increased "conduction loss" term as well (3.39 Ohms instead
        % of 2.35 from FET data sheet- see FET section)
    case 2
        % system has the original high frequency transformer PWB1010-1L:
        Lim.ChXfmrTempRise = 30; % 120 degree absolute limit, with 90 degree PCB temperature
        Rth = 1200; % Thermal resistance of the Transformer, in degrees C. per Watt
        tau = 50e3; % 50 msec thermal time constant scaled to usec for transformer
        windingR = 0.4; % total transformer winding resistance in Ohms
        % coefficients for Transmit FET switching loss
        % (see definitions in the FET section)
        Cswl = [0.1875, 0, 0, 0.01, 0.03, 0];
    case 3
        % system has the low frequency transformer MSD7342-105ML_1000 :
        Lim.ChXfmrTempRise = 30; % 120 degree absolute limit, with 90 degree PCB temperature
        Rth = 130; % Thermal resistance of the Transformer, in degrees C. per Watt
        tau = 50e3; % 50 msec thermal time constant scaled to usec for transformer
        windingR = 15.6; % total transformer winding resistance in Ohms from the transformer data sheet
        % coefficients for Transmit FET switching loss
        % (see definitions in the FET section)
        Cswl = [0.0695, 0, 0, 0, 0.03, 0]; % These values are TBD,
        % currently just a copy of the standard frequency TXindes = 1
        % values
    case 4
        % system has the newer high frequency transformer PWB1010L:
        Lim.ChXfmrTempRise = 30; % 120 degree absolute limit, with 90 degree PCB temperature
        Rth = 750; % Thermal resistance of the Transformer, in degrees C. per Watt
        tau = 50e3; % 50 msec thermal time constant scaled to usec for transformer
        windingR = 0.64; % total transformer winding resistance in Ohms
        % coefficients for Transmit FET switching loss
        % (see definitions in the FET section)
        Cswl = [0.1875, 0, 0, 0.01, 0.03, 0];
    case 5
        % Active clamp system with the standard frequency transformer DA2319:
        Lim.ChXfmrTempRise = 30; % 120 degree absolute limit, with 90 degree PCB temperature
        Rth = 130; % Thermal resistance of the Transformer, in degrees C. per Watt
        tau = 50e3; % 50 msec thermal time constant scaled to usec for transformer
        windingR = 1.45; % total transformer winding resistance in Ohms from the DA2319 data sheet
        % coefficients for Transmit FET switching loss
        % (see definitions in the FET section)
        Cswl = [0.0695, 0, 0, 0, 0.03, 0]; % These values from Larry's analysis Aug. 3,
        % and increased "conduction loss" term as well (3.39 Ohms instead
        % of 2.35 from FET data sheet- see FET section)
    case 6
        % Active clamp system with high frequency transformer PWB1010L:
        Lim.ChXfmrTempRise = 30; % 120 degree absolute limit, with 90 degree PCB temperature
        Rth = 750; % Thermal resistance of the Transformer, in degrees C. per Watt
        tau = 50e3; % 50 msec thermal time constant scaled to usec for transformer
        windingR = 0.64; % total transformer winding resistance in Ohms
        % coefficients for Transmit FET switching loss
        % (see definitions in the FET section)
        Cswl = [0.1875, 0, 0, 0.01, 0.03, 0];
    case 7
        % Active clamp system with low frequency transformer MSD7342-105ML_1000 :
        Lim.ChXfmrTempRise = 30; % 120 degree absolute limit, with 90 degree PCB temperature
        Rth = 130; % Thermal resistance of the Transformer, in degrees C. per Watt
        tau = 50e3; % 50 msec thermal time constant scaled to usec for transformer
        windingR = 15.6; % total transformer winding resistance in Ohms from the transformer data sheet
        % coefficients for Transmit FET switching loss
        % (see definitions in the FET section)
        Cswl = [0.0695, 0, 0, 0, 0.03, 0]; % These values are TBD,
        % currently just a copy of the standard frequency TXindes = 1
        % values
    otherwise
        % unrecognized TXindex value
        fprintf(2,'ERROR: Unrecognized value of Resource.SysConfig.TXindex. \n');
        error('TXEventCheck: Exiting due to error condition identified above.');
end
iterate = 1;
% if a TPC highVoltageLimit adjustment was needed, repeat the temperature
% rise evaluation since different channels may exceed the limits while
% using different TPC profiles
while iterate
    iterate = 0; % don't repeat if no adjustments were needed
    prevValue = 0; % stores the value found in previous event; initialized to zero
    maxValue = 0; % highest value encountered so far
    oldmaxValue = 0; % stores the maxValue from previous pass throgh while loop, to check if we have reached steady state
    rising = 1; % controls when to exit the while loop
    while rising
        for TEnum = 1:numTXEvents
            HVLim = TPC(TXEvent(TEnum).profile).highVoltageLimit - clampV; % HVLim for profile in this event
            TXthermal(TEnum).XfmrP = windingR * HVLim * HVLim * 0.5 * TXEvent(TEnum).ChIpk1V.*TXEvent(TEnum).ChIpk1V;
            TXthermal(TEnum).ChXfmrTempRise  = TXthermal(TEnum).XfmrP*Rth  -  (TXthermal(TEnum).XfmrP*Rth - prevValue)*exp(-Bdur(TEnum)/tau);
            maxValue = max(maxValue, max(TXthermal(TEnum).ChXfmrTempRise));
            % this gives us temperature at the end of the current event's transmit
            % burst, based on prevTemp from the previous TXEvent
            prevValue = TXthermal(TEnum).ChXfmrTempRise*exp(-idleT(TEnum)/tau);
%             if maxValue > Lim.ChXfmrTempRise
%                 % note which TPC profile triggered overlimit condition
%                 limProfile = TXEvent(TEnum).profile;
%             end
        end
        if maxValue <= oldmaxValue || (maxValue - oldmaxValue)/maxValue < .01
            rising = 0;
        end
        oldmaxValue = maxValue;
    end

    if maxValue > Lim.ChXfmrTempRise
        HVLim = TPC(limProfiles(1)).highVoltageLimit - clampV;
        newHVLim = HVLim * sqrt(Lim.ChXfmrTempRise/maxValue); % new HVLim value needed to stay at the limit
        newHVLim = max(floor(newHVLim + clampV), 2.5) - clampV; % round hvLimit down to integer value
        % we need to reduce HVLim; but only warn the user if verbose is enabled
        if verbose > 0 % display warning only if enabled
            disp(' ')
            disp(['Transformer internal temperature rise of ',num2str(maxValue,'%.0f'), ' Deg. C. exceeds ',num2str(Lim.ChXfmrTempRise,'%.0f'), ' Degree limit.']);
            disp(['Reducing TPC Profile hvLimit to ', num2str(newHVLim+clampV, '%.0f'), ' Volts to conform to the limit.']);
            disp(' ')
        end
        if newHVLim >= HVLim
            % last iteration did not produce a new value, so stop iterations
            iterate = 0;
        else
            HVLim = newHVLim; % Update HVLim with the new limiting value
            for pnum = 1:length(limProfiles)
                TPC(limProfiles(pnum)).highVoltageLimit = HVLim + clampV;
            end
            for TEnum = 1:numTXEvents
                TXthermal(TEnum).ChXfmrTempRise = TXthermal(TEnum).ChXfmrTempRise * Lim.ChXfmrTempRise/maxValue;
            end
            iterate = 1; % repeat temp rise loop since an adjustment was needed
        end
    end
end


% ***************** TRANSMIT FET DISSIPATION LIMIT *********************
Lim.ChFETTempRise = 60; % 150 degree absolute die limit, with 90 degree PCB temperature
if TPC(5).inUse == 2
    % HIFU configuration with heat sinks on the acquisition boards
    Rth = 64; % Thermal resistance of the FET device package, in degrees C. per Watt
else
    % Extended Transmit option without added heat sinks
    Rth = 100; % Higher resistance for no heat sink configuration
end
% Thermal resistance derived from an operating state assumed to be at the
% dissipation limit: 1 MHz continuous transmit at 33 V. and 0.4 A. RMS load
% current, leading to 0.94 W. FET package dissipation (0.56 switching +
% 0.38 conduction).
Rdiss = 3.39; % Effective "on resistance" from Larry Aug. 3
% (actually just the current squared term from Larry's model derived from
% his simulation and analysis)

% % Rdiss = 2.35; % effective on-resistance of FET's that leads to the conduction loss dissipation, 2.35 Ohms from data sheet
tau = 100e3; % 100 msec thermal time constant scaled to usec for FET package


% This while iterate loop will evaluate FET temperature rise from the
% thermal model, and then check against the limit.  If limit is exceeded
% HVLim will be reduced but since dissipation depends on a combination of
% HVLim (switching loss) and HVLim squared (conduction loss), the reduction
% may not have been enough so the while iterate loop will repeat the entire
% process until we reach a stable state where HVLim was not reduced on the
% last pass.
iterate = 1;
reduceHVstring = []; % if it stays empty there are no reductions to report
while iterate
    iterate = 0;
    prevValue = 0; % stores the value found in previous event; initialized to zero
    maxValue = 0; % highest value encountered so far
    oldmaxValue = 0; % stores the maxValue from previous pass throgh while loop, to check if we have reached steady state
    rising = 1; % controls when to exit the while loop
    for TEnum = 1:numTXEvents
        HVLim = TPC(TXEvent(TEnum).profile).highVoltageLimit - clampV; % HVLim for profile in this event
        % estimate switching loss for current value of HV:
        % the following expression estimates switching loss using the
        % Cswl coefficient array defined above (see "Transformer Dissipation limit"
        % section for values and meaning of each term)

        TXthermal(TEnum).FETSWloss = Numpulses(TEnum, :)/(2*Bdur(TEnum)).*  ...
            ((Cswl(1) + Cswl(2)*HVLim + Cswl(3) * HVLim * HVLim +...
            Cswl(4) * HVLim * HVLim .* TXEvent(TEnum).ChIpk1V) +...
            (Cswl(5) * HVLim .* max((HVLim .* TXEvent(TEnum).ChIpk1V - 1.5), 0)));


        % Now find Isquared R conduction loss dissipation level for this
        % burst (divide peak value by two to get RMS current level squared
        % during the burst, for the sinusoidal approximation of actual
        % current waveform)
        TXthermal(TEnum).FETcondloss = Rdiss * HVLim^2*TXEvent(TEnum).ChIpk1V.*TXEvent(TEnum).ChIpk1V / 2;
    end
    while rising
        for TEnum = 1:numTXEvents
            % Add the gate drive loss, switching losses, and conduction loss to get total FET
            % power dissipation
            chFETdiss = TXthermal(TEnum).FETSWloss + TXthermal(TEnum).FETcondloss;
            TXthermal(TEnum).ChFETTempRise  = chFETdiss*Rth  -  (chFETdiss*Rth - prevValue)*exp(-Bdur(TEnum)/tau);
            maxValue = max(maxValue, max(TXthermal(TEnum).ChFETTempRise));
            % this gives us temperature at the end of the current event's transmit
            % burst, based on prevTemp from the previous TXEvent
            prevValue = TXthermal(TEnum).ChFETTempRise*exp(-idleT(TEnum)/tau);
        end
% % %         if maxValue > Lim.ChFETTempRise
% % %             % note which TPC profile triggered overlimit condition
% % %             limProfile = TXEvent(TEnum).profile;
% % %         end
        if maxValue <= oldmaxValue || (maxValue - oldmaxValue)/maxValue < .01
            rising = 0;
        end
        oldmaxValue = maxValue;
    end


    if maxValue > Lim.ChFETTempRise
        HVLim = TPC(limProfiles(1)).highVoltageLimit - clampV;
        newHVLim = HVLim * sqrt(Lim.ChFETTempRise/maxValue); % new HVLim value needed to stay at the limit
        newHVLim = max(floor(newHVLim + clampV), 2.5) - clampV; % round hvLimit down to integer value
        if newHVLim == HVLim
            % if newly reduced HVLim is not lower than HVLim from previous
            % iteration, that means we've already gone as low as we can go
            % by reducing HV so exit with an error condition
            fprintf(2, 'ERROR: TX FET temperature rise exceeds maximum limit due to excessive switching losses.\n');
            fprintf(2, 'TX burst frequency and/or ratio of burst duration to PRI must be reduced.\n');
            error('TXEventCheck: Exiting due to error condition identified above.');
        end
        % Create error string on the first iteration
        if isempty(reduceHVstring)
            reduceHVstring = ['TX FET estimated temperature rise of ',num2str(maxValue,'%.1f'), ' Degrees C. exceeds ', num2str(Lim.ChFETTempRise,'%.1f'), ' Degree maximum limit.'];
        end
        HVLim = newHVLim; % Update HVLim with the new limiting value
        for pnum = 1:length(limProfiles)
            TPC(limProfiles(pnum)).highVoltageLimit = HVLim + clampV;
        end
        iterate = 1; % go back and recalculate FET dissipation at new HVLim
%         disp(['Maximum FET  temp rise ',num2str(maxValue,'%.1f'), ' Deg. for this iteration.']); % uncomment this line to monitor iterations
    end
end % end of the "while iterate" loop

if ~isempty(reduceHVstring) && verbose > 1
    disp(' ')
    disp(reduceHVstring);
    disp(['Reducing TPC Profile hvLimit to ', num2str(TPC(limProfiles(1)).highVoltageLimit, '%.0f'), ' Volts to conform to the limit.']);
    disp(' ')
end





%% Estimate output power, system dissipation, and HV supply current

% these tests only apply if profile 5 is in use
% VTS-2540 bug fix from VTS-2509: TPC(%).inUse could be either 1 or 2
% (Extended Transmit or HIFU) so test for any non-zero value, not just 1.
if TPC(5).inUse
    HVLim = TPC(5).highVoltageLimit - clampV;
    % create a while loop in which we first estimate output power, system
    % dissipation, and HV supply input power and then compare those levels to
    % the push cap and HV supply current limits.  If a limit is exceeded a new
    % HVLim value will be needed, and then another pass through the while loop
    % to re-evaluate the power levels at the new HVLim.


    % Set current-carrying limit for the transmit HV energy storage capacitor
    % provided within the system.  Note that at present the Vantage system only
    % supports one type of capacitor; there is no provision for supporting more
    % than this one type: a 15 mF 100 V. capacitor rated at 30 A. RMS,
    % 90 A. peak (Epcos B41560A9159M000  15 mF +/-20% 100 V. 2.5" dia. 4.1" long)
    Lim.pushCapIpk = 90; % push capacitor peak current limit in Amps
    Cpush = 15e3; % push capacitor is 15 mF, or 15,000 uF (uF units are used here since time is in usec)

    % Set max limit on per-board hv supply current through backplane TPC5 switch
    if P5HVperBdindex
        % limit is programmable, but we always enforce a maximum limit of 90
        % Amps per board
        Lim.BdHvSupplyIout = 90;
        bdLimScale = 1.5; % scale the per-board estimated current by this factor to determine the per-board limit
        bdLimMinOffset = 5; % add this offset to the estimated current, and use the larger of this result or the scaled result
    else
        % Limit is not programmable and is fixed at 40 Amps.  So the max limit
        % allowed by TxEventCheck is set to 36 Amps to avoid error conditions
        % due to component tolerances and inaccuracy in the current estimate
        % when there is no actual fault
        Lim.BdHvSupplyIout = 36;
    end

    % Stet the thresholds (half-cycle period in usec) beyond which the current
    % at the Push Cap and TX bus distribution to individual acq boards needs to
    % be set at the burst average power level (for burst frequencies high
    % enough that the hv decoupling capacitors on the acq boards absorb the
    % load current variation doring a cycle of the transmit waveform), or at
    % the peak power level within the waveform (for frequencies low enough that
    % the current drawn from the push capacitor will follow the waveform).  In
    % between these two thresholds, linear interpolation between the limits
    % will be used.
    hvIpkTh = 5; % 5 usec (100 KHz burst frequency) beyond which peak current must be used
    hvIavgTh = 1; % 1 usec (500 KHz burst frequency) beyond which average current can be used

    % first, compute per-channel power levels over all channels.  Then
    % aggregate into per-board totals and system total, for evaluating the HV
    % supply limits.

    iterate = 1;
    reduceHVstring = []; % if it stays empty there are no reductions to report
    % each time a new limit restriction is found, it will overwrite
    % reduceHVstring so the string actually displayed will be the last (and
    % thus most restrictive) limit found
    while iterate
        iterate = 0;
        % reinitialize everything for each iteration (i.e. start from scratch,
        % but with a new HVLim value)
        HVsupplyIoutMax = 0;
        cumHVsupplyQout = 0;
        BdTotHVPwr = zeros(numTXEvents, numBoards);
        BdHvSupplyIout = zeros(numTXEvents, numBoards);
        for TEnum = 1:numTXEvents
            scalePeriod = min(hvIpkTh, max(hvIavgTh, 1/2*TW(TXEvent(TEnum).tw).estimatedAvgFreq)); % period to use in finding hvSupplyIscale
            % The scaling coefficient hvSupplyIscale ranges from 1 to 2, based
            % on the assumption that the transmit load current waveform is a
            % sinusoid (fundamental term of the Fourier series for trilevel
            % waveform at estimatedAvgFreq).  Thus the average supply current is
            % simply the RMS value of total load current, and the peak current
            % is twice that amount.
            hvSupplyIscale = 1 + (scalePeriod - hvIavgTh)/(hvIpkTh - hvIavgTh);
            ChIsqRMS = HVLim^2 * TXEvent(TEnum).ChIpk1V.*TXEvent(TEnum).ChIpk1V / 2;
            TXEvent(TEnum).ChPwrOut = ChIsqRMS * abs(Zload(TEnum)) * cos(angle(Zload(TEnum)));
            CondLoss = ChIsqRMS * real(Zsource(TEnum));
            ClampLoss = HVLim * TXEvent(TEnum).ChIpk1V * clampV * 2/pi; % average waveform current times fixed clamp drop
            % sum conduction loss, clamp loss, and typical FET switching loss
            % to estimate total HV-driven dissipation in the system.
            ChSysDissHV = CondLoss + ClampLoss  + Cswl(6) * Numpulses(TEnum, :)/(2*Bdur(TEnum));
            for bdnum = 1:numBoards
                chindex = (64*(bdnum-1) + (1:64));
                BdTotHVPwr(TEnum, bdnum) = sum(ChSysDissHV(chindex)) + sum(TXEvent(TEnum).ChPwrOut(chindex));
            end
            BdHvSupplyIout(TEnum, :) = hvSupplyIscale * BdTotHVPwr(TEnum, :)/HVLim; % per-board supply current for this event
            BdHvSupplyIoutMax = max(BdHvSupplyIout(TEnum, :)); % max per-board supply current for this event

            if BdHvSupplyIoutMax > Lim.BdHvSupplyIout % maximum per-board load current
                newHVLim = HVLim * Lim.BdHvSupplyIout/BdHvSupplyIoutMax; % new HVLim value needed to stay at the limit
                newHVLim = max(floor(newHVLim + clampV), 2) - clampV; % round hvLimit down to integer value
                % Create error string
                reduceHVstring = ['Per-board Peak HV supply current of ',num2str(BdHvSupplyIoutMax,'%.0f'), ' Amps exceeds ', ...
                    num2str(Lim.BdHvSupplyIout,'%.0f'), ' Amp maximum limit.'];
                if HVLim ~= newHVLim
                    % force another iteration, and zero the power estimates for this
                    % iteration so following limit checks in this iteration
                    % won't see artifically high values
                    iterate = 1;
                    BdTotHVPwr = zeros(numTXEvents, numBoards);
                    BdHvSupplyIout = zeros(numTXEvents, numBoards);
                    TXEvent(TEnum).ChPwrOut = 0*TXEvent(TEnum).ChPwrOut;
                end
                HVLim = newHVLim; % Update HVLim with the new limiting value
            end

            TXEvent(TEnum).BdHvSupplyIout = BdHvSupplyIout(TEnum, :);
            TXEvent(TEnum).TotalHVPowerIn = sum(BdTotHVPwr(TEnum, :));
            TXEvent(TEnum).TotalPowerOut = sum(TXEvent(TEnum).ChPwrOut);
            TXEvent(TEnum).HVLim = HVLim;
            HVsupplyIoutMax = max(HVsupplyIoutMax, hvSupplyIscale*TXEvent(TEnum).TotalHVPowerIn/HVLim);
            cumHVsupplyQout = cumHVsupplyQout + Bdur(TEnum) * TXEvent(TEnum).TotalHVPowerIn/HVLim; % Amp-usec cumulative output from HV supply
        end
        if iterate
            % skip the rest of this iteration and start over with a new one, so
            % all individual TXEvents will have consistent values from the same
            % HVLim
            continue
        end
        % now check pushcap limit against total HV supply current. For external
        % HIFU, supply current should be zero at start of burst and thus load
        % current will initially all come from the push capacitor, setting the
        % peak current value.  For Extended Transmit systems, the Aux. supply
        % output current is so small that nearly all of the transmit burst
        % current comes from the capacitor.
        if HVsupplyIoutMax > Lim.pushCapIpk % maximum channel peak current
    %         disp(['HVsupplyIoutMax ',num2str(HVsupplyIoutMax,'%.0f'), ' A. for this iteration.']);  % uncomment to display and see if we are iterating
            newHVLim = HVLim * Lim.pushCapIpk/HVsupplyIoutMax; % new HVLim value needed to stay at the limit
            newHVLim = max(floor(newHVLim + clampV), 2) - clampV; % round hvLimit down to integer value
            % Create error string
            reduceHVstring = ['Push Capacitor Peak current of ',num2str(HVsupplyIoutMax,'%.0f'), ...
                ' Amps exceeds ',num2str(Lim.pushCapIpk,'%.0f'), ' Amp maximum limit.'];
            if HVLim ~= newHVLim
                iterate = 1;
                HVLim = newHVLim; % Update HVLim with the new limiting value
                continue % start a new iteration so all per-event values will be updated
            end
        end

        %% ************* Transducer max average power limit
        % If this limit exists in the Trans structure, calculate the total
        % transmit outut power level averaged over the entire event sequence,
        % and compare it to the limit.  Rather than reducing HVLim to stay at
        % the limit, notify the user with an error message and exit if it is
        % exceeded.
        TxPwrOutAvg = 0;
        cumPRI = 0;
        for TEnum = 1:numTXEvents
            TxPwrOutAvg = TxPwrOutAvg + TXEvent(TEnum).TotalPowerOut * TW(TXEvent(TEnum).tw).Bdur;
            cumPRI = cumPRI + TXEvent(TEnum).cumPRI;
        end
        TxPwrOutAvg = TxPwrOutAvg / cumPRI; % Long-term average TX output power in Watts
        TX_Limits.TxPwrOutAvg = TxPwrOutAvg;

        %% ************* PER-BOARD HV DISTRIBUTION DISSIPATION LIMIT ***************
        % Thermal model for dissipation in PCB traces and connector pins, to
        % enforce an RMS current limit of 15 Amps max.
        Lim.BdTxBusTempRise = 35; % 125 degree absolute temperature limit, with 90 degree PCB temperature
        Rth = 3.1; % Thermal resistance of per-board TX bus distribution, in degrees C. per Watt
        Rdiss = 0.05; % effective series resistance of per-board TX bus distribution
        tau = 100e3; % thermal time constant in usec for TX Bus PCB trace
        prevValue = 0; % stores the value found in previous event; initialized to zero
        maxValue = 0.1; % highest value encountered so far
        oldmaxValue = 0; % stores the maxValue from previous pass throgh while loop, to check if we have reached steady state
        rising = 1; % controls when to exit the while loop
        while rising
            for TEnum = 1:numTXEvents
                % find Isquared R dissipation level for this burst
                TxBusdiss = Rdiss * BdTotHVPwr(TEnum, :).*BdTotHVPwr(TEnum, :)/(HVLim^2);
                TXthermal(TEnum).BdTxBusTempRise  = TxBusdiss*Rth  -  (TxBusdiss*Rth - prevValue)*exp(-Bdur(TEnum)/tau);
                maxValue = max(maxValue, max(TXthermal(TEnum).BdTxBusTempRise));
                % this gives us temperature at the end of the current event's transmit
                % burst, based on prevTemp from the previous TXEvent
                prevValue = TXthermal(TEnum).BdTxBusTempRise*exp(-idleT(TEnum)/tau);
            end
            if maxValue <= oldmaxValue || (maxValue - oldmaxValue)/maxValue < .01 || maxValue > Lim.BdTxBusTempRise
                rising = 0;
            end
            oldmaxValue = maxValue;
        end

        if maxValue > Lim.BdTxBusTempRise
            newHVLim = HVLim * sqrt(Lim.BdTxBusTempRise/maxValue); % new HVLim value needed to stay at the limit
            newHVLim = max(floor(newHVLim + clampV), 2) - clampV; % round hvLimit down to integer value
            % Create error string
            reduceHVstring = ['TX Distribution temperature rise of ',num2str(maxValue,'%.1f'), ...
                ' Degrees C. exceeds ', num2str(Lim.BdTxBusTempRise,'%.1f'), ' Degree maximum limit.'];
            if HVLim ~= newHVLim
                % only go back and recalculate only if HVLim is actually changing
                HVLim = newHVLim; % Update HVLim with the new limiting value
                iterate = 1; % go back and recalculate all power levels at new HVLim
            end
        end
    end


    TX_Limits.TXEvent = TXEvent; % save the system state output from TXEventCheck
    TX_Limits.TXthermal = TXthermal; % and the thermal temperature rise results
    assignin('base', 'TX_Limits', TX_Limits);
    if isfield(Trans, 'maxAvgPower') && ~isempty(Trans.maxAvgPower)
        if TxPwrOutAvg > Trans.maxAvgPower
            error('TXEventCHeck: TX average output power of %d Watts exceeds transducer limit of %d Watts.\n', ...
                round(TxPwrOutAvg), round(Trans.maxAvgPower));
        end
    end

    if ~isempty(reduceHVstring) && verbose > 1
        disp(' ')
        disp(reduceHVstring);
        disp(['Reducing hvLimit to ', num2str(HVLim+clampV, '%.0f'), ' Volts to conform to the limit.']);
        disp(' ')
    end

    if P5HVperBdindex
        % reprogram the per-board P5 HV current limits
        PerSlotMaxIout = max(BdHvSupplyIout, [], 1); % maximum over all events for each board
        BdLim = max((bdLimScale * PerSlotMaxIout), (bdLimMinOffset + PerSlotMaxIout));
        % Note empty slots have a limit of 4 Amps (a limit less than 2 Amps may
        % result in a spurious fault condition, due to offset errors and
        % component tolerances in the current monitoring circuitry)
        switch numBoards
            case 1
                result = com.verasonics.hal.tpc.Tpc.setAllAcqBoardsProfile5BusAmperageLimits([BdLim(1), 4, 4, 4]);
            case 2
                result = com.verasonics.hal.tpc.Tpc.setAllAcqBoardsProfile5BusAmperageLimits([BdLim(1), 4, 4, BdLim(2)]);
            case 4
                result = com.verasonics.hal.tpc.Tpc.setAllAcqBoardsProfile5BusAmperageLimits([BdLim(1), BdLim(2), BdLim(3), BdLim(4)]);
        end
        if ~result
        %if ~strcmpi(result, 'Success') && ~strcmpi(result, 'Hardware Not Open')
            % ERROR!  Failed to set per-board.
            fprintf(2, 'TXEventCheck: ERROR!  Failed to set Backplane TPC 5 current limits.');
            error(' ');
        end
    end



    % ********** FIND PUSH CAP DROOP FOR EXTENDED TRANSMIT SYSTEMS ***********
    % Push capacitor will not droop when using external HIFU supply, so
    % don't bother to evaluate the droop at all

    auxHVPsImax = 0.5; % maximum output current in Amps for aux HV supply from TPC


    % First check for total transmit HV supply load (in Amp-usec) over all
    % events for this event sequence.  If that exceeds
    % the aux HV supply output capacity, droop will
    % increase with each pass through event sequence until the HV voltage
    % drops to a sustainable level, with the aux supply operating as a
    % current source.  This is not regarded as an error condition, but we
    % cannot estimate per-event droop and will provide a status message to
    % the user.
    if cumHVsupplyQout > auxHVPsImax*totalPRI
        % The HVLim voltage setting cannot be sustained.
        if verbose > 1
            disp(' ');
            disp(['NOTICE: Estimated average HV transmit supply load of ', num2str(cumHVsupplyQout/totalPRI,'%.1f'),...
                ' Amps exceeds power supply capacity of ', num2str(auxHVPsImax,'%.1f'), ' Amp.']);
            disp('Actual HV transmit level may drop below the selected value.');
        end
    else
        % now that we know droop will reach a stable value, loop through event
        % sequence to find estimated droop for each event.
        prevVdroop = 0;
        maxVdroop = 0;
        oldmaxVdroop = 0;
        rising = 1; % controls when to exit the while loop
        while rising
            for TEnum = 1:numTXEvents
                TXEvent(TEnum).pushcapVdroop = prevVdroop + Bdur(TEnum) * max(0, TXEvent(TEnum).TotalHVPowerIn/HVLim - auxHVPsImax) / Cpush;
                maxVdroop = max(maxVdroop, TXEvent(TEnum).pushcapVdroop);
                prevVdroop = max(0, TXEvent(TEnum).pushcapVdroop - auxHVPsImax*idleT(TEnum)/Cpush);
            end
            if maxVdroop<=oldmaxVdroop
                rising = 0;
            end
            oldmaxVdroop = maxVdroop;
        end
    end
    TPC(5).highVoltageLimit = HVLim + clampV; % update TPC structure

    % Set TPC current limit 25% above estimated current
    supplyupdate = 0;
    newCurrentLimit = 1.25 * HVsupplyIoutMax;
    % don't do the current limit update if TPC(5) not in use (VTS-2509)
    if TPC(5).currentLimit ~= newCurrentLimit
        TPC(5).currentLimit = newCurrentLimit;
        supplyupdate = 1;
    end
end

TX_Limits.TXEvent = TXEvent; % save the system state output from TXEventCheck


%% Save results, update HV supply if needed

% send new current limit value to the base workspace so the call to
% external power supply control function below will be able to see it
assignin('base','TPC', TPC);
if TPC(5).inUse
    % check and adjust actual power supply settings, but only if HW is
    % present
    hv = TPC(5).hv; % get the existing high voltage setting
    if hv>TPC(5).highVoltageLimit || supplyupdate % if current limit changed or hv is greater than new limit, send commands to update hv
        hv = min(TPC(5).highVoltageLimit, hv);
        if evalin('base', 'VDAS') && TPC(5).inUse == 2
            result = setTpcProfileHighVoltage(hv, 5, isVsxInit); % VTS-1796 pass VSX initialization status to setTpcProfileHighVoltage
            if ~strcmpi(result, 'Success')
                fprintf(2, 'VSX: Communication error with external power supply.');
                error(' ');
            end
        end
        % note these set commands don't generate errors if the tag is not
        % found:
        set(findobj('Tag','hv2Sldr'),'Value',hv);
        set(findobj('Tag','hv2Value'),'String',num2str(hv,'%.1f'));
        if  voltageTrackP5 > 0
            % if voltage tracking enabled, also set the profile that is to track P5
            TPC(voltageTrackP5).hv = hv;
            if voltageTrackP5 == 1
                % HV slider needs to be updated to show new voltage setting
                set(findobj('Tag','hv1Sldr'),'Value',hv);
                set(findobj('Tag','hv1Value'),'String',num2str(hv,'%.1f'));
            end
        end
        TPC(5).hv = hv;
    end
end
% save results variables back to the workspace
assignin('base','TPC', TPC);
assignin('base','TX',TX); % save the stuff added in TX structure
% NOTE: make sure each of the arrays listed below have been added into
% TX_Limits after they have been updated; if that is the case we don't need
% to do it again here.
% % TX_Limits.Event = Event; % save stuff added to Event structure
% % TX_Limits.TXEvent = TXEvent; % save the system state output from TXEventCheck
% % TX_Limits.TXthermal = TXthermal; % and the thermal temperature rise results
assignin('base','TX_Limits', TX_Limits); % push it all back to the base workspace

% For TXEventCheck we are now finished, except for one last thing: create
% update&Run command for TPC, if any highVoltageLimit values have changed
updateTPC = 0; % flag to track whether TPC update is needed
for i = 1:5
    if TPC(i).inUse
        % note unused profiles may have empty values and we will leave them
        % that way.
        if PreviousHvLimit(i) ~= TPC(i).highVoltageLimit
            hv = min(TPC(i).highVoltageLimit, TPC(i).hv);
            TPC(i).hv = hv; % update to new value
            % a limit has been changed, so we need to do update%Run to TPC.
            % Note this will also make TPC change hv, if we have reduced it
            % to stay under new highVoltageLimit
            updateTPC = 1; % create update&Run command below
        end
    end
end

if updateTPC && evalin('base', 'exist(''VSX_Control'', ''var'')')
    % if VSX_Control doesn't exist this was an initialization
    % call to TXEventCheck and we don't need to do anything
    % since initialization will be loading everything
    VSX_Control = evalin('base', 'VSX_Control');
    n=length(VSX_Control)+1;
    VSX_Control(n).Command = 'update&Run';
    VSX_Control(n).Parameters = {'TPC'};
    assignin('base', 'VSX_Control', VSX_Control);
    assignin('base', 'TPC', TPC);
end

return

