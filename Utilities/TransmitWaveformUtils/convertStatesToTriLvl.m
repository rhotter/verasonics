function [TriLvlWvfm, Result] = convertStatesToTriLvl(States, repeatLim)
% convertStatesToTriLvl: utility function to create a TriLvl waveform
% vector that is the equivalent of the States array input argument.
%
%   The second input argument "repeatLim" is optional and can be used to
%   limit repeat loops to no more than the repeat count specified by
%   repeatLim; a repeatLim value of 0 means leave the loop counts
%   unchanged.  If repeatLim is not specified, a default value of 20 will be
%   assigned.
%
%   The output argument TriLvlWvfm will be a column vector of waveform
%   samples at the system clock rate, with each sample having a value of 1,
%   0, or -1. If input States array has multiple channels, the output array
%   will have a separate column for each channel; additional zeros will be
%   added at end of the channel's tri level waveform if it is shorter than
%   the longest one.
%   Output argument "Result" will be a string set to "Success" if the
%   conversion was completed successfullly or a string identifying the
%   problem if an error occurred.  If there was an error, TriLvlWvfm will
%   be set to empty.

% Revision History
% April 4, 2020 initial version for 4.3.0 release (VTS-1365)

% initialize return values
TriLvlWvfm = [];
Result = 'Error: Unidentified Fault Condition';

if nargin == 1
    repeatLim = 20; % default value if not specified
else
    repeatLim = round(repeatLim);
    if repeatLim < 0
        Result = 'Error: repeatLim must be a non-negative integer';
        return
    end  
end

% confirm correct format for States input array
if size(States, 2) ~= 2
    Result = 'Error: States array must have two columns';
    return
elseif size(States, 3) ~= 1
    Result = 'Error: Multiple States input arrays not supported';
    return
end

numCH = size(States, 3); % number of channels in the input waveform

States = round(States); % all values must be integers

for chnum = 1:numCH
    StatesCh = States(:, :, chnum); % States array for this channel
    loopLevel = 4; % max supported loop level
    ChTriLvlWvfm = []; % tri level waveform to be created for this channel
    while loopLevel > 0
        % This while loop will find all commands in the StatesCh array, and
        % unwrap the highest loop level that is found.  After unwrapping the
        % repeating waveform segment, the associated loop start and loop end
        % commands are deleted.  The while loop will repeat that process until
        % all commands are deleted and we have a "flat" StatesCh array.
        foundWvEnd = 0; % flag to indicate we've processed waveform end command
        % find all the commands in states array
        CmdAdr = find(abs(StatesCh(:, 1)) > 1);
        if isempty(CmdAdr)
            % if there are no commands left to process, we are done and ready
            % to convert flattened states array to trilevel.  Set loopLevel to
            % zero and continue, to break out of the while loop
            loopLevel = 0;
            continue
        end
        CmdID = abs(StatesCh(CmdAdr, 1));
        CmdArg = StatesCh(CmdAdr, 2);
        CmdLoopLvl = CmdArg;
        % step through all commands to check for correct order and loop level
        SequenceLpLvl = 0; % always start at zero
        for cmdnum = 1:length(CmdAdr)
            if CmdID(cmdnum) == 30
                % Waveform end command; waveform terminates here so delete the
                % rest of the waveform unless we are in the middle of a loop
                % (an error condition). Note that the StatesCh array could
                % contain loops that are beyond the waveform end command, so we
                % check for and process waveform end first.  It will delete all
                % rows of the StatesCh array beyond the Waveform end command.
                if SequenceLpLvl > 0
                    Result = 'Error: StatesCh has waveform end command within a loop';
                    TriLvlWvfm = [];
                    return
                else
                    % valid waveform end; truncate states waveform here and
                    % continue (delete the waveform end command and everything
                    % beyond it)
                    stsEnd = CmdAdr(cmdnum) - 1;
                    if stsEnd == 0
                        % this is a null waveform so create States
                        % array with one-clock zero State
                        StatesCh = [0 1];
                    else
                        newStates = StatesCh(1:stsEnd, :);
                        StatesCh = newStates;
                    end
                    foundWvEnd = 1;
                    break % break out of cmdnum loop
                end
            elseif CmdID(cmdnum) == 10
                % loop start command; increment loop level and place it in
                % CmdLoopLvl but first check for valid loop count > 1
                if CmdArg(cmdnum) < 2
                    Result = 'Error: Loop start command with loop count < 2';
                    TriLvlWvfm = [];
                    return
                end
                SequenceLpLvl = SequenceLpLvl + 1;
                CmdLoopLvl(cmdnum) = SequenceLpLvl;
            elseif CmdID(cmdnum) == 20
                % loop end command; check loop level and report error if
                % mismatch otherwise proceed
                if CmdLoopLvl(cmdnum) ~= SequenceLpLvl
                    Result = 'Error: StatesCh has loop end command without matching loop start';
                    TriLvlWvfm = [];
                    return
                end
                SequenceLpLvl = SequenceLpLvl - 1;
            else
                % unrecognized command
                Result = 'Error: Unrecognized Command ID- not 10, 20 or 30';
                TriLvlWvfm = [];
                return
            end
        end
        if foundWvEnd
            % jump back to start of while loop to reprocess StatesCh array after
            % truncating at waveform end
            continue
        end
        % ready to process a command; first set loop level down to max that we
        % actually have.  Note that if no commands were found, CmdLoopLvl will
        % be empty and thus loopLevel will be set to empty as well- but this
        % was alreay detected above, and triggered an exit from the while loop
        loopLevel = max(CmdLoopLvl);

        % find the first command at that loop level; because we already
        % confirmed command order, it will be loop start
        lpStCmdindx = find( CmdLoopLvl == loopLevel, 1);
        if CmdID(lpStCmdindx) ~= 10
            % If this is not a loop start, there may be a bug in decoding logic
            % above
            Result = 'Error: Potential logic bug- expected loop start command is something else';
            TriLvlWvfm = [];
            return
        end
        preSegEnd = CmdAdr(lpStCmdindx) - 1;
        lpSegStrt = CmdAdr(lpStCmdindx) + 1; % start of repeating segment
        % next command will be loop end
        if CmdID(lpStCmdindx+1) ~= 20
            % If this is not a loop send, there may be a bug in decoding logic
            % above
            Result = 'Error: Potential logic bug- expected loop end command is something else';
            TriLvlWvfm = [];
            return
        end
        lpSegEnd = CmdAdr(lpStCmdindx+1) - 1; % end of repeating segment
        postSegStrt = CmdAdr(lpStCmdindx+1) + 1; % start of post segment
        postSegEnd = size(StatesCh, 1); % end of post segment
        % make sure repeating segment is not empty
        if lpSegEnd < lpSegStrt
            Result = 'Error: Repeat Loop commands with no waveform segment to repeat';
            TriLvlWvfm = [];
            return
        end
        % find the loop count to be used
        if repeatLim == 0
            repeatCount = CmdArg(lpStCmdindx);
        else
            repeatCount = min(repeatLim, CmdArg(lpStCmdindx));
        end
        % Create the repeated segment
        RepeatSeg = repmat(StatesCh(lpSegStrt:lpSegEnd, :), repeatCount, 1);
        % rebuild the StatesCh array with repeated segment and associated loop
        % commands removed
        newStates = [StatesCh(1:preSegEnd, :); RepeatSeg; StatesCh(postSegStrt:postSegEnd, :)];
        StatesCh = newStates;
    end % end of the while loop to process commands and unwrap loops


    % All done processing commands, and can now convert StatesCh
    % to TriLevel
    if max(StatesCh(:, 1)) > 1 || min(StatesCh(:, 1)) < -1
        % this might mean there is a bug in the command processing above...
        Result = 'Error: Potential logic bug- flattened StatesCh array has illegal value';
        return
    end

    for rownum = 1:size(StatesCh, 1)
        newState = StatesCh(rownum, 1) * ones(StatesCh(rownum, 2), 1);
        ChTriLvlWvfm = [ChTriLvlWvfm; newState];
    end
    
    chLgth = length(ChTriLvlWvfm); % length of waveform for this channel
    if chnum == 1
        % create output array for all channels, filled with zeros
        TriLvlWvfm = zeros(chLgth, numCH); % Initialize output to null waveform for all channels
        maxChLgth = chLgth; % will track max length over all channels
    elseif chLgth > maxChLgth
        % Need to pad output array with more zeros since new waveform is
        % longer than previous max
        TriLvlWvfm = [TriLvlWvfm; zeros(chLgth - maxChLgth, numCH)];
        maxChLgth = chLgth; % we have a new max length
    end
    
    % now paste in the new channel waveform
    TriLvlWvfm(1:chLgth, chnum) = ChTriLvlWvfm;

end % end of for loop over all channels


Result = 'Success';
    
end

