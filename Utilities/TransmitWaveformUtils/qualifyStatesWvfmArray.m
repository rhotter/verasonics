function [StatesOut, result] = qualifyStatesWvfmArray(StatesIn)
% qualifyStatesWvfmArray: utility function to consolidate unused or
% redundant entries in a States array and check for errors or unreognized
% values
%   The input argument StatesIn can be any States array, either for a
%   single waveform (size N X 2) or an array of per-channel waveforms (size
%   N X 2 X numCH).
%
%   The return argument "result" will be a string set to "Success" if the
%   qualify processing was completed with no errors, or a string
%   identifying the problem if the input StatesIn array is not valid.
%   In that case, the returned StatesOut array will be empty.
%
% The following conditions will be checked for in this function:
%   * ensure all entries are integers
%   * report an error if any value is invalid or unrecognized
%   * find contiguous states at the same level, and combine them into a
%       single state with the same total duration
%   * find states of zero duration and eliminate them
%   * check for proper nesting and termination of loopStart - loopEnd
%       commands
%   * check for repeating loop segment duration of at least 2 clocks
%   * check for presence of waveformEnd command, and add one if needed
%
% Note that the returned StatesOut array may be shorter than StatesIn, if
% any consolidation was done.  If StatesIn has multiple waveforms,
% StatesOut will match the length of the longest individual output array.
% If other arrays are shorter, they will be padded at the end with
% additional waveformEnd commands.

% Created April 12, 2020 for release 4.3.0 (VTS-1365)

% initialize output variables
StatesOut = []; % this will be returned if there was an error
result = 'Success'; % this will be overwritten if there is an error

if isempty(StatesIn)
    result = 'StatesIn input array is empty';
    return
else
    % check size of input array
    [numRowsIn, numcol, numCH] = size(StatesIn);
end
if numcol ~= 2
    result = 'StatesIn input array does not have two columns';
    return
end
StatesIn = round(StatesIn); % make sure all input values are integers

if any(StatesIn(:, 2, :) < 0)
    result = 'illegal negative value found in column 2 of StatesIn input array.';
    return
end

Sout = zeros(numRowsIn, 2, numCH); % create placeholder output array same size as StatesIn
Sout(:, 1, :) = 30; % fill Sout with waveformEnd commands

numRowsOut = 1; % To keep track of overall length as each CH is prcessed

for chnum = 1:numCH
    % process one waveform from StatesIn
    looplvl = 0; % to track loop level throughout the array
    stlvl = 7; % track most recent state level; start with dummy value so first state will be a new level
    abscmds = abs(StatesIn(:, 1, chnum)); % absolute value of command column
    outrow = 0; % tracks last output row created for this channel
    rptDur = zeros(1, 4); % tracks duration of repeat segments
    wvend = 0; % indicates we never found an end command
    for rownum = 1:numRowsIn
        % process all rows for this channel's input array
        if abscmds(rownum) < 2
            % this is a state, not a command
            if StatesIn(rownum, 2, chnum) == 0
                % this state has zero duration so just ignore it and go to
                % next row
                continue
            else
                % state with nonzero duration
                if looplvl > 0
                    % in a repeating segment so add to rptDur; if this is a
                    % nested loop the added duration counts for lower level
                    % loops as well
                    rptDur(1:looplvl) = rptDur(1:looplvl) + StatesIn(rownum, 2, chnum);
                end
                if StatesIn(rownum, 1, chnum) == stlvl
                    % state is at same level as previous state so add
                    % duration on to it
                    Sout(outrow, 2, chnum) = Sout(outrow, 2, chnum) + StatesIn(rownum, 2, chnum);
                else
                    % new state at a new level so add a new row
                    outrow = outrow + 1;
                    Sout(outrow, :, chnum) = StatesIn(rownum, :, chnum); % copy over the level and duration
                    stlvl = StatesIn(rownum, 1, chnum); % this is current level for checking next row from StatesIn
                end
            end
        elseif abscmds(rownum) == 10
            % loopStart command
            if StatesIn(rownum, 2, chnum) < 2
                % illegal repeat count value
                result = ['loopStart with repeat count less than 2 found in channel ', num2str(chnum), ' of StatesIn input array.'];
                return
            else
                % legitimate loopStart so copy it over and update looplvl,
                % rptDur, and stlvl
                outrow = outrow + 1;
                Sout(outrow, :, chnum) = StatesIn(rownum, :, chnum); % copy over the whle row
                looplvl = looplvl + 1; % entering new loop level
                rptDur(looplvl) = 0; % reset loop segment duration to zero for this new loop
                stlvl = 7; % next state row has to be copied over- not a continuation
            end
        elseif abscmds(rownum) == 20
            % loopEnd command
            if StatesIn(rownum, 2, chnum) ~= looplvl
                % illegal loop level
                result = ['loopEnd without matching loopStart in channel ', num2str(chnum), ' of StatesIn input array.'];
                return
            elseif rptDur(looplvl) < 2
                % inadequate repeating segment duration
                result = ['repeating loop segment has duration less than 2 in channel ', num2str(chnum), ' of StatesIn input array.'];
                return
            else
                % legitimate loopEnd so copy it over and update looplvl,
                % rptDur, and stlvl
                outrow = outrow + 1;
                Sout(outrow, :, chnum) = StatesIn(rownum, :, chnum); % copy over the whle row
                rptDur(looplvl) = 0; % reset loop segment duration to zero since we are ending this loop
                looplvl = looplvl - 1; % going down to lower loop level
                stlvl = 7; % next state row has to be copied over- not a continuation
            end
        elseif abscmds(rownum) == 30
            % waveformEnd command
            if StatesIn(rownum, 2, chnum) ~= 0
                % illegal waveformEnd argument
                result = ['waveformEnd command with nonzero argument in channel ', num2str(chnum), ' of StatesIn input array.'];
                return
            elseif looplvl > 0 || max(rptDur) > 0
                % waveform end with unfinished loop
                result = ['waveform end with unterminated repeat loop in channel ', num2str(chnum), ' of StatesIn input array.'];
                return
            else
                % legitimate loopEnd; don't need to copy it over since Sout
                % was filled with loopEnd commands initially, but do need
                % to increment outrow
                outrow = outrow + 1;
                wvend = 1; % indicates we actually found an end command
                break % ignore the rest of this waveform since we hit the end command
            end
        else
            % unrecognized command value
            result = ['unrecognized command value found in channel ', num2str(chnum), ' of StatesIn input array.'];
            return
        end
    end
    if ~wvend
        % never encountered an end command so we need to add one
        if looplvl > 0 || max(rptDur) > 0
            % waveform end with unfinished loop
            result = ['waveform end with unterminated repeat loop in channel ', num2str(chnum), ' of StatesIn input array.'];
            return
        else
            % legitimate waveform end so increment outrow
            outrow = outrow + 1;
            if outrow > numRowsIn
                % StatesIn wasn't big enough so add another row of end
                % commands to Sout
                endrow = zeros(1, 2, numCH);
                endrow(:, 2, :) = 30; % fill it with end commands
                Sout = [Sout; endrow];
            end
        end
    end
    % done with this waveform; update numRowsOut
    numRowsOut = max(numRowsOut, outrow);
end

% done processing all waveforms; copy active rows from Sout into StatesOut
StatesOut = Sout(1:numRowsOut, :, :);
end

