function [ States, result ] = convertTriLvlToStates( TriLvlWvfm )
% Convert a Tri-level waveform vector array into a States array equivalent
%
% The input TriLvlWvfm array must be made up of 1-D Matlab column vectors
% for each individual waveform, either N X 1 for a single waveform or N X
% numCH for an array of per-channel waveforms, where N is the number of
% samples (at the 250 MHz system clock rate) in the waveforms. The sample
% values in each waveform represent the output level for that sample,
% restricted to +1, 0, or -1).  In a multichannel TriLvlWvfm array where
% all waveforms are not the same length, the shorter ones should be padded
% with zero samples at the end to match the length of the longest waveform.
%
% The return argument "States" will be a representation of the input
% waveform in the States array format.  No attempt will be made to dectect
% repeating patterns and turn them into a loop in the States array.  If
% more than one channel waveform is provided in the input TRriLvlWvfm
% array, the States ouput array will have a third dimension to match the
% number of waveforms provided in the input TriLvlWvfm.  If a waveform ends
% with a zero state of duration greater than one, that final state will be
% truncated to a duration of one clock.
%
% The return argument "result" will be a string set to "Success" if the
% conversion was completed with no errors, or a string identifying the
% problem if the input TriLvlWvfm array was not valid.  In that case, the
% returned States array will be empty.
%
% Created April 14, 2020 for release 4.3.0

% initialize output variables
States = []; % this will be returned if there was an error
result = 'Success'; % this will be overwritten if there is an error

if isempty(TriLvlWvfm)
    result = 'TriLvlWvfm input array is empty';
    return
else
    numSamples = size(TriLvlWvfm, 1); % number of samples in the longest dimension
    numCH = size(TriLvlWvfm, 2); % number of waveforms to be converted
end


Tlvl = sign(TriLvlWvfm); % convert all entries to +1, 0, or -1
if ~isequal(Tlvl, TriLvlWvfm)
    result = 'TriLvlWvfm input array has values that are not +1, 0, or -1';
    return
end

for chnum = 1:numCH
    % convert one channel, into single-channel StatesCH array
    StatesCH = [0 0]; % initialize StatesCH to N X 2 dimensions
    stRow = 0; % last completed row in States array
    currentLvl = -2; % tracks level we are currently working on

    for sampnum = 1:numSamples
        % append each entry in TriLvlWvfm to StatesCH
        if Tlvl(sampnum, chnum) ~= currentLvl
            % this is a new level, so start a new state
            currentLvl = Tlvl(sampnum, chnum);
            stRow = stRow + 1;
            StatesCH(stRow, :) = [currentLvl  1 ]; % new state with duration 1
        else
            % this is another sample at same level, so increment current state
            % duration
            StatesCH(stRow, 2) = StatesCH(stRow, 2) + 1;
        end
    end
    % end of the loop for converting one channel from TriLvlWvfm into StatesCH
    if StatesCH(stRow, 1) == 0
        % If last state is at zero level, clip its duration to one clock
        StatesCH(stRow, 2) = 1;
    end
    StatesCH = [StatesCH; [30 0]]; % add end command at the end of the waveform
    stRow = stRow + 1; % increment to match overall length with end command

    % Now copy the completed StatesCH waveform into corresponding channel of
    % output States array, after padding States array as needed to match
    % length of the new array
    if chnum == 1
        % This is the first waveform; create States output array of correct
        % size and fill it with end commands
        States = zeros(stRow, 2, numCH);
        States(:, 1, :) = 30; % All entries are now waveformEnd commands
        maxRows = stRow; % to track longest States array created
    elseif stRow > maxRows
        % new StatesCH waveform is longer than States array, so Pad it on
        % all channels with more end commands
        Padstates = zeros(stRow-maxRows, 2, numCH);
        Padstates(:, 1, :) = 30; % All entries are now waveformEnd commands
        States = [States; Padstates];
        maxRows = stRow; % this is longest waveform so far
    end
    States(1:stRow, :, chnum) = StatesCH;
end

end

