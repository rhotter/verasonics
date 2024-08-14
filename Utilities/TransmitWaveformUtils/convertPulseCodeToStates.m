function [ States, result ] = convertPulseCodeToStates( PulseCode)
% Translate a legacy PulseCode array into its States array equivalent
%
% The input Pulse code array can be either a single waveform (an array of
% size (numPCrows X 5), or an array of waveforms (an array of size
% numPCrows X 5 X numCh).  The output States array will be of size
% (numStRows X 2 X numCh), where numStRows will be the length of longest
% States waveform array.  If there are channels with a shorter States array
% than numStRows, they will be terminated with an end command and the
% remaining rows will be filled with additional end commands.
%
% The return argument "result" will be a string set to "Success" if the
% conversion was completed with no errors, or a string identifying the
% problem if the input PulseCode array was not valid.  In that case, the
% returned States array will be empty.
%
% Created Mar. 30, 2020 for release 4.3.0

% initialize output variables
States = []; % this will be returned if there was an error
result = 'Success'; % this will be overwritten if there is an error

if size(PulseCode, 2) ~= 5
    result = 'PulseCode array does not have five columns';
    return
end

numCh = size(PulseCode, 3); % number of channels to be converted
numPCrows = size(PulseCode, 1); % number of rows in PulseCode array to be converted

% Create StatesOut array with initial length equal to numPCrows and fill it
% with waveform end commands; length will be increased as needed while each
% channel is converted.  If conversion is completed successfully, StatesOut
% will be copied into the returned States Array.
numStRows = numPCrows; % Initial length of StatesOut
StatesOut = zeros(numStRows, 2, numCh); % Create with all zeros
StatesOut(:, 1, :) = 30; % turn every row into an end command


for chnum = 1:numCh
    % Convert one channel from PulseCode to States
    stRow = 0; % last completed row in States array for this channel
    % Convert one row at a time from the PulseCode array
    for pcRow = 1:numPCrows
        % Check value of repeat count for this row
        if PulseCode(pcRow, 5, chnum) == 0
            % This is PulseCode waveform end; we are done with this
            % channel.  Ignore the first four entries of this PulseCode row
            % and break out of pcRow for-loop
            break
        elseif PulseCode(pcRow, 5, chnum) < 0
            result = 'Unrecognized Repeat count in PulseCode array';
            return
        elseif PulseCode(pcRow, 5, chnum) > 1
            % This PulseCode row has a repeat count, so add loopStart
            % command prior to the states for this row
            stRow = stRow + 1;
            if stRow > numStRows
                % Add one more row to StatesOut array and fill it with
                % WaveformEnd commands for all channels
                numStRows = numStRows + 1;
                StatesOut(numStRows, :, :) = 0;
                StatesOut(numStRows, 1, :) = 30; % end command in added row for every channel
            end
            % Add the loopStart command with repeat count
            StatesOut(stRow, :, chnum) = [10 PulseCode(pcRow, 5, chnum)];
        end
        % now ready to process all four states from this row
        for stcol = 1:4
            if PulseCode(pcRow, stcol, chnum) == 0
                % zero duration so skip it
                continue
            end
            if stcol == 1 || stcol == 3
                % this is a zero state
                lvl = 0;
            else
                % this is an active pulse
                lvl = sign(PulseCode(pcRow, stcol, chnum));
            end
            dur = abs(PulseCode(pcRow, stcol, chnum));
            stRow = stRow + 1;
            if stRow > numStRows
                % Add one more row to StatesOut array and fill it with
                % WaveformEnd commands for all channels
                numStRows = numStRows + 1;
                StatesOut(numStRows, :, :) = 0;
                StatesOut(numStRows, 1, :) = 30; % end command in added row for every channel
            end
            % Add the new state to StatesOut
            StatesOut(stRow, :, chnum) = [lvl dur];
        end
        % Done processing the state entries for this row from PulseCode
        if PulseCode(pcRow, 5, chnum) > 1
            % this was a repeating row so add the loopEnd command
            stRow = stRow + 1;
            if stRow > numStRows
                % Add one more row to StatesOut array and fill it with
                % WaveformEnd commands for all channels
                numStRows = numStRows + 1;
                StatesOut(numStRows, :, :) = 0;
                StatesOut(numStRows, 1, :) = 30; % end command in added row for every channel
            end
            % Add the loopEnd command
            StatesOut(stRow, :, chnum) = [20 1];
        end
    end % end of the loop for converting all PulseCode rows for this channel
    if stRow == numStRows
        % New States array for this channel has filled entire StatesOut
        % array length, so add one more row to StatesOut for the
        % WaveformEnd command
        numStRows = numStRows + 1;
        StatesOut(numStRows, :, :) = 0;
        StatesOut(numStRows, 1, :) = 30; % end command in added row for every channel
    end
end % end of the loop for converting all channels of PulseCode into StatesOut

% if we reached this point, the entire conversion has been successful. Copy
% StatesOut into States and return
States = StatesOut;
end













% 
% for Chnum = 1:numCh
%     stRow = 0; % last completed row in States array
%     for rownum = 1:numRows
%         rpt = PulseCode(rownum, 5, Chnum); % repeat count from PulseCode
%         if rpt <= 0
%             % waveform termination; the first four entries of this row and
%             % all subsequent PulseCode rows are ignored so break out of
%             % loop
%             break
%         end
%         rpt2 = 0; % to keep track of whether a nested repeat was used
%         if rpt > 1
%             lplvl = 1;
%             % we need to add a repeat command for this segment
%             if rpt > 4096
%                 lplvl = 2;
%                 % create a nested repeat for large repeat count
%                 exp = ceil(log2(rpt)) - 12;
%                 if exp > 12
%                     error('PulseCode2WvDef: Maximum segment nested repeat count of 2^24 has been exceeded.');
%                 end
%                 rpte = 2^exp;
%                 rpt = round(rpt/rpte)-1;
%                 rpt2 = rpte;
%                 % nested repeat for this segment
%                 stRow = stRow + 1;
%                 States(stRow, :, Chnum) = [ 10 rpt2 ];
%                 RSeg = zeros(1, 2, numCh); % this will be the duplicate copy of repeating segment
%                 rSegRow = 0; % pointer to last row used
%             end
%             stRow = stRow + 1;
%             States(stRow, :, Chnum) = [ 10, rpt ];
%             % 10 identifies the repeat command; second entry is nested loop
%             % level
%         end
%         % repeat comand is complete; now add up to four States rows for the
%         % four states in this row of PulseCode array.  Only those entries
%         % with nonzero duration are copied into the States array:
%         for i = 1:4
%             % step through the four entries in this row
%             if PulseCode(rownum, i, Chnum) ~= 0 % skip this entry if duration is zero
%                 if i == 1 || i == 3
%                     lvl = 0;
%                 else
%                     lvl = sign(PulseCode(rownum, i, Chnum));
%                 end
%                 dur = abs(PulseCode(rownum, i, Chnum));
%                 if stRow > 0 && States(stRow, 1, Chnum) == lvl
%                     % level matches current state so add to the duration
%                     States(stRow, 2, Chnum) = States(stRow, 2, Chnum) + dur;
%                     if rpt2 > 1
%                         RSeg(rSegRow, 2, Chnum) = RSeg(rSegRow, 2, Chnum) + dur;
%                     end
%                 else
%                     stRow = stRow + 1;
%                     States(stRow, :, Chnum) = [lvl dur];
%                     if rpt2 > 1
%                         rSegRow = rSegRow + 1;
%                         RSeg(rSegRow, :, Chnum) = [lvl dur];
%                     end
%                 end
%             end
%         end
%         if rpt > 1
%             if rpt2 > 1
%                 % insert end of nested repeat loop command
%                 stRow = stRow + 1;
%                 States(stRow, :, Chnum) = [ 20 2 ];
%                 % add the repeated segment for nested loop
%                 States = vertcat(States, RSeg);
%                 stRow = stRow + rSegRow;
%             end
%             % insert end of repeat loop command
%             stRow = stRow + 1;
%             States(stRow, :, Chnum) = [ 20 1 ];
%         end
%     end
% 
%     % Add [0 0] at the end of States
%     if isequal(size(States,1),stRow) && ~isequal(States(stRow,:,Chnum),[0 0])
%         States(stRow+1,:,Chnum) = [0 0];
%     end
% 
%     % check for null waveform, and replace with end command
%     if isequal(States(stRow+1, :, Chnum), [0 0])
%         States(stRow+1, :, Chnum) = [30 0]; % waveform end command
%     end
% end
% 
% 
% end

