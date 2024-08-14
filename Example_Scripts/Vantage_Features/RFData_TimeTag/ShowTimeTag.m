% ShowTimeTag

% Read RcvData buffer 1, and extract the time tag value from the first two
% samples from channel 1 in each frame (note this would be the time tag
% from only the first acquisition of a frame that included multiple
% acquisitions).  Create an array of time stamp values in seconds for all
% frames in the receive buffer

numFrames = size(RcvData{1}, 3);  % number of frames in Rcv buffer 1 in Matlab Workspace

% create a column vector for time tag value from first acquisition of each
% frame
TimeTagValues = zeros(numFrames, 1);

W = zeros(2, 1); % W will be the two time tag words after conversion to unsigned integer values
for n = 1:numFrames
    for i=1:2
        % read first and second sample from column 1 of the frame, and
        % convert to double
        W(i) = double(RcvData{1}(i, 1, n));
        if W(i) < 0
            % translate 2's complement negative values to their unsigned integer
            % equivalents
            W(i) = W(i) + 65536;
        end
    end
    % first sample is 16 LSB's of time tag value and second sample is 16
    % MSb's so scale and sum to get actual time tag value
    TimeTagValues(n, 1) = W(1) + 65536 * W(2);
end

% the 32 bit time tag counter increments every 25 usec, so we have to scale
% by 25 * 1e-6 to convert to a value in seconds

TimeTagValues = TimeTagValues/4e4;
% This array will be left in workspace so you can plot it, list the values,
% etc.
