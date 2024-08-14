% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: SetUpL11_5vFlashHFR_saveRFperFrame.m "High Frame Rate" acquisition
% for ultrafast imaging with saving RF data per frame in realtime. This script
% shows how to use external function to save the RF data of each frame in
% realtime. Here, only non-zeros channel data will be saved to reduce the saving
% time. The hardware and software are synchronized with a sync command after
% data acquisition of each 'super' frame to prevent the data from being
% overwritten. Therefore, please note that the interval between the last
% sub-frame of the current super frame and the first sub-frame of the next super
% frame will be longer than the timeToNextAcq defined in the SeqControl(2). In
% addition, the default frame limit is set to 100 "superframes". Please change
% to a bigger number for more superframes.
%
% To support higher acquisition frame rates with reduced DMA overhead, this
% script acquires a large set of T/R acquisitions into each RcvBuffer 'super'
% frame, performing a transferToHost only after each group of "numAcqs"
% acquisitions.
%
% Last update:
% 12/28/2018 - test with SW 4.0

clear all
filename = 'L11-5vFlashHFR_saveRFperFrame'; % used to launch VSX automatically

% --- Frequently Modified Parameters ------------------------------------------
P.startDepth = 5;   % Acquisition depth in wavelengths
P.endDepth = 192;   % This should preferrably be a multiple of 128 samples.

P.numAcqs = 100;      % no. of Acquisitions in a Receive frame (this is a "superframe")
P.numFrames = 4;      % no. of Receive frames (real-time images are produced 1 per frame)
simulateMode = 0;   % set to acquire data using Vantage 128 hardware

% -----------------------------------------------------------------------------

% Define system parameters.
Resource.Parameters.numTransmit = 128;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;   % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = simulateMode;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Specify Trans structure array.
Trans.name = 'L11-5v';
Trans.units = 'wavelengths'; % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);    % L11-5v transducer is 'known' transducer so we can use computeTrans.
Trans.maxHighVoltage = 50;      % set maximum high voltage limit for pulser supply.

% Specify PData structure array.
PData(1).PDelta = [Trans.spacing, 0, 0.5];
PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3)); % startDepth, endDepth and pdelta set PData(1).Size.
PData(1).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(3) = 1;      % single image page
PData(1).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of upper lft crnr.
% No PData.Region specified, so a default Region for the entire PData array will be created by computeRegions.

% Specify Media object. 'pt1.m' script defines array of point targets.
pt1;
Media.attenuation = -0.5;
Media.function = 'movePoints';

% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 4096*P.numAcqs;   % this size allows for maximum range
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = P.numFrames;       % number of 'super frames'

% Specify Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,.67,2,1];

% Specify TX structure array.
TX(1).waveform = 1;            % use 1st TW structure.
TX(1).Origin = [0.0,0.0,0.0];  % flash transmit origin at (0,0,0).
TX(1).focus = 0;
TX(1).Steer = [0.0,0.0];       % theta, alpha = 0.
TX(1).Apod = ones(1,Trans.numelements);
TX(1).Delay = computeTXDelays(TX(1));

% Specify TGC Waveform structure.
TGC.CntrlPts = [0,298,395,489,618,727,921,1023];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Receive structure arrays -
%   endDepth - add additional acquisition depth to account for some channels
%              having longer path lengths.
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength,...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 1), 1,P.numFrames*P.numAcqs);  % movepoints EVERY acquisition to illustrate superframe concept
                                                                    % real-time images will look "jerky" but using the reconstructAll script,
                                                                    % playback process all acquisitions and shows smooth displacement

% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
%     Receive(P.numAcqs*(i-1) + 1).callMediaFunc = 1;  % move points only once per frame
    for j = 1:P.numAcqs
        % -- Acquisitions for 'super' frame.
        rcvNum = P.numAcqs*(i-1) + j;
        Receive(rcvNum).Apod(:)=1;
        Receive(rcvNum).framenum = i;
        Receive(rcvNum).acqNum = j;
    end
end


% Save Data Process
for N = 1: P.numFrames  % no. of Receive frames (real-time images are produced 1 per frame)
    Process(N).classname = 'External';
    Process(N).method = 'saveRFperFrame'; % calls the 'saveRFperFrame' function
    Process(N).Parameters = {'srcbuffer','receive',... % name of buffer to process.
        'srcbufnum',1,...
        'srcframenum',N,... % process the most recent frame.
        'dstbuffer','none'};
end

% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; % jump back to start
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(2).argument = 300;  % usec
SeqControl(3).command = 'sync';

nsc = 4; % nsc is count of SeqControl objects

% Acquire all frames defined in RcvBuffer
n = 1; % n is count of Events
for i = 1:Resource.RcvBuffer(1).numFrames
    for j = 1:P.numAcqs
        Event(n).info = 'Acquire RF';
        Event(n).tx = 1;
        Event(n).rcv = P.numAcqs*(i-1) + j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;
    end
    % Set last acquisitions SeqControl for transferToHost.
    Event(n-1).seqControl = [3,nsc];
        SeqControl(nsc).command = 'transferToHost'; % transfer all acqs in one super frame
        nsc = nsc + 1;

    Event(n).info = 'Save RF data';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = i;
    Event(n).seqControl = 0;
    n=n+1;

end

Event(n).info = 'Jump back to first event';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 1;

% User specified UI Control Elements
EF(1).Function = text2cell('%saveRFperFrame');

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = P.numAcqs;

% Save all the structures to a .mat file.
save(['MatFiles/',filename]);

return

%saveRFperFrame - save RF
saveRFperFrame(RData)

persistent frameNum folder

if isempty(folder)
    folder = fullfile('UserData', ['RcvData_PerFrame_',datestr(now,'dd-mmmm-yyyy_HH-MM-SS')]);
    mkdir(folder);
end

% file size can be reduced by elimating all zeros
TXApod = evalin('base','TX.Apod');
endSample = evalin('base','Receive(end).endSample');

frameNumLimit = 99;
numDigits = floor(log10(abs(frameNumLimit)))+1;

if isempty(frameNum)
    frameNum = 1;
else
    frameNum = frameNum + 1;
end

% Note: The Value of Version is set to V6 for saving smaller data size without
% compression. Set it to V7.3 for larger files. More information, see:
% https://www.mathworks.com/help/matlab/ref/save.html
if frameNum <= frameNumLimit
    tic
    RFfilename = fullfile( folder, ['Frame', sprintf('_%0*i', numDigits, frameNum) '.mat'] );
    RFdata = RData(1:endSample,find(TXApod),:);
    save(RFfilename,'RFdata','-v6');
    fprintf('Saving time for frame %g = %g s \n',frameNum,toc);
else
    error(sprintf('Exceeded the limit of frame number, %g. Please increase the frameNumLimit \n',frameNumLimit));
end
return
%saveRFperFrame

