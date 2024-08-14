% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: SetUpL11_4vAcquireRF_64LE.m - Example of RF data acquisition
%
% Description:
%   Sequence programming file for L11-4v Linear array, acquiring RF data of
%   a single plane wave transmit with 2-1 synthetic aperture (two
%   acquisitions per line to acquire data from all 128 receive elements).
%   All 128 transmit and 64 receive channels are active for each synthetic
%   aperture acquisition. External processing is used asynchronous with
%   respect to acquisition.
%
% Last update:
% 12/14/2015 - modified for SW 3.0
% 05/06/2020 - Update to SW 4.3 format for new UIControls and function definitions (VTS 1691). 
%   More info:(.../Example_Scripts/Vantage_Features/New UI Scheme/SetUpL11_5vFlash_NewUI)

clear all

% Specify system parameters
Resource.Parameters.numTransmit = 128;      % no. of transmit channels (2 brds).
Resource.Parameters.numRcvChannels = 64;    % no. of receive channels (2 brds).
Resource.Parameters.speedOfSound = 1540;    % speed of sound in m/sec
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;       % runs script in simulate mode
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Specify media points
Media.MP(1,:) = [0,0,100,1.0]; % [x, y, z, reflectivity]
Media.function = 'movePoints';

% Specify Trans structure array.
Trans.name = 'L11-4v';
Trans.units = 'wavelengths'; % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);  % L11-4v transducer is 'known' transducer.

% Specify Resource buffers.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 4096;  % this allows for 1/4 maximum range
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numTransmit;
Resource.RcvBuffer(1).numFrames = 10;       % allocate 10 frames.

% Specify Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,.67,2,1];

% Specify TX structure array.
TX(1).waveform = 1;            % use 1st TW structure.
TX(1).focus = 0;
TX(1).Apod = ones(1,Trans.numelements);
TX(1).Delay = computeTXDelays(TX(1));

% Specify TGC Waveform structure.
TGC(1).CntrlPts = [500,590,650,710,770,830,890,950];
TGC(1).rangeMax = 200;
TGC(1).Waveform = computeTGCWaveform(TGC);

% Specify Receive structure array -
Receive = repmat(struct(...
                'Apod', zeros(1,Trans.numelements), ...
                'startDepth', 0, ...
                'endDepth', 200, ...
                'TGC', 1, ...
                'mode', 0, ...
                'bufnum', 1, ...
                'framenum', 1, ...
                'acqNum', 1, ...
                'sampleMode', 'NS200BW'), ...
                1,2*Resource.RcvBuffer(1).numFrames);

% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    % -- 1st synthetic aperture acquisition.
    Receive(2*i-1).Apod(1:Resource.Parameters.numRcvChannels) = 1;
    Receive(2*i-1).framenum = i;
    Receive(2*i-1).acqNum = 1;
    Receive(2*i-1).callMediaFunc = 1;
    % -- 2nd synthetic aperture acquisition.
    Receive(2*i).Apod = zeros(1, 128);
    Receive(2*i).Apod(65:128) = 1;
    Receive(2*i).framenum = i;
    Receive(2*i).acqNum = 2;   % two acquisitions per frame
end

% Specify an external processing event.
Process(1).classname = 'External';
Process(1).method = 'myProcFunction';
Process(1).Parameters = {'srcbuffer','receive',... % name of buffer to process.
                         'srcbufnum',1,...
                         'srcframenum',-1,... % process the most recent frame.
                         'dstbuffer','none'};

% Specify sequence events.
SeqControl(1).command = 'timeToNextAcq';
SeqControl(1).argument = 200;
SeqControl(2).command = 'timeToNextAcq';
SeqControl(2).argument = 9800;          % time to next frame is 10ms
SeqControl(3).command = 'jump';
SeqControl(3).argument = 1;
nsc = 4; % start index for new SeqControl

n = 1;   % start index for Events
for i = 1:Resource.RcvBuffer(1).numFrames

	Event(n).info = 'Aquisition 1st half syn. aperture.';
	Event(n).tx = 1;
	Event(n).rcv = 2*i-1;
	Event(n).recon = 0;
	Event(n).process = 0;
	Event(n).seqControl = 1;
    n = n+1;

	Event(n).info = 'Aquisition 2nd half syn. aperture.';
	Event(n).tx = 1;
	Event(n).rcv = 2*i;
	Event(n).recon = 0;
	Event(n).process = 0;
	Event(n).seqControl = [2,nsc];
	   SeqControl(nsc).command = 'transferToHost';
	   nsc = nsc + 1;
    n = n+1;

	Event(n).info = 'Call external Processing function.';
	Event(n).tx = 0;
	Event(n).rcv = 0;
	Event(n).recon = 0;
	Event(n).process = 1;
	Event(n).seqControl = 0;
	n = n+1;
end
Event(n).info = 'Jump back to Event 1.';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 3;

%% User specified UI Control Elements
nr = Resource.Parameters.numRcvChannels;
UI(1).Control = vsv.seq.uicontrol.VsSliderControl('LocationCode','UserB1','Label','Plot Channel',...
                   'SliderMinMaxVal',[1,128,64],...
                   'SliderStep',[1/(2*nr),8/nr],'ValueFormat', '%3.0f');
UI(1).Callback = @(src,evt, UIValue)assignin('base', 'myPlotChnl', round(UIValue));

% External function delimiter definition
EF(1).Function = vsv.seq.function.ExFunctionDef('myProcFunction',@myProcFunction);

% Save all the structures to a .mat file.
save('MatFiles/L11-4vAcquireRF_64LE');
return

%% **** Callback routines used by External function definition (EF) ****

function myProcFunction(RData)
    persistent myHandle
    Receive = evalin('base','Receive');
    % If myPlotChnl exists, read it for the channel to plot.
    if evalin('base','exist(''myPlotChnl'',''var'')')
        channel = evalin('base','myPlotChnl');
    else
        channel = 64;  % Channel no. to plot
    end
    % Create the figure if it doesn't exist.
    if isempty(myHandle)||~ishandle(myHandle)
        figure('name','Receive Signal','NumberTitle','off');
        myHandle = axes('XLim',[0,Receive(1).endSample],'YLim',[-2048 2048], ...
                        'NextPlot','replacechildren');
    end
    % Plot the RF data.
    if channel < 65, plot(myHandle,RData(1:Receive(1).endSample,channel));
    else plot(myHandle,RData(Receive(2).startSample:Receive(2).endSample,channel)); end
    drawnow limitrate
end

