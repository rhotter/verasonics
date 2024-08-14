% ErrorScript1 has a bug about an empty field in the event sequence.
% EventAnalysisTool shows a warning message if any of event has an empty field

clear all

P.startDepth = 5;   % Acquisition depth in wavelengths
P.endDepth = 192;   % This should preferrably be a multiple of 128 samples.

% Define system parameters.
Resource.Parameters.numTransmit = 128;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;    % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Specify Trans structure array.
Trans.name = 'L7-4';
Trans.units = 'mm'; % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);  % L7-4 transducer is 'known' transducer so we can use computeTrans.
Trans.maxHighVoltage = 50;  % set maximum high voltage limit for pulser supply.

% Specify PData structure array.
PData(1).PDelta = [Trans.spacing, 0, 0.5];
PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3)); % startDepth, endDepth and pdelta set PData(1).Size.
PData(1).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(3) = 1;      % single image page
PData(1).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of upper lft crnr.

% Specify Media object. 'pt1.m' script defines array of point targets.
pt1;
Media.function = 'movePoints';

% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 4096;   % this size allows for maximum range
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 1;       % 10 frames used for RF cineloop.
Resource.InterBuffer(1).numFrames = 1;   % one intermediate buffer needed.
Resource.ImageBuffer(1).numFrames = 10;
Resource.DisplayWindow.Title = 'Image Display';
Resource.DisplayWindow.pdelta = 0.3;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = gray(256);

% Specify Transmit waveform structure.
TW.type = 'parametric';
TW.Parameters = [Trans.frequency,.67,2,1];   % A, B, C, D

% Specify TX structure array.
TX(1).waveform = 1;            % use 1st TW structure.
TX(1).Origin = [0.0,0.0,0.0];  % flash transmit origin at (0,0,0).
TX(1).focus = 0;
TX(1).Steer = [0.0,0.0];       % theta, alpha = 0.
TX(1).Apod = ones(1,Trans.numelements);
TX(1).Delay = computeTXDelays(TX(1));

% Specify TGC Waveform structure.
TGC.CntrlPts = [0,200,344,452,606,747,870,920];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Receive structure arrays -
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
wlsPer128 = 128/(4*2); % wavelengths in 128 samples for 4 samplesPerWave
maxAcqLength = wlsPer128*ceil(maxAcqLength/wlsPer128);
Receive(1).Apod = ones(1,Trans.numelements);
Receive(1).startDepth = P.startDepth;
Receive(1).endDepth = maxAcqLength;
Receive(1).TGC = 1;
Receive(1).bufnum = 1;
Receive(1).framenum = 1;
Receive(1).acqNum = 1;
Receive(1).sampleMode = 'NS200BW';
Receive(1).mode = 0;

% Specify Recon structure arrays.
Recon(1).senscutoff = 0.6;
Recon(1).pdatanum = 1;
Recon(1).rcvBufFrame = 1;
Recon(1).IntBufDest = [1,1];
Recon(1).ImgBufDest = [1,-1];
Recon(1).RINums = 1;

% Define ReconInfo structures.
ReconInfo(1).mode = 'replaceIntensity';
ReconInfo(1).txnum = 1;
ReconInfo(1).rcvnum = 1;
ReconInfo(1).regionnum = 1;

% external processing
Process(1).classname = 'External';
Process(1).method = 'TestDisplay';
Process(1).Parameters = {'srcbuffer','receive',...
                         'srcbufnum',1,...
                         'srcframenum',1,...
                         'dstbuffer','none'};

% scan sequence
n = 1;
nsc = 1;

Event(n).info = 'acquisition';
Event(n).tx = 1;         % use 1st TX structure.
Event(n).rcv = 1;      % use 1st Rcv structure.
Event(n).recon = 0;      % no reconstruction.
Event(n).process = 0;    % no processing
Event(n).seqControl = nsc;
        SeqControl(nsc).command = 'noop';
        SeqControl(nsc).argument = 500;
        nsc = nsc + 1;
n = n+1;

Event(n).info = 'transfer to host';
Event(n).tx = 0;         % use 1st TX structure.
Event(n).rcv = 0;      % use 2nd Rcv structure.
Event(n).recon = 0;      % no reconstruction.
Event(n).process = 0;    % no processing
Event(n).seqControl = nsc; % time between frames, SeqControl struct defined below.
        SeqControl(nsc).command = 'transferToHost';
        nsc = nsc + 1;
n = n+1;

Event(n).info = 'Reconstruct';
Event(n).tx = 0;         % no transmit
Event(n).rcv = 0;        % no rcv
Event(n).recon = 1;      % reconstruction
Event(n).process = 0;    % processing
n = n+1;

Event(n).info = 'external processing';
Event(n).tx = 0;        % no TX
Event(n).rcv = 0;       % no Rcv
Event(n).recon = 0;     % no Recon
Event(n).process = 1;
Event(n).seqControl = nsc; % jump command
        SeqControl(nsc).command = 'jump';
        SeqControl(nsc).argument = 1;
        nsc = nsc + 1;
n = n+1;

Event(n).info = 'Jump back to first event';
Event(n).tx = 0;        % no TX
Event(n).rcv = 0;       % no Rcv
Event(n).recon = 0;     % no Recon
Event(n).process = 0;
Event(n).seqControl = nsc; % jump command
        SeqControl(nsc).command = 'jump';
        SeqControl(nsc).argument = 1;
        nsc = nsc + 1;
n = n+1;

filename = 'MatFiles/Error1';
save(filename);
% VSX;
