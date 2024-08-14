% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name SetUpL11_5vWideBeamSC_64.m - Example of wide beam scan sequence with spatial compounding.
%
% Description:
%   Generate sequence program for L11_5v linear array using wide beam transmit scanning with
%   spatial compounding.  Three wide beam scans are performed using P.numRays scan lines in each
%   frame, with different steering directions (steered left, no steering, steered right).  The
%   frames are averaged after image processing using a running average that updates each frame.
%   The transmit aperture size can be set by the variable - numTx. The receive aperture always
%   covers the full 64 element aperture.

% In the following aperture examples, each space represents 2 elements.
% Aperture 1 (P.numTx=32):
%   tttttttttttttttt------------------------------------------------
%   rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr--------------------------------
% Aperture for almost middle wide beam 24 (P.numTx=32, P.numRays=24):
%   ------------------------tttttttttttttttt------------------------
%   ----------------rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr----------------
% Aperture for last wide beam (P.numTx=32, P.numRays=48):
%   ------------------------------------------------tttttttttttttttt
%   --------------------------------rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
%
% Last update:
% 08/16/2018 - script improvement

clear all

P.numTx = 33;   % no. of elements in TX aperture. Should be odd, max is 63.
if (P.numTx/2 == floor(P.numTx/2)), P.numTx = P.numTx-1; end  % force odd
P.numRays = 48; % no. of rays in frame
P.txFocusMm = 300;   % transmit focal pt in wavelengths
P.dtheta = 10*(pi/180);  % angle delta for beams steered left or right
P.startDepthMm = 2;  % startDepth in mm
P.endDepthMm = 40;  % endDepth in mm
P.maxDepthMm = 70;  % maxDepth for RangeChange and RcvBuffer

% Define system parameters.
Resource.Parameters.speedOfSound = 1540;  % speed of sound in m/s
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.
Resource.System.SoftwareVersion = [3 5 0]; % Minimum software release for this script.
Resource.System.UTA = '260-MUX'; % This script requires the 260-MUX UTA.
Resource.System.transmitChannels = 64;
Resource.System.receiveChannels = 64;

% Specify Trans structure array.
Trans.name = 'L11-5v';
Trans.units = 'wavelengths';
Trans = computeTrans(Trans);
Trans = computeUTAMux64(Trans); % Add HVMux field for use with UTA 260-Mux

% Convert mm to wavelength
demodFreq = Trans.frequency; % demodulation frequency
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);
P.txFocus = P.txFocusMm*scaleToWvl;
P.startDepth = P.startDepthMm*scaleToWvl;  % startDepth in wavelength
P.endDepth = P.endDepthMm*scaleToWvl;
P.maxDepth = P.maxDepthMm*scaleToWvl;
maxBufLength = ceil(sqrt(P.maxDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
maxBufSizePerAcq = 128*ceil(maxBufLength*8*(demodFreq/Trans.frequency)/128);

% Specify PData structure arrays.
PData.PDelta = [Trans.spacing/2, 0, 0.5];
PData.Size(1,1) = ceil((P.endDepth-P.startDepth)/PData.PDelta(3)); % range, origin and pdelta set PData(1).Size.
PData.Size(1,2) = ceil((Trans.numelements*Trans.spacing)/PData.PDelta(1));
PData.Size(1,3) = 1;      % single image page
PData.Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of uppr lft crnr.
% Predefine all the Region structures, then modify them below.
PData.Region = repmat(struct('Shape',struct('Name','Parallelogram',...
                                               'Position',[0,0,P.startDepth],...
                                               'width',6*Trans.spacing*Trans.numelements/(P.numRays),...
                                               'height',P.endDepth-P.startDepth,...
                                               'angle', 0.0)),1,3*P.numRays); % default is no steering
% Compute the x coords of the TX beam centers
TxOrgX = (-63.5*Trans.spacing):(127*Trans.spacing/(P.numRays-1)):(63.5*Trans.spacing);
% Specify P.numRays rectangular regions centered on TX beam origins (use default angle of 0.0).
for n = 1:P.numRays, PData.Region(n).Shape.Position(1) = TxOrgX(n); end
m = P.numRays;
% Define numRays steered left parallelogram regions, centered on TX beam origins. Adjust the angle
%   so that the steering goes to zero over 8 beams at the left and right edge.
for n = 1:P.numRays
    if n<=8
        angle = -((n-1)/8)*P.dtheta;
    elseif n>(P.numRays-8)
        angle = -((P.numRays-n)/8)*P.dtheta;
    else
        angle = -P.dtheta;
    end
    PData.Region(n+m).Shape.Position(1) = TxOrgX(n);
    PData.Region(n+m).Shape.height = (P.endDepth-P.startDepth)/cos(angle);
    PData.Region(n+m).Shape.angle = angle;
end
m = m + P.numRays;
% Define numRays steered right parallelogram regions, centered on TX beam origins. Adjust the angle
%   so that the steering goes to zero over 8 beams at the left and right edge.
for n = 1:P.numRays
    if n<=8
        angle = ((n-1)/8)*P.dtheta;
    elseif n>(P.numRays-8)
        angle = ((P.numRays-n)/8)*P.dtheta;
    else
        angle = P.dtheta;
    end
    PData.Region(n+m).Shape.Position(1) = TxOrgX(n);
    PData.Region(n+m).Shape.height = (P.endDepth-P.startDepth)/cos(angle);
    PData.Region(n+m).Shape.angle = angle;
end
PData.Region = computeRegions(PData);

% Specify Media object.
pt1;
Media.attenuation = -0.5;
Media.function = 'movePoints';

% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = P.numRays*maxBufSizePerAcq; % this size allows for all rays, with range of up to 400 wvlngths
Resource.RcvBuffer(1).colsPerFrame = Resource.System.receiveChannels;
Resource.RcvBuffer(1).numFrames = 12;
Resource.RcvBuffer(2).datatype = 'int16';
Resource.RcvBuffer(2).rowsPerFrame = P.numRays*maxBufSizePerAcq; % this size allows for all rays, with range of up to 400 wvlngths
Resource.RcvBuffer(2).colsPerFrame = Resource.System.receiveChannels;
Resource.RcvBuffer(2).numFrames = 12;  % all RcvBuffers should have same no. of frames
Resource.RcvBuffer(3).datatype = 'int16';
Resource.RcvBuffer(3).rowsPerFrame = P.numRays*maxBufSizePerAcq; % this size allows for all rays, with range of up to 400 wvlngths
Resource.RcvBuffer(3).colsPerFrame = Resource.System.receiveChannels;
Resource.RcvBuffer(3).numFrames = 12;
Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).numFrames = 1;  % one intermediate buffer needed.
Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).numFrames = 10;
Resource.DisplayWindow(1).Title = 'L11-5vWideBeamSC';
Resource.DisplayWindow(1).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData.Size(2)*PData.PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData.Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData.Origin(1),0,PData.Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).numFrames = 50;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = gray(256);

% Specify Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,0.67,2,1];

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'Apod', zeros(1,Resource.System.transmitChannels), ...
                   'focus', P.txFocus, ...
                   'Steer', [0.0,0.0], ...
                   'Delay', zeros(1,Resource.System.transmitChannels), ...
                   'TXPD', [], ...
                   'peakCutOff', 0.3, ...
                   'peakBLMax', 4.0), 1, 3*P.numRays);

% - Set event specific TX attributes.  We will specify a MUX aperture the same size
%   as the number of TX & Receive channels, even though the number of active TX channels
%   set by P.numTx is usually smaller.  This will allow the same MUX aperture to be used
%   for receive, eliminating the time for the MUX to switch to a different aperture.
Ce = fix(1:(127/(P.numRays-1)):128); % calculate center element numbers.
ApodFunction = kaiser(P.numTx,1); % Apodization function for TX
for n = 1:P.numRays  % specify P.numRays transmit events
    TX(n).Origin(1) = TxOrgX(n);

    apodStartIndx = 1; % default start index of apodization function.
    lft = Ce(n) - floor(P.numTx/2);
    if lft < 1
        apodStartIndx = -lft+2;
        lft = 1;
    end
    apodEndIndx = P.numTx;
    rt = Ce(n) + floor(P.numTx/2);
    if rt > Trans.numelements
        apodEndIndx = P.numTx - (rt - Trans.numelements);
        rt = Trans.numelements;
    end
    TX(n).aperture = max(1,rt-63);
    TX(n).Apod((lft:rt)-TX(n).aperture+1) = 1;
    % Apply apodization function.
    [~,CIndices,~] = find(TX(n).Apod);
    V = ApodFunction(apodStartIndx:apodEndIndx);
    TX(n).Apod(CIndices) = V;
    % Compute transmit delays
    TX(n).Delay = computeTXDelays(TX(n));
    % Compute transmit pixel data
end
m = P.numRays;
for n = 1:P.numRays
    TX(n+m).Origin(1) = TX(n).Origin(1);
    TX(n+m).aperture = TX(n).aperture;
    TX(n+m).Apod = TX(n).Apod;
    if n<=8
        TX(n+m).Steer = [-((n-1)/8)*P.dtheta,0.0];
    elseif n>(P.numRays-8)
        TX(n+m).Steer = [-((P.numRays-n)/8)*P.dtheta,0.0];
    else
        TX(n+m).Steer = [-P.dtheta,0.0];
    end
    % Compute transmit delays
    TX(n+m).Delay = computeTXDelays(TX(n+m));
end
m = m + P.numRays;
for n = 1:P.numRays
    TX(n+m).Origin(1) = TX(n).Origin(1);
    TX(n+m).aperture = TX(n).aperture;
    TX(n+m).Apod = TX(n).Apod;
    if n<=8
        TX(n+m).Steer = [((n-1)/8)*P.dtheta,0.0];
    elseif n>(P.numRays-8)
        TX(n+m).Steer = [((P.numRays-n)/8)*P.dtheta,0.0];
    else
        TX(n+m).Steer = [P.dtheta,0.0];
    end
    % Compute transmit delays
    TX(n+m).Delay = computeTXDelays(TX(n+m));
end

% calculate transmit TXPD
h = waitbar(0,'Program TX parameters, please wait!');
steps = 3*P.numRays;
for n = 1:steps
    TX(n).TXPD = computeTXPD(TX(n),PData);
    waitbar(n/steps)
end
close(h)

% Specify Receive structure arrays.
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
Receive = repmat(struct('Apod', ones(1,Resource.System.receiveChannels), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, 3*P.numRays*Resource.RcvBuffer(1).numFrames);
% - Set event specific Receive attributes for each frame.
Receive(1).callMediaFunc = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    k = 3*P.numRays*(i-1);
    Receive(k+1).callMediaFunc = 1;
    for j = 1:P.numRays
        Receive(k+j).aperture = TX(j).aperture;
        Receive(k+j).framenum = i;
        Receive(k+j).acqNum = j;
        Receive(k+P.numRays+j).aperture = TX(j).aperture;
        Receive(k+P.numRays+j).bufnum = 2;
        Receive(k+P.numRays+j).framenum = i;
        Receive(k+P.numRays+j).acqNum = j;
        Receive(k+2*P.numRays+j).aperture = TX(j).aperture;
        Receive(k+2*P.numRays+j).bufnum = 3;
        Receive(k+2*P.numRays+j).framenum = i;
        Receive(k+2*P.numRays+j).acqNum = j;
    end
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [0,225,410,475,549,659,736,791];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Recon structure arrays.
% - We need three Recon structures, one for each steering direction.
Recon = repmat(struct('senscutoff', 0.5, ...
               'pdatanum', 1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...
               'rcvBufFrame',-1, ...
               'RINums',1:P.numRays), 1, 3);

% - Set specific Recon attributes.
Recon(2).RINums = (P.numRays+1):(2*P.numRays);
Recon(3).RINums = (2*P.numRays+1):(3*P.numRays);

% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 'accumIQ', ...  % default is to accumulate IQ data.
                   'Pre',[],...
                   'Post',[],...
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'scaleFactor', 0.5, ...
                   'regionnum', 0), 1, 3*P.numRays);

% - Set specific ReconInfo attributes.
ReconInfo(1).Pre = 'clearInterBuf';
for j = 1:P.numRays
    ReconInfo(j).txnum = j;
    ReconInfo(j).rcvnum = j;
    ReconInfo(j).regionnum = j;
end
ReconInfo(P.numRays).Post = 'IQ2IntensityImageBuf';
k = P.numRays;
ReconInfo(k+1).Pre = 'clearInterBuf';
for j = 1:P.numRays
    ReconInfo(j+k).txnum = j+k;
    ReconInfo(j+k).rcvnum = j+k;
    ReconInfo(j+k).regionnum = j+k;
end
ReconInfo(k+P.numRays).Post = 'IQ2IntensityImageBuf';
k = k + P.numRays;
ReconInfo(k+1).Pre = 'clearInterBuf';
for j = 1:P.numRays
    ReconInfo(j+k).txnum = j+k;
    ReconInfo(j+k).rcvnum = j+k;
    ReconInfo(j+k).regionnum = j+k;
end
ReconInfo(k+P.numRays).Post = 'IQ2IntensityImageBuf';

% Specify Process structure array.
pers = 20;
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',1.0,...            % pgain is image processing gain
                         'reject',2,...      % reject level
                         'persistMethod','simple',...
                         'persistLevel',pers,...
                         'interpMethod','4pt',...
                         'grainRemoval','low',...
                         'processMethod','reduceSpeckle2',...
                         'averageMethod','runAverage3',...
                         'compressMethod','power',...
                         'compressFactor',45,...
                         'mappingMethod','full',...
                         'display',1,...      % display image after processing
                         'displayWindow',1};

% Specify SeqControl structure arrays.
%  - Jump back to start.
SeqControl(1).command = 'jump';
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(2).argument = 220;  % 220 usec
SeqControl(3).command = 'timeToNextAcq';  % optional time between frames (not needed for synchronous operation)
SeqControl(3).argument = 30000;  % 30000 usec = 30msec time between frames
SeqControl(4).command = 'returnToMatlab';
nsc = 5; % nsc is count of SeqControl objects

% Specify Event structure arrays.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    % Acquire frame with straight ahead wide beams
    k = 3*P.numRays*(i-1);
    for j = 1:P.numRays        % Acquire frame
        Event(n).info = 'acquire aperture.';
        Event(n).tx = j;
        Event(n).rcv = k+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;
    end
    Event(n-1).seqControl = [3,nsc]; % modify last acquisition Event's seqControl
      SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
      nsc = nsc+1;

    Event(n).info = 'recon and process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 1;
    Event(n).seqControl = 0;
    n = n+1;

    % Acquire frame with steered left wide beams
    for j = 1:P.numRays        % Acquire frame
        Event(n).info = 'acquire aperture.';
        Event(n).tx = j+P.numRays;
        Event(n).rcv = k+P.numRays+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;
    end
    Event(n-1).seqControl = [3,nsc]; % modify last acquisition Event's seqControl
      SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
      nsc = nsc+1;

    Event(n).info = 'recon and process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 2;
    Event(n).process = 1;
    Event(n).seqControl = 0;
    n = n+1;

    % Acquire frame with steered right wide beams
    for j = 1:P.numRays        % Acquire frame
        Event(n).info = 'acquire aperture.';
        Event(n).tx = j+2*P.numRays;
        Event(n).rcv = k+2*P.numRays+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;
    end
    Event(n-1).seqControl = [3,nsc]; % modify last acquisition Event's seqControl
      SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
      nsc = nsc+1;

    Event(n).info = 'recon and process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 3;
    Event(n).process = 1;
    Event(n).seqControl = 0;
    if i ~= Resource.RcvBuffer(1).numFrames  % Exit to Matlab every 3rd frame
        Event(n).seqControl = 4;
    end
    n = n+1;
end

Event(n).info = 'Jump back';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 1;


% User specified UI Control Elements
% - Sensitivity Cutoff
UI(1).Control =  {'UserB7','Style','VsSlider','Label','Sens. Cutoff',...
                  'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                  'SliderStep',[0.025,0.1],'ValueFormat','%1.3f'};
UI(1).Callback = text2cell('%SensCutOffCallback');

% - Range Change
scaleToWvl = 1;
AxesUnit = 'mm';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if ~strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
        AxesUnit = 'wls';
        scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);
    end
end
UI(2).Control = {'UserA1','Style','VsSlider','Label',['Range (',AxesUnit,')'],...
                 'SliderMinMaxVal',[2,P.maxDepthMm,P.endDepthMm]*scaleToWvl,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f'};
UI(2).Callback = text2cell('%RangeChangeCallback');

% - Transmit focus change
UI(3).Control = {'UserB4','Style','VsSlider','Label',['TX Focus (',AxesUnit,')'],...
                 'SliderMinMaxVal',[10,600,P.txFocusMm]*scaleToWvl,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f'};
UI(3).Callback = text2cell('%TxFocusCallback');

% - TX Aperture change
UI(4).Control = {'UserB3','Style','VsSlider','Label','TX Aper',...
                 'SliderMinMaxVal',[1,128,P.numTx],'SliderStep',[1/128,1/32],'ValueFormat','%3.0f'};
UI(4).Callback = text2cell('%ApertureCallback');

% - Peak CutOff
UI(5).Control = {'UserB2','Style','VsSlider','Label','Peak Cutoff',...
                  'SliderMinMaxVal',[0,20.0,TX(1).peakCutOff],...
                  'SliderStep',[0.005,0.020],'ValueFormat','%1.3f'};
UI(5).Callback = text2cell('%PeakCutOffCallback');

% - Max. Burst Length
UI(6).Control = {'UserB1','Style','VsSlider','Label','Max. BL',...
                  'SliderMinMaxVal',[0,20.0,TX(1).peakBLMax],...
                  'SliderStep',[0.005,0.020],'ValueFormat','%1.3f'};
UI(6).Callback = text2cell('%MaxBLCallback');

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 3;

% Save all the structures to a .mat file.
filename = 'MatFiles/L11-5vWideBeamSC_64';
save(filename)

return

% **** Callback routines to be converted by text2cell function. ****
%SensCutOffCallback - Sensitivity cutoff change
ReconL = evalin('base', 'Recon');
for i = 1:size(ReconL,2)
    ReconL(i).senscutoff = UIValue;
end
assignin('base','Recon',ReconL);
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'Recon'};
assignin('base','Control', Control);
return
%SensCutOffCallback

%RangeChangeCallback - Range change
Trans = evalin('base','Trans');
Resource = evalin('base','Resource');
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);

P = evalin('base','P');
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
        P.endDepthMm = UIValue;
        P.endDepth = UIValue*scaleToWvl;
    else
        P.endDepth = UIValue;
        P.endDepthMm = P.endDepth/scaleToWvl;
    end
end
assignin('base','P',P);

% No range change if in simulate mode 2.
if Resource.Parameters.simulateMode == 2
    return
end

% Modify PData and Regions for new depth.
PData = evalin('base','PData');
PData.Size(1) = ceil((P.endDepth-P.startDepth)/PData.PDelta(3)); % P.startDepth, P.endDepth and pdelta set PData.Size(1).
PData.Region = repmat(struct('Shape',struct( ...
                    'Name','Parallelogram',...
                    'Position',[0,0,P.startDepth],...
                    'width',5*Trans.spacing*Trans.numelements/(P.numRays),...
                    'height',P.endDepth-P.startDepth,...
                    'angle',0.0)),1,3*P.numRays);
TxOrgX = (-63.5*Trans.spacing):(127*Trans.spacing/(P.numRays-1)):(63.5*Trans.spacing);
for n = 1:P.numRays, PData.Region(n).Shape.Position(1) = TxOrgX(n); end
m = P.numRays;
for n = 1:P.numRays
    if n<=8
        angle = -((n-1)/8)*P.dtheta;
    elseif n>(P.numRays-8)
        angle = -((P.numRays-n)/8)*P.dtheta;
    else
        angle = -P.dtheta;
    end
    PData.Region(n+m).Shape.Position(1) = TxOrgX(n);
    PData.Region(n+m).Shape.height = (P.endDepth-P.startDepth)/cos(angle);
    PData.Region(n+m).Shape.angle = angle;
end
m = m + P.numRays;
% Define numRays steered right parallelogram regions, centered on TX beam origins. Adjust the angle
%   so that the steering goes to zero over 8 beams at the left and right edge.
for n = 1:P.numRays
    if n<=8
        angle = ((n-1)/8)*P.dtheta;
    elseif n>(P.numRays-8)
        angle = ((P.numRays-n)/8)*P.dtheta;
    else
        angle = P.dtheta;
    end
    PData.Region(n+m).Shape.Position(1) = TxOrgX(n);
    PData.Region(n+m).Shape.height = (P.endDepth-P.startDepth)/cos(angle);
    PData.Region(n+m).Shape.angle = angle;
end
PData.Region = computeRegions(PData);
assignin('base','PData',PData);
evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData.Size(1)*PData.PDelta(3)/Resource.DisplayWindow(1).pdelta);');

% Update TXPD data of TX structures.
TX = evalin('base','TX');
h = waitbar(0,'Program TX parameters, please wait!');
steps = 3*P.numRays;
for ind = 1:steps
    TX(ind).TXPD = computeTXPD(TX(ind),PData);
    waitbar(ind/steps)
end
close(h)
assignin('base','TX',TX);

% Update Receive and TGC
Receive = evalin('base', 'Receive');
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
for i = 1:size(Receive,2)
    Receive(i).endDepth = maxAcqLength;
end
assignin('base','Receive',Receive);
evalin('base','TGC.rangeMax = P.endDepth;');
evalin('base','TGC.Waveform = computeTGCWaveform(TGC);');

Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'PData','InterBuffer','ImageBuffer','DisplayWindow','TX','Receive','TGC','Recon'};
assignin('base','Control', Control);
assignin('base', 'action', 'displayChange');
return
%RangeChangeCallback

%TxFocusCallback
simMode = evalin('base','Resource.Parameters.simulateMode');
% No focus change if in simulate mode 2.
if simMode == 2
    set(hObject,'Value',evalin('base','P.txFocus'));
    return
end
Trans = evalin('base','Trans');
Resource = evalin('base','Resource');
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);

P = evalin('base','P');
P.txFocus = UIValue;
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
        P.txFocus = UIValue*scaleToWvl;
    end
end
assignin('base','P',P);

% - Redefine event specific TX attributes for the new focus.
TX = evalin('base', 'TX');
PData = evalin('base','PData');
h = waitbar(0,'Program TX parameters, please wait!');
steps = 3*P.numRays;
for ind = 1:steps
    % write new focus value to TX
    TX(ind).focus = P.txFocus;
    TX(ind).Delay = computeTXDelays(TX(ind));
    TX(ind).TXPD = computeTXPD(TX(ind),PData);
    waitbar(ind/steps)
end
close(h)
assignin('base','TX',TX);
% Set Control command to update TX
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'TX'};
assignin('base','Control', Control);
return
%TxFocusCallback

%ApertureCallback
simMode = evalin('base','Resource.Parameters.simulateMode');
% No focus change if in simulate mode 2.
if simMode == 2
    set(hObject,'Value',evalin('base','P.numTx'));
    return
end
numTx = UIValue;
if (numTx/2 == floor(numTx/2)),numTx = numTx-1; end
evalin('base',['P.numTx = ',int2str(numTx),';']);
Trans = evalin('base', 'Trans');
numRays = evalin('base', 'P.numRays');
dtheta = evalin('base','P.dtheta');
TX = evalin('base', 'TX');
PData = evalin('base','PData');
scaleToWvl = evalin('base','scaleToWvl');
TxOrgX = (-63.5*Trans.spacing):(127*Trans.spacing/(numRays-1)):(63.5*Trans.spacing);
Ce = fix(1:(127/(numRays-1)):128);
ApodFunction = kaiser(numTx,1); % Apodization function for TX
% - Redefine event specific TX attributes for the new aperture.
h = waitbar(0,'Program TX parameters, please wait!');
steps = 3*numRays;
for n = 1:numRays  % specify P.numRays transmit events
    TX(n).Origin(1) = TxOrgX(n);

    apodStartIndx = 1; % default start index of apodization function.
    lft = Ce(n) - floor(numTx/2);
    if lft < 1
        apodStartIndx = -lft+2;
        lft = 1;
    end
    apodEndIndx = numTx;
    rt = Ce(n) + floor(numTx/2);
    if rt > Trans.numelements
        apodEndIndx = numTx - (rt - Trans.numelements);
        rt = Trans.numelements;
    end
    TX(n).aperture = max(1,rt-63);
    TX(n).Apod((lft:rt)-TX(n).aperture+1) = 1;
    % Apply apodization function.
    [~,CIndices,~] = find(TX(n).Apod);
    V = ApodFunction(apodStartIndx:apodEndIndx);
    TX(n).Apod(CIndices) = V;
    % Compute transmit delays
    TX(n).Delay = computeTXDelays(TX(n));
    % Compute transmit pixel data
    TX(n).TXPD = computeTXPD(TX(n),PData(1));
%     waitbar(n/steps);
end
m = numRays;
for n = 1:numRays
    TX(n+m).Origin(1) = TX(n).Origin(1);
    TX(n+m).aperture = TX(n).aperture;
    TX(n+m).Apod = TX(n).Apod;
    if n<=8
        TX(n+m).Steer = [-((n-1)/8)*dtheta,0.0];
    elseif n>(numRays-8)
        TX(n+m).Steer = [-((numRays-n)/8)*dtheta,0.0];
    else
        TX(n+m).Steer = [-dtheta,0.0];
    end
    % Compute transmit delays
    TX(n+m).Delay = computeTXDelays(TX(n+m));
    % Compute transmit pixel data
    TX(n+m).TXPD = computeTXPD(TX(n+m),PData);
    waitbar((n+numRays)/steps);
end
m = m + numRays;
for n = 1:numRays
    TX(n+m).Origin(1) = TX(n).Origin(1);
    TX(n+m).aperture = TX(n).aperture;
    TX(n+m).Apod = TX(n).Apod;
    if n<=8
        TX(n+m).Steer = [((n-1)/8)*dtheta,0.0];
    elseif n>(numRays-8)
        TX(n+m).Steer = [((numRays-n)/8)*dtheta,0.0];
    else
        TX(n+m).Steer = [dtheta,0.0];
    end
    % Compute transmit delays
    TX(n+m).Delay = computeTXDelays(TX(n+m));
    % Compute transmit pixel data
    TX(n+m).TXPD = computeTXPD(TX(n+m),PData);
    waitbar((n+2*numRays)/steps);
end
close(h)
assignin('base','TX', TX);
assignin('base','Trans', Trans);
% Set Control command to update TX
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'TX','Trans'};
assignin('base','Control', Control);
return
%ApertureCallback

%PeakCutOffCallback
TX = evalin('base', 'TX');
for i=1:size(TX,2)
    TX(i).peakCutOff = UIValue;
end
assignin('base','TX',TX);
% Set Control.Command to set TX
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'TX'};
assignin('base','Control', Control);
%PeakCutOffCallback

%MaxBLCallback
TX = evalin('base', 'TX');
for i=1:size(TX,2)
    TX(i).peakBLMax = UIValue;
end
assignin('base','TX',TX);
% Set Control.Command to set TX
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'TX'};
assignin('base','Control', Control);
%MaxBLCallback
