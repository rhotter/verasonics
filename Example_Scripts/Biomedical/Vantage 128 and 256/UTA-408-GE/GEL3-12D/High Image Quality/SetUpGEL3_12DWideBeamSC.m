% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: SetUpGEL3_12DWideBeamSC.m - Example of wide beam imaging
%                                        with deep focal zone and spatial
%                                        compounding
% Description:
%    Generate .mat file for GEL3-12D linear array for Vantage system.
%    Transmit/Receive is performed with a 128 element mux aperture that is
%    moved across the 256 element aperture. The Tx aperture is moved with each
%    ray line, but the mux aperture is moved only for ray lines in the
%    central portion of the transducer aperture where the full 128 channel
%    aperture can be centered around the Tx aperture.
%
%    In the following aperture examples, each space represents 4 elements.
%    Aperture 1 (P.numTx=12):
%      ttt-------------------------------------------------------------
%      rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr--------------------------------
%    Aperture for middle beam 128 (TX.aperture=65, P.numTx=24, numRay=128):
%      -----------------------------tttttt-----------------------------
%      ----------------rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr----------------
%    Aperture for wide beam (TX.aperture=129, P.numTx=12, numRay=256):
%      -------------------------------------------------------------ttt
%      --------------------------------rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
%
%   The receive data from each aperture are stored under different acqNums in the Receive buffer.
%   This version does asynchronous acquisition and processing.
%
% Last update:
% 11/30/2016 - test with SW 3.2

clear all

P.numTx = 60;    % number of active transmitters in TX aperture.
P.numRays = 64;
P.startDepthMm = 0;  % startDepth in mm
P.endDepthMm = 50;  % endDepth in mm
P.maxDepthMm = 100;  % maxDepth for RangeChange and RcvBuffer
P.txFocusMm = 300;   % transmit focal pt in wavelengths
P.dtheta = 10*(pi/180);  % angle delta for beams steered left or right

% Define system parameters.
Resource.Parameters.numTransmit = 128;  % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;  % number of receive channels.
Resource.Parameters.speedOfSound = 1540;  % speed of sound in m/s
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Specify Trans structure array.
Trans.name = 'GEL3-12D';
Trans.units = 'wavelengths'; % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);

% Convert mm to wavelength
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);
P.txFocus = P.txFocusMm*scaleToWvl;
P.startDepth = P.startDepthMm*scaleToWvl;  % startDepth in wavelength
P.endDepth = P.endDepthMm*scaleToWvl;
P.maxDepth = P.maxDepthMm*scaleToWvl;

% Specify PData structure array.
PData.PDelta = [Trans.spacing, 0, 0.5];
PData.Size(1) = ceil((P.endDepth-P.startDepth)/PData.PDelta(3)); % P.startDepth, P.endDepth and pdelta set PData.Size.
PData.Size(2) = ceil((Trans.numelements*Trans.spacing)/PData.PDelta(1));
PData.Size(3) = 1;      % single image page
PData.Origin = [-((Trans.numelements-1)/2)*Trans.spacing,0,P.startDepth]; % x,y,z of upper lft crnr.

% Define P.numRays rectangular regions, positioned around the origin of the wide beams.
TxOrgX = (-127.5*Trans.spacing):(255*Trans.spacing/(P.numRays-1)):(127.5*Trans.spacing); % x coords of beam centers
for n = 1:P.numRays
    PData.Region(n) = struct('Shape',struct('Name','Rectangle','Position',[TxOrgX(n),0,0],'width',32,'height',P.endDepth));
end
m = P.numRays;
% Define P.numRays steered left parallelogram regions, positioned around the origin of the wide beams.
for n = 1:P.numRays
    if n<=10
        angle = -((n-1)/10)*P.dtheta;
    elseif n>(P.numRays-10)
        angle = -((P.numRays-n)/10)*P.dtheta;
    else
        angle = -P.dtheta;
    end
    PData.Region(n+m) = struct('Shape',struct('Name','Parallelogram','Position',[TxOrgX(n),0,0],'width',32,'height',P.endDepth/cos(P.dtheta),'angle',angle));
end
m = m + P.numRays;
% Define P.numRays steered right parallelogram regions, positioned around the origin of the wide beams.
for n = 1:P.numRays
    if n<=10
        angle = ((n-1)/10)*P.dtheta;
    elseif n>(P.numRays-10)
        angle = ((P.numRays-n)/10)*P.dtheta;
    else
        angle = P.dtheta;
    end
    PData.Region(n+m) = struct('Shape',struct('Name','Parallelogram','Position',[TxOrgX(n),0,0],'width',32,'height',P.endDepth/cos(P.dtheta),'angle',angle));
end
PData.Region = computeRegions(PData);

% Specify Media object. 'pt1.m' script defines array of point targets.
pt1;
Media.attenuation = -0.5;
Media.function = 'movePoints';

demodFreq = Trans.frequency;
maxBufLength = ceil(sqrt(P.maxDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
maxBufSizePerAcq = 128*ceil(maxBufLength*8*(demodFreq/Trans.frequency)/128);

% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = P.numRays*maxBufSizePerAcq; % this size allows for all rays, with max range
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 12;    % should be multiple of 3
Resource.RcvBuffer(2).datatype = 'int16';
Resource.RcvBuffer(2).rowsPerFrame = P.numRays*maxBufSizePerAcq; % this size allows for all rays, with max range
Resource.RcvBuffer(2).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(2).numFrames = 12;    % should be multiple of 3
Resource.RcvBuffer(3).datatype = 'int16';
Resource.RcvBuffer(3).rowsPerFrame = P.numRays*maxBufSizePerAcq; % this size allows for all rays, with max range
Resource.RcvBuffer(3).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(3).numFrames = 12;    % should be multiple of 3
Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).numFrames = 1;  % one intermediate buffer needed.
Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).numFrames = 10;
Resource.DisplayWindow(1).Title = 'GEL3-12DWideBeamSC';
Resource.DisplayWindow(1).pdelta = 0.40;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
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
                   'aperture', 1, ...
                   'Apod', zeros(1,Resource.Parameters.numTransmit), ...
                   'focus', P.txFocus, ...
                   'Steer', [0.0,0.0], ...
                   'Delay', zeros(1,Resource.Parameters.numTransmit), ...
                   'TXPD', [], ...
                   'peakCutOff', 0.5, ...
                   'peakBLMax', 15.0), 1, 3*P.numRays);

if strcmp(Trans.units, 'wavelengths')
    scaleToWvl = 1;
end

% - Set event specific TX attributes.
for n = 1:P.numRays  % specify P.numRays transmit events
    % Set transmit Origins.
    TX(n).Origin(1) = TxOrgX(n);
    % Compute available transmit mux aperture
    [Dummy,ce] = min(abs(scaleToWvl*Trans.ElementPos(:,1)-TxOrgX(n))); % ce is closest element to center of aperture.
    lft = round(ce - 64);
    if lft < 1, lft = 1; end
    if lft > 129, lft = 129; end
    TX(n).aperture = lft;
    % Compute TX.Apod within mux aperture.
    lft = round(ce - P.numTx/2);
    if lft < 1, lft = 1; end
    rt = round(ce + P.numTx/2);
    if rt > 256, rt = 256; end
    TX(n).Apod((lft-(TX(n).aperture-1)):(rt-(TX(n).aperture-1))) = 1;
    % Apply apodization function.
    [RIndices,CIndices,V] = find(TX(n).Apod);
    V = kaiser(size(V,2),1);
    TX(n).Apod(CIndices) = V;
end
m = P.numRays;
for n = 1:P.numRays
    TX(n+m).Origin(1) = TX(n).Origin(1);
    TX(n+m).aperture = TX(n).aperture;
    TX(n+m).Apod = TX(n).Apod;
    if n<=10
        TX(n+m).Steer = [-((n-1)/10)*P.dtheta,0.0];
    elseif n>(P.numRays-10)
        TX(n+m).Steer = [-((P.numRays-n)/10)*P.dtheta,0.0];
    else
        TX(n+m).Steer = [-P.dtheta,0.0];
    end
end
m = m + P.numRays;
for n = 1:P.numRays
    TX(n+m).Origin(1) = TX(n).Origin(1);
    TX(n+m).aperture = TX(n).aperture;
    TX(n+m).Apod = TX(n).Apod;
    if n<=10
        TX(n+m).Steer = [((n-1)/10)*P.dtheta,0.0];
    elseif n>(P.numRays-10)
        TX(n+m).Steer = [((P.numRays-n)/10)*P.dtheta,0.0];
    else
        TX(n+m).Steer = [P.dtheta,0.0];
    end
end

% calculate TXPD and TX delay
h = waitbar(0,'Program TX parameters, please wait!');
steps = 3*P.numRays;
for ind = 1:steps
    TX(ind).Delay = computeTXDelays(TX(ind));
    TX(ind).TXPD = computeTXPD(TX(ind),PData);
    waitbar(ind/steps)
end
close(h)

% Specify TGC Waveform structure.
TGC.CntrlPts = [0,271,498,617,767,903,1000,1023];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Receive structure arrays

maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
Receive = repmat(struct('aperture', 1, ...
                        'Apod', zeros(1,128), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 0),1,3*P.numRays*Resource.RcvBuffer(1).numFrames);

% - Set event specific Receive attributes for each frame.
for i = 1:Resource.RcvBuffer(1).numFrames
    k = 3*P.numRays*(i-1);
    Receive(k+1).callMediaFunc = 1;
    for j = 1:P.numRays
        Receive(k+j).aperture = TX(j).aperture; % mux aperture same as transmit
        Receive(k+j).Apod(1:128) = 1.0;
        Receive(k+j).framenum = i;
        Receive(k+j).acqNum = j;
        Receive(k+P.numRays+j).aperture = TX(j).aperture; % mux aperture same as transmit
        Receive(k+P.numRays+j).Apod(1:128) = 1.0;
        Receive(k+P.numRays+j).bufnum = 2;
        Receive(k+P.numRays+j).framenum = i;
        Receive(k+P.numRays+j).acqNum = j;
        Receive(k+2*P.numRays+j).aperture = TX(j).aperture; % mux aperture same as transmit
        Receive(k+2*P.numRays+j).Apod(1:128) = 1.0;
        Receive(k+2*P.numRays+j).bufnum = 3;
        Receive(k+2*P.numRays+j).framenum = i;
        Receive(k+2*P.numRays+j).acqNum = j;
    end
end

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
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'scaleFactor', 1.0, ...
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
                         'pgain',1.0,...     % pgain is image processing gain
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
SeqControl(2).command = 'timeToNextAcq';  % time between ray line acquisitions
SeqControl(2).argument = 200;  % 200 usec
SeqControl(3).command = 'timeToNextAcq';
SeqControl(3).argument = 20000; %20000 usec = 20msec time between frames
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
    if floor(i/3) == i/3     % Exit to Matlab every 3rd frame
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
import vsv.seq.uicontrol.VsSliderControl

% - Sensitivity Cutoff
UI(1).Control = VsSliderControl('LocationCode','UserB7',...
    'Label','Sens. Cutoff',...
    'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
    'SliderStep',[0.025,0.1],'ValueFormat','%1.3f',...
    'Callback',@SensCutOffCallback);

% - Range Change
scaleToWvl = 1;
AxesUnit = 'mm';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if ~strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
        AxesUnit = 'wls';
        scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);
    end
end
UI(2).Control = VsSliderControl('LocationCode','UserA1',...
    'Label',['Range (',AxesUnit,')'],'SliderMinMaxVal',[5,P.maxDepthMm,P.endDepthMm]*scaleToWvl,...
    'SliderStep',[0.1,0.2],'ValueFormat','%3.0f',...
    'Callback',@RangeChangeCallback);

% - Transmit focus change
UI(3).Control = VsSliderControl('LocationCode','UserB4',...
    'Label',['TX Focus (',AxesUnit,')'],...
    'SliderMinMaxVal',[20,400,P.txFocusMm]*scaleToWvl,...
    'SliderStep',[0.1,0.2],'ValueFormat','%3.0f',...
    'Callback',@TxFocusCallback);

% - TX Aperture change
UI(4).Control = VsSliderControl('LocationCode','UserB3',...
    'Label','TX Aper',...
    'SliderMinMaxVal',[1,128,P.numTx],...
    'SliderStep',[1/128,1/32],'ValueFormat','%3.0f',...
    'Callback',@ApertureCallback);

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 3;
% Save all the structures to a .mat file.
save('MatFiles/GEL3-12DWideBeamSC');

%% Automatic VSX Execution:
% Uncomment the following line to automatically run VSX every time you run
% this SetUp script (note that if VSX finds the variable 'filename' in the
% Matlab workspace, it will load and run that file instead of prompting the
% user for the file to be used):

% filename = 'GEL3-12DWideBeamSC';  VSX;

% **** Callback routines used by UIControls ****
function SensCutOffCallback(~,~,UIValue)
    ReconL = evalin('base', 'Recon');
    for i = 1:size(ReconL,2)
        ReconL(i).senscutoff = UIValue;
    end
    assignin('base','Recon',ReconL);
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'Recon'};
    assignin('base','Control', Control);
end
 
function RangeChangeCallback(hObject,~,UIValue)
    simMode = evalin('base','Resource.Parameters.simulateMode');
    %No range change if in simulate mode 2.
    if simMode == 2
        set(hObject,'Value',evalin('base','P.endDepth'));
        return
    end
    Trans = evalin('base','Trans');
    Resource = evalin('base','Resource');
    scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);

    P = evalin('base','P');
    P.endDepth = UIValue;
    if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
        if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
            P.endDepth = UIValue*scaleToWvl;
        end
    end
    assignin('base','P',P);

    PData = evalin('base','PData');
    PData.Size(1) = ceil((P.endDepth-P.startDepth)/PData.PDelta(3)); % P.startDepth, P.endDepth and pdelta set PData.Size(1).
    for i = 1:size(PData.Region,2)
        PData.Region(i).Shape.height = P.endDepth;
    end
    PData.Region = computeRegions(PData);
    assignin('base','PData',PData);
    %Update TXPD data of TX structures.
    TX = evalin('base','TX');
    h = waitbar(0,'Program TX parameters, please wait!');
    steps = 3*P.numRays;
    for ind = 1:steps
        TX(ind).TXPD = computeTXPD(TX(ind),PData);
        waitbar(ind/steps)
    end
    close(h)
    assignin('base','TX',TX);
    evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData.Size(1)*PData.PDelta(3)/Resource.DisplayWindow(1).pdelta);');
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
end

function TxFocusCallback(hObject,~,UIValue)
    simMode = evalin('base','Resource.Parameters.simulateMode');
    %No focus change if in simulate mode 2.
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

    %- Redefine event specific TX attributes for the new focus.
    TX = evalin('base', 'TX');
    PData = evalin('base','PData');
    h = waitbar(0,'Program TX parameters, please wait!');
    steps = 3*P.numRays;
    for ind = 1:steps
        %write new focus value to TX
        TX(ind).focus = P.txFocus;
        TX(ind).Delay = computeTXDelays(TX(ind));
        TX(ind).TXPD = computeTXPD(TX(ind),PData);
        waitbar(ind/steps)
    end
    close(h)
    assignin('base','TX',TX);
    %Set Control command to update TX
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'TX'};
    assignin('base','Control', Control);
end

function ApertureCallback(hObject,~,UIValue)
    simMode = evalin('base','Resource.Parameters.simulateMode');
    %No focus change if in simulate mode 2.
    if simMode == 2
        set(hObject,'Value',evalin('base','P.numTx'));
        return
    end
    Trans = evalin('base', 'Trans');
    P = evalin('base','P');
    P.numTx = UIValue;
    assignin('base','P',P);

    TX = evalin('base', 'TX');
    scaleToWvl = evalin('base','scaleToWvl');
    TxOrgX = evalin('base','TxOrgX');
    PData = evalin('base','PData');
    %- Redefine event specific TX attributes for the new aperture.
    h = waitbar(0,'Program TX parameters, please wait!');
    steps = 3*P.numRays;
    for n = 1:P.numRays  % specify P.numRays transmit events
        TX(n).Origin(1) = TxOrgX(n);
        %Compute available transmit mux aperture
        [~,ce] = min(abs(scaleToWvl*Trans.ElementPos(:,1)-TxOrgX(n))); % ce is closest element to center of aperture.
        lft = round(ce - 64);
        if lft < 1, lft = 1; end
        if lft > 129, lft = 129; end
        TX(n).aperture = lft;
        %Compute TX.Apod within mux aperture.
        lft = round(ce - P.numTx/2);
        if lft < 1, lft = 1; end
        rt = round(ce + P.numTx/2);
        if rt > 256, rt = 256; end
        TX(n).Apod((lft-(TX(n).aperture-1)):(rt-(TX(n).aperture-1))) = 1;
        %Apply apodization function.
        [~,CIndices,V] = find(TX(n).Apod);
        V = kaiser(size(V,2),1);
        TX(n).Apod(CIndices) = V;
        %Compute transmit delays
        TX(n).Delay = computeTXDelays(TX(n));
        %Compute transmit pixel data
        TX(n).TXPD = computeTXPD(TX(n),PData);
        waitbar((n+P.numRays)/steps);
    end
    m = P.numRays;
    for n = 1:P.numRays
        TX(n+m).Origin(1) = TX(n).Origin(1);
        TX(n+m).aperture = TX(n).aperture;
        TX(n+m).Apod = TX(n).Apod;
        if n<=10
            TX(n+m).Steer = [-((n-1)/10)*P.dtheta,0.0];
        elseif n>(P.numRays-10)
            TX(n+m).Steer = [-((P.numRays-n)/10)*P.dtheta,0.0];
        else
            TX(n+m).Steer = [-P.dtheta,0.0];
        end
        %Compute transmit delays
        TX(n+m).Delay = computeTXDelays(TX(n+m));
        %Compute transmit pixel data
        TX(n+m).TXPD = computeTXPD(TX(n+m),PData);
        waitbar((n+P.numRays)/steps);
    end
    m = m + P.numRays;
    for n = 1:P.numRays
        TX(n+m).Origin(1) = TX(n).Origin(1);
        TX(n+m).aperture = TX(n).aperture;
        TX(n+m).Apod = TX(n).Apod;
        if n<=10
            TX(n+m).Steer = [((n-1)/10)*P.dtheta,0.0];
        elseif n>(P.numRays-10)
            TX(n+m).Steer = [((P.numRays-n)/10)*P.dtheta,0.0];
        else
            TX(n+m).Steer = [P.dtheta,0.0];
        end
        %Compute transmit delays
        TX(n+m).Delay = computeTXDelays(TX(n+m));
        %Compute transmit pixel data
        TX(n+m).TXPD = computeTXPD(TX(n+m),PData);
        waitbar((n+2*P.numRays)/steps);
    end
    close(h)
    assignin('base','TX', TX);
    %Set Control command to update TX
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'TX'};
    assignin('base','Control', Control);
end

