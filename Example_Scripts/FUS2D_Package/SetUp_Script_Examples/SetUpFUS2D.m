% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use.
%
%   Currently, only one focus is predefined and different UI controls will
%   be generated based on the pairing number. If the HIFU array is annular,
%   i.e., FUS 2D pairing 01 to 03, only one slider with Z direction will
%   appear. If the FUS 2D pairing is 04 to 06, three sliders with Z, R
%   and Rotation Angle will be generated for focal postion adjustment.
%
%   Definition of three sliders:
%   Z is the axial location of the focus. The default value, 0, is the natual
%   focus. Positive value moves the focus deeper and negative value moves the
%   focus toward to the transducer. R is the lateral distance between the focus
%   and the axial axis. Rotation Angle is the angle between the focus and the
%   Bmode plane, which is also the XZ plane along the element 1 of the HIFU
%   array
%
% File name: SetUpFUS2D.m

clear all
delete(instrfind)

HIFU.PRF = 1;  % Hz
HIFU.Duty = 0.6; % DutyCycle

%% System Check
[pairingNum,TransDual,V128,VDAS,extPwrParam] = fus2DSystemCheck();

if ~isempty(extPwrParam)
    Resource.HIFU = extPwrParam;
end

if ismember(pairingNum,[1,2,3])
    colsPerFrame = 128;
    elementsUsedForImg = 64;
    ImgApod = ones(1,64);
elseif V128 % pairing 04 on V128 system only
    colsPerFrame = 128;
    elementsUsedForImg = 64;
    ImgApod = [ones(1,64),zeros(1,64)]; % Used for TXPD calculation
else % pairing 4-6 with V256
    colsPerFrame = 256;
    elementsUsedForImg = 128;
    ImgApod = ones(1,128);
end


% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 4;

numBmodeFrames = 2;
RcvBufferFrames = numBmodeFrames + mod(numBmodeFrames,2); % can't be an odd number

%% Define system parameters.
Resource.Parameters.Connector = [1,2]; % two connectors
% Note: Since SW 4.0, numTransmit and numRcvChannels are no longer required
Resource.Parameters.speedOfSound = 1540;
Resource.Parameters.speedCorrectionFactor = 1.0;
Resource.Parameters.verbose = 1;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

%% Specify Imaging Parameters
Trans.name = TransDual.Img;
Trans.units = 'wavelengths'; % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);

TransImg = Trans;
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);

% PData For imaging array
P(1).numRays = 24; % no. of raylines to program
P(1).startDepthMm = 0;
if ismember(pairingNum,[4,5,6])
    isAnnular = 0; ZSldrMinMax = [-20 20]; RSldrMinMax = [0 8];
    P(1).endDepthMm = 220;  % endDepth in mm
elseif ismember(pairingNum,[1,2,3])
    isAnnular = 1; ZSldrMinMax = [-25 25];
    P(1).endDepthMm = 100;  % endDepth in mm
end

P(1).startDepth = P(1).startDepthMm*scaleToWvl;
P(1).endDepth = P(1).endDepthMm*scaleToWvl;   % Acquisition depth in wavelengths
P(1).txFocus = -900;

scanangle = pi/4;
P(1).aperture = elementsUsedForImg*Trans.spacing; % aperture based on Trans.numelements
P(1).theta = -scanangle/2;
P(1).radius = (P(1).aperture/2)/tan(-P(1).theta); % dist. to virt. apex
P(1).rayDelta = scanangle/(P(1).numRays-1);

beamWidth = 0.13; % approximate beam width in radians
Angle = (P(1).theta-beamWidth):(2*(-P(1).theta+beamWidth)/(P(1).numRays-1)):(-P(1).theta+beamWidth);
pDataShift = 0;
if V128 && isequal(pairingNum,4)
    pDataShift = P(1).aperture/2;
end
PData(1).PDelta = [1.0, 0, 0.5];
PData(1).Size(1) = 10 + ceil((P(1).endDepth-P(1).startDepth)/PData(1).PDelta(3));
PData(1).Size(2) = 10 + ceil(2*(P(1).endDepth + P(1).radius)*sin(-P(1).theta)/PData(1).PDelta(1));
PData(1).Size(3) = 1;
PData(1).Origin = [-(PData(1).Size(2)/2+pDataShift)*PData(1).PDelta(1),0,P(1).startDepth];
PData(1).Region = repmat(struct(...
    'Shape',struct('Name','SectorFT', ...
    'Position',[-pDataShift,0,-P(1).radius], ...
    'z',P(1).startDepth, ...
    'r',P(1).radius+P(1).endDepth, ...
    'angle',P(1).rayDelta*10, ...
    'steer',0)),1,P(1).numRays+1);

% Define the first PData region as the entire scan region to use as a mask.
PData(1).Region(1).Shape.angle = scanangle;
for n = 2:P(1).numRays+1
    PData(1).Region(n).Shape.steer = Angle(n-1);%P(1).theta + (n-1)*P(1).rayDelta;
    PData(1).Region(n).Shape.andWithPrev = 1;
end
PData(1).Region = computeRegions(PData(1));

% Imaging pulse
TW(1).type = 'parametric';
TW(1).Parameters = [TransImg.frequency,.67,2,1];

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
    'Origin', [0.0,0.0,0.0], ...
    'focus', P(1).txFocus, ...
    'FocalPtMm', [],...
    'Steer', [0.0,0.0], ...
    'Apod', ImgApod, ...
    'Delay', zeros(1,TransImg.numelements), ...
    'TXPD', [], ...
    'peakCutOff', 1, ...
    'peakBLMax', 4), 1, P(1).numRays) ;

% - Set event specific TX attributes.
beamWidth = 0.15; % approximate beam width in radians
Angle = (P(1).theta-beamWidth):(2*(-P(1).theta+beamWidth)/(P(1).numRays-1)):(-P(1).theta+beamWidth);
h = waitbar(0,'Program TX parameters for Imaging pulse, please wait!');
for n = 1:P(1).numRays
    TX(n).Origin = [P(1).radius*tan(Angle(n)), 0.0, 0.0];
    TX(n).Steer = [Angle(n),0.0];
    TX(n).Delay = computeTXDelays(TX(n));
    TX(n).TXPD = computeTXPD(TX(n),PData(1),'Attenuation',-0.4);
    waitbar(n/P(1).numRays);
end
close(h)
TXImg = TX;

%% HIFU Parameters. Note that the frequency of imaging arrary should be used
clear Trans
Trans.name = TransDual.HIFU;
Trans.units = 'wavelengths';

% The delay calculation should use the wavelength of the imaging transducer
Trans.frequency = TransImg.frequency;
Trans = computeTrans(Trans);
TransHIFU = Trans;

% Get the center frequency for trnasmission
TransHIFU1 = rmfield(TransHIFU, 'frequency');
TransHIFU1 = computeTrans(TransHIFU1);
TransHIFU.frequency = TransHIFU1.frequency;

geoFocusMm = TransHIFU.radiusMm; % mm
FociPosMm = [0,0,geoFocusMm];
numFocus = size(FociPosMm,1);

% Define Rotate for slider callbacks
Rotate.z = geoFocusMm;
Rotate.r = 0;
Rotate.theta = 0;

% HIFU pulse
TW(2).type = 'parametric';
HIFU.PD = 1e3*HIFU.Duty/HIFU.PRF;% ms
C = HIFU.PD*1e3*TransHIFU.frequency*2;
TW(2).Parameters = [TransHIFU.frequency,.67,C,1];

% Specify TX structure array for flash transmit.
TXHIFU = repmat(struct('waveform', 2, ...
    'Origin', [0.0,0.0,0.0], ...
    'focus', 0, ...
    'Steer', [0.0,0.0], ...
    'Apod', ones(1,TransHIFU.numelements), ...
    'Delay', zeros(1,TransHIFU.numelements),...
    'TXPD', [],...
    'peakCutOff',0.5,...
    'peakBLMax',4), 1, numFocus);

for n = 1:numFocus
    TXHIFU(n).FocalPtMm = FociPosMm(n,:);
    TXHIFU(n).Delay = computeTXDelays(TXHIFU(n));
end
TX = orderfields(TX);
TX(P(1).numRays+1:P(1).numRays+numFocus) = orderfields(TXHIFU);

%% Combine Img and HIFU parameters, only few fields need to be changed
clear Trans
% Start from
Trans = TransHIFU;

% Use Image array for common fields
f = fieldnames(TransImg);
for i = 1:length(f)
    Trans.(f{i}) = TransImg.(f{i});
end

% Fields need to be set for combination of two transducers
Trans.name = 'custom';      % Must be 'custom' to prevent confusion from the two unique transducer ID's that will actually be connected
Trans.id = hex2dec('0000'); % Dummy ID to be used with the 'fake scanhead' feature
Trans.impedance = TransHIFU.impedance;  % Power estimation requires the impedance of the HIFU transducer
Trans.ElementPos = [TransImg.ElementPos(1:elementsUsedForImg,:); TransHIFU.ElementPos];
Trans.numelements = elementsUsedForImg+TransHIFU.numelements;    % total over both connectors

% Trans.Connector should be corrected based on the system channels
if V128
    if isequal(pairingNum,4)
        Trans.ConnectorES = [TransImg.ConnectorES(1:64); TransHIFU.ConnectorES];
        Trans.numelements = 128;
    else
        Trans.ConnectorES = [repmat(TransImg.ConnectorES,1,size(TransHIFU.ConnectorES,2)); TransHIFU.ConnectorES];
    end
else
    Trans.ConnectorES = [repmat(TransImg.ConnectorES,1,size(TransHIFU.ConnectorES,2)); TransHIFU.ConnectorES + 128];
end

% Set TX parameters correctly for 2 connectors
for n = 1:P(1).numRays
    TX(n).Apod  = [TXImg(n).Apod(1:elementsUsedForImg), zeros(1,TransHIFU.numelements)];
    TX(n).Delay = [TXImg(n).Delay(1:elementsUsedForImg), zeros(1,TransHIFU.numelements)];
end

for n = P(1).numRays+1:P(1).numRays+numFocus
    TX(n).Apod  = [zeros(1,elementsUsedForImg),TXHIFU(n-P(1).numRays).Apod];
    TX(n).Delay = [zeros(1,elementsUsedForImg),TXHIFU(n-P(1).numRays).Delay];
end

% Specify Media object.
pt1;
Media.MP(:,3) = Media.MP(:,3)+100;
Media.attenuation = -0.5;
Media.function = 'movePoints';

%% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = P(1).numRays*4096*2;
Resource.RcvBuffer(1).colsPerFrame = colsPerFrame; % either 128 or 256
Resource.RcvBuffer(1).numFrames = RcvBufferFrames;
Resource.InterBuffer(1).numFrames = 1;
Resource.ImageBuffer(1).numFrames = 10;
Resource.DisplayWindow(1).Title = [TransDual.Img,'WideBeam'];
Resource.DisplayWindow(1).pdelta = 0.6;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
    DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).numFrames = 20;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow.Colormap = gray(256);

% TPC
TPC(1).maxHighVoltage = Trans.maxHighVoltage;
TPC(5).maxHighVoltage = TransHIFU.maxHighVoltage;
TPC(5).hv = 4;

%% Specify Receive structure arrays.
maxAcqLength = 1.5*ceil(sqrt(P(1).aperture^2 + P(1).endDepth^2 - 2*P(1).aperture*P(1).endDepth*cos(P(1).theta-pi/2)) - P(1).startDepth);
wlsPer128 = 128/(4*2); % wavelengths in 128 samples for 4 samplesPerWave

Receive = repmat(struct('Apod', [ones(1,elementsUsedForImg),zeros(1,TransHIFU.numelements)], ...
    'startDepth', P(1).startDepth, ...
    'endDepth', P(1).startDepth + wlsPer128*ceil(maxAcqLength/wlsPer128), ...
    'TGC', 1, ...
    'bufnum', 1, ...
    'framenum', 1, ...
    'acqNum', 1, ...
    'sampleMode', 'NS200BW', ...
    'mode', 0, ...
    'callMediaFunc', 0),1,P(1).numRays*Resource.RcvBuffer(1).numFrames);
% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    Receive(P(1).numRays*(i-1)+1).callMediaFunc = 1;
    for j = 1:P(1).numRays
        Receive(P(1).numRays*(i-1)+j).framenum = i;
        Receive(P(1).numRays*(i-1)+j).acqNum = j;
    end
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [443,600,740,811,822,875,971,1023];
TGC.rangeMax = P(1).endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

%% Specify Recon and Process structure arrays.
Recon = struct('senscutoff', 0.5, ...
    'newFrameTimeout', 1e5,...
    'pdatanum', 1, ...
    'rcvBufFrame', -1, ...
    'IntBufDest', [1,1], ...
    'ImgBufDest', [1,-1], ...
    'RINums', 1:P(1).numRays);

% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 'accumIQ', ...  % accumulate IQ
    'Pre',[],...
    'Post',[],...
    'txnum', 1, ...
    'rcvnum', 1, ...
    'scaleFactor', 0.2, ...
    'threadSync', 1,...
    'regionnum', 0), 1, P(1).numRays);

% - Set specific ReconInfo attributes.
ReconInfo(1).Pre = 'clearInterBuf';     % clear entire Interbuffer frame
for j = 1:P(1).numRays
    ReconInfo(j).txnum = j;
    ReconInfo(j).rcvnum = j;
    ReconInfo(j).regionnum = j+1;
end
ReconInfo(P(1).numRays).Post = 'IQ2IntensityImageBuf';  % detect IQ data

% Specify Process structure array.
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
    'framenum',-1,...   % (-1 => lastFrame)
    'pdatanum',1,...    % number of PData structure to use
    'pgain',5,...            % pgain is image processing gain
    'reject',2,...      % reject level
    'persistMethod','simple',...
    'persistLevel',30,...
    'interpMethod','4pt',...
    'grainRemoval','none',...
    'processMethod','none',...
    'averageMethod','none',...
    'compressMethod','power',...
    'compressFactor',40,...
    'mappingMethod','full',...
    'display',1,...      % display image after processing
    'displayWindow',1};

%% Specify SeqControl structure arrays.  Missing fields are set to NULL.
TTNA1 = 1;
SeqControl(TTNA1).command = 'timeToNextAcq';
SeqControl(TTNA1).argument = 500;

TTNA2 = 2;
SeqControl(TTNA2).command = 'timeToNextAcq';
SeqControl(TTNA2).argument = 20000 - (P(1).numRays-1)*SeqControl(TTNA1).argument;

TTNEB = 3;
TTNAHIFU = 1e6/HIFU.PRF;
SeqControl(TTNEB).command = 'timeToNextEB';
SeqControl(TTNEB).argument = TTNAHIFU;

RTNMatlab = 4;
SeqControl(RTNMatlab).command = 'returnToMatlab';

SYNC = 5;
SeqControl(SYNC).command = 'sync';

JUMP = 6;
SeqControl(JUMP).command = 'jump';
SeqControl(JUMP).argument = 1;
nsc = 7; % nsc is count of SeqControl objects

%% Event Sequence
n = 1;
% Bmode Frame
for tn = 1:numFocus
    rcvNum = 1;
    for fn = 1:numBmodeFrames
        % Initial TPC profile 5
        Event(n).info = 'select TPC profile 1 to start';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = nsc; % set TPC profile command.
        n = n+1;
        SeqControl(nsc).command = 'setTPCProfile';
        SeqControl(nsc).argument = 1;
        SeqControl(nsc).condition = 'immediate';
        nsc = nsc + 1;

        for j = 1:P(1).numRays                 % Acquire frame
            Event(n).info = 'Acquire ray line.';
            Event(n).tx = j;
            Event(n).rcv = rcvNum;
            Event(n).recon = 0;
            Event(n).process = 0;
            Event(n).seqControl = TTNA1;
            n = n+1;
            rcvNum = rcvNum +1;
        end
        Event(n-1).seqControl = nsc; % modify last acquisition Event's seqControl
        SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
        nsc = nsc+1;

        Event(n).info = 'Reconstruct and process';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 1;
        Event(n).process = 1;
        Event(n).seqControl = 0;
        if floor(i/frameRateFactor) == i/frameRateFactor
            Event(n).seqControl = RTNMatlab;
        end
        n = n+1;

        % Initial TPC profile 5
        Event(n).info = 'select TPC profile 5';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = nsc; % set TPC profile command.
        n = n+1;
        SeqControl(nsc).command = 'setTPCProfile';
        SeqControl(nsc).argument = 5;
        SeqControl(nsc).condition = 'immediate';
        nsc = nsc + 1;

        Event(n).info = 'noop delay';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = nsc; % noop to allow time for TPC profile transition.
        n = n+1;
        SeqControl(nsc).command = 'noop';
        SeqControl(nsc).argument = 5e3;
        nsc = nsc + 1;

        Event(n).info = 'HIFU Pulse';
        Event(n).tx = tn+P(1).numRays;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = [TTNEB,RTNMatlab];
        n = n+1;
    end
end

Event(n).info = 'Jump back to first event';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = JUMP;

% User specified UI Control Elements
% - Sensitivity Cutoff
UI(1).Control =  {'UserB7','Style','VsSlider','Label','Sens. Cutoff',...
    'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
    'SliderStep',[0.025,0.1],'ValueFormat','%1.3f'};
UI(1).Callback = text2cell('%SensCutOffCallback');

% - Z
UI(2).Control =  {'UserB4','Style','VsSlider','Label','Z (mm)',...
    'SliderMinMaxVal',[ZSldrMinMax,0],...
    'SliderStep',[1/diff(ZSldrMinMax),5/diff(ZSldrMinMax)],'ValueFormat','%2.1f'};
UI(2).Callback = text2cell('%ZCallback');

if ~isAnnular
    % - R
    UI(3).Control =  {'UserB3','Style','VsSlider','Label','R (mm)',...
        'SliderMinMaxVal',[RSldrMinMax,0],...
        'SliderStep',[1/diff(RSldrMinMax),2/diff(RSldrMinMax)],'ValueFormat','%2.1f'};
    UI(3).Callback = text2cell('%RCallback');

    % - Rotation Angle
    UI(4).Control =  {'UserB2','Style','VsSlider','Label','Rotation Angle',...
        'SliderMinMaxVal',[-180,180,0],...
        'SliderStep',[1/360,5/360],'ValueFormat','%3.0f'};
    UI(4).Callback = text2cell('%RotationAngleCallback');
end

% Save all the structures to a .mat file.
matFilename = ['FUS2D_Pairing0',num2str(pairingNum)];
save(fullfile('MatFiles', matFilename));

return

% **** Callback routines to be converted by text2cell function. ****s
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

%ZCallback
Rotate = evalin('base','Rotate');
geoFocusMm = evalin('base','geoFocusMm');
Rotate.z = geoFocusMm + UIValue;
r = Rotate.r;
TX = evalin('base','TX');
TXHIFU = evalin('base','TXHIFU');

TransFinal = evalin('base','Trans');
TransImg = evalin('base','TransImg');
elementsUsedForImg = evalin('base','elementsUsedForImg');
% Delay calculatioin needs corret Trans from HIFU array but with imaging
% array's center frequency
evalin('base','Trans = TransHIFU;');
evalin('base','Trans.frequency = TransImg.frequency;');
TXHIFU.FocalPtMm = [r*cos(Rotate.theta), r*sin(Rotate.theta), Rotate.z];
TXHIFU.Delay = computeTXDelays(TXHIFU);
TX(end).Delay(elementsUsedForImg+1:end) = TXHIFU.Delay;

assignin('base','Trans',TransFinal);
assignin('base','Rotate',Rotate);
assignin('base','TX',TX);
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'TX'};
assignin('base','Control', Control);

return
%ZCallback

%RCallback
Rotate = evalin('base','Rotate');
Rotate.r = UIValue;
r = Rotate.r;
TX = evalin('base','TX');
TXHIFU = evalin('base','TXHIFU');

TransFinal = evalin('base','Trans');
TransImg = evalin('base','TransImg');
elementsUsedForImg = evalin('base','elementsUsedForImg');

% Delay calculatioin needs corret Trans from HIFU array but with imaging
% array's center frequency
evalin('base','Trans = TransHIFU;');
evalin('base','Trans.frequency = TransImg.frequency;');
TXHIFU.FocalPtMm = [r*cos(Rotate.theta), r*sin(Rotate.theta), Rotate.z];
TXHIFU.Delay = computeTXDelays(TXHIFU);
TX(end).Delay(elementsUsedForImg+1:end) = TXHIFU.Delay;

assignin('base','Trans',TransFinal);
assignin('base','Rotate',Rotate);
assignin('base','TX',TX);
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'TX'};
assignin('base','Control', Control);

return
%RCallback

%RotationAngleCallback
Rotate = evalin('base','Rotate');
Rotate.theta = UIValue*pi/180;
r = Rotate.r;
TX = evalin('base','TX');
TXHIFU = evalin('base','TXHIFU');

TransFinal = evalin('base','Trans');
TransImg = evalin('base','TransImg');
elementsUsedForImg = evalin('base','elementsUsedForImg');

% Delay calculatioin needs corret Trans from HIFU array but with imaging
% array's center frequency
evalin('base','Trans = TransHIFU;');
evalin('base','Trans.frequency = TransImg.frequency;');
TXHIFU.FocalPtMm = [r*cos(Rotate.theta), r*sin(Rotate.theta), Rotate.z];
TXHIFU.Delay = computeTXDelays(TXHIFU);
TX(end).Delay(elementsUsedForImg+1:end) = TXHIFU.Delay;

assignin('base','Trans',TransFinal);
assignin('base','Rotate',Rotate);
assignin('base','TX',TX);
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'TX'};
assignin('base','Control', Control);

return
%RotationAngleCallback
