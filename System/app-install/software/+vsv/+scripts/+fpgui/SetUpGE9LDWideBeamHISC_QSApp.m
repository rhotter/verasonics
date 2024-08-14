% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use.
%
% File Name: SetUpGE9LDvWideBeamHISC.m - Wide beam, spatial compounding with
%   harmonic imaging.
%
% Description:
%   Sequence programming file for GE9LD curved array using wide beam acquisitions
%   with pulse inversion transmits for each beam position. On transmit, only P.numTx
%   elements are used, with P.numTx/2 transmitters on each side of the center element
%   (where possible). The wide beam transmits are steered between successive frames -
%   left, straight ahead, right.
%   All 192 elements are used on receive, although the element sensitivity cutoff will
%   limit the useful aperture. The receive acquisitions use 200% bandwidth.
%
%   The acquisition is asynchronous with processing. Frames are acquired from each steering
%   direction sequentially, with each frame stored in a separate RcvBuffer. The processing uses
%   the most recent frame from each RcvBuffer to output successive frames from each steering
%   direction.  Three frames are then averaged with a sliding window to produce the display image.
%
%   This script requires a Vantage256 system with the GE Universal Transducer Adapter.
%
% Last update:
%   10/04/2017 - test with SW 3.3

%clear[ 	]+all

demodFreq = 6.9444;  % Closest realizable demodFrequency

P.numTx = 58;   % no. of elements in TX aperture. was 140
P.numRays = 64; % no. of rays in frame was 120
P.txFocusMm = 100; % focus in mm
P.startDepthMm = 2;  % startDepth in mm
P.endDepthMm = 50;
P.maxDepthMm = 90;  % maxDepth for RangeChange and RcvBuffer
P.dtheta = 8*(pi/180); % steering angle for beams

% Specify system parameters.
Resource.Parameters.speedOfSound = 1540;
Resource.Parameters.numTransmit = 256;  % number of transmit channels.
Resource.Parameters.numRcvChannels = 256;  % number of receive channels.
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Specify Trans structure array.
Trans.name = 'GE9LD';
Trans.units = 'wavelengths'; % required in Gen3 to prevent default to mm units
Trans = computeTrans(Trans);  % GEC1-6D transducer is 'known' transducer so we can use computeTrans.

% Convert mm to wavelength
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);
P.txFocus = P.txFocusMm*scaleToWvl;
P.startDepth = P.startDepthMm*scaleToWvl;  % startDepth in wavelength
P.endDepth = P.endDepthMm*scaleToWvl;
P.maxDepth = P.maxDepthMm*scaleToWvl;
maxBufLength = ceil(sqrt(P.maxDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
maxBufSizePerAcq = 128*ceil(maxBufLength*8*(demodFreq/Trans.frequency)/128); % demodulation frequency

% Specify PData structure array.
PData.PDelta = [1.0, 0, 0.5];
PData.Size(1) = ceil((P.endDepth-P.startDepth)/PData.PDelta(3)); % P.startDepth, P.endDepth and pdelta set PData.Size.
PData.Size(2) = ceil((Trans.numelements*Trans.spacing)/PData.PDelta(1));
PData.Size(3) = 1;      % single image page
PData.Origin = [-Trans.spacing*(Trans.numelements/2),0,P.startDepth]; % x,y,z of upper lft crnr.
% Define numRays rectangular regions, positioned around the origin of the wide beams.
TxOrgX = (-95.5*Trans.spacing):(191*Trans.spacing/(P.numRays-1)):(95.5*Trans.spacing); % x coords of beam centers
for n = 1:P.numRays
    PData.Region(n) = struct('Shape',struct('Name','Rectangle','Position',[TxOrgX(n),0,0],'width',24,'height',P.endDepth));
end
m = P.numRays;
% Define numRays steered left parallelogram regions, positioned around the origin of the wide beams.
for n = 1:P.numRays
    if n<=10
        angle = -((n-1)/10)*P.dtheta;
    elseif n>(P.numRays-10)
        angle = -((P.numRays-n)/10)*P.dtheta;
    else
        angle = -P.dtheta;
    end
    PData.Region(n+m) = struct('Shape',struct('Name','Parallelogram','Position',[TxOrgX(n),0,0],'width',24,'height',P.endDepth/cos(P.dtheta),'angle',angle));
end
m = m + P.numRays;
% Define numRays steered right parallelogram regions, positioned around the origin of the wide beams.
for n = 1:P.numRays
    if n<=10
        angle = ((n-1)/10)*P.dtheta;
    elseif n>(P.numRays-10)
        angle = ((P.numRays-n)/10)*P.dtheta;
    else
        angle = P.dtheta;
    end
    PData.Region(n+m) = struct('Shape',struct('Name','Parallelogram','Position',[TxOrgX(n),0,0],'width',24,'height',P.endDepth/cos(P.dtheta),'angle',angle));
end
PData.Region = computeRegions(PData);

% Specify Media object. 'pt1.m' script defines array of point targets.
pt1;
Media.attenuation = -0.5;
Media.function = 'movePoints';

% Specify Resources.
rowsPerFrame = P.numRays*maxBufSizePerAcq;
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = rowsPerFrame; % this size allows for all rays, max range
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = numRcvFrames;
Resource.RcvBuffer(2).datatype = 'int16';
Resource.RcvBuffer(2).rowsPerFrame = rowsPerFrame; % this size allows for all rays, with max range
Resource.RcvBuffer(2).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(2).numFrames = numRcvFrames;  % all RcvBuffers should have same no. of frames
Resource.RcvBuffer(3).datatype = 'int16';
Resource.RcvBuffer(3).rowsPerFrame = rowsPerFrame; % this size allows for all rays, with max range
Resource.RcvBuffer(3).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(3).numFrames = numRcvFrames;
Resource.InterBuffer.datatype = 'complex double';
Resource.InterBuffer.numFrames = 1;
Resource.ImageBuffer.datatype = 'double';
Resource.ImageBuffer.numFrames = 6;
Resource.DisplayWindow(1).Title = 'GE9LDWideBeamHISC';
Resource.DisplayWindow(1).pdelta = 0.30;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).numFrames = numImageFrames;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow.Colormap = gray(256);

% Specify Transmit waveform structure.
TW(1) = struct('type','parametric','Parameters',[3.5,.67,2,-1]);
TW(2) = struct('type','parametric','Parameters',[3.5,.67,2,1]);

% Specify TX structure array.
% - We need 6*P.numRays transmit specifications.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'Apod', zeros(1,Trans.numelements), ...
                   'focus', P.txFocus, ...
                   'Steer', [0.0,0.0], ...
                   'TXPD', [], ...
                   'peakCutOff', 1.5, ...
                   'peakBLMax', 5.0, ...
                   'Delay', zeros(1,Trans.numelements)), 1, 6*P.numRays);

% - Set event specific TX attributes.
%    P.numTx/2 is the number of elements to include on each side of the
%    center element, for the specified focus and sensitivity cutoff.
%    Thus the full transmit aperture will be P.numTx + 1 elements.
for n = 1:P.numRays   % 2*P.numRays transmit events
    % Set transmit Origins.
    TX(n).Origin = TxOrgX(n);
    ce = round(1+(n-1)*191/(P.numRays));
    lft = round(ce - P.numTx/2);
    if lft < 1, lft = 1; end
    rt = round(ce + P.numTx/2);
    if rt > Trans.numelements, rt = Trans.numelements; end
    TX(n).Apod(lft:rt) = 1;
    % Apply apodization function.
    [RIndices,CIndices,V] = find(TX(n).Apod);
    V = kaiser(size(V,2),1);
    TX(n).Apod(CIndices) = V;
    % Compute transmit delays
    TX(n).Delay = computeTXDelays(TX(n));
    % Inverted pulse TXs
    TX(n+P.numRays).waveform = 2;
    TX(n+P.numRays).Origin = TX(n).Origin;
    TX(n+P.numRays).Apod = TX(n).Apod;
    TX(n+P.numRays).Delay = TX(n).Delay;
end
m = 2*P.numRays;
% Steered left beams
for n = 1:P.numRays   % 2*P.numRays transmit events
    TX(n+m).Origin = TX(n).Origin;
    TX(n+m).Apod = TX(n).Apod;
    if n<=10
        TX(n+m).Steer = [-((n-1)/10)*P.dtheta,0.0];
    elseif n>(P.numRays-10)
        TX(n+m).Steer = [-((P.numRays-n)/10)*P.dtheta,0.0];
    else
        TX(n+m).Steer = [-P.dtheta,0.0];
    end
    TX(n+m).Delay = computeTXDelays(TX(n+m));
    % Inverted pulse TXs
    TX(n+m+P.numRays).waveform = 2;
    TX(n+m+P.numRays).Origin = TX(n+m).Origin;
    TX(n+m+P.numRays).Apod = TX(n+m).Apod;
    TX(n+m+P.numRays).Delay = TX(n+m).Delay;
end
m = m + 2*P.numRays;
% Steered right beams
for n = 1:P.numRays   % 2*P.numRays transmit events
    TX(n+m).Origin = TX(n).Origin;
    TX(n+m).Apod = TX(n).Apod;
    if n<=10
        TX(n+m).Steer = [((n-1)/10)*P.dtheta,0.0];
    elseif n>(P.numRays-10)
        TX(n+m).Steer = [((P.numRays-n)/10)*P.dtheta,0.0];
    else
        TX(n+m).Steer = [P.dtheta,0.0];
    end
    TX(n+m).Delay = computeTXDelays(TX(n+m));
    % Inverted pulse TXs
    TX(n+m+P.numRays).waveform = 2;
    TX(n+m+P.numRays).Origin = TX(n+m).Origin;
    TX(n+m+P.numRays).Apod = TX(n+m).Apod;
    TX(n+m+P.numRays).Delay = TX(n+m).Delay;
end

% Calculating TXPD
h = waitbar(0,'Program TX parameters, please wait!');
steps = 3*P.numRays;
for i = 1:3
    for j = 1:P.numRays
        ind = P.numRays*2*(i-1)+j;
        TX(ind).TXPD = computeTXPD(TX(ind),PData);
        TX(ind+P.numRays).TXPD = TX(ind).TXPD;
        waitbar((P.numRays*(i-1)+j)/steps)
    end
end
close(h)

LowPassCoef = ...
[-0.00156 +0.01019 -0.00931 -0.00894 +0.02786 -0.01614...
 -0.03119 +0.06653 -0.02164 -0.11957 +0.28329 +0.64093];

InputFilter = ...
[+0.00073 -0.00092 -0.00214 +0.00183 +0.00247 -0.00058 +0.00299 ...
 -0.00427 -0.01605 +0.00958 +0.02676 -0.00839 -0.01273 -0.00259 ...
 -0.04535 +0.01663 +0.14044 -0.02106 -0.23227 +0.00961 +0.27057];

% Specify Receive structure arrays.
% -- Compute the maximum receive path length, using the law of cosines.
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
% We need 2*3*numRays*Frames receives for pulse inversion, steering
% directions, number of rays and number of frames.
Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'demodFrequency', demodFreq, ...
                        'mode', 0, ...
                        'LowPassCoef', LowPassCoef,...
                        'InputFilter', InputFilter,...
                        'callMediaFunc', 0), 1, 6*P.numRays*Resource.RcvBuffer(1).numFrames);

m = P.numRays;
% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    k = 6*m*(i-1);
    Receive(k+1).callMediaFunc = 1;
    for j = 1:m
        % Buffer 1
        Receive(j+k).framenum = i;
        Receive(j+k).acqNum = j;
        Receive(j+k+m).framenum = i;
        Receive(j+k+m).acqNum = j;
        Receive(j+k+m).mode = 1;
        % Buffer 2
        Receive(j+k+2*m).bufnum = 2;
        Receive(j+k+2*m).framenum = i;
        Receive(j+k+2*m).acqNum = j;
        Receive(j+k+3*m).bufnum = 2;
        Receive(j+k+3*m).framenum = i;
        Receive(j+k+3*m).acqNum = j;
        Receive(j+k+3*m).mode = 1;
        % Buffer 3
        Receive(j+k+4*m).bufnum = 3;
        Receive(j+k+4*m).framenum = i;
        Receive(j+k+4*m).acqNum = j;
        Receive(j+k+5*m).bufnum = 3;
        Receive(j+k+5*m).framenum = i;
        Receive(j+k+5*m).acqNum = j;
        Receive(j+k+5*m).mode = 1;
    end
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [0,250,590,710,770,830,890,950];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Recon structure arrays.
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
    ReconInfo(j+k).txnum = 2*P.numRays+j;
    ReconInfo(j+k).rcvnum = 2*P.numRays+j;
    ReconInfo(j+k).regionnum = j+k;
end
ReconInfo(k+P.numRays).Post = 'IQ2IntensityImageBuf';
k = k + P.numRays;
ReconInfo(k+1).Pre = 'clearInterBuf';
for j = 1:P.numRays
    ReconInfo(j+k).txnum = 4*P.numRays+j;
    ReconInfo(j+k).rcvnum = 4*P.numRays+j;
    ReconInfo(j+k).regionnum = j+k;
end
ReconInfo(k+P.numRays).Post = 'IQ2IntensityImageBuf';


% Specify Process structure array.
%   First processing structure specifies how to display image.
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,... % frame number in src buffer (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure (defines output figure).
                         'norm',1,...        % normalization method(1 means fixed)
                         'pgain',1,...       % pgain is image processing gain
                         'reject',2,...
                         'persistMethod','simple',...
                         'persistLevel',20,...
                         'interpMethod','4pt',...
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','runAverage3',...
                         'compressMethod','power',...      % X^0.5 normalized to output word size
                         'compressFactor',45,...
                         'mappingMethod','full',...
                         'display',1,...      % display image after processing
                         'displayWindow',1};

% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; %  jump back to start.
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq'; % time between rays
SeqControl(2).argument = 100; % 100usec
SeqControl(3).command = 'timeToNextAcq';  % time between frames
SeqControl(3).argument = 1000;  %1ms
SeqControl(4).command = 'returnToMatlab';
nsc = 5;

% Specify Event structure arrays.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    % Acquire frame with straight ahead wide beams
    k = 6*m*(i-1);
    for j = 1:m     % Acquire frame
        Event(n).info = 'Acquire beam, normal TX.';
        Event(n).tx = j;
        Event(n).rcv = j+k;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;

        Event(n).info = 'Acquire beam, inverted TX.';
        Event(n).tx = j+m;
        Event(n).rcv = j+k+m;
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
    for j = 1:m        % Acquire frame
        Event(n).info = 'Acquire beam, normal TX.';
        Event(n).tx = j+2*m;
        Event(n).rcv = j+k+2*m;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;

        Event(n).info = 'Acquire beam, inverted TX.';
        Event(n).tx = j+3*m;
        Event(n).rcv = j+k+3*m;
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
    for j = 1:m        % Acquire frame
        Event(n).info = 'Acquire beam, normal TX.';
        Event(n).tx = j+4*m;
        Event(n).rcv = j+k+4*m;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;

        Event(n).info = 'Acquire beam, inverted TX.';
        Event(n).tx = j+5*m;
        Event(n).rcv = j+k+5*m;
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
    Event(n).seqControl = 4;
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
MinMaxMm = [20,P.maxDepthMm]; % min max in mm
if isfield(Resource.DisplayWindow(1),'AxesUnits') && strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
    AxesUnit = 'mm';
    MinMaxVal = [MinMaxMm, P.endDepth/scaleToWvl];
else
    AxesUnit = 'wls';
    MinMaxVal = [MinMaxMm*scaleToWvl, P.endDepth];
end
stepBase = MinMaxVal(2)-MinMaxVal(1);
UI(2).Control = VsSliderControl('LocationCode','UserA1',...
    'Label',['Range (',AxesUnit,')'],'SliderMinMaxVal',MinMaxVal,...
    'SliderStep',[5/stepBase,10/stepBase],'ValueFormat','%3.0f',...
    'Callback',@RangeChangeCallback);

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 1;

% Save all the structures to a .mat file.
filename = fullfile( vsv.file.getVSXDir, 'MatFiles/SetUpGE9LDWideBeamHISC_QSApp.mat');

save(filename);


% **** Callback routines used by UIControls ****
function SensCutOffCallback(~, ~, UIValue)
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

function RangeChangeCallback(hObject, ~, UIValue)
    simMode = evalin('base','Resource.Parameters.simulateMode');
    % No range change if in simulate mode 2.
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

    % Modify PData for new range
    PData = evalin('base','PData');
    PData.Size(1) = ceil((P.endDepth-P.startDepth)/PData.PDelta(3)); % P.startDepth, P.endDepth and pdelta set PData.Size.

    for n = 1:P.numRays
        PData.Region(n).Shape.height = P.endDepth;
    end
    m = P.numRays;
    for n = 1:P.numRays
        if n<=10
            angle = -((n-1)/10)*P.dtheta;
        elseif n>(P.numRays-10)
            angle = -((P.numRays-n)/10)*P.dtheta;
        else
            angle = -P.dtheta;
        end
        PData.Region(n+m).Shape.height = P.endDepth/cos(angle);
    end
    m = m + P.numRays;
    for n = 1:P.numRays
        if n<=10
            angle = ((n-1)/10)*P.dtheta;
        elseif n>(P.numRays-10)
            angle = ((P.numRays-n)/10)*P.dtheta;
        else
            angle = P.dtheta;
        end
        PData.Region(n+m).Shape.height = P.endDepth/cos(angle);
    end

    PData.Region = computeRegions(PData);
    assignin('base','PData',PData);
    evalin('base','PData(1).Region = computeRegions(PData(1));');
    evalin('base','Resource.DisplayWindow(1).Position(3) = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);');
    evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);');
    evalin('base','Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];');

    % Update TXPD data of TX structures.
    TX = evalin('base','TX');
    PData = evalin('base','PData');
    h = waitbar(0,'Program TX parameters, please wait!');
    steps = 3*P.numRays;
    for i = 1:3
        for j = 1:P.numRays
            ind = P.numRays*2*(i-1)+j;
            TX(ind).TXPD = computeTXPD(TX(ind),PData);
            TX(ind+P.numRays).TXPD = TX(ind).TXPD;
            waitbar((P.numRays*(i-1)+j)/steps)
        end
    end
    close(h)
    assignin('base','TX',TX);

    % Update Receive structures
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
    Control.Parameters = {'PData','InterBuffer','ImageBuffer','TX','Receive','TGC','Recon','DisplayWindow','TGC'};
    assignin('base','Control', Control);
    assignin('base', 'action', 'displayChange');
end

