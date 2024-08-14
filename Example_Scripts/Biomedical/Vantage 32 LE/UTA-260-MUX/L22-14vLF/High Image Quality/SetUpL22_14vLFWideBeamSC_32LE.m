% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use.
%
% File name SetUpL22_14vLFWideBeamSC_32LE.m: example of linear array wide beam
% transmit with spatial compounding
%
% Generate .mat Sequence Object file for L22-14vLF Linear array wide beam transmit.
% The transmit aperture size can be set by the variable - numTx. The transmit aperture is
% translated across the array with the number of rays set by the variable - numRays.
% The receive aperture always covers 32 element aperture.
%
% In the following aperture examples, each space represents 2 elements.
% Aperture 1 (P.numTx=32):
%   tttttttttttttttt------------------------------------------------
%   rrrrrrrrrrrrrrrr------------------------------------------------
% Aperture for almost middle wide beam 24 (P.numTx=32, P.numRays=48):
%   ------------------------tttttttttttttttt------------------------
%   ------------------------rrrrrrrrrrrrrrrr------------------------
% Aperture for last wide beam (P.numTx=32, P.numRays=48):
%   ------------------------------------------------tttttttttttttttt
%   ------------------------------------------------rrrrrrrrrrrrrrrr
%
% Last update:
% 09-03-2018 modified for 32LE system.

clear all

P.startDepth = 5;
P.endDepth = 160;
P.numRays = 48;     % no. of rays in frame
P.txFocus = 3*P.endDepth;
P.numTx = 32;       % no. of elements in TX aperture (maximum of 64).
P.dtheta = 10*(pi/180);  % angle delta for beams steered left or right

% Define system parameters.
Resource.Parameters.speedOfSound = 1540;  % speed of sound in m/s
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.
% Resource.System.Product = 'Vantage64';
Resource.System.SoftwareVersion = [3 5 0]; % Minimum software release for this script.
Resource.System.UTA = '260-MUX'; % This script requires the 260-MUX UTA.
Resource.System.transmitChannels = 64;
Resource.System.receiveChannels = 32;

% Specify Trans structure array.
Trans.name = 'L22-14v LF';
Trans.units = 'wavelengths';
Trans = computeTrans(Trans);
Trans.maxHighVoltage = 30;  % set maximum high voltage limit for pulser supply.
Trans = computeUTAMux64(Trans); % Add HVMux field for use with UTA 260-Mux

scaleToWvl = 1;
if strcmp(Trans.units, 'mm')
    scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);
end

% Specify PData structure arrays.
PData.PDelta = [0.75*Trans.spacing, 0, 0.4];
PData.Size(1,1) = ceil((P.endDepth-P.startDepth)/PData.PDelta(3)); % range, origin and pdelta set PData.Size.
PData.Size(1,2) = ceil((Trans.numelements*Trans.spacing)/PData.PDelta(1));
PData.Size(1,3) = 1;      % single image page
PData.Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of uppr lft crnr.
% Predefine all the Region structures, then modify them below.
PData.Region = repmat(struct('Shape',struct('Name','Parallelogram',...
                                               'Position',[0,0,P.startDepth],...
                                               'width',Trans.spacing*(0.55*P.numTx),...
                                               'height',P.endDepth-P.startDepth,...
                                               'angle',0.0)),1,3*P.numRays);
% Compute the x coords of the TX beam centers
TxOrgX = (-63.5*Trans.spacing):(127*Trans.spacing/(P.numRays-1)):(63.5*Trans.spacing);
for n = 1:P.numRays
    PData.Region(n).Shape.Position(1) = TxOrgX(n);
    PData.Region(P.numRays+n).Shape.Position(1) = TxOrgX(n);
    PData.Region(2*P.numRays+n).Shape.Position(1) = TxOrgX(n);
end
% Define P.numRays steered left parallelogram regions, centered on TX beam origins. Adjust the angle
%   so that the steering goes to zero over 8 beams at the left and right edge.
m = P.numRays;
for n = 1:P.numRays
    if n<=8
        angle = -((n-1)/8)*P.dtheta;
    elseif n>(P.numRays-8)
        angle = -((P.numRays-n)/8)*P.dtheta;
    else
        angle = -P.dtheta;
    end
    PData.Region(n+m).Shape.height = (P.endDepth-P.startDepth)/cos(angle);
    PData.Region(n+m).Shape.angle = angle;
end
% Define numRays steered right parallelogram regions, centered on TX beam origins. Adjust the angle
%   so that the steering goes to zero over 8 beams at the left and right edge.
m = m + P.numRays;
for n = 1:P.numRays
    if n<=8
        angle = ((n-1)/8)*P.dtheta;
    elseif n>(P.numRays-8)
        angle = ((P.numRays-n)/8)*P.dtheta;
    else
        angle = P.dtheta;
    end
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
Resource.RcvBuffer(1).rowsPerFrame = P.numRays*4096;
Resource.RcvBuffer(1).colsPerFrame = Resource.System.receiveChannels;
Resource.RcvBuffer(1).numFrames = 20;    % 20 frames stored in RcvBuffer.
Resource.RcvBuffer(2).datatype = 'int16';
Resource.RcvBuffer(2).rowsPerFrame = P.numRays*4096;
Resource.RcvBuffer(2).colsPerFrame = Resource.System.receiveChannels;
Resource.RcvBuffer(2).numFrames = 20;    % 20 frames stored in RcvBuffer.
Resource.RcvBuffer(3).datatype = 'int16';
Resource.RcvBuffer(3).rowsPerFrame = P.numRays*4096;
Resource.RcvBuffer(3).colsPerFrame = Resource.System.receiveChannels;
Resource.RcvBuffer(3).numFrames = 20;    % 20 frames stored in RcvBuffer.
Resource.InterBuffer(1).numFrames = 1;   % one intermediate buffer needed.
Resource.ImageBuffer(1).numFrames = 10;
Resource.DisplayWindow(1).Title = 'L22-14vLFWideBeamSC_32LE';
Resource.DisplayWindow(1).pdelta = 0.25;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData.Size(2)*PData.PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData.Size(1)*PData.PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData.Origin(1),0,PData.Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).numFrames = 20;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = gray(256);

% Specify Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [18,0.67,3,1];

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'Apod', zeros(1,Resource.System.transmitChannels), ...
                   'focus', P.txFocus, ...
                   'Steer', [0.0,0.0], ...
                   'Delay', zeros(1,Resource.System.transmitChannels), ...
                   'TXPD', [], ...
                   'peakCutOff', 0.8, ...
                   'peakBLMax', 4.0), 1, 3*P.numRays);

% - Set event specific TX attributes.  We will specify a MUX aperture the same size
%   as the number of TX & Receive channels, even though the number of active TX channels
%   set by P.numTx is usually smaller.  This will allow the same MUX aperture to be used
%   for receive, eliminating the time for the MUX to switch to a different aperture.
Ce = fix(1:(127/(P.numRays-1)):128); % calculate center element numbers.
for n = 1:P.numRays
    TX(n).Origin(1) = TxOrgX(n);
    % Now compute TX(n).Apod based on P.numTX for straight ahead beams
    lft = Ce(n) - floor(P.numTx/2);
    if lft < 1, lft = 1; end
    rt = Ce(n) + floor(P.numTx/2);
    if rt > Trans.numelements, rt = Trans.numelements; end
    TX(n).aperture = min(lft,65);
    TX(n).Apod((lft:rt)-TX(n).aperture+1) = 1;
    % - Apply apodization function.
    [RIndices,CIndices,V] = find(TX(n).Apod);
    V = kaiser(size(V,2),1);
    TX(n).Apod(CIndices) = V;
    % - Compute transmit delays
    TX(n).Delay = computeTXDelays(TX(n));
    % Compute TX for steered left beams
    TX(P.numRays+n).Origin(1) = TX(n).Origin(1);
    TX(P.numRays+n).Apod = TX(n).Apod;
    TX(P.numRays+n).aperture = TX(n).aperture;
    if n<=8
        TX(P.numRays+n).Steer = [-((n-1)/8)*P.dtheta,0.0];
    elseif n>(P.numRays-8)
        TX(P.numRays+n).Steer = [-((P.numRays-n)/8)*P.dtheta,0.0];
    else
        TX(P.numRays+n).Steer = [-P.dtheta,0.0];
    end
    % - Compute transmit delays
    TX(P.numRays+n).Delay = computeTXDelays(TX(P.numRays+n));
    % Compute TX for steered right beams
    TX(2*P.numRays+n).Origin(1) = TX(n).Origin(1);
    TX(2*P.numRays+n).Apod = TX(n).Apod;
    TX(2*P.numRays+n).aperture = TX(n).aperture;
    if n<=8
        TX(2*P.numRays+n).Steer = [((n-1)/8)*P.dtheta,0.0];
    elseif n>(P.numRays-8)
        TX(2*P.numRays+n).Steer = [((P.numRays-n)/8)*P.dtheta,0.0];
    else
        TX(2*P.numRays+n).Steer = [P.dtheta,0.0];
    end
    % - Compute transmit delays
    TX(2*P.numRays+n).Delay = computeTXDelays(TX(2*P.numRays+n));
end

% Compute transmit pixel data
h = waitbar(0,'Program TX parameters, please wait!');
steps = 3*P.numRays;
for i = 1:steps
    TX(i).TXPD = computeTXPD(TX(i),PData);
    waitbar(i/steps)
end
close(h)

% Specify TGC Waveform structure.
TGC.CntrlPts = [60,425,677,777,890,956,1023,1023];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Receive structure arrays.
% Sampling center frequency is 15.625, but we want the bandpass filter centered on
% the actual transducer center frequency of 18 MHz with 67% bandwidth, or 12 to 24
% MHz.  Coefficients below were set using "G3_BPFdevelopment" with normalized cf=1.15
% (18 MHz), bw=0.85, xsn wdth=0.41 resulting in -3 dB 0.71 to 1.6 (11.1 to 25 MHz),
% and -20 dB 0.57 to 1.74 (8.9 to 27.2 MHz)
%
BPF1 = [ -0.00009 -0.00128 +0.00104 +0.00085 +0.00159 +0.00244 -0.00955 ...
         +0.00079 -0.00476 +0.01108 +0.02103 -0.01892 +0.00281 -0.05206 ...
         +0.01358 +0.06165 +0.00735 +0.09698 -0.27612 -0.10144 +0.48608 ];

maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
Receive = repmat(struct('Apod', [zeros(1,16), ones(1,32),zeros(1,16)], ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'InputFilter', BPF1, ...
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, 3*P.numRays*Resource.RcvBuffer(1).numFrames);
% - Set event specific Receive attributes for each frame.
for i = 1:Resource.RcvBuffer(1).numFrames
    k = 3*P.numRays*(i-1);
    Receive(k+1).callMediaFunc = 1;
    for j = 1:P.numRays
        Receive(k+j).aperture = TX(j).aperture; % mux aperture same as transmit
        Receive(k+j).framenum = i;
        Receive(k+j).acqNum = j;
        % Receives for steered left beams
        Receive(k+P.numRays+j).aperture = Receive(k+j).aperture;
        Receive(k+P.numRays+j).bufnum = 2;
        Receive(k+P.numRays+j).framenum = i;
        Receive(k+P.numRays+j).acqNum = j;
        % Receives for steered right beams
        Receive(k+2*P.numRays+j).aperture = Receive(k+j).aperture;
        Receive(k+2*P.numRays+j).bufnum = 3;
        Receive(k+2*P.numRays+j).framenum = i;
        Receive(k+2*P.numRays+j).acqNum = j;
    end
end

% Specify Recon structure arrays.
% - We need one Recon structure which will be used for each frame.
Recon = repmat(struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame',-1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...
               'RINums',(1:P.numRays)), 1, 3);
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
                         'pgain',1.7,...            % pgain is image processing gain
                         'reject',2,...      % reject level
                         'persistMethod','simple',...
                         'persistLevel',pers,...
                         'interpMethod','4pt',...
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','runAverage3',...
                         'compressMethod','power',...
                         'compressFactor',45,...
                         'mappingMethod','full',...
                         'display',1,...      % display image after processing
                         'displayWindow',1};

% Specify SeqControl structure arrays.
TTNA = round(2*384*(1/Trans.frequency))+20; % acq. time in usec for max depth, plus 20 usec for overhead and HVMux settling

%  - Jump back to start.
SeqControl(1).command = 'jump';
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(2).argument = TTNA;  % 220 usec
SeqControl(3).command = 'timeToNextAcq';
SeqControl(3).argument = 12000;  % 10 ms between frames
SeqControl(4).command = 'returnToMatlab';
nsc = 5; % nsc is count of SeqControl objects

% Specify Event structure arrays.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    % Acquire frame with straight ahead wide beams
    k = 3*P.numRays*(i-1);
    for j = 1:P.numRays        % Acquire frame
        Event(n).info = 'acquire widebeam';
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
    % Reconstruct straight ahead frame
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
    % Reconstruct steered left frame
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
    % Reconstruct steered right frame
    Event(n).info = 'recon and process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 3;
    Event(n).process = 1;
    Event(n).seqControl = 0;
    if i ~= Resource.RcvBuffer(1).numFrames  % Exit to Matlab after 3 frames
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
MinMaxVal = [64,250,P.endDepth]; % default unit is wavelength
AxesUnit = 'wls';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
        AxesUnit = 'mm';
        MinMaxVal = MinMaxVal * (Resource.Parameters.speedOfSound/1000/Trans.frequency);
    end
end
UI(2).Control = VsSliderControl('LocationCode','UserA1',...
    'Label',['Range (',AxesUnit,')'],'SliderMinMaxVal',MinMaxVal,...
    'SliderStep',[0.1,0.2],'ValueFormat','%3.0f',...
    'Callback',@RangeChangeCallback);

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 3;

% Save all the structures to a .mat file.
save('MatFiles/L22-14vLFWideBeamSC_32LE');


% **** Callback functions used by UI controls. ****
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

    PData = evalin('base','PData');
    PData.Size(1) = ceil((P.endDepth-P.startDepth)/PData.PDelta(3));
    PData.Region = repmat(struct('Shape',struct( ...
                        'Name','Parallelogram',...
                        'Position',[0,0,P.startDepth],...
                        'width',Trans.spacing*(0.55*P.numTx),...
                        'height',P.endDepth-P.startDepth,...
                        'angle',0.0)),1,3*P.numRays);
    TxOrgX = (-63.5*Trans.spacing):(127*Trans.spacing/(P.numRays-1)):(63.5*Trans.spacing);
    for n = 1:P.numRays
        PData.Region(n).Shape.Position(1) = TxOrgX(n);
        PData.Region(P.numRays+n).Shape.Position(1) = TxOrgX(n);
        PData.Region(2*P.numRays+n).Shape.Position(1) = TxOrgX(n);
    end
    % Define P.numRays steered left parallelogram regions, centered on TX beam origins. Adjust the angle
    %   so that the steering goes to zero over 8 beams at the left and right edge.
    m = P.numRays;
    for n = 1:P.numRays
        if n<=8
            angle = -((n-1)/8)*P.dtheta;
        elseif n>(P.numRays-8)
            angle = -((P.numRays-n)/8)*P.dtheta;
        else
            angle = -P.dtheta;
        end
        PData.Region(n+m).Shape.height = (P.endDepth-P.startDepth)/cos(angle);
        PData.Region(n+m).Shape.angle = angle;
    end
    % Define numRays steered right parallelogram regions, centered on TX beam origins. Adjust the angle
    %   so that the steering goes to zero over 8 beams at the left and right edge.
    m = m + P.numRays;
    for n = 1:P.numRays
        if n<=8
            angle = ((n-1)/8)*P.dtheta;
        elseif n>(P.numRays-8)
            angle = ((P.numRays-n)/8)*P.dtheta;
        else
            angle = P.dtheta;
        end
        PData.Region(n+m).Shape.height = (P.endDepth-P.startDepth)/cos(angle);
        PData.Region(n+m).Shape.angle = angle;
    end
    PData.Region = computeRegions(PData);
    assignin('base','PData',PData);
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
    % Update Receive structures.
    Receive = evalin('base', 'Receive');
    maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
    for i = 1:size(Receive,2)
        Receive(i).endDepth = maxAcqLength;
    end
    assignin('base','Receive',Receive);
    evalin('base','TGC.rangeMax = P.endDepth;');
    evalin('base','TGC.Waveform = computeTGCWaveform(TGC);');
    evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData.Size(1)*PData.PDelta(3)/Resource.DisplayWindow(1).pdelta);');
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'PData','InterBuffer','ImageBuffer','TX','Receive','TGC','Recon','DisplayWindow'};
    assignin('base','Control', Control);
    assignin('base', 'action', 'displayChange');
end
