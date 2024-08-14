% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: SetUpL11_5gHFlashAnglesSC_64LE.m - Plane wave transmits with steering
%                                         angles and spatial compounding
% Description:
%   Sequence programming file for L11-5gH Linear array, using plane wave
%   transmits with multiple steering angles and spatile compouding over 3
%   angle ranges. Acquisition is over na uniformly spaced angles over +/- 20
%   degrees. Reconstruction is on na/3 angles at a time, with the detected
%   result averaged multiplicatively (corresponded to documentation page 79).
%   All 128 transmit and 64 receive channels are active for 2-1 synthetic
%   aperture acquisition. Processing is asynchronous with respect to acquisition.
%
%  09/20/21 modified for L11-5gH

clear all

na = 15;      % Set na = number of angles (na should be a multiple of 3).
if (na > 1)
  dtheta = (40*pi/180)/(na-1);
  startAngle = -40*pi/180/2; % set dtheta to range over +/- 20 degrees.
else
  error('na should be a multiple of 3.\n');
end

% Specify system parameters.
Resource.Parameters.numTransmit = 128;  % number of transmit channels.
Resource.Parameters.numRcvChannels = 64;  % number of receive channels.
Resource.Parameters.speedOfSound = 1540;
Resource.Parameters.simulateMode = 0; % 0 means no simulation, if hardware is present.
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.
Resource.Parameters.verbose = 2;
% Resource.Parameters.initializeOnly = 0;

% Specify Trans structure array.
Trans.name = 'L11-5gH';
Trans.units = 'mm'; % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);  % L11-5gH transducer is 'known' transducer so we can use computeTrans.
Trans.maxHighVoltage = 50;  % set maximum high voltage limit for pulser supply.

% Specify P structure array.
P.radius = 0;              % ROC for curved lin. or dist. to virt. apex
P.theta = 0;
P.numRays = 1;      % no. of Rays (1 for Flat Focus)
P.rayDelta = 128*Trans.spacing;  % spacing in radians(sector) or dist. between rays (wvlnghts)
P.startDepth = 5;   % Acquisition depth in wavelengths
P.endDepth = 192;   % This should preferrably be a multiple of 128 samples.

% Specify PData structure array.
PData(1).PDelta = [ 1.0, 0, 0.5];
PData(1).Size(1,1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3)); % range, origin and pdelta set PData(1).Size.
PData(1).Size(1,2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(1,3) = 1;      % single image page
PData(1).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of uppr lft crnr.

% Specify Media object.
pt1;
Media.attenuation = -0.5;
Media.function = 'movePoints';

% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = na*4096*2; % this size allows for maximum range
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numTransmit;
Resource.RcvBuffer(1).numFrames = 30;    % 30 frames stored in RcvBuffer.
Resource.InterBuffer(1).numFrames = 1;   % one intermediate buffer needed.
Resource.ImageBuffer(1).numFrames = 10;
Resource.DisplayWindow(1).Title = 'L11-5gHFlashAnglesSC_64LE';
Resource.DisplayWindow(1).pdelta = 0.35;
Resource.DisplayWindow(1).Position = [250,150, ...
    ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta), ... % width
    ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta)];    % height
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)]; % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).numFrames = 20;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = gray(256);

% Specify Transmit waveform structure.
TW = struct('type','parametric', ...
            'Parameters',[Trans.frequency,.67,2,1]); % center frequency of 6.43

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'focus', 0.0, ...
                   'Steer', [0.0,0.0], ...
                   'Apod', kaiser(Resource.Parameters.numTransmit,1)', ...
                   'Delay', zeros(1,Trans.numelements)), 1, na);
% - Set event specific TX attributes.
for n = 1:na   % na transmit events
    TX(n).Steer = [(startAngle+(n-1)*dtheta),0.0];
    TX(n).Delay = computeTXDelays(TX(n));
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [0,297,424,515,627,764,922,1000];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Receive structure arrays.
% - We need 2*na Receives for every frame.
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
Receive = repmat(struct('Apod', zeros(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength,...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, 2*na*Resource.RcvBuffer(1).numFrames);
% - Set event specific Receive attributes for each frame.
for i = 1:Resource.RcvBuffer(1).numFrames
    Receive(2*na*(i-1)+1).callMediaFunc = 1;
    for j = 1:2:2*na
        Receive(2*na*(i-1)+j).Apod(1:Resource.Parameters.numRcvChannels) = 1.0;
        Receive(2*na*(i-1)+j).framenum = i;
        Receive(2*na*(i-1)+j).acqNum = j;      % two acquisitions per frame
        Receive(2*na*(i-1)+j+1).Apod((Resource.Parameters.numRcvChannels+1):Trans.numelements) = 1.0;
        Receive(2*na*(i-1)+j+1).framenum = i;
        Receive(2*na*(i-1)+j+1).acqNum = j+1;  % two acquisitions per frame
    end
end

% Specify Recon structure arrays.
% - We need one Recon structure which will be used for each frame.
Recon = struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame',-1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...
               'RINums',(1:2*na));

% Define ReconInfo structures.
% We need 2*na ReconInfo structures for na steering angles andn the 2-1 synthetic aperture acquisition.
ReconInfo = repmat(struct('mode', 'accumIQ', ...  % default is to accumulate IQ data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'normPower',0, ...
                   'regionnum', 1), 1, 2*na);
% - Set specific ReconInfo attributes.
for j = 1:2:2*na  % For each row in the column
    ReconInfo(j).txnum = (j+1)/2;
    ReconInfo(j).rcvnum = j;
    ReconInfo(j+1).txnum = (j+1)/2;
    ReconInfo(j+1).rcvnum = j+1;
end
ReconInfo(1).mode = 'replaceIQ';  % replace IQ data
j = round(2*na/3);
ReconInfo(j).mode = 'accumIQ_replaceIntensity';  % accum and detect
ReconInfo(j).normPower = 0.33;
ReconInfo(j+1).mode = 'replaceIQ';
j = round(4*na/3);
ReconInfo(j).mode = 'accumIQ_multiplyIntensity'; % accum, detect and multiply-accum
ReconInfo(j).normPower = 0.33;
ReconInfo(j+1).mode = 'replaceIQ';
ReconInfo(2*na).mode = 'accumIQ_multiplyIntensity'; % accum, detect and multiply-accum
ReconInfo(2*na).normPower = 0.33;

% Specify Process structure array.
pers = 20;
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1, ...   % number of buffer to process.
                         'framenum',-1, ...   % frame number in src buffer (-1 => lastFrame)
                         'pdatanum',1, ...    % number of PData structure (defines output figure).
                         'pgain',3.0, ...            % pgain is image processing gain
                         'reject',2,...       % reject level
                         'persistMethod','simple', ...
                         'persistLevel',pers, ...
                         'interpMethod','4pt',...
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',40,...
                         'mappingMethod','full',...
                         'display',1, ...     % display image after processing
                         'displayWindow',1};

% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; % jump back to start
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(2).argument = 160;  % 160 usec
SeqControl(3).command = 'timeToNextAcq';  % time between frames
SeqControl(3).argument = 20000 - (na-1)*160;  % 20 msec
SeqControl(4).command = 'returnToMatlab';
nsc = 5; % nsc is count of SeqControl objects

% Specify Event structure arrays.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    for j = 1:2:2*na                      % Acquire frame
        Event(n).info = '1st half of aperture.';
        Event(n).tx = (j+1)/2;
        Event(n).rcv = 2*na*(i-1)+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;

        Event(n).info = '2nd half of aperture.';
        Event(n).tx = (j+1)/2;
        Event(n).rcv = 2*na*(i-1)+j+1;
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

    if floor(i/2) == i/2     % Exit to Matlab every 2nd frame
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
    'Callback',@SensCutoffCallback);

% - Range Change
MinMaxVal = [64,300,P.endDepth]; % default unit is wavelength
AxesUnit = 'wls';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
        AxesUnit = 'mm';
        MinMaxVal = MinMaxVal * (Resource.Parameters.speedOfSound/1000/Trans.frequency);
    end
end
UI(2).Control = VsSliderControl('LocationCode','UserA1',...
    'Label',['Range (',AxesUnit,')'],...
    'SliderMinMaxVal',MinMaxVal,'SliderStep',[0.1,0.2],...
    'ValueFormat','%3.0f',...
    'Callback',@RangeChangeCallback);

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 2;

% Save all the structures to a .mat file.
save('MatFiles/L11-5gHFlashAnglesSC_64LE');

% **** Callback functions to be used by UI Controls. ****
function SensCutoffCallback(~,~,UIValue)
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

    evalin('base','PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));');
    evalin('base','PData(1).Region = computeRegions(PData(1));');
    evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);');
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
    Control.Parameters = {'PData','InterBuffer','ImageBuffer','DisplayWindow','Receive','TGC','Recon'};
    assignin('base','Control', Control);
    assignin('base', 'action', 'displayChange');
end
