% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use.
%
% File name: SetUpC5_2gHWideBeamHI.m - Example of curved array imaging with
%                                   wide beam transmits and harmonic imaging.
% Description:
%   Sequence programming file for C5-2gH curved array using widebeam pulse
%   inversion transmits and receive acquisitions. Of the 128 transmit channels, only
%   numTx are used, with numTx/2 transmitters on each side of the center element
%   (where possible). All 128 receive channels are used, although the
%   element sensitivity cutoff will limit the useful aperture. The receive acquisitions
%   use 200% bandwidth. Reconstruction and Processing are asynchronous with
%   respect to acquisition.
%
%   The PData regions are generated to provide an overlap of 5 regions at every pixel.
%   The transmit focus is set to below the bottom of the PData region to generate a wide
%   beam that is as wide at the bottom of a single region.
%
% Last update:
% 05/06/2020 - Update to SW 4.3 format for new user UIControls and External function definitions (VTS 1691).
%   More info:(?/Example_Scripts/Vantage_Features/New UI Scheme/SetUpL11_5vFlash_NewUI)
% 10/01/2021 - modified for C5-2gH

clear all

P.numTx = 32;   % no. of elements in TX aperture.
P.numRays = 60; % no. of rays in frame
P.txFocus = 480;
P.startDepth = 5;  % startDepth in wavelength
P.endDepth = 192;

% Specify system parameters.
Resource.Parameters.speedOfSound = 1540;
Resource.Parameters.numTransmit = 128;  % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;  % number of receive channels.
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Specify Trans structure array.
Trans.name = 'C5-2gH';
Trans.units = 'wavelengths'; % required in Gen3 to prevent default to mm units
Trans = computeTrans(Trans);  % C5-2gH transducer is 'known' transducer so we can use computeTrans.

radius = Trans.radius;
scanangle = Trans.numelements*Trans.spacing/radius;
dtheta = scanangle/P.numRays; % angle between rays
theta = -(scanangle/2) + 0.5*dtheta; % angle to left edge from centerline
Angle = theta:dtheta:(-theta);

% Specify PData structure array.
PData(1).PDelta = [1.0, 0, 0.5];  % x, y and z pdeltas
sizeRows = 10 + ceil((P.endDepth + radius - (radius * cos(scanangle/2)))/PData(1).PDelta(3));
sizeCols = 10 + ceil(2*(P.endDepth + radius)*sin(scanangle/2)/PData(1).PDelta(1));
PData(1).Size = [sizeRows,sizeCols,1];     % size, origin and pdelta set region of interest.
PData(1).Origin(1,1) = (P.endDepth+radius)*sin(-scanangle/2) - 5;
PData(1).Origin(1,2) = 0;
PData(1).Origin(1,3) = ceil(radius * cos(scanangle/2)) - radius - 5;
% Define PData Regions for numRays scanlines
for n = 1:P.numRays
    PData(1).Region(n) = struct(...
        'Shape',struct('Name','Sector',...
                       'Position',[0,0,-radius],...
                       'r1',radius+P.startDepth,...
                       'r2',radius+P.endDepth,...
                       'angle',dtheta*5,...
                       'steer',Angle(n)));
end
PData(1).Region = computeRegions(PData(1));

%  Media points for curved array.
% - Uncomment for speckle
% Media.numPoints = (20000);
% Media.MP = rand(Media.numPoints,4);
% Media.MP(:,2) = 0;
% Media.MP(:,4) = 0.01 + 0.04*Media.MP(:,4);  % Random low amplitude
% RandR = P.endDepth*Media.MP(:,1)+radius;
% RandTheta = scanangle*(Media.MP(:,3)-0.5);
% Media.MP(:,1) = RandR.*sin(RandTheta);
% Media.MP(:,3) = RandR.*cos(RandTheta)-radius;
% - Define points
%Media.MP(1,:) = [0,0,radius+70,1.0];
Media.MP(1,:) = [0,0,10,1.0];
Media.MP(2,:) = [(radius+10)*sin(-0.2608),0,(radius+10)*cos(-0.2608)-radius,1.0];
Media.MP(3,:) = [(radius+10)*sin(0.2608),0,(radius+10)*cos(0.2608)-radius,1.0];
Media.MP(4,:) = [(radius+10)*sin(-0.5267),0,(radius+10)*cos(-0.5267)-radius,1.0];
Media.MP(5,:) = [(radius+10)*sin(0.5267),0,(radius+10)*cos(0.5267)-radius,1.0];
Media.MP(6,:) = [0,0,40,1.0];
Media.MP(7,:) = [0,0,70,1.0];
Media.MP(8,:) = [(radius+70)*sin(-0.2608),0,(radius+70)*cos(-0.2608)-radius,1.0];
Media.MP(9,:) = [(radius+70)*sin(0.2608),0,(radius+70)*cos(0.2608)-radius,1.0];
Media.MP(10,:) = [(radius+70)*sin(-0.5267),0,(radius+70)*cos(-0.5267)-radius,1.0];
Media.MP(11,:) = [(radius+70)*sin(0.5267),0,(radius+70)*cos(0.5267)-radius,1.0];
Media.MP(12,:) = [0,0,100,1.0];
Media.MP(13,:) = [0,0,130,1.0];
Media.MP(14,:) = [(radius+130)*sin(-0.2608),0,(radius+130)*cos(-0.2608)-radius,1.0];
Media.MP(15,:) = [(radius+130)*sin(0.2608),0,(radius+130)*cos(0.2608)-radius,1.0];
Media.MP(16,:) = [(radius+130)*sin(-0.5267),0,(radius+130)*cos(-0.5267)-radius,1.0];
Media.MP(17,:) = [(radius+130)*sin(0.5267),0,(radius+130)*cos(0.5267)-radius,1.0];
Media.MP(18,:) = [0,0,160,1.0];
Media.MP(19,:) = [0,0,190,1.0];
Media.MP(20,:) = [(radius+190)*sin(-0.2608),0,(radius+190)*cos(-0.2608)-radius,1.0];
Media.MP(21,:) = [(radius+190)*sin(0.2608),0,(radius+190)*cos(0.2608)-radius,1.0];
Media.MP(22,:) = [(radius+190)*sin(-0.5267),0,(radius+190)*cos(-0.5267)-radius,1.0];
Media.MP(23,:) = [(radius+190)*sin(0.5267),0,(radius+190)*cos(0.5267)-radius,1.0];
Media.function = 'movePoints';

% Specify Resources.
Resource.RcvBuffer.datatype = 'int16';
Resource.RcvBuffer.rowsPerFrame = P.numRays*4480;    
Resource.RcvBuffer.colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer.numFrames = 10;
Resource.InterBuffer.datatype = 'complex double';
Resource.InterBuffer.numFrames = 1;
Resource.ImageBuffer.datatype = 'double';
Resource.ImageBuffer.numFrames = 10;
Resource.DisplayWindow(1).Title = 'C5-2gHWideBeamHI';
Resource.DisplayWindow(1).pdelta = 0.45;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).numFrames = 20;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow.Colormap = gray(256);

% Specify Transmit waveform structure.
TW(1) = struct('type','parametric','Parameters',[2.25,.67,2,-1]);
TW(2) = struct('type','parametric','Parameters',[2.25,.67,2,1]);

% Specify TX structure array.
% - We need 128 transmit specifications.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'focus', P.txFocus, ...
                   'Steer', [0.0,0.0], ...
                   'Apod', zeros(1,Trans.numelements), ...
                   'Delay', zeros(1,Trans.numelements), ...
                   'TXPD', [], ...
                   'peakCutOff', 1.0, ...
                   'peakBLMax', 4.0), 1, 2*P.numRays);

% - Set event specific TX attributes.
%    numTx/2 is the number of elements to include on each side of the
%    center element, for the specified focus and sensitivity cutoff.
for n = 1:P.numRays   % numRays transmit events
    TX(n + P.numRays).waveform = 2;
    % Set transmit Origins.
    TX(n).Origin = [radius*sin(Angle(n)), 0.0, radius*cos(Angle(n))-radius];
    TX(n+P.numRays).Origin = TX(n).Origin;
    ce = round(1+127*(Angle(n) - theta)/(-2*theta));
    % Set transmit Apodization so that a maximum of numTx + 1 transmitters are active.
    lft = round(ce - P.numTx/2);
    if lft < 1, lft = 1; end
    rt = round(ce + P.numTx/2);
    if rt > Trans.numelements, rt = Trans.numelements; end
    TX(n).Apod(lft:rt) = 1.0;
    % Apply apodization function.
    [RIndices,CIndices,V] = find(TX(n).Apod);
    V = kaiser(size(V,2),1);
    TX(n).Apod(CIndices) = V;
    TX(n+P.numRays).Apod = TX(n).Apod;
    TX(n).Delay = computeTXDelays(TX(n));
    TX(n+P.numRays).Delay = TX(n).Delay;
    TX(n).TXPD = computeTXPD(TX(n),PData(1));
    TX(n+P.numRays).TXPD = TX(n).TXPD;
end

% Specify Receive structure arrays.
% -- Compute the maximum receive path length, using the law of cosines.
maxAcqLength = ceil(sqrt((P.endDepth+radius)^2 + radius^2 - ...
                     2*(P.endDepth+radius)*radius*cos(scanangle)));
Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'demodFrequency',4.5, ...
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, 2*P.numRays*Resource.RcvBuffer(1).numFrames);
% - Set event specific Receive attributes.
m = P.numRays;
for i = 1:Resource.RcvBuffer(1).numFrames
    Receive(m*(i-1)+1).callMediaFunc = 1;
    for j = 1:m
        Receive(2*m*(i-1)+j).framenum = i;
        Receive(2*m*(i-1)+j).acqNum = j;
    end
    for j = 1:m
        Receive(2*m*(i-1)+j+m).framenum = i;
        Receive(2*m*(i-1)+j+m).mode = 1;
        Receive(2*m*(i-1)+j+m).acqNum = j;
    end
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [0,438,600,658,770,830,941,1023];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Recon structure arrays.
Recon = struct('senscutoff', 0.55, ...
               'pdatanum', 1, ...
               'rcvBufFrame',-1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...
               'RINums', 1:P.numRays);

% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 'accumIQ', ...  % accumulate IQ
                   'Pre',[],...
                   'Post',[],...
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'scaleFactor', 0.2, ...
                   'regionnum', 0), 1, P.numRays);
% - Set specific ReconInfo attributes.
ReconInfo(1).Pre = 'clearInterBuf';     % clear entire Interbuffer frame
for j = 1:P.numRays
    ReconInfo(j).txnum = j;
    ReconInfo(j).rcvnum = j;
    ReconInfo(j).regionnum = j;
end
ReconInfo(P.numRays).Post = 'IQ2IntensityImageBuf';  % detect IQ data

% Specify Process structure array.
%   First processing structure specifies how to display image.
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,... % frame number in src buffer (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure (defines output figure).
                         'norm',1,...        % normalization method (1 means fixed)
                         'pgain',1.0,...            % pgain is image processing gain
                         'reject',5,...
                         'persistMethod','simple',...
                         'persistLevel',20,...
                         'interpMethod','4pt',...
                         'compressMethod','power',...      % X^0.5 normalized to output word size
                         'compressFactor',40,...
                         'averageMethod','none',...
                         'mappingMethod','full',...
                         'display',1,...      % display image after processing
                         'displayWindow',1};

% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; %  jump back to start.
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq'; % time between rays
SeqControl(2).argument = 300; % 300usec
SeqControl(3).command = 'timeToNextAcq'; % time between frames
SeqControl(3).argument = 10000; % 10msec
SeqControl(4).command = 'returnToMatlab';
nsc = 5;

% Specify Event structure arrays.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    for j = 1:m                 % Acquire rays
        Event(n).info = 'Acquire ray line, normal TX';
        Event(n).tx = j;
        Event(n).rcv = 2*m*(i-1)+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;

        Event(n).info = 'Acquire ray line, inverted TX';
        Event(n).tx = j+m;
        Event(n).rcv = 2*m*(i-1)+j+m;
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
    if (floor(i/4) == i/4)&&(i < Resource.RcvBuffer(1).numFrames)     % Exit to Matlab every 4th frame reconstructed
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
wls2mm = 1;
AxesUnit = 'wls';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
        AxesUnit = 'mm';
        wls2mm = Resource.Parameters.speedOfSound/1000/Trans.frequency;
    end
end
UI(2).Control = VsSliderControl('LocationCode','UserA1',...
    'Label',['Range (',AxesUnit,')'],...
    'SliderMinMaxVal',[64,300,P.endDepth]*wls2mm,...
    'SliderStep',[0.1,0.2],'ValueFormat','%3.0f',...
    'Callback',@RangeChangeCallback);

% - Transmit focus change
UI(3).Control = VsSliderControl('LocationCode','UserB4',...
    'Label',['TX Focus (',AxesUnit,')'],...
    'SliderMinMaxVal',[100,600,P.txFocus]*wls2mm,...
    'SliderStep',[0.1,0.2],'ValueFormat','%3.0f',...;
    'Callback',@TxFocusCallback);

% - numTx
UI(4).Control = VsSliderControl('LocationCode','UserB3',...
    'Label','numTx','SliderMinMaxVal',[10,64,P.numTx],...
    'SliderStep',[4/54,8/54],'ValueFormat','%3.0f',...
    'Callback',@numTxCallback);

% - Peak CutOff
UI(5).Control = VsSliderControl('LocationCode','UserB2',...
    'Label','Peak Cutoff',...
    'SliderMinMaxVal',[0,20.0,TX(1).peakCutOff],...
    'SliderStep',[0.005,0.020],'ValueFormat','%1.3f',...
    'Callback',@PeakCutOffCallback);

% - Max. Burst Length
UI(6).Control = VsSliderControl('LocationCode','UserB1',...
    'Label','Max. BL',...
    'SliderMinMaxVal',[0,20.0,TX(1).peakBLMax],...
    'SliderStep',[0.005,0.020],'ValueFormat','%1.3f',...
    'Callback',@MaxBLCallback);

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 4;

% Save all the structures to a .mat file.
save('MatFiles/C5-2gHWideBeamHI');

%% Automatic VSX Execution:
% Uncomment the following line to automatically run VSX every time you run
% this SetUp script (note that if VSX finds the variable 'filename' in the
% Matlab workspace, it will load and run that file instead of prompting the
% user for the file to be used):

% filename = 'C5-2gHWideBeamHI';  VSX;


%% **** Callback routines used by UIControls (UI) ****

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

    scanangle = evalin('base','scanangle');
    radius = evalin('base','radius');
    theta = evalin('base','theta');
    dtheta = -2*theta/P.numRays;
    Angle = theta:dtheta:(-theta);
    height = P.endDepth + radius - (radius*cos(scanangle/2));
    % Modify PData for new range
    PData = evalin('base','PData');
    PData(1).Size(1) = 10+ceil(height/PData(1).PDelta(3));
    for j = 1:P.numRays
        PData(1).Region(j).Shape = struct( ...
                           'Name','Sector',...
                           'Position',[0,0,-radius],...
                           'r1',radius+P.startDepth,...
                           'r2',radius+P.endDepth,...
                           'angle',dtheta*5,...
                           'steer',Angle(j));
        PData(1).Region(j).PixelsLA = [];
        PData(1).Region(j).numPixels = 0;
    end
    assignin('base','PData',PData);
    evalin('base','PData(1).Region = computeRegions(PData(1));');
    evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);');
    % Update TXPD data of TX structures.
    TX = evalin('base','TX');
    PData = evalin('base','PData');
    for i = 1:P.numRays
        TX(i).TXPD = computeTXPD(TX(i),PData(1));
        TX(i+P.numRays).TXPD = TX(i).TXPD;
    end
    assignin('base','TX',TX);
    % Update Receive structures
    Receive = evalin('base', 'Receive');
    maxAcqLength = ceil(sqrt((P.endDepth+radius)^2 + radius^2 - ...
                         2*(P.endDepth+radius)*radius*cos(scanangle)));
    for i = 1:size(Receive,2)
        Receive(i).endDepth = maxAcqLength;
    end
    assignin('base','Receive',Receive);
    evalin('base','TGC.rangeMax = P.endDepth;');
    evalin('base','TGC.Waveform = computeTGCWaveform(TGC);');
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'PData','InterBuffer','ImageBuffer','TX','Receive','TGC','Recon','DisplayWindow'};
    assignin('base','Control', Control);
    assignin('base', 'action', 'displayChange');
end

function TxFocusCallback(hObject, ~, UIValue)
    simMode = evalin('base','Resource.Parameters.simulateMode');
    % No focus change if in simulate mode 2.
    if simMode == 2
        set(hObject,'Value',evalin('base','P.txFocus'));
        return
    end
    P = evalin('base','P');
    P.txFocus = UIValue;
    assignin('base','P',P);
    PData = evalin('base','PData');
    TX = evalin('base', 'TX');
    for n = 1:P.numRays
        TX(n).focus = P.txFocus;
        TX(n+P.numRays).focus = TX(n).focus;
        TX(n).Delay = computeTXDelays(TX(n));
        TX(n+P.numRays).Delay = TX(n).Delay;
        TX(n).TXPD = computeTXPD(TX(n),PData(1));
        TX(n+P.numRays).TXPD = TX(n).TXPD;
    end
    assignin('base','TX', TX);
    % Set Control command to update TX
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'TX'};
    assignin('base','Control', Control);
end

function numTxCallback(hObject, ~, UIValue)
% number of transmitters in TX aperture
    simMode = evalin('base','Resource.Parameters.simulateMode');
    % No numTx change if in simulate mode 2.
    if simMode == 2
        set(hObject,'Value',evalin('base','P.numTx'));
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
    TX = evalin('base', 'TX');
    Angle = evalin('base','Angle');
    Trans = evalin('base','Trans');
    radius = Trans.radius;
    theta = evalin('base','theta');
    for n = 1:P.numRays
        TX(n).Origin = [radius*sin(Angle(n)), 0.0, radius*cos(Angle(n))-radius];
        TX(n+P.numRays).Origin = TX(n).Origin;
        ce = round(1+127*(Angle(n) - theta)/(-2*theta));
        % Set transmit Apodization so that a maximum of numTx + 1 transmitters are active.
        lft = round(ce - P.numTx/2);
        if lft < 1, lft = 1; end
        rt = round(ce + P.numTx/2);
        if rt > Trans.numelements, rt = Trans.numelements; end
        TX(n).Apod = zeros(1,Trans.numelements);
        TX(n).Apod(lft:rt) = 1.0;
        % Apply apodization function.
        [~,CIndices,V] = find(TX(n).Apod);
        V = kaiser(size(V,2),1);
        TX(n).Apod(CIndices) = V;
        TX(n+P.numRays).Apod = TX(n).Apod;
        TX(n).Delay = computeTXDelays(TX(n));
        TX(n+P.numRays).Delay = TX(n).Delay;
        TX(n).TXPD = computeTXPD(TX(n),PData(1));
        TX(n+P.numRays).TXPD = TX(n).TXPD;
    end
    assignin('base','TX', TX);
    % Set Control command to update TX
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'TX'};
    assignin('base','Control', Control);
end

function PeakCutOffCallback(~, ~, UIValue)
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
end

function MaxBLCallback(~, ~, UIValue)
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
end