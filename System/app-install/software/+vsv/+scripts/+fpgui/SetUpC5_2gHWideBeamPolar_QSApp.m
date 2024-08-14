% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: SetUpC5_2gHWideBeamPolar.m - Example of curved array imaging with
%                             wide beam transmits and polar coordinate PData.
% Description:
%   Sequence programming file for C5-2gH curved array using P.numRays wide beam
%   transmits and receive acquisitions. 
%   Processing is asynchronous with respect to acquisition.
%
% Last update:
%   09/28/2021 Adapted script from SetUpC5_2vWideBeamPolar.m

%clear[ 	]+all

P.startDepth = 5; % startDepth and endDepth in wavelengths
P.endDepth = 300;
P.numRays = 64;   % no. of transmit beams to program
P.numTx = 60;     % no. of elements in transmit aperture - 1, should not larger than 62
P.txFocus = 4*P.endDepth;  % focal point in wavelengths

% Specify system parameters.
Resource.Parameters.numTransmit = 128;
Resource.Parameters.speedOfSound = 1540;
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

Resource.System.SoftwareVersion = [4 0 0]; % Minimum software release for this script.

% Specify Trans structure array.
Trans.name = 'C5-2gH';
Trans.units = 'wavelengths'; % required in Gen3 to prevent default to mm units
Trans = computeTrans(Trans);  % C5-2gH transducer is 'known' transducer so we can use computeTrans.


% useful parameters
radius = Trans.radius;  % radius of curved array in wavelengths
scanangle = (Trans.numelements-1)*Trans.spacing/radius;
dtheta = scanangle/(P.numRays-1);
theta = -(scanangle/2); % angle to left edge from centerline
Angle = theta:dtheta:(-theta);

% Specify PData structure array.
PData.Coord = 'polar';
PData.PDelta = [scanangle/127, 0.5, 0]; % [theta,r,z] theta increment is equivalent to 128 ray line scan
sizeRows = 10 + ceil((P.endDepth + radius)/PData.PDelta(2));
sizeCols = 128;
PData.Size = [sizeRows,sizeCols,1];     % size, origin and pdelta set region of interest.
PData.Origin = [0,0,-radius];   % Origin is in rectangular coords.
% Define PData Regions for numRays scanlines
for j = 1:P.numRays
    PData.Region(j) = struct('Shape',struct( ...
                       'Name','Sector',...
                       'Position',PData.Origin,...
                       'r1',radius+P.startDepth,...
                       'r2',radius+P.endDepth,...
                       'angle',dtheta*7,...
                       'steer',Angle(j)));
end
PData.Region = computeRegions(PData);

%  Media points for curved array.
% - Uncomment for speckle
% Media.numPoints = (20000);
% Media.MP = rand(Media.numPoints,4);
% Media.MP(:,2) = 0;
% Media.MP(:,4) = 0.01 + 0.04*Media.MP(:,4);  % Random low amplitude
% RandR = endDepthMm *Media.MP(:,1)+radius;
% RandTheta = scanangle*(Media.MP(:,3)-0.5);
% Media.MP(:,1) = RandR.*sin(RandTheta);
% Media.MP(:,3) = RandR.*cos(RandTheta)-radius;
% - Define points
%Media.MP(1,:) = [0,0,70,1.0];
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
Media.MP(24,:) = [P.startDepth*sin(-0.2608),0,P.startDepth*cos(-0.2608)-radius,1.0];
Media.MP(25,:) = [P.startDepth*sin(0.2608),0,P.startDepth*cos(0.2608)-radius,1.0];
Media.MP(26,:) = [P.startDepth*sin(-0.5267),0,P.startDepth*cos(-0.5267)-radius,1.0];
Media.MP(27,:) = [P.startDepth*sin(0.5267),0,P.startDepth*cos(0.5267)-radius,1.0];
Media.function = 'movePoints';

% Specify Resources.
Resource.RcvBuffer.datatype = 'int16';
Resource.RcvBuffer.rowsPerFrame = P.numRays*4096;    % numRays segments of 4096 samples
Resource.RcvBuffer.colsPerFrame = Resource.Parameters.numTransmit;
Resource.RcvBuffer.numFrames = 10;
Resource.InterBuffer.datatype = 'complex double';
Resource.InterBuffer.numFrames = 1;
Resource.ImageBuffer.datatype = 'double';
Resource.ImageBuffer.numFrames = 10;
Resource.DisplayWindow(1).Title = 'C5-2gHWideBeamPolar';
Resource.DisplayWindow(1).pdelta = 0.45;
ScrnSize = get(0,'ScreenSize');
pdeltaR = PData.PDelta(2);
pdeltaT = PData.PDelta(1);
DwWidth = ceil((2*PData.Size(1)*pdeltaR*sin(pdeltaT*(PData.Size(2)-1)/2))/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil((P.endDepth + radius - (radius*cos(scanangle/2)))/Resource.DisplayWindow(1).pdelta)+10;
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [-DwWidth/2*Resource.DisplayWindow(1).pdelta,0, ...
                                         -(radius-radius * cos(scanangle/2))];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).numFrames = numImageFrames;
Resource.DisplayWindow.Colormap = gray(256);

% Specify Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,.67,2,1];

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'focus', P.txFocus, ...
                   'Steer', [0.0,0.0], ...
                   'Apod', zeros(1,Trans.numelements), ...
                   'Delay', zeros(1,Trans.numelements), ...
                   'TXPD', [], ...
                   'peakCutOff', 1.0, ...
                   'peakBLMax', 4.0), 1, P.numRays);

% - Set event specific TX attributes.
%    floor(P.numTx/2) is the number of elements to include on each side of the
%    transmit origin center element, for the specified focus and sensitivity cutoff.
if floor(P.numTx/2) == P.numTx/2, L = P.numTx+1; else, L = P.numTx; end
W = kaiser(L,2.0);
Ce = zeros(1,P.numRays);
for n = 1:P.numRays   % numRays transmit events
    % Set transmit Origins.
    TX(n).Origin = [radius*sin(Angle(n)), 0.0, radius*cos(Angle(n))-radius];
    Ce(n) = round(1+127*(Angle(n) - theta)/scanangle);
    % Set transmit Apodization so that a maximum of numTx transmitters are active.
    w1 = 1;
    w2 = L;
    lft = Ce(n) - floor(P.numTx/2);
    if lft < 1, w1 = -lft+2; lft = 1; end
    rt = Ce(n) + floor(P.numTx/2);
    if rt > Trans.numelements, w2 = L - (rt - Trans.numelements); rt = Trans.numelements; end
    TX(n).Apod(lft:rt) = W(w1:w2);
    TX(n).Delay = computeTXDelays(TX(n));
    TX(n).TXPD = computeTXPD(TX(n),PData);
end

% Specify Receive structure arrays.
% - We need numRays Receives for each frame.
% -- Compute the maximum receive path length, using the law of cosines.
maxAcqLength = ceil(sqrt((P.endDepth+radius)^2 + radius^2 - 2*(P.endDepth+radius)*radius*cos(scanangle)));
Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength+50, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode','NS200BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, P.numRays*Resource.RcvBuffer(1).numFrames);
% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    Receive(P.numRays*(i-1)+1).callMediaFunc = 1;
    for j = 1:P.numRays
        Receive(P.numRays*(i-1)+j).framenum = i;
        Receive(P.numRays*(i-1)+j).acqNum = j;
    end
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [153,308,410,520,605,665,705,760];
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
                   'Pre',[], ...
                   'Post',[], ...
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'scaleFactor', 0.25, ...
                   'regionnum', 0), 1, P.numRays);
% - Set specific ReconInfo attributes.
ReconInfo(1).Pre = 'clearInterBuf';     % clear entire Interbuffer frame
for i = 1:P.numRays
    ReconInfo(i).txnum = i;
    ReconInfo(i).rcvnum = i;
    ReconInfo(i).regionnum = i;
end
ReconInfo(P.numRays).Post = 'IQ2IntensityImageBuf';  % detect IQ data

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
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',45,...
                         'display',1,...      % display image after processing
                         'displayWindow',1};

% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; %  jump back to start.
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';
SeqControl(2).argument = 250; % 250 usec/ray * 128 = 32msec
SeqControl(3).command = 'timeToNextAcq';
SeqControl(3).argument = 1000; % 1 msec
SeqControl(4).command = 'returnToMatlab';
nsc = 5;

% Specify Event structure arrays.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    for j = 1:P.numRays                      % Acquire rays
        Event(n).info = 'Acquire ray line';
        Event(n).tx = j;
        Event(n).rcv = P.numRays*(i-1)+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;
    end
    Event(n-1).seqControl = 3; % time between frames

    Event(n).info = 'Transfer frame to host.';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = nsc;
       SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
       nsc = nsc+1;
    n = n+1;

    Event(n).info = 'recon and process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 1;
    Event(n).seqControl = 0;
    if (floor(i/4) == i/4)&&(i ~= Resource.RcvBuffer(1).numFrames)  % Exit to Matlab every 4th frame
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
    'SliderStep',[0.025,0.1],'ValueFormat','%1.3f',....
    'Callback',@SensCutoffCallback);

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

% - numTx
UI(3).Control = VsSliderControl('LocationCode','UserB3',...
    'Label','numTx',...
    'SliderMinMaxVal',[10,62,P.numTx],...
    'SliderStep',[4/54,8/54],'ValueFormat','%3.0f',...
    'Callback',@numTxCallback);

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 5;

% Save all the structures to a .mat file.
filename = fullfile( vsv.file.getVSXDir, 'MatFiles/SetUpC5_2gHWideBeamPolar_QSApp.mat');

save(filename);

%% Automatic VSX Execution:
% Uncomment the following line to automatically run VSX every time you run
% this SetUp script (note that if VSX finds the variable 'filename' in the
% Matlab workspace, it will load and run that file instead of prompting the
% user for the file to be used):

% filename = 'C5-2gHWideBeamPolar';  VSX;


%% **** Callback routines used by UIControls (UI) ****

function SensCutoffCallback(~, ~, UIValue)
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
    theta = -(scanangle/2);
    dtheta = scanangle/(P.numRays-1);
    Angle = theta:dtheta:(-theta);
    % Modify PData for new range
    PData = evalin('base','PData');
    PData(1).Size(1) = 10+ceil((P.endDepth + radius)/PData.PDelta(2));
    for j = 1:P.numRays
        PData(1).Region(j).Shape = struct( ...
                           'Name','Sector',...
                           'Position',PData.Origin,...
                           'r1',radius+P.startDepth,...
                           'r2',radius+P.endDepth,...
                           'angle',dtheta*7,...
                           'steer',Angle(j));
        PData(1).Region(j).PixelsLA = [];
        PData(1).Region(j).numPixels = 0;
    end
    assignin('base','PData',PData);
    evalin('base','PData.Region = computeRegions(PData);');
    DwHeight = ceil((P.endDepth + radius - (radius*cos(scanangle/2)))/Resource.DisplayWindow(1).pdelta)+10;
    assignin('base','DwHeight',DwHeight);
    evalin('base','Resource.DisplayWindow(1).Position(4) = DwHeight;');
    % Update TXPD data of TX structures.
    TX = evalin('base','TX');
    PData = evalin('base','PData');
    for i = 1:size(TX,2)
        TX(i).TXPD = computeTXPD(TX(i),PData);
    end
    assignin('base','TX',TX);
    % Update Receive structures
    Receive = evalin('base', 'Receive');
    maxAcqLength = ceil(sqrt((P.endDepth+radius)^2 + radius^2 - 2*(P.endDepth+radius)*radius*cos(scanangle)));
    for i = 1:size(Receive,2)
        Receive(i).endDepth = maxAcqLength+50;
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

function numTxCallback(hObject, ~, UIValue)
    simMode = evalin('base','Resource.Parameters.simulateMode');
    % No numTx change if in simulate mode 2.
    if simMode == 2
        set(hObject,'Value',evalin('base','P.numTX'));
        return
    end
    P = evalin('base','P');
    P.numTx = UIValue;
    assignin('base','P',P);
    PData = evalin('base','PData');
    TX = evalin('base', 'TX');
    Angle = evalin('base','Angle');
    Trans = evalin('base','Trans');
    radius = Trans.radius;
    theta = evalin('base','theta');
    for n = 1:P.numRays
        TX(n).Origin = [radius*sin(Angle(n)), 0.0, radius*cos(Angle(n))-radius];
        ce = round(1+127*(Angle(n) - theta)/(-2*theta));
        % Set transmit Apodization so that a maximum of numTx + 1 transmitters are active.
        lft = round(ce - P.numTx/2);
        if lft < 1, lft = 1; end
        rt = round(ce + P.numTx/2);
        if rt > Trans.numelements, rt = Trans.numelements; end
        TX(n).Apod = zeros(1,Trans.numelements);
        TX(n).Apod(lft:rt) = 1.0;
        TX(n).Delay = computeTXDelays(TX(n));
        TX(n).TXPD = computeTXPD(TX(n),PData(1));
    end
    assignin('base','TX', TX);
    % Set Control command to update TX
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'TX'};
    assignin('base','Control', Control);
end