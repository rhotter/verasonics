% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: SetUpC5_2vWideBeamPolarSC.m - Example of curved array imaging with wide
%                beam transmits and spatial compounding using polar coordinate PData.
% Description:
%   Sequence programming file for C5-2v curved array using wide beam transmits
%   and r/theta scanning with spatial compounding. Of the 128 transmit channels, only
%   P.numTx are used, with P.numTx transmitters on each side of the center element
%   (where possible). All 128 receive channels are used, although the
%   element sensitivity cutoff will limit the useful aperture. The receive acquisitions
%   use 200% bandwidth. Reconstruction and Processing are asynchronous with
%   respect to acquisition.
%
%   Three wide beam scans are performed using P.numRays scan lines in each frame,
%   with different transmit steering directions (steered left, no steering, steered right).  The
%   frames are averaged after image processing using a running average that updates every frame.
%
% Last update:
% 09/11/2020 - Update to SW 4.3 format for new user UIControls and External function definitions (VTS-1691).
%   More info:(?/Example_Scripts/Vantage_Features/New UI Scheme/SetUpL11_5vFlash_NewUI)

clear all
P.startDepthMm = 2; % startDepth in mm
P.endDepthMm = 180; % endDepth in mm
P.numRays = 64;   % no. of transmit beams to program
P.txFocusMm = 600; % focus in mm - set to larger than endDepth
P.numTx = 30;     % no. of elements in transmit aperture - 1

% Specify system parameters.
Resource.Parameters.numTransmit = 128;  % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;  % number of receive channels.
Resource.Parameters.speedOfSound = 1540;
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Specify Trans structure array.
Trans.name = 'C5-2v';
Trans.units = 'wavelengths'; % required in Gen3 to prevent default to mm units
Trans = computeTrans(Trans);  % C5-2v transducer is 'known' transducer so we can use computeTrans.
Trans.maxHighVoltage = 50;  % set maximum high voltage limit for pulser supply.
radius = Trans.radius;  % radius of curved array in wavelengths
scanangle = (Trans.numelements-16)*Trans.spacing/radius; % set scanangle smaller than array angle
dtheta = scanangle/(P.numRays-1);
theta = -(scanangle/2); % angle to left edge from centerline
Angle = theta:dtheta:(-theta);
mm2wl = 1000*Trans.frequency/Resource.Parameters.speedOfSound;

startDepth = P.startDepthMm*mm2wl;
endDepth = P.endDepthMm*mm2wl;
txFocus = P.txFocusMm*mm2wl;

% Specify PData structure array.
PData.Coord = 'polar';
PData.PDelta = [scanangle/127, 0.5, 0]; % [theta,r,z] theta increment is equivalent to 128 ray line scan
sizeRows = 10 + ceil((endDepth + radius)/PData.PDelta(2));
sizeCols = 128;
PData.Size = [sizeRows,sizeCols,1];     % size, origin and pdelta set region of interest.
PData.Origin = [0,0,-radius];   % Origin is in rectangular coords.
% Define PData Regions for numRays scanlines
for j = 1:P.numRays
    PData.Region(j) = struct('Shape',struct( ...
                       'Name','Sector',...
                       'Position',PData.Origin,...
                       'r1',radius+startDepth,...
                       'r2',radius+endDepth,...
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
Media.MP(24,:) = [0,0,220,1.0];
Media.MP(25,:) = [0,0,250,1.0];
Media.MP(26,:) = [(radius+250)*sin(-0.2608),0,(radius+250)*cos(-0.2608)-radius,1.0];
Media.MP(27,:) = [(radius+250)*sin(0.2608),0,(radius+250)*cos(0.2608)-radius,1.0];
Media.MP(28,:) = [(radius+250)*sin(-0.5267),0,(radius+250)*cos(-0.5267)-radius,1.0];
Media.MP(29,:) = [(radius+250)*sin(0.5267),0,(radius+250)*cos(0.5267)-radius,1.0];
Media.MP(30,:) = [startDepth*sin(-0.2608),0,startDepth*cos(-0.2608)-radius,1.0];
Media.MP(31,:) = [startDepth*sin(0.2608),0,startDepth*cos(0.2608)-radius,1.0];
Media.MP(32,:) = [startDepth*sin(-0.5267),0,startDepth*cos(-0.5267)-radius,1.0];
Media.MP(33,:) = [startDepth*sin(0.5267),0,startDepth*cos(0.5267)-radius,1.0];
Media.function = 'movePoints';

% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 3*P.numRays*4096; % size for 3 angles, max depth
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 10;
Resource.RcvBuffer(2).datatype = 'int16';
Resource.RcvBuffer(2).rowsPerFrame = 3*P.numRays*4096; % size for 3 angles, max depth
Resource.RcvBuffer(2).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(2).numFrames = 10;  % all RcvBuffers should have same no. of frames
Resource.RcvBuffer(3).datatype = 'int16';
Resource.RcvBuffer(3).rowsPerFrame = 3*P.numRays*4096; % size for 3 angles, max depth
Resource.RcvBuffer(3).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(3).numFrames = 10;
Resource.InterBuffer.datatype = 'complex double';
Resource.InterBuffer.numFrames = 1;
Resource.ImageBuffer.datatype = 'double';
Resource.ImageBuffer.numFrames = 20;
Resource.DisplayWindow(1).Title = 'C5-2vWideBeamPolarSC';
Resource.DisplayWindow(1).pdelta = 0.45;
ScrnSize = get(0,'ScreenSize');
pdeltaR = PData.PDelta(2);
pdeltaT = PData.PDelta(1);
DwWidth = ceil((2*PData.Size(1)*pdeltaR*sin(pdeltaT*(PData.Size(2)-1)/2))/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil((endDepth + radius - (radius*cos(scanangle/2)))/Resource.DisplayWindow(1).pdelta)+10;
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [-DwWidth/2*Resource.DisplayWindow(1).pdelta,0, ...
                                         -(radius-radius * cos(scanangle/2))];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).numFrames = 20;
Resource.DisplayWindow.Colormap = gray(256);

% Specify Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,.67,3,1];

% Specify TX structure array.
% - We need 128 transmit specifications.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'focus', txFocus, ...
                   'Steer', [0.0,0.0], ...
                   'Apod', zeros(1,Trans.numelements), ...
                   'Delay', zeros(1,Trans.numelements), ...
                   'TXPD', [], ...
                   'peakCutOff', 0.5, ...
                   'peakBLMax', 4.0), 1, 3*P.numRays);
% - Set event specific TX attributes.
%    numTx/2 is the number of elements to include on each side of the
%    center element, for the specified focus and sensitivity cutoff.
%    Thus the full transmit aperture will be numTx + 1 elements.
for n = 1:P.numRays   % numRays transmit events
    TX(n).waveform = 1;  % Set transmit waveform
    % Set transmit Origins.
    TX(n).Origin = [radius*sin(Angle(n)), 0.0, radius*cos(Angle(n))-radius];
    ce = 8 + round(112*(Angle(n) - theta)/scanangle);
    % Set transmit Apodization so that a maximum of numTx + 1 transmitters are active.
    lft = round(ce - P.numTx/2);
    if lft < 1, lft = 1; end
    rt = round(ce + P.numTx/2);
    if rt > Trans.numelements, rt = Trans.numelements; end
    TX(n).Apod(lft:rt) = 1.0;
    [RIndices,CIndices,V] = find(TX(n).Apod);
    V = kaiser(size(V,2),1);
    TX(n).Apod(CIndices) = V;
end
m = P.numRays;
for n = 1:P.numRays   % P.numRays steered left transmit events
    TX(n+m).Origin = TX(n).Origin;
    TX(n+m).Apod = TX(n).Apod;
    if n<=10
        TX(n+m).Steer = [-((n-1)/10)*dtheta,0.0];
    elseif n>(P.numRays-10)
        TX(n+m).Steer = [-((P.numRays-n)/10)*dtheta,0.0];
    else
        TX(n+m).Steer = [-dtheta,0.0];
    end
end
m = m + P.numRays;
for n = 1:P.numRays   % P.numRays steered right transmit events
    TX(n+m).Origin = TX(n).Origin;
    TX(n+m).Apod = TX(n).Apod;
    if n<=10
        TX(n+m).Steer = [((n-1)/10)*dtheta,0.0];
    elseif n>(P.numRays-10)
        TX(n+m).Steer = [((P.numRays-n)/10)*dtheta,0.0];
    else
        TX(n+m).Steer = [dtheta,0.0];
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

% Specify Receive structure arrays.
% - We need numRays Receives for each frame.
% -- Compute the maximum receive path length, using the law of cosines.
maxAcqLength = ceil(sqrt((endDepth+radius)^2 + radius^2 - 2*(endDepth+radius)*radius*cos(scanangle)));
Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', startDepth, ...
                        'endDepth', maxAcqLength+20, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode','NS200BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, 3*P.numRays*Resource.RcvBuffer(1).numFrames);
% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    k = 3*P.numRays*(i-1);
    Receive(k+1).callMediaFunc = 1;
    for j = 1:P.numRays
        Receive(k+j).framenum = i;
        Receive(k+j).acqNum = j;
        Receive(k+P.numRays+j).bufnum = 2;
        Receive(k+P.numRays+j).framenum = i;
        Receive(k+P.numRays+j).acqNum = j;
        Receive(k+2*P.numRays+j).bufnum = 3;
        Receive(k+2*P.numRays+j).framenum = i;
        Receive(k+2*P.numRays+j).acqNum = j;
    end
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [153,308,410,520,605,665,705,760];
TGC.rangeMax = endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

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
                   'scaleFactor', 0.33, ...
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
    ReconInfo(j+k).regionnum = j;
end
ReconInfo(k+P.numRays).Post = 'IQ2IntensityImageBuf';
k = k + P.numRays;
ReconInfo(k+1).Pre = 'clearInterBuf';
for j = 1:P.numRays
    ReconInfo(j+k).txnum = j+k;
    ReconInfo(j+k).rcvnum = j+k;
    ReconInfo(j+k).regionnum = j;
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
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','runAverage3',...
                         'compressMethod','power',...
                         'compressFactor',40,...
                         'display',1,...      % display image after processing
                         'displayWindow',1};

% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; %  jump back to start.
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';
SeqControl(2).argument = 320; % 320 usec/ray * 64 * 3 = 61.4msec
SeqControl(3).command = 'timeToNextAcq';
SeqControl(3).argument = 1000; % 1 msec
SeqControl(4).command = 'returnToMatlab';
nsc = 5;

% Specify Event structure arrays.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    % Acquire frame with straight ahead wide beams
    k = 3*P.numRays*(i-1);
    for j = 1:P.numRays        % Acquire frame
        Event(n).info = 'acquire aperture.';
        Event(n).tx = j;
        Event(n).rcv = j+k;
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
        Event(n).rcv = j+k+P.numRays;
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
        Event(n).rcv = j+k+2*P.numRays;
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
UI(1).Control =  VsSliderControl('LocationCode','UserB7',...
    'Label','Sens. Cutoff',...
    'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
    'SliderStep',[0.025,0.1],'ValueFormat','%1.3f',...
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
    'SliderMinMaxVal',[64,425,endDepth]*wls2mm,...
    'SliderStep',[0.1,0.2],'ValueFormat','%3.0f',...
    'Callback',@RangeChangeCallback);

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 5;

% Save all the structures to a .mat file.
save('MatFiles/C5-2vWideBeamPolarSC');

%% Automatic VSX Execution:
% Uncomment the following line to automatically run VSX every time you run
% this SetUp script (note that if VSX finds the variable 'filename' in the
% Matlab workspace, it will load and run that file instead of prompting the
% user for the file to be used):

% filename = 'C5-2vWideBeamPolarSC';  VSX;


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
        set(hObject,'Value',evalin('base','endDepth'));
        return
    end
    Trans = evalin('base','Trans');
    Resource = evalin('base','Resource');
    mm2wl = 1000*Trans.frequency/Resource.Parameters.speedOfSound;

    P = evalin('base','P');
    endDepth = UIValue;
    if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
        if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
            endDepth = endDepth*mm2wl;
        end
    end
    startDepth = P.startDepthMm*mm2wl;
    P.endDepthMm = endDepth/mm2wl;
    assignin('base','P',P);

    scanangle = evalin('base','scanangle');
    radius = evalin('base','radius');
    theta = -(scanangle/2);
    dtheta = scanangle/(P.numRays-1);
    Angle = theta:dtheta:(-theta);
    % Modify PData for new range
    PData = evalin('base','PData');
    PData(1).Size(1) = 10+ceil((endDepth + radius)/PData.PDelta(2));
    for j = 1:P.numRays
        PData(1).Region(j).Shape = struct( ...
                           'Name','Sector',...
                           'Position',PData.Origin,...
                           'r1',radius+startDepth,...
                           'r2',radius+endDepth,...
                           'angle',dtheta*7,...
                           'steer',Angle(j));
        PData(1).Region(j).PixelsLA = [];
        PData(1).Region(j).numPixels = 0;
    end
    assignin('base','PData',PData);
    evalin('base','PData.Region = computeRegions(PData);');
    DwHeight = ceil((endDepth + radius - (radius*cos(scanangle/2)))/Resource.DisplayWindow(1).pdelta)+10;
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
    maxAcqLength = ceil(sqrt((endDepth+radius)^2 + radius^2 - 2*(endDepth+radius)*radius*cos(scanangle)));
    for i = 1:size(Receive,2)
        Receive(i).endDepth = maxAcqLength+20;
    end
    assignin('base','Receive',Receive);
    evalin('base','TGC.rangeMax = endDepth;');
    evalin('base','TGC.Waveform = computeTGCWaveform(TGC);');
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'PData','InterBuffer','ImageBuffer','TX','Receive','TGC','Recon','DisplayWindow'};
    assignin('base','Control', Control);
    assignin('base', 'action', 'displayChange');
end