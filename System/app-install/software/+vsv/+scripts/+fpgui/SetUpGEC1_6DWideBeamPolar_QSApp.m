% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use.
%
% File name: SetUpGEC1_6DWideBeamPolar.m - Example of curved array imaging with
%            wide beam transmits, and polar coordinate PData
%
% Description:
%   Sequence programming file for GEC1-6D curved array using m wide beam
%   transmits and receive acquisitions. Of the 192 transmit channels, only
%   P.numTx are used, with P.numTx/2 transmitters on each side of the center element
%   (where possible). All 192 receive channels are used, although the element
%   sensitivity cutoff will limit the useful aperture. The receive acquisitions
%   use 200% bandwidth. Reconstruction and Processing are asynchronous with
%   respect to acquisition.
%
%   The PData regions are generated in polar coordinates to provide an overlap of at least
%   9 regions at every pixel. The transmit focus is set to below the bottom of the PData region.
%
% Last update:
% 7/23/2019 - update and test with SW 4.1.0

%clear[ 	]+all

mm2wl = 3.9/1.54;  % scale factor to convert mm to wavelengths for Trans.frequency = 3.9 MHz.
P.startDepth = 2*mm2wl;  % startDepth (2 mm)
P.endDepth = 120*mm2wl;  % endDepth (10 cm)
P.numTx = 40;    % no. of elements in TX aperture (-1).
P.numRays = 96; % no. of rays in frame
P.txFocus = 2*P.endDepth; % focus in wavelengths
P.maxDepthMm = 150;  % maxDepth for RangeChange and RcvBuffer
P.maxDepth = P.maxDepthMm*mm2wl;   % maxDepth for RangeChange and RcvBuffer

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
Trans.name = 'GEC1-6D';
Trans.units = 'wavelengths'; % required in Gen3 to prevent default to mm units
Trans = computeTrans(Trans); % GEC1-6D transducer is 'known' transducer so we can use computeTrans.
radius = Trans.radius;
scanangle = (Trans.numelements-1)*Trans.spacing/radius;
theta = -(scanangle/2);      % angle to left edge from centerline
dtheta = scanangle/(P.numRays-1);
Angle = theta:dtheta:(-theta);

% Specify PData structure array.
PData.Coord = 'polar';
PData.PDelta = [scanangle/192, 0.5, 0]; % [theta,r,z] theta increment is equivalent to 192 ray line scan
sizeRows = 10 + ceil((P.endDepth + radius)/PData.PDelta(2));
sizeCols = 192;
PData.Size = [sizeRows,sizeCols,1];     % size, origin and pdelta set region of interest.
PData.Origin = [0,0,-radius];   % Origin is in rectangular coords.
% Define PData Regions for numRays scanlines
for j = 1:P.numRays
    PData.Region(j) = struct('Shape',struct( ...
                       'Name','Sector',...
                       'Position',PData.Origin,...
                       'r1',radius+P.startDepth,...
                       'r2',radius+P.endDepth,...
                       'angle',dtheta*9,...
                       'steer',Angle(j)));
end
PData.Region = computeRegions(PData);

%  Media points for curved array.
% - Uncomment for speckle
% Media.numPoints = (20000);
% Media.MP = rand(Media.numPoints,4);
% Media.MP(:,2) = 0;
% Media.MP(:,4) = 0.01 + 0.04*Media.MP(:,4);  % Random low amplitude
% RandR = P.endDepth*scaleToWvl*Media.MP(:,1)+radius;
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
Media.MP(24,:) = [P.startDepth*sin(-0.2608),0,P.startDepth*cos(-0.2608)-radius,1.0];
Media.MP(25,:) = [P.startDepth*sin(0.2608),0,P.startDepth*cos(0.2608)-radius,1.0];
Media.MP(26,:) = [P.startDepth*sin(-0.5267),0,P.startDepth*cos(-0.5267)-radius,1.0];
Media.MP(27,:) = [P.startDepth*sin(0.5267),0,P.startDepth*cos(0.5267)-radius,1.0];
Media.MP(28,:) = [0,0,220,1.0];
Media.numPoints = 28;
Media.function = 'movePoints';
Media.attenuation = -0.5;

% Specify Resources.
maxBufLength = ceil(sqrt((P.maxDepth+radius)^2 + radius^2 - ...
                     2*(P.maxDepth+radius)*radius*cos(scanangle)));
maxBufSizePerAcq = 128*ceil(maxBufLength*8*(Trans.frequency/Trans.frequency)/128); % 8 to match the receive sample mode
Resource.RcvBuffer.datatype = 'int16';
Resource.RcvBuffer.rowsPerFrame = P.numRays*maxBufSizePerAcq; % allows for max range
Resource.RcvBuffer.colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer.numFrames = 12;
Resource.InterBuffer.datatype = 'complex double';
Resource.InterBuffer.numFrames = 1;
Resource.ImageBuffer.datatype = 'double';
Resource.ImageBuffer.numFrames = 12;
Resource.DisplayWindow(1).Title = 'GEC1-6DWideBeamPolar';
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
Resource.DisplayWindow(1).numFrames = numImageFrames;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow.Colormap = gray(256);

% Specify Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,.67,2,1];

% Specify TX structure array.
% - We need P.numRays transmit specifications.
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
%    P.numTx/2 is the number of elements to include on each side of the
%    center element, for the specified focus and sensitivity cutoff.
%    Thus the full transmit aperture will be P.numTx + 1 elements.
% - Straight ahead beams
for n = 1:P.numRays   % numRays transmit events
    % Set transmit Origins.
    TX(n).Origin = [radius*sin(Angle(n)), 0.0, radius*cos(Angle(n))-radius];
    ce = round(1+Trans.numelements*(Angle(n) - theta)/scanangle);
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

% calculate TXPD and TX delay
h = waitbar(0,'Program TX parameters, please wait!');
steps = P.numRays;
for j = 1:steps
    TX(j).Delay = computeTXDelays(TX(j));
    TX(j).TXPD = computeTXPD(TX(j),PData);
    waitbar(j/steps)
end
close(h)

% Specify Receive structure arrays.
% -- Compute the maximum receive path length.
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
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, P.numRays*Resource.RcvBuffer(1).numFrames);
% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    k = P.numRays*(i-1);
    Receive(k+1).callMediaFunc = 1;
    for j = 1:P.numRays
        Receive(k+j).framenum = i;
        Receive(k+j).acqNum = j;
    end
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [0,250,410,520,605,665,725,800];
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
for j = 1:P.numRays
    ReconInfo(j).txnum = j;  % beams with no steering
    ReconInfo(j).rcvnum = j;
    ReconInfo(j).regionnum = j;
end
ReconInfo(P.numRays).Post = 'IQ2IntensityImageBuf';  % detect IQ data

% Specify Process structure array.
%   Processing structure specifies how to display image.
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,... % frame number in src buffer (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure (defines output figure).
                         'norm',1,...        % normalization method(1 means fixed)
                         'pgain',3,...            % pgain is image processing gain
                         'reject',2,...
                         'persistMethod','simple',...
                         'persistLevel',10,...
                         'interpMethod','4pt',...
                         'averageMethod','none',...
                         'compressMethod','power',...      % X^0.5 normalized to output word size
                         'compressFactor',50,...
                         'mappingMethod','full',...
                         'display',1,...      % display image after processing
                         'displayWindow',1};

% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; %  jump back to start.
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq'; % time between rays
SeqControl(2).argument = 200; % 200usec
SeqControl(3).command = 'timeToNextAcq';  % time between frames
SeqControl(3).argument = 3050;  % 3000 usec = 3msec time between frames
SeqControl(4).command = 'returnToMatlab';
nsc = 5;

% Specify Event structure arrays.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    k = P.numRays*(i-1);
    for j = 1:P.numRays          % Acquire frame
        Event(n).info = 'Acquire ray lines';
        Event(n).tx = j;
        Event(n).rcv = k+j;
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
    MinMaxVal = [MinMaxMm, P.endDepth/mm2wl];
else
    AxesUnit = 'wls';
    MinMaxVal = [MinMaxMm*mm2wl, P.endDepth];
end

UI(2).Control = VsSliderControl('LocationCode','UserA1',...
    'Label',['Range (',AxesUnit,')'],'SliderMinMaxVal',MinMaxVal,...
    'SliderStep',[0.1,0.2],'ValueFormat','%3.0f',...
    'Callback',@RangeChangeCallback);

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 12;

% Save all the structures to a .mat file.
filename = fullfile( vsv.file.getVSXDir, 'MatFiles/SetUpGEC1_6DWideBeamPolar_QSApp.mat');

save(filename);

% Automatic VSX Execution:
% Uncomment the following line to automatically run VSX every time you run
% this SetUp script (note that if VSX finds the variable 'filename' in the
% Matlab workspace, it will load and run that file instead of prompting the
% user for the file to be used):

% filename = 'GEC1-6DWideBeamPolar';  VSX;


% **** Callback functions that used by UI Controls. ****
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
    P.txFocus = 2*P.endDepth;
    assignin('base','P',P);

    scanangle = evalin('base','scanangle');
    radius = Trans.radius;
    theta = -(scanangle/2);
    dtheta = scanangle/(P.numRays-1);
    Angle = theta:dtheta:(-theta);
    % Modify PData for new range
    PData = evalin('base','PData');
    PData.Size(1) = 10+ceil((P.endDepth + radius)/PData.PDelta(2));
    for j = 1:P.numRays
        PData.Region(j).Shape = struct( ...
                           'Name','Sector',...
                           'Position',PData.Origin,...
                           'r1',radius+P.startDepth,...
                           'r2',radius+P.endDepth,...
                           'angle',dtheta*9,...
                           'steer',Angle(j));
        PData.Region(j).PixelsLA = [];
        PData.Region(j).numPixels = 0;
    end
    PData.Region = computeRegions(PData);
    % Update TXPD data of TX structures.
    TX = evalin('base','TX');
    % calculate TXPD and TX delay
    h = waitbar(0,'Programming TX parameters, please wait!');
    steps = P.numRays;
    for ind = 1:steps
        TX(ind).focus = P.txFocus;
        TX(ind).Delay = computeTXDelays(TX(ind));
        TX(ind).TXPD = computeTXPD(TX(ind),PData);
        waitbar(ind/steps)
    end
    close(h)
    assignin('base','TX',TX);
    assignin('base','PData',PData);
    % Update Receive structures
    Receive = evalin('base', 'Receive');
    maxAcqLength = ceil(sqrt((P.endDepth+radius)^2 + radius^2 - 2*(P.endDepth+radius)*radius*cos(scanangle)));
    for i = 1:size(Receive,2)
        Receive(i).endDepth = maxAcqLength;
    end
    assignin('base','Receive',Receive);
    evalin('base','TGC.rangeMax = P.endDepth;');
    evalin('base','TGC.Waveform = computeTGCWaveform(TGC);');
    % Update DisplayWindow
    DwHeight = ceil((P.endDepth + radius - (radius*cos(scanangle/2)))/Resource.DisplayWindow(1).pdelta)+10;
    assignin('base','DwHeight',DwHeight);
    evalin('base','Resource.DisplayWindow(1).Position(4) = DwHeight;');
    % Set Control structure
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'PData','InterBuffer','ImageBuffer','TX','Receive','TGC','Recon','DisplayWindow'};
    assignin('base','Control', Control);
    assignin('base', 'action', 'displayChange');
end
