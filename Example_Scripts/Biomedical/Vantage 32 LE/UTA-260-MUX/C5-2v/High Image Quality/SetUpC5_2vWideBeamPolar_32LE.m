% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: SetUpC5_2vWideBeamPolar_32LE.m - Example of curved array imaging with
%                             wide beam transmits and polar coordinate PData.
% Description:
%   Sequence programming file for C5-2v curved array using P.numRays wide beam
%   transmits and receive acquisitions. Of the 64 transmit channels,
%   only P.numTx are used, with P.numTx/2 transmitters on each side of the center
%   element (where possible). The receive aperture always covers a 32 element
%   aperture centered about the beam origin (where possible).
%   Processing is asynchronous with respect to acquisition.
%
% Last update:
%  03/18/2019 - updated for use with 4.0

clear all

P.startDepth = 5; % startDepth and endDepth in wavelengths
P.endDepth = 280;
P.numRays = 64;   % no. of transmit beams to program
P.numTx = 36;     % no. of elements in transmit aperture - 1, should not larger than 62
P.txFocus = 4*P.endDepth;  % focal point in wavelengths

% Specify system parameters.
Resource.Parameters.speedOfSound = 1540;
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.
% Resource.System.Product = 'Vantage64';
Resource.System.SoftwareVersion = [4 0 0]; % Minimum software release for this script.
Resource.System.UTA = '260-MUX'; % This script requires the 260-MUX UTA.
Resource.Parameters.numTransmit = 64;

% Specify Trans structure array.
Trans.name = 'C5-2v';
Trans.units = 'wavelengths'; % required in Gen3 to prevent default to mm units
Trans = computeTrans(Trans);  % C5-2v transducer is 'known' transducer so we can use computeTrans.
Trans.maxHighVoltage = 50;  % set maximum high voltage limit for pulser supply.
Trans = computeUTAMux64(Trans); % Add HVMux field for use with UTA 260-Mux
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
Resource.DisplayWindow(1).Title = 'C5-2vWideBeamPolar_32LE';
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
Resource.DisplayWindow(1).numFrames = 20;
Resource.DisplayWindow.Colormap = gray(256);

% Specify Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,.67,2,1];

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
                   'peakBLMax', 4.0), 1, P.numRays);
% - Set event specific TX attributes.
%    numTx/2 is the number of elements to include on each side of the
%    center element, for the specified focus and sensitivity cutoff.
%    Thus the full transmit aperture will be numTx + 1 elements.
%   We will specify a MUX aperture the same size
%   as the number of TX & Receive channels, eliminating the time for the MUX to switch to a different aperture.
for n = 1:P.numRays   % numRays transmit events
    TX(n).waveform = 1;  % Set transmit waveform
    % Set transmit Origins.
    TX(n).Origin = [radius*sin(Angle(n)), 0.0, radius*cos(Angle(n))-radius];
    ce = round(1+127*(Angle(n) - theta)/scanangle);
    % Compute 64 channel transmit aperture centered around center element, ce
    ActiveElements = zeros(1,Trans.numelements);
    lft = ce - 31;
    if lft < 1, lft = 1; end
    rt = ce + 32;
    if rt > Trans.numelements, rt = Trans.numelements; end
    ActiveElements(lft:rt) = 1;
    % Compute MUX aperture for all active elements.
    TX(n).aperture = computeMuxAperture(ActiveElements, Trans);
    % Set transmit Apodization so that a maximum of numTx + 1 transmitters are active.
    lft = round(ce - P.numTx/2);
    if lft < 1, lft = 1; end
    rt = round(ce + P.numTx/2);
    if rt > Trans.numelements, rt = Trans.numelements; end
    TX(n).Apod(lft:rt) = 1.0;
    [RIndices,CIndices,V] = find(TX(n).Apod);
    V = kaiser(size(V,2),1);
    TX(n).Apod(CIndices) = V;
    TX(n).Delay = computeTXDelays(TX(n));
    TX(n).TXPD = computeTXPD(TX(n),PData);
end

% Specify Receive structure arrays.
% - We need numRays Receives for each frame.
% -- Compute the maximum receive path length, using the law of cosines.
maxAcqLength = ceil(sqrt((P.endDepth+radius)^2 + radius^2 - 2*(P.endDepth+radius)*radius*cos(scanangle)));
Receive = repmat(struct('Apod', zeros(1,Trans.numelements), ...
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
        Receive(P.numRays*(i-1)+j).aperture = TX(j).aperture; % mux aperture same as transmit
        ce = round(1+127*(Angle(j) - theta)/scanangle);
        % Compute 32 channel transmit aperture centered around center element, ce
        lft = ce - 15;
        if lft < 1, lft = 1; end
        rt = ce + 16;
        if rt > Trans.numelements, rt = Trans.numelements; end
        Receive(P.numRays*(i-1)+j).Apod(lft:rt) = 1;
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
% - Sensitivity Cutoff
UI(1).Control =  {'UserB7','Style','VsSlider','Label','Sens. Cutoff',...
                  'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                  'SliderStep',[0.025,0.1],'ValueFormat','%1.3f'};
UI(1).Callback = text2cell('%SensCutoffCallback');

% - Range Change
wls2mm = 1;
AxesUnit = 'wls';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
        AxesUnit = 'mm';
        wls2mm = Resource.Parameters.speedOfSound/1000/Trans.frequency;
    end
end
UI(2).Control = {'UserA1','Style','VsSlider','Label',['Range (',AxesUnit,')'],...
                 'SliderMinMaxVal',[64,300,P.endDepth]*wls2mm,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f'};
UI(2).Callback = text2cell('%RangeChangeCallback');

% - Transmit focus change
UI(3).Control = {'UserB4','Style','VsSlider','Label',['TX Focus (',AxesUnit,')'],...
                 'SliderMinMaxVal',[P.endDepth,2*P.txFocus,P.txFocus],'SliderStep',[10/500,20/500],'ValueFormat','%3.0f'};
UI(3).Callback = text2cell('%TxFocusCallback');

% - numTx
UI(4).Control = {'UserB3','Style','VsSlider','Label','numTx',...
                  'SliderMinMaxVal',[10,62,P.numTx],...
                  'SliderStep',[4/54,8/54],'ValueFormat','%3.0f'};
UI(4).Callback = text2cell('%numTxCallback');

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 5;

% Save all the structures to a .mat file.
save('MatFiles/C5-2vWideBeamPolar_32LE');
return


% **** Callback routines to be converted by text2cell function. ****
%SensCutoffCallback - Sensitivity cutoff change
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
%SensCutoffCallback

%RangeChangeCallback - Range change
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
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
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
return
%RangeChangeCallback

%TxFocusCallback - TX focus changel
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
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
        P.txFocus = UIValue*scaleToWvl;
    end
end
assignin('base','P',P);
PData = evalin('base','PData');
TX = evalin('base', 'TX');
for n = 1:size(TX,2)   % for all transmit events
    TX(n).focus = P.txFocus;
    TX(n).Delay = computeTXDelays(TX(n));
    TX(n).TXPD = computeTXPD(TX(n),PData);
end
assignin('base','TX', TX);
% Set Control command to update TX
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'TX'};
assignin('base','Control', Control);
return
%TxFocusCallback

%numTxCallback - number of transmitters in TX aperture
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
scanangle = evalin('base','scanangle');
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
return
%numTxCallback
