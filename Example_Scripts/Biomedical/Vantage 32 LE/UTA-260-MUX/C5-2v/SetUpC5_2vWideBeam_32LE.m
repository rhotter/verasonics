% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use.
%
% File name SetUpC5_2vWideBeam_32LE.m: Example of curved array imaging with
%                                   wide beam transmits
% Generate .mat Sequence Object file for C5_2v Linear array wide beam transmit.
% The transmit aperture size can be set by the variable - numTx. The transmit aperture is
% translated across the array with the number of rays set by the variable - numRays.
% The receive aperture always covers the 32 element aperture.
%
%
% In the following aperture examples, each space represents 2 elements.
% Aperture 1 (P.numTx=32):
%   tttttttttttttttt------------------------------------------------
%   rrrrrrrrrrrrrrrr------------------------------------------------
% Aperture for almost middle wide beam 24 (P.numTx=32, P.numRays=24):
%   ------------------------tttttttttttttttt------------------------
%   ------------------------rrrrrrrrrrrrrrrr------------------------
% Aperture for last wide beam (P.numTx=32, P.numRays=48):
%   ------------------------------------------------tttttttttttttttt
%   ------------------------------------------------rrrrrrrrrrrrrrrr
%
% Last update:
% 07/07/2020 VTS-1755 Fix bugs found in script
% 05/18/2018 Modified for dynamic MUX programming

clear all

P.numTx = 49;   % no. of elements in TX aperture.
P.numRays = 60; % no. of rays in frame
P.txFocus = 600;
P.startDepth = 5;  % startDepth in wavelength
P.endDepth = 240;

% Define system parameters.
Resource.Parameters.speedOfSound = 1540;  % speed of sound in m/s
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.
Resource.System.SoftwareVersion = [3 5 0]; % Minimum software release for this script.
Resource.System.UTA = '260-MUX'; % This script requires the 260-MUX UTA.
Resource.System.transmitChannels = 64;
Resource.System.receiveChannels = 32;

% Specify Trans structure array.
Trans.name = 'C5-2v';
Trans.units = 'wavelengths';
Trans = computeTrans(Trans);
Trans.maxHighVoltage = 50;  % set maximum high voltage limit for pulser supply.
Trans = computeUTAMux64(Trans); % Add HVMux field for use with UTA 260-Mux

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
                       'angle',dtheta*7,...
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
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = P.numRays*4096;
Resource.RcvBuffer(1).colsPerFrame = Resource.System.receiveChannels;
Resource.RcvBuffer(1).numFrames = 30;    % 30 frames stored in RcvBuffer.
Resource.InterBuffer(1).numFrames = 1;   % one intermediate buffer needed.
Resource.ImageBuffer(1).numFrames = 10;
Resource.DisplayWindow(1).Title = 'C5-2vWideBeam_32LE';
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
Resource.DisplayWindow(1).Colormap = gray(256);

% Specify Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,0.67,1,1];

scaleToWvl = 1;
if strcmp(Trans.units, 'mm')
    scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);
end

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'Apod', zeros(1,Resource.System.transmitChannels), ...
                   'focus', P.txFocus, ...
                   'Steer', [0.0,0.0], ...
                   'Delay', zeros(1,Resource.System.transmitChannels), ...
                   'TXPD', [], ...
                   'peakCutOff', 0.9, ...
                   'peakBLMax', 4.0), 1, P.numRays);

% - Set event specific TX attributes.
%   We will specify a MUX aperture the same size
%   as the number of TX & Receive channels, even though the number of active TX channels
%   set by P.numTx is usually smaller.  This will allow the same MUX aperture to be used
%   for receive, eliminating the time for the MUX to switch to a different aperture.
Ce = zeros(1,P.numRays); % Ce gets the nearest element no. to the aperture center.
for n = 1:P.numRays
    TX(n).Origin = [radius*sin(Angle(n)), 0.0, radius*cos(Angle(n))-radius];
    Ce(n) = round(1+127*(Angle(n) - theta)/(-2*theta));
    % Now compute TX(n).Apod based on P.numTX
    L = P.numTx; % length of apodization function
    W = kaiser(L,1.5);
    w1 = 1; % w1 and w2 are the start and end indices of the apod function
    w2 = L;
    lft = Ce(n) - floor(P.numTx/2);
    if lft < 1, w1 = -lft+2; lft = 1; end
    rt = Ce(n) + floor(P.numTx/2);
    if rt > Trans.numelements, w2 = L - (rt - Trans.numelements); rt = Trans.numelements; end
    TX(n).aperture = min(lft,65);
    TX(n).Apod((lft:rt)-TX(n).aperture+1) = W(w1:w2);
    % Compute transmit delays
    TX(n).Delay = computeTXDelays(TX(n));
    % Compute transmit pixel data
    TX(n).TXPD = computeTXPD(TX(n),PData(1));
end

% Specify Receive structure arrays.
maxAcqLength = ceil(sqrt((P.endDepth+radius)^2 + radius^2 - ...
                     2*(P.endDepth+radius)*radius*cos(scanangle)));
Receive = repmat(struct('Apod', [zeros(1,16),ones(1,32),zeros(1,16)], ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, P.numRays*Resource.RcvBuffer(1).numFrames);
% - Set event specific Receive attributes for each frame.
for i = 1:Resource.RcvBuffer(1).numFrames
    Receive(P.numRays*(i-1)+1).callMediaFunc = 1;
    for j = 1:P.numRays
        Receive(P.numRays*(i-1)+j).aperture = TX(j).aperture; % mux aperture same as transmit
        Receive(P.numRays*(i-1)+j).framenum = i;
        Receive(P.numRays*(i-1)+j).acqNum = j;
    end
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [51,245,410,475,549,659,736,791];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Recon structure arrays.
% - We need one Recon structure which will be used for each frame.
Recon = struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame',-1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...
               'RINums',(1:P.numRays));

% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 'accumIQ', ...  % default is to accumulate IQ data.
                   'Pre',[],...
                   'Post',[],...
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'scaleFactor', 0.2, ...
                   'regionnum', 0, ...
                   'threadSync', 1), 1, P.numRays);
% - Set specific ReconInfo attributes.
ReconInfo(1).Pre = 'clearInterBuf';
for j = 1:P.numRays
    ReconInfo(j).txnum = j;
    ReconInfo(j).rcvnum = j;
    ReconInfo(j).regionnum = j;
end
ReconInfo(P.numRays).Post = 'IQ2IntensityImageBuf';

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
                         'compressFactor',40,...
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
SeqControl(3).argument = 20000-(P.numRays-1)*TTNA;  % 20 ms between frames
SeqControl(4).command = 'returnToMatlab';
nsc = 5; % nsc is count of SeqControl objects

% Specify Event structure arrays.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    for j = 1:P.numRays        % Acquire frame
        Event(n).info = 'acquire widebeam';
        Event(n).tx = j;
        Event(n).rcv = P.numRays*(i-1)+j;
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
    if floor(i/3) == i/3     % Exit to Matlab every 3rd frame reconstructed
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
UI(1).Callback = text2cell('%SensCutOffCallback');

% - Range Change
MinMaxVal = [64,300,P.endDepth]; % default unit is wavelength
AxesUnit = 'wls';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
        AxesUnit = 'mm';
        MinMaxVal = MinMaxVal * (Resource.Parameters.speedOfSound/1000/Trans.frequency);
    end
end
UI(2).Control = {'UserA1','Style','VsSlider','Label',['Range (',AxesUnit,')'],...
                 'SliderMinMaxVal',MinMaxVal,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f'};
UI(2).Callback = text2cell('%RangeChangeCallback');

% - Peak CutOff
UI(3).Control = {'UserB2','Style','VsSlider','Label','Peak Cutoff',...
                  'SliderMinMaxVal',[0,20.0,TX(1).peakCutOff],...
                  'SliderStep',[0.005,0.020],'ValueFormat','%1.3f'};
UI(3).Callback = text2cell('%PeakCutOffCallback');

% - Max. Burst Length
UI(4).Control = {'UserB1','Style','VsSlider','Label','Max. BL',...
                  'SliderMinMaxVal',[0,20.0,TX(1).peakBLMax],...
                  'SliderStep',[0.005,0.020],'ValueFormat','%1.3f'};
UI(4).Callback = text2cell('%MaxBLCallback');

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 3;

% Save all the structures to a .mat file.
save('MatFiles/C5-2vWideBeam_32LE');
return

% **** Callback routines to be converted by text2cell function. ****
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
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
        P.endDepth = UIValue*scaleToWvl;
    end
end
assignin('base','P',P);

scanangle = evalin('base','scanangle');
radius = evalin('base','radius');
theta = evalin('base','theta');
dtheta = scanangle/P.numRays; % angle between rays
Angle = theta:dtheta:(-theta);
height = P.endDepth + radius - (radius*cos(scanangle/2));
% Modify PData for new range
PData = evalin('base','PData');
PData.Size(1) = 10+ceil(height/PData(1).PDelta(3));
for n = 1:P.numRays
    PData(1).Region(n).Shape = struct(...
                       'Name','Sector',...
                       'Position',[0,0,-radius],...
                       'r1',radius+P.startDepth,...
                       'r2',radius+P.endDepth,...
                       'angle',dtheta*7,...
                       'steer',Angle(n));
    PData(1).Region(n).PixelsLA = [];
    PData(1).Region(n).numPixels = 0;
end
assignin('base','PData',PData);
evalin('base','PData(1).Region = computeRegions(PData(1));');

% Update TXPD data of TX structures.
TX = evalin('base','TX');
PData = evalin('base','PData');
for i = 1:size(TX,2)
    TX(i).TXPD = computeTXPD(TX(i),PData);
end
assignin('base','TX',TX);
% Update Receive structures.
Receive = evalin('base', 'Receive');
maxAcqLength = ceil(sqrt((P.endDepth+radius)^2 + radius^2 - ...
                     2*(P.endDepth+radius)*radius*cos(scanangle)));
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
return
%RangeChangeCallback

%PeakCutOffCallback
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
%PeakCutOffCallback

%MaxBLCallback
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
%MaxBLCallback
