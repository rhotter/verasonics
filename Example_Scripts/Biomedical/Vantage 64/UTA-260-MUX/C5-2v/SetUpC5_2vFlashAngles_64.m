% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: SetUpC5_2vFlashAngles_64.m - Example of curved array flash imaging
%                                      with steering angle transmits
% Description:
%   Sequence programming file for C5-2v curved array with multiple steering angles and 3-1 synthetic aperture
%   (three acquisitions per line to acquire data from all 128 receive
%   elements). 64 transmit channels and 42 or 44 receive channels
%   are active and positioned as follows (each char represents 2 elements)
%   for each of the 3 synthetic apertures.
%
%   Element Nos.                                                               1
%                               3   4           6         8    9               2
%               1               3   2           5         6    6               8
%   Aperture 1: |               |   |           |         |    |               |
%               tttttttttttttttttttttttttttttttt--------------------------------
%               rrrrrrrrrrrrrrrrrrrrr-------------------------------------------
%               |               |    |          |         |    |               |
%   Aperture 2: |               |    |          |         |    |               |
%               ----------------tttttttttttttttttttttttttttttttt----------------
%               ---------------------rrrrrrrrrrrrrrrrrrrrrr---------------------
%               |               |    |          |         |    |               |
%   Aperture 3: |               |    |          |         |    |               |
%               --------------------------------tttttttttttttttttttttttttttttttt
%               -------------------------------------------rrrrrrrrrrrrrrrrrrrrr
%               |               |    |          |         |    |               |
%
%   The receive data from each of these apertures are stored under
%   different acqNums in the Receive buffer. The reconstruction sums the
%   IQ data from the 3 aquisitions and multiple angles and computes intensity values
%   to produce the full frame. Processing is asynchronous with respect to acquisition.
%
%
% 02/06/2019 updated for use with 4.0 releases

clear all
P.startDepth = 5;  % P.startDepth and P.endDepth are in wavelength
P.endDepth = 240;

na = 7;      % Set na = number of angles.
if (na > 1)
    dtheta = (24*pi/180)/(na-1);
    startAngle = -12*pi/180;
else
    dtheta = 0;
    startAngle=0;
end % set dtheta to range over +/- 12 degrees.

% Specify system parameters.
Resource.Parameters.numTransmit = 64;
Resource.Parameters.speedOfSound = 1540;
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.
Resource.System.SoftwareVersion = [4 0 0]; % Minimum software release for this script.
Resource.System.UTA = '260-MUX'; % This script requires the 260-MUX UTA.

% Specify Trans structure array.
Trans.name = 'C5-2v';
Trans.units = 'wavelengths'; % required in Gen3 to prevent default to mm units
Trans = computeTrans(Trans);  % C5-2v transducer is 'known' transducer so we can use computeTrans.
Trans.maxHighVoltage = 50;  % set maximum high voltage limit for pulser supply.
Trans = computeUTAMux64(Trans); % Add HVMux field for use with UTA 260-Mux
radius = Trans.radius;
scanangle = Trans.numelements*Trans.spacing/radius;
theta = -(scanangle/2); % angle to left edge from centerline

% Specify PData structure array.
PData(1).PDelta = [1.0, 0, 0.5];  % x, y and z pdeltas
sizeRows = 10 + ceil((P.endDepth + radius - (radius * cos(scanangle/2)))/PData(1).PDelta(3));
sizeCols = 10 + ceil(2*(P.endDepth + radius)*sin(scanangle/2)/PData(1).PDelta(1));
PData(1).Size = [sizeRows,sizeCols,1];     % size, origin and pdelta set region of interest.
PData(1).Origin(1,1) = (P.endDepth+radius)*sin(-scanangle/2) - 5;
PData(1).Origin(1,2) = 0;
PData(1).Origin(1,3) = ceil(radius * cos(scanangle/2)) - radius - 5;
PData(1).Region = struct(...
    'Shape',struct('Name','Sector',...
                   'Position',[0,0,-radius],...
                   'r1',radius+P.startDepth,...
                   'r2',radius+P.endDepth,...
                   'angle',scanangle));
PData(1).Region = computeRegions(PData(1));

%  Media points for curved array.
% - Uncomment for speckle
% Media.numPoints = (20000);
% Media.MP = rand(Media.numPoints,4);
% Media.MP(:,2) = 0;
% Media.MP(:,4) = 0.01 + 0.04*Media.MP(:,4);  % Random low amplitude
% RandR = P.endDepth *Media.MP(:,1)+radius;
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
Media.function = 'movePoints';

% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 3*na*4096;
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numTransmit;
Resource.RcvBuffer(1).numFrames = 20;     % 20 frames for rcvDataLoop buffer.
Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).numFrames = 1;  % one intermediate buffer needed.
Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).numFrames = 10;
Resource.DisplayWindow(1).Title = 'C5-2vFlashAngles_64';
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
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,.67,2,1];

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'focus', 0.0, ...
                   'Steer', [0.0,0.0], ...
                   'Apod', zeros(1,Trans.numelements), ...
                   'Delay', zeros(1,Trans.numelements)), 1, 3*na);
% - Set event specific TX attributes.
angle = startAngle;
for j = 1:3:3*na   % na transmit events
    TX(j).Apod(1:64) = 1;
    TX(j).aperture = computeMuxAperture(TX(j).Apod, Trans);
    TX(j).Steer = [angle,0.0];
    TX(j).Delay = computeTXDelays(TX(j));
    TX(j+1).Apod(33:96) = 1;
    TX(j+1).aperture = computeMuxAperture(TX(j+1).Apod, Trans);
    TX(j+1).Steer = [angle,0.0];
    TX(j+1).Delay = computeTXDelays(TX(j+1));
    TX(j+2).Apod(65:128) = 1;
    TX(j+2).aperture = computeMuxAperture(TX(j+2).Apod, Trans);
    TX(j+2).Steer = [angle,0.0];
    TX(j+2).Delay = computeTXDelays(TX(j+2));
    angle = angle + dtheta;
end

% Specify Receive structure arrays.
% - We need na Receives for every frame.
% -- Compute the maximum receive path length, using the law of cosines.
maxAcqLength = ceil(sqrt((P.endDepth+radius)^2 + radius^2 - ...
                     2*(P.endDepth+radius)*radius*cos(scanangle)));
Receive = repmat(struct('Apod', zeros(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, 3*na*Resource.RcvBuffer(1).numFrames);
% - Set event specific Receive attributes for each frame.
for i = 1:Resource.RcvBuffer(1).numFrames
    k = 3*na*(i-1);
    Receive(k+1).callMediaFunc = 1;
    for j = 1:3:3*na
        % -- 1st synthetic aperture acquisition for aperture 1.
        Receive(k+j).Apod(1:42) = 1.0;
        Receive(k+j).aperture = computeMuxAperture(Receive(k+j).Apod, Trans);
        Receive(k+j).framenum = i;
        Receive(k+j).acqNum = j;      % three acquisitions per frame
        % -- 2nd synthetic aperture acquisition for aperture 33.
        Receive(k+j+1).Apod(43:86) = 1.0;
        Receive(k+j+1).aperture = computeMuxAperture(Receive(k+j+1).Apod, Trans');
        Receive(k+j+1).framenum = i;
        Receive(k+j+1).acqNum = j+1;  % three acquisitions per frame
        % -- 3rd synthetic aperture acquisition for aperture 65.
        Receive(k+j+2).Apod(87:128) = 1.0;
        Receive(k+j+2).aperture = computeMuxAperture(Receive(k+j+2).Apod, Trans);
        Receive(k+j+2).framenum = i;
        Receive(k+j+2).acqNum = j+2;  % three acquisitions per frame
    end
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [153,308,410,520,605,665,705,760];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Recon structure arrays.
% - We need one Recon structure which will be used for each frame.
Recon = struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame',-1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...
               'RINums',1:3*na);

% Define ReconInfo structures.
% We need 3*na ReconInfo structures for na steering angles.
ReconInfo = repmat(struct('mode', 4, ...  % accumulate IQ data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 1), 1, 3*na);
% - Set specific ReconInfo attributes.
if na>1
    ReconInfo(1).mode = 'replaceIQ';
    for j = 1:3:3*na  % For each row in the column
        ReconInfo(j).txnum = j;
        ReconInfo(j).rcvnum = j;
        ReconInfo(j+1).txnum = j+1;
        ReconInfo(j+1).rcvnum = j+1;
        ReconInfo(j+2).txnum = j+2;
        ReconInfo(j+2).rcvnum = j+2;
    end
    ReconInfo(3*na).mode = 'accumIQ_replaceIntensity';  % accumulate and detect
else
    ReconInfo(1).mode = 'replaceIntensity';
end

% Specify Process structure array.
pers = 20;
cmpFactor = 40;
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',1.0,...     % pgain is image processing gain
                         'reject',2,...
                         'grainRemoval','none',...
                         'persistMethod','none',...
                         'persistLevel',pers,...
                         'interpMethod','4pt',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',cmpFactor,...
                         'mappingMethod','full',...
                         'display',1,...      % display image after processing
                         'displayWindow',1};

% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; % jump back to start
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(2).argument = 290;  % 290 usec
SeqControl(3).command = 'timeToNextAcq';  % time between frames
SeqControl(3).argument = 10000;  % 10 msec
SeqControl(4).command = 'returnToMatlab';
nsc = 5; % nsc is count of SeqControl objects

% Specify Event structure arrays.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    k = 3*na*(i-1);
    for j = 1:3:3*na                      % Acquire frame
        Event(n).info = '1st aperture.';
        Event(n).tx = j;
        Event(n).rcv = k+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;

        Event(n).info = '2nd aperture.';
        Event(n).tx = j+1;
        Event(n).rcv = k+j+1;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;

        Event(n).info = '3rd aperture.';
        Event(n).tx = j+2;
        Event(n).rcv = k+j+2;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n=n+1;
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
    if (floor(i/5) == i/5)&&(i ~= Resource.RcvBuffer(1).numFrames)  % Exit to Matlab every 5th frame
        Event(n).seqControl = 4;
    else
        Event(n).seqControl = 0;
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

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 5;

% Save all the structures to a .mat file.
save('MatFiles/C5-2vFlashAngles_64');
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

radius = evalin('base','radius');
scanangle = evalin('base','scanangle');
height = P.endDepth + radius - (radius*cos(scanangle/2));
% Modify PData for new range
PData = evalin('base','PData');
PData(1).Size(1) = 10+ceil(height/PData(1).PDelta(3));
PData(1).Region = struct(...
    'Shape',struct('Name','Sector','Position',[0,0,-radius],'r1',radius+P.startDepth,'r2',radius+P.endDepth,'angle',scanangle));
assignin('base','PData',PData);
evalin('base','PData(1).Region = computeRegions(PData(1));');
evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);');
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
Control.Parameters = {'PData','InterBuffer','ImageBuffer','Receive','TGC','Recon','DisplayWindow'};
assignin('base','Control', Control);
assignin('base', 'action', 'displayChange');
return
%RangeChangeCallback
