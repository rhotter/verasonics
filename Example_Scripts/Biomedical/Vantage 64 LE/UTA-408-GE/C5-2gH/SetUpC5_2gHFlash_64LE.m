% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: SetUpC5_2gHFlash_64LE.m - Example of curved array imaging with flash
%                                transmit
% Description:
%   Sequence programming file for C5-2gH curved array, using a flash transmit
%   and 2-1 synthetic aperture receive acquisition. All 128 transmit and 64
%   receive channels are active for each acquisition. Processing is
%   asynchronous with respect to acquisition.
%
% Last update:
%   12/14/2015 Modified for SW 3.0
%   10/04/2021 modified for C5-2gH

clear all
P.startDepth = 5;  % P.startDepth is in wavelength
P.endDepth = 192;  % P.endDepth is in wavelength

% Specify system parameters.
Resource.Parameters.numTransmit = 128;  % number of transmit channels.
Resource.Parameters.numRcvChannels = 64;  % number of receive channels.
Resource.Parameters.speedOfSound = 1540;
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Specify Trans structure array.
Trans.name = 'C5-2gH';
Trans.units = 'wavelengths'; % required in Gen3 to prevent default to mm units
Trans = computeTrans(Trans);  % C5-2 transducer is 'known' transducer so we can use computeTrans.

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
    'Shape',struct('Name','Sector','Position',[0,0,-radius],'r1',radius+P.startDepth,'r2',radius+P.endDepth,'angle',scanangle));
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
Resource.RcvBuffer(1).rowsPerFrame = 4096*2;
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numTransmit;
Resource.RcvBuffer(1).numFrames = 20;     % 20 frames for rcvDataLoop buffer.
Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).numFrames = 1;  % one intermediate buffer needed.
Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).numFrames = 10;
Resource.DisplayWindow(1).Title = 'C5-2gHFlash_64LE';
Resource.DisplayWindow(1).numFrames = 20;
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
TX(1) = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,radius], ...
                   'focus', 0.0, ...
                   'Steer', [0.0,0.0], ...
                   'Apod', ones(1,Trans.numelements), ...
                   'Delay', zeros(1,Trans.numelements)), 1, 1);
% - Set event specific TX attributes.
TX(1).Delay = computeTXDelays(TX(1));

% Specify Receive structure arrays.
% -- Compute the maximum receive path length, using the law of cosines.
maxAcqLength = ceil(sqrt((P.endDepth+radius)^2 + radius^2 - ...
                     2*(P.endDepth+radius)*radius*cos(scanangle)));
Receive = repmat(struct('Apod', zeros(1,Trans.numelements), ...
                        'mode', 0, ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'callMediaFunc', 0), 1, 2*Resource.RcvBuffer(1).numFrames);

% - Set event specific Receive attributes for each frame.
for i = 1:Resource.RcvBuffer(1).numFrames
    % -- 1st synthetic aperture acquisition for full frame.
    Receive(2*i-1).Apod(1:Resource.Parameters.numRcvChannels) = 1.0;
    Receive(2*i-1).framenum = i;
    Receive(2*i-1).acqNum = 1;
    Receive(2*i-1).callMediaFunc = 1;
    % -- 2nd synthetic aperture acquisition for full frame.
    Receive(2*i).Apod((Resource.Parameters.numRcvChannels+1):Trans.numelements) = 1.0;
    Receive(2*i).framenum = i;
    Receive(2*i).acqNum = 2;   % two acquisitions per frame
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [0,296,433,548,631,737,816,924];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Recon structure arrays.
% - We need one Recon structure which will be used for each frame.
Recon = struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame',-1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...
               'RINums',1:2);

% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 'replaceIQ', ...
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 1),1,2);
ReconInfo(2).mode = 'accumIQ_replaceIntensity';
ReconInfo(2).rcvnum = 2;

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
SeqControl(1).command = 'jump'; % jump back to start
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(2).argument = 290;  % 290 usec
SeqControl(3).command = 'timeToNextAcq';  % time between frames
SeqControl(3).argument = 5000;  % 5 msec
SeqControl(4).command = 'returnToMatlab';
nsc = 5; % nsc is count of SeqControl objects

% Specify Event structure arrays.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames                     % Acquire frame
    % Acquire 1st half of frame
    Event(n).info = '1st half of Aperture.';
    Event(n).tx = 1;
    Event(n).rcv = 2*i-1;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = 2;
    n = n+1;

    % Acquire 2nd half of frame
    Event(n).info = '2nd half of Aperture.';
    Event(n).tx = 1;
    Event(n).rcv = 2*i;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = [3,nsc];
    SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
    nsc = nsc+1;
    n = n+1;

    Event(n).info = 'recon and process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 1;
    Event(n).seqControl = 0;
    if (floor(i/5) == i/5)&&(i~=Resource.RcvBuffer(1).numFrames) % Exit to Matlab every 5th frame
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
import vsv.seq.uicontrol.VsSliderControl

% - Sensitivity Cutoff
UI(1).Control = VsSliderControl('LocationCode','UserB7',...
    'Label','Sens. Cutoff',...
    'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
    'SliderStep',[0.025,0.1],'ValueFormat','%1.3f',...
    'Callback',@SensCutOffCallback);

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
    'SliderMinMaxVal',MinMaxVal,...
    'SliderStep',[0.1,0.2],'ValueFormat','%3.0f',...
    'Callback',@RangeChangeCallback);

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 5;

% Save all the structures to a .mat file.
save('MatFiles/C5-2gHFlash_64LE');


% **** Callback routines used by UI Controls. ****
function SensCutOffCallback(~,~,UIValue)
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
end
