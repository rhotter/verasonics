% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use.
%
% File name: SetUpP4_2gHWideBeamPolar.m - Example of phased array imaging with
%                            wide beam transmits using polar coordinate PData
% Description:
%   Sequence programming file for P4-2gH phased array in virtual apex format,
%   using m wide beam tranmits and receive acquisitions. All 96 transmit and
%   receive channels are active for each acquisition. Processing is
%   asynchronous with respect to acquisition. Note: The P4-2gH is a 96 element
%   probe that is wired to the scanhead connector with the elements 1-96
%   connected to I/O channels 1-96.  We need a Trans.Connector array
%   to specify the connector channels used, which will be defined
%   in the Trans structure.
%
% Last update:
%   10-30-2020 Adapted from SetUpP4_2vWideBeamPolar.m script.
%   09-23-2021 modified for P4-2gH

%clear[ 	]+all

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
Trans.name = 'P4-2gH';
Trans.units = 'wavelengths'; % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);

mm2wl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);

% Define parameters for settings
P.numRays = 64; % no. of raylines to program
P.startDepth = 0;
P.endDepth = 120*mm2wl;  % Acquisition depth in wavelengths
P.txFocus = 400*mm2wl;
P.theta = 75*pi/180;   % 75 degree angle of sector
P.rayDelta = P.theta/(P.numRays-1);
aperture = 96*Trans.spacing; % aperture based on 96 elements
dapex = (aperture/4)/tan(P.theta/2); % dist. to virt. apex

% Set up PData structure.
PData.Coord = 'polar';
PData.PDelta = [P.theta/127,0.5,0]; % [dtheta,dr,dz]  define 128 ray lines for recon
PData.Origin = [0,0,-dapex];
PData.Size(1) = ceil((P.endDepth - PData.Origin(3))/PData.PDelta(2)) + 10; % rows
PData.Size(2) = ceil(P.theta/PData.PDelta(1) + 1); % cols
PData.Size(3) = 1; % pages
PData.Region = repmat(struct(...
            'Shape',struct('Name','SectorFT', ...
            'Position',PData.Origin, ...
            'z',P.startDepth, ...
            'r',dapex+P.endDepth, ...
            'angle',P.rayDelta*10, ...
            'steer',0)),1,P.numRays);
% - set position of regions to correspond to beam spacing.
Angle = (-P.theta/2):P.rayDelta:(P.theta/2);
for j = 1:P.numRays
    PData.Region(j).Shape.steer = Angle(j);
end
PData.Region = computeRegions(PData);

% Specify Media.  Use point targets in middle of PData.
% Set up Media points
% - Uncomment for speckle
% Media.MP = rand(40000,4);
% Media.MP(:,2) = 0;
% Media.MP(:,4) = 0.04*Media.MP(:,4) + 0.04;  % Random amplitude
% Media.MP(:,1) = 2*halfwidth*(Media.MP(:,1)-0.5);
% Media.MP(:,3) = P.acqDepth*Media.MP(:,3);
Media.MP(1,:) = [-45,0,30,1.0];
Media.MP(2,:) = [-15,0,30,1.0];
Media.MP(3,:) = [15,0,30,1.0];
Media.MP(4,:) = [45,0,30,1.0];
Media.MP(5,:) = [-15,0,60,1.0];
Media.MP(6,:) = [-15,0,90,1.0];
Media.MP(7,:) = [-15,0,120,1.0];
Media.MP(8,:) = [-15,0,150,1.0];
Media.MP(9,:) = [-45,0,120,1.0];
Media.MP(10,:) = [15,0,120,1.0];
Media.MP(11,:) = [45,0,120,1.0];
Media.MP(12,:) = [-10,0,69,1.0];
Media.MP(13,:) = [-5,0,75,1.0];
Media.MP(14,:) = [0,0,78,1.0];
Media.MP(15,:) = [5,0,80,1.0];
Media.MP(16,:) = [10,0,81,1.0];
Media.MP(17,:) = [-75,0,120,1.0];
Media.MP(18,:) = [75,0,120,1.0];
Media.MP(19,:) = [-15,0,180,1.0];
Media.numPoints = 19;
Media.attenuation = -0.5;
Media.function = 'movePoints';

% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = P.numRays*4096;
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = numRcvFrames;     % 10 frames used for RF cineloop.
Resource.InterBuffer(1).numFrames = 1;
Resource.ImageBuffer(1).numFrames = 10;
Resource.DisplayWindow(1).Title = 'P4-2gHWideBeamPolar';
Resource.DisplayWindow(1).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
pdeltaR = PData.PDelta(2);
pdeltaT = PData.PDelta(1);
DwWidth = ceil((2*PData.Size(1)*pdeltaR*sin(pdeltaT*(PData.Size(2)-1)/2))/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil((PData.Size(1)*pdeltaR)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [-DwWidth/2*Resource.DisplayWindow(1).pdelta,0,0];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).numFrames = numImageFrames;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow.Colormap = gray(256);

% Specify TW structure array.
TW.type = 'parametric';
TW.Parameters = [Trans.frequency,0.67,1,1];

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'focus', P.txFocus, ...
                   'Steer', [0.0,0.0], ...
                   'Apod', zeros(1,96), ...  % set TX.Apod for 96 elements
                   'Delay', zeros(1,96), ...
                   'TXPD', [], ...
                   'peakCutOff', 1.5, ...
                   'peakBLMax', 3.5), 1, P.numRays);
% - Set event specific TX attributes.
m = 200; % numelements for hann window
W = hann(m);
TXorgs = dapex*tan(Angle);
for n = 1:P.numRays   % P.numRays transmit events
    TX(n).Origin = [TXorgs(n),0.0,0.0];
    ce = round(Trans.numelements*(TXorgs(n) - Trans.ElementPos(1,1))/aperture);
    TX(n).Apod(:) = W((m/2-(ce-1)):(m/2+(96-ce)));
    TX(n).Steer = [Angle(n),0.0];
    TX(n).Delay = computeTXDelays(TX(n));
    TX(n).TXPD = computeTXPD(TX(n),PData);
end

% Specify Receive structure arrays.
maxAcqLength = ceil(P.endDepth + aperture); % sufficient for 90 degree scan
Receive = repmat(struct('Apod', ones(1,96), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 0),1,P.numRays*Resource.RcvBuffer(1).numFrames);
% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    Receive(P.numRays*(i-1)+1).callMediaFunc = 1;
    for j = 1:P.numRays
        Receive(P.numRays*(i-1)+j).framenum = i;
        Receive(P.numRays*(i-1)+j).acqNum = j;
    end
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [0,292,470,643,733,858,940,972];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Recon structure arrays.
Recon = struct('senscutoff', 0.5, ...
               'pdatanum', 1, ...
               'rcvBufFrame', -1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...
               'RINums', 1:P.numRays);

% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 'accumIQ', ...  % accumulate IQ
                   'Pre',[],...
                   'Post',[],...
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'scaleFactor', 0.5, ...
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
%  - Jump back to start.
SeqControl(1).command = 'jump';
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq'; % time between rays
SeqControl(2).argument = 260;
SeqControl(3).command = 'timeToNextAcq'; % time between frames
SeqControl(3).argument = 25000 - (P.numRays-1)*260;
SeqControl(4).command = 'returnToMatlab';
nsc = 5;

% Specify Event structure arrays.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    for j = 1:P.numRays                 % Acquire frame
        Event(n).info = 'Acquire ray line.';
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
frameRateFactor = 3;

% Save all the structures to a .mat file.
filename = fullfile( vsv.file.getVSXDir, 'MatFiles/SetUpP4_2gHWideBeamPolar_QSApp.mat');

save(filename);

% **** Callback routines to be used by UI Controls ****
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

function RangeChangeCallback(~,~,UIValue)
    AxesUnits = 'wl';
    mm2wl = evalin('base','mm2wl');
    P = evalin('base','P');
    Resource = evalin('base','Resource');
    if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
        if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
            AxesUnits = 'mm';
        end
    end
    simMode = evalin('base','Resource.Parameters.simulateMode');
    % No range change if in simulate mode 2.
    if simMode == 2
        if strcmp(AxesUnits,'mm')
            set(hObject,'Value',P.endDepth/mm2wl); % convert to mm
        else
            set(hObject,'Value',P.endDepth);
        end
        return
    end
    Trans = evalin('base','Trans');
    P.endDepth = UIValue;
    if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
        if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
            P.endDepth = UIValue*mm2wl;
        end
    end
    assignin('base','P',P);
    aperture = 96*Trans.spacing; % aperture based on 96 elements
    dapex = (aperture/4)/tan(P.theta/2); % dist. to virt. apex
    PData = evalin('base','PData');
    PData.Size(1) = ceil((P.endDepth - PData.Origin(3))/PData.PDelta(2)) + 10; % rows
    PData.Size(2) = ceil(P.theta/PData.PDelta(1) + 1); % cols
    PData.Size(3) = 1; % pages
    PData.Region = repmat(struct(...
                'Shape',struct('Name','SectorFT', ...
                'Position',PData.Origin, ...
                'z',P.startDepth, ...
                'r',dapex+P.endDepth, ...
                'angle',P.rayDelta*10, ...
                'steer',0)),1,P.numRays);
    % - set position of regions to correspond to beam spacing.
    Angle = (-P.theta/2):P.rayDelta:(P.theta/2);
    for j = 1:P.numRays
        PData.Region(j).Shape.steer = Angle(j);
    end
    PData.Region = computeRegions(PData);
    assignin('base','PData',PData);
    evalin('base','Resource.DisplayWindow(1).Position(4) = ceil((PData.Size(1)*PData.PDelta(2))/Resource.DisplayWindow(1).pdelta);');
    % Update TXPD data of TX structures.
    TX = evalin('base','TX');
    for i = 1:size(TX,2)
        TX(i).TXPD = computeTXPD(TX(i),PData);
    end
    assignin('base','TX',TX);
    % Update Receive structures.
    Receive = evalin('base', 'Receive');
    maxAcqLength = ceil(P.endDepth + aperture); % sufficient for 90 degree scan
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
