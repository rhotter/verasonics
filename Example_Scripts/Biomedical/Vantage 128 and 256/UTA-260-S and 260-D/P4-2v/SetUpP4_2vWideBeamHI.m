% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use.
%
% File name: SetUpP4_2vWideBeamHI.m - Example of phased array pulse inversion
%                                     imaging with wide beam transmits
% Description:
%   Sequence programming file for P4-2v phased array in virtual apex format,
%   using m wide beam pulse inversion transmits and receive acquisitions.
%   All 64 transmit and receive channels are active for each acquisition.
%   Processing is asynchronous with respect to acquisition. Note: The P4-2v
%   is a 64 element probe that is wired to the scanhead connector with
%   elements 1-64 connected to inputs 33-96. We therefore need a Trans.Connector
%   array to specify the connector channels used, which will be defined by
%   the computeTrans function.
%
% Last update:
% 05/18/2016

clear all

P.numRays = 64; % no. of raylines to program
P.startDepth = 0;
P.endDepth = 320;   % Acquisition depth in wavelengths
P.txFocus = -700;

% Specify system parameters.
Resource.Parameters.numTransmit = 128;  % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;  % number of receive channels.
Resource.Parameters.speedOfSound = 1540;
Resource.Parameters.speedCorrectionFactor = 1.0;
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Specify Trans structure array.
Trans.name = 'P4-2v';
Trans.units = 'wavelengths'; % Explicit declaration avoids warning message when selected by default
% note nominal center frequency in computeTrans is 2.976 MHz
Trans = computeTrans(Trans);
Trans.maxHighVoltage = 50;  % set maximum high voltage limit for pulser supply.

P.theta = (pi/180)*75;  % 75 degree sector
P.aperture = 64*Trans.spacing; % aperture based on 64 elements
P.dApex = (P.aperture/4)/tan(P.theta/2); % dist. to virt. apex
P.rayDelta = P.theta/(P.numRays-1); % angle increment between rays

% Set up PData structure. The number of Regions matches the number of rays, but the Regions
% are larger in angle.  The apex of the SectorFT region is positioned above the virtual apex of the
% scan, to obtain a wider sector at the transducer face.
PData(1).PDelta = [0.75, 0, 0.5];
PData(1).Size(1) = 10 + ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));
PData(1).Size(2) = 10 + ceil(2*(P.endDepth + P.dApex)*sin(P.theta/2)/PData(1).PDelta(1));
PData(1).Size(3) = 1;
PData(1).Origin = [-(PData(1).Size(2)/2)*PData(1).PDelta(1),0,P.startDepth];
PData(1).Region = repmat(struct(...
            'Shape',struct('Name','SectorFT', ...
            'Position',[0,0,-P.dApex], ...
            'z',P.startDepth, ...
            'r',P.dApex+P.endDepth, ...
            'angle',P.rayDelta*5.5, ...
            'steer',0, ...
            'andWithPrev',1)),1,P.numRays+1);
% First Region covers entire field of view
PData(1).Region(1).Shape.andWithPrev = 0;
PData(1).Region(1).Shape.angle = P.theta;
% Set ray line Regions Position and steering angles
Angle = (-P.theta/2):P.rayDelta:(P.theta/2);
dz = P.dApex/2;  % distance from virtual apex to SectorFT apex.
for n = 2:P.numRays+1
    PData(1).Region(n).Position = [-dz*tan(Angle(n-1)),0,(P.dApex+dz)];
    PData(1).Region(n).Shape.steer = Angle(n-1);
end
PData(1).Region = computeRegions(PData(1));

% Specify Media.  Use point targets in middle of PData.
% Set up Media points
% - Uncomment for speckle
% Media.MP = rand(10000,4);
% Media.MP(:,2) = 0; % no y values
% Media.MP(:,4) = 0.04*Media.MP(:,4) + 0.04;  % Random amplitude
% Media.MP(:,3) = P.endDepth*Media.MP(:,3);   % random z values
% Media.MP(:,1) = (P.dApex+Media.MP(:,3)).*tan(P.theta*(Media.MP(:,1)-0.5));
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
Resource.RcvBuffer(1).rowsPerFrame = P.numRays*4096*2;
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 20;     % 20 frames used for RF cineloop.
Resource.InterBuffer(1).numFrames = 1;
Resource.ImageBuffer(1).numFrames = 10;
Resource.DisplayWindow(1).Title = 'P4-2vWideBeamHI';
Resource.DisplayWindow(1).pdelta = 0.45;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).numFrames = 50;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow.Colormap = gray(256);

% Specify TW structure array.
% Specify Transmit waveform structure.  These structures are persistent and we
%   only need to specify what changes in subsequent structures.
TW(1).type = 'parametric';
TW(1).Parameters = [2.25,0.67,1,1];
TW(2).type = 'parametric';
TW(2).Parameters = [2.25,0.67,1,-1];

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'focus', P.txFocus, ...
                   'Steer', [0.0,0.0], ...
                   'Apod', ones(1,64), ...  % set TX.Apod for 64 elements
                   'Delay', zeros(1,64), ...
                   'TXPD', [], ...
                   'peakCutOff', 1.8, ...
                   'peakBLMax', 10.0), 1, 2*P.numRays);

% - Set event specific TX attributes.
m = 100;
W = hann(m);
TXorgs = P.dApex*tan(Angle);
h = waitbar(0,'Program TX parameters, please wait!');
for n = 1:P.numRays
    TX(n).Origin = [TXorgs(n), 0.0, 0.0];
    ce = round(Trans.numelements*(TXorgs(n) - Trans.ElementPos(1,1))/P.aperture);
    TX(n).Apod(ce) = W(m/2);
    for j=1:64
        lft=ce-j;
        if lft<1, lft=1; elseif lft<(ce-m/2), lft = ce-m/2; end
        rt=ce+j;
        if rt>64, rt=64; elseif rt>(ce+m/2), rt = ce+m/2; end
    end
    TX(n).Apod(lft:rt) = W((m/2-(ce-lft)):(m/2+(rt-ce)));
    TX(n).Steer = [Angle(n),0.0];
    TX(n).Delay = computeTXDelays(TX(n));
    TX(n).TXPD = computeTXPD(TX(n),PData);
    TX(n+P.numRays).waveform = 2;
    TX(n+P.numRays).Origin = TX(n).Origin;
    TX(n+P.numRays).Apod = TX(n).Apod;
    TX(n+P.numRays).Steer = TX(n).Steer;
    TX(n+P.numRays).Delay = TX(n).Delay;
    TX(n+P.numRays).TXPD = TX(n).TXPD;
    waitbar(n/P.numRays)
end
close(h)

% Specify Receive structure arrays.
maxAcqLength = sqrt(P.aperture^2 + P.endDepth^2 - 2*P.aperture*P.endDepth*cos(P.theta/2+pi/2)) - P.startDepth;
Receive = repmat(struct('Apod', ones(1,64), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'InputFilter',[-0.00110 +0.00000 +0.00238 +0.00000 +0.00174 +0.00000 -0.01175 ...
                                       +0.00000 +0.01129 +0.00000 +0.01758 +0.00000 -0.05240 +0.00000 ...
                                       +0.02454 +0.00000 +0.10825 +0.00000 -0.28284 +0.00000 +0.36462], ...
                        'sampleMode', 'NS200BW', ...
                        'demodFrequency',4.1667, ...
                        'mode', 0, ...
                        'callMediaFunc', 0),1,2*P.numRays*Resource.RcvBuffer(1).numFrames);
% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
%     Receive(2*P.numRays*(i-1)+1).callMediaFunc = 1;
    for j = 1:P.numRays
        Receive(2*P.numRays*(i-1)+j).framenum = i;
        Receive(2*P.numRays*(i-1)+j).acqNum = j;
        Receive(2*P.numRays*(i-1)+j+P.numRays).mode = 1;
        Receive(2*P.numRays*(i-1)+j+P.numRays).framenum = i;
        Receive(2*P.numRays*(i-1)+j+P.numRays).acqNum = j;
    end
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [200,400,590,709,769,830,890,950];
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
                   'threadSync', 1,...
                   'regionnum', 0), 1, P.numRays);

% - Set specific ReconInfo attributes.
ReconInfo(1).Pre = 'clearInterBuf';     % clear entire Interbuffer frame
for j = 1:P.numRays
    ReconInfo(j).txnum = j;
    ReconInfo(j).rcvnum = j;
    ReconInfo(j).regionnum = j+1;
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
                         'display',1,...      % display image after processing
                         'displayWindow',1};

% Specify SeqControl structure arrays.
%  - Jump back to start.
SeqControl(1).command = 'jump';
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq'; % time between rays
SeqControl(2).argument = 250;
SeqControl(3).command = 'timeToNextAcq'; % time between frames
SeqControl(3).argument = 40000 - (2*P.numRays-1)*250;
SeqControl(4).command = 'returnToMatlab';
nsc = 5;

% Specify Event structure arrays.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    for j = 1:P.numRays                 % Acquire frame
        Event(n).info = 'Acquire full aperture.';
        Event(n).tx = j;
        Event(n).rcv = 2*P.numRays*(i-1)+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;

        Event(n).info = 'Acquire full aperture, inverted transmit';
        Event(n).tx = j+P.numRays;
        Event(n).rcv = 2*P.numRays*(i-1)+j+P.numRays;
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
MinMaxVal = [64,360,P.endDepth]; % default unit is wavelength
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

% - Wide Beam Region width
UI(5).Control = {'UserC2','Style','VsSlider','Label','Reg. Width',...
                  'SliderMinMaxVal',[1,9,PData.Region(2).Shape.angle/P.rayDelta],...
                  'SliderStep',[0.0625,0.125],'ValueFormat','%1.3f'};
UI(5).Callback = text2cell('%RegionWidthCallback');


% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 3;

% Save all the structures to a .mat file.
save('MatFiles/P4-2vWideBeamHI');

%% Automatic VSX Execution:
% Uncomment the following line to automatically run VSX every time you run
% this SetUp script (note that if VSX finds the variable 'filename' in the
% Matlab workspace, it will load and run that file instead of prompting the
% user for the file to be used):

% filename = 'P4-2vWideBeamHI';  VSX;

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
mm2wl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);

P = evalin('base','P');
P.endDepth = UIValue;
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
        P.endDepth = UIValue*mm2wl;
    end
end
assignin('base','P',P);

PData = evalin('base','PData');
if P.endDepth <= 192
   PData(1).PDelta(3) = 0.5;
elseif P.endDepth <= 256
   PData(1).PDelta(3) = 0.75;
elseif P.endDepth <= 320
   PData(1).PDelta(3) = 1.0;
end
PData(1).Size(1) = 10 + ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));
PData(1).Size(2) = 10 + ceil(2*(P.endDepth + P.dApex)*sin(P.theta/2)/PData(1).PDelta(1));
PData(1).Size(3) = 1;
PData(1).Origin = [-(PData(1).Size(2)/2)*PData(1).PDelta(1),0,P.startDepth];
PData(1).Region = repmat(struct(...
            'Shape',struct('Name','SectorFT', ...
            'Position',[0,0,-P.dApex], ...
            'z',P.startDepth, ...
            'r',P.dApex+P.endDepth, ...
            'angle',P.rayDelta*5, ...
            'steer',0, ...
            'andWithPrev',1)),1,P.numRays+1);
% First Region covers entire field of view
PData(1).Region(1).Shape.andWithPrev = 0;
PData(1).Region(1).Shape.angle = P.theta;
% - set steering of Regions to correspond to beam spacing.
Angle = (-P.theta/2):P.rayDelta:(P.theta/2);
for j = 2:P.numRays+1
    PData(1).Region(j).Shape.steer = Angle(j-1);
end
PData(1).Region = computeRegions(PData(1));
assignin('base','PData',PData);
evalin('base','Resource.DisplayWindow(1).Position(3) = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);');
evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);');
evalin('base','Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];');
% Update TXPD data of TX structures.
TX = evalin('base','TX');
h = waitbar(0,'Program TX parameters, please wait!');
for j = 1:P.numRays
    TX(j).TXPD = computeTXPD(TX(j),PData(1));
    TX(j+P.numRays).TXPD = TX(j).TXPD;
    waitbar(j/P.numRays);
end
close(h)
assignin('base','TX',TX);
% Update Receive structures.
Receive = evalin('base', 'Receive');
maxAcqLength = sqrt(P.aperture^2 + P.endDepth^2 - 2*P.aperture*P.endDepth*cos(P.theta/2+pi/2)) - P.startDepth;
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
Control.Parameters = {'TX','Recon'};
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
Control.Parameters = {'TX','Recon'};
assignin('base','Control', Control);
%MaxBLCallback

%RegionWidthCallback
P = evalin('base','P');
PData = evalin('base','PData');
PData.Region = repmat(struct(...
            'Shape',struct('Name','SectorFT', ...
            'Position',[0,0,-P.dApex], ...
            'z',P.startDepth, ...
            'r',P.dApex+P.endDepth, ...
            'angle',P.rayDelta*UIValue, ...
            'steer',0, ...
            'andWithPrev',1)),1,P.numRays+1);
% First Region covers entire field of view
PData.Region(1).Shape.andWithPrev = 0;
PData.Region(1).Shape.angle = P.theta;
% - set steering of Regions to correspond to beam spacing.
Angle = (-P.theta/2):P.rayDelta:(P.theta/2);
for j = 2:P.numRays+1
    PData.Region(j).Shape.steer = Angle(j-1);
end
PData.Region = computeRegions(PData);
assignin('base','PData',PData);
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'PData','InterBuffer','ImageBuffer','Recon'};
assignin('base','Control', Control);
%RegionWidthCallback
