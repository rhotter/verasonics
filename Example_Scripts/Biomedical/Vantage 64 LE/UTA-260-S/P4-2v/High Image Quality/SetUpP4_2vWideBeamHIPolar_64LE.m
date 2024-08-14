% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: SetUpP4_2vWideBeamPolarHI_64LE.m - Example of phased array imaging with
%                                        wide beam transmits and harmonic imaging.
% Description:
%   Sequence programming file for P4-2v phased array in virtual apex format,
%   using m wide beam transmit pairs with pulse inversion. All 64
%   transmit and receive channels are active for each acquisition.
%   The wide beams are reconstucted on a polar coordinate pixel grid containing 128
%   angles spaced equally over the scan angle.
%   Processing is asynchronous with respect to acquisition.
%
% Last update:
% 03/08/2019 updated for 4.0.1 release

clear all

% Define system parameters.
Resource.Parameters.numTransmit = 128;
Resource.Parameters.numRcvChannels = 64;  % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.
Resource.System.SoftwareVersion = [4 0 0]; % Minimum software release for this script.

% Specify Trans structure array.
Trans.name = 'P4-2v';
Trans.units = 'wavelengths'; % required in Gen3 to prevent default to mm units
Trans = computeTrans(Trans);
Trans.maxHighVoltage = 50;  % set maximum high voltage limit for pulser supply.
mm2wl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);

% Define parameters for settings
P.numRays = 64; % no. of raylines to program
P.startDepth = 0;
P.endDepth = 140*mm2wl;  % Acquisition depth in wavelengths
P.txFocus = -360*mm2wl;
P.theta = 75*pi/180;   % 75 degree angle of sector
P.rayDelta = P.theta/(P.numRays-1);
aperture = 64*Trans.spacing; % aperture based on 64 elements
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
            'angle',P.rayDelta*6, ...
            'steer',0, ...
            'andWithPrev',1)),1,P.numRays+1);
% First Region covers entire field of view
PData(1).Region(1).Shape.andWithPrev = 0;
PData(1).Region(1).Shape.angle = P.theta;
% - set position of regions to correspond to beam spacing.
Angle = (-P.theta/2):P.rayDelta:(P.theta/2);
for j = 1:P.numRays
    PData.Region(j+1).Shape.steer = Angle(j);
end
PData.Region = computeRegions(PData);

% Specify Media.  Use point targets in middle of PData.
% Set up Media points
% - Uncomment for speckle
% Media.MP = rand(40000,4);
% Media.MP(:,2) = 0;
% Media.MP(:,4) = 0.04*Media.MP(:,4) + 0.04;  % Random amplitude
% Media.MP(:,1) = 2*halfwidth*(Media.MP(:,1)-0.5);
% Media.MP(:,3) = SFormat.acqDepth*Media.MP(:,3);
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
Resource.RcvBuffer(1).rowsPerFrame = 4224*P.numRays; % This is for the max range on range slider.
Resource.RcvBuffer(1).colsPerFrame = 128; % 128 columns in RcvBuffer but only 64 used
Resource.RcvBuffer(1).numFrames = 10;
Resource.InterBuffer(1).numFrames = 1;  % 1 frame defined but no intermediate buffer needed.
Resource.ImageBuffer(1).numFrames = 10;
Resource.DisplayWindow(1).Title = 'P4-2vWideBeamPolarHI_64';
Resource.DisplayWindow(1).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
pdeltaR = PData(1).PDelta(2);
pdeltaT = PData(1).PDelta(1);
DwWidth = ceil((2*PData.Size(1)*pdeltaR*sin(pdeltaT*(PData.Size(2)-1)/2))/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(P.endDepth/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [-DwWidth/2*Resource.DisplayWindow(1).pdelta,0,0];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).numFrames = 20;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow.Colormap = gray(256);

% Specify TW structure array.
% - Note nominal center frequency in computeTrans is 2.976 MHz
TW(1).type = 'parametric';
TW(1).Parameters = [2.083,0.67,1,1];
TW(2).type = 'parametric';
TW(2).Parameters = [2.083,0.67,1,-1];

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'focus', P.txFocus, ...
                   'Steer', [0.0,0.0], ...
                   'Apod', zeros(1,Trans.numelements), ...  % set TX.Apod for 64 elements
                   'Delay', zeros(1,Trans.numelements), ...
                   'TXPD', [], ...
                   'peakCutOff', 1.5, ...
                   'peakBLMax', 3.5), 1, 2*P.numRays);
% - Set event specific TX attributes.
m = 200; % numelements for hann window
W = hann(m);
TXorgs = dapex*tan(Angle);
for n = 1:P.numRays   % P.numRays transmit events
    TX(n).Origin = [TXorgs(n),0.0,0.0];
    ce = round(64*(TXorgs(n) - Trans.ElementPos(1,1))/aperture);
%     ce = round((Trans.numelements-1) * (TX(n).Origin(1)-TX(1).Origin(1))/(aperture/2)) + 1;
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
    k = n+P.numRays;
    TX(k).waveform = 2;
    TX(k).Origin = TX(n).Origin;
    TX(k).Apod = TX(n).Apod;
    TX(k).Steer = TX(n).Steer;
    TX(k).Delay = TX(n).Delay;
    TX(k).TXPD = TX(n).TXPD;
end

% Specify Receive structure arrays.
% - We need P.numRays Receives for each frame.
% -- Compute the maximum receive path length, using the law of cosines.
maxAcqLength = ceil(sqrt((aperture)^2 + P.endDepth^2 - 2*(aperture)*P.endDepth*cos(P.theta/2+pi/2)) - P.startDepth);
Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
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
    k = 2*P.numRays*(i-1);
    Receive(k+1).callMediaFunc = 1;
    for j = 1:P.numRays
        Receive(k+j).framenum = i;
        Receive(k+j).acqNum = j;
        Receive(k+j+P.numRays).mode = 1;
        Receive(k+j+P.numRays).framenum = i;
        Receive(k+j+P.numRays).acqNum = j;
    end
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [0,314,558,706,830,879,934,968];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Recon structure array.
Recon = struct('senscutoff', 0.4, ...
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
                   'scaleFactor', 1.0, ...
                   'threadSync', 1,...
                   'regionnum', 0), 1, P.numRays);

% - Set specific ReconInfo attributes.
ReconInfo(1).Pre = 'clearInterBuf';     % clear entire Interbuffer frame
for j = 1:P.numRays
    ReconInfo(j).txnum = j;
    ReconInfo(j).rcvnum = j;
    ReconInfo(j).regionnum = j+1; % add 1 to skip first region
end
ReconInfo(P.numRays).Post = 'IQ2IntensityImageBuf';  % detect IQ data

% Specify Process structure array.
pers = 20;
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',5.0,...            % pgain is image processing gain
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
SeqControl(1).command = 'jump';
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq'; % time between rays
SeqControl(2).argument = 280;
SeqControl(3).command = 'timeToNextAcq'; % time between frames
SeqControl(3).argument = 25000 - (P.numRays-1)*280;
SeqControl(4).command = 'returnToMatlab';
nsc = 5;

% Specify Event structure arrays.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    k = 2*P.numRays*(i-1);
    for j = 1:P.numRays                      % Acquire frame
        Event(n).info = 'Acquire non-inverted transmit ray line';
        Event(n).tx = j;
        Event(n).rcv = k+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;

        Event(n).info = 'Acquire inverted transmit ray line';
        Event(n).tx = j+P.numRays;
        Event(n).rcv = k+j+P.numRays;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;
    end
    Event(n-1).seqControl = [3,nsc]; % Replace last event's seqControl for frame time and transferToHost.
       SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
       nsc = nsc+1;

    Event(n).info = 'recon and process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 1;
    Event(n).seqControl = 0;
    % Exit to Matlab every 3rd frame reconstructed
    if (floor(i/3) == i/3)&&(i<Resource.RcvBuffer(1).numFrames)
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
MinMaxVal = [64,300,P.endDepth]; % default unit is wavelength
AxesUnit = 'wls';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
        AxesUnit = 'mm';
        wls2mm = Resource.Parameters.speedOfSound/1000/Trans.frequency;
        MinMaxVal = MinMaxVal * wls2mm;
    end
end
UI(2).Control = {'UserA1','Style','VsSlider','Label',['Range (',AxesUnit,')'],...
                 'SliderMinMaxVal',MinMaxVal,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f'};
UI(2).Callback = text2cell('%RangeChangeCallback');

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 3;

% Save all the structures to a .mat file.
save('MatFiles/P4-2vWideBeamPolarHI_64LE');
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
mm2wl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);
P = evalin('base','P');
P.endDepth = UIValue;
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
        P.endDepth = UIValue*mm2wl;
    end
end
assignin('base','P',P);
aperture = 64*Trans.spacing; % aperture based on 64 elements
dapex = (aperture/4)/tan(P.theta/2); % dist. to virt. apex

% Modify PData for new range
PData = evalin('base','PData');
if P.endDepth <= 192
   PData.PDelta(2) = 0.5;
elseif P.endDepth <= 256
   PData.PDelta(2) = 0.7;
elseif P.endDepth <= 320
   PData.PDelta(2) = 0.9;
end
PData.Size(1) = ceil((P.endDepth - PData.Origin(3))/PData.PDelta(2)) + 10; % rows
PData.Region = repmat(struct(...
            'Shape',struct('Name','SectorFT', ...
            'Position',PData.Origin, ...
            'z',P.startDepth, ...
            'r',dapex+P.endDepth, ...
            'angle',P.rayDelta*6, ...
            'steer',0, ...
            'andWithPrev',1)),1,P.numRays+1);
% First Region covers entire field of view
PData(1).Region(1).Shape.andWithPrev = 0;
PData(1).Region(1).Shape.angle = P.theta;
% - set steering of Regions to correspond to beam spacing.
Angle = (-P.theta/2):P.rayDelta:(P.theta/2);
for j = 1:P.numRays
    PData.Region(j+1).Shape.steer = Angle(j);
end
PData.Region = computeRegions(PData);
assignin('base','PData',PData);
evalin('base','pdeltaR = PData.PDelta(2);');
evalin('base','DwWidth = ceil((2*PData.Size(1)*pdeltaR*sin(pdeltaT*(PData.Size(2)-1)/2))/Resource.DisplayWindow(1).pdelta);');
evalin('base','DwHeight = ceil(P.endDepth/Resource.DisplayWindow(1).pdelta);');
evalin('base','Resource.DisplayWindow(1).Position(3) = DwWidth;');
evalin('base','Resource.DisplayWindow(1).Position(4) = DwHeight;');
evalin('base','Resource.DisplayWindow(1).ReferencePt(1) = -DwWidth/2*Resource.DisplayWindow(1).pdelta;');
% Update TXPD data of TX structures.
TX = evalin('base','TX');
for j = 1:P.numRays
    TX(j).TXPD = computeTXPD(TX(j),PData);
    TX(j+P.numRays).TXPD = TX(j).TXPD;
end
assignin('base','TX',TX);
% Update Receive structures.
Receive = evalin('base', 'Receive');
maxAcqLength = sqrt(aperture^2 + P.endDepth^2 - 2*aperture*P.endDepth*cos(P.theta/2+pi/2)) - P.startDepth;
for i = 1:size(Receive,2)
    Receive(i).endDepth = ceil(maxAcqLength);
end
assignin('base','Receive',Receive);
evalin('base','TGC.rangeMax = P.endDepth;');
evalin('base','TGC.Waveform = computeTGCWaveform(TGC);');
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'PData','InterBuffer','ImageBuffer','TX','Receive','TGC','Recon','DisplayWindow','ImageBuffer'};
assignin('base','Control', Control);
assignin('base', 'action', 'displayChange');
return
%RangeChangeCallback
