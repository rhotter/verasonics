% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use.
%
% File name: SetUp_NDT_LineMode_64_SI.m - Example of imaging using conventioanl line
% by line scan
%
% Description:
%   Sequence programming file for Olympus 5MHz Phase array, using line by line focused wave scan.
%   For each acquisition, all channels are active for transmit and receive.
%   Reconstruction and processing can be set to synchronised with respect
%   to acquisition with Resource.Parameters.waitForProcessing = 1.
%
% Last update:
% 06/23/2020 - Update to SW 4.3 format for new user UIControls and External function definitions (VTS 1691).
%   More info:(?/Example_Scripts/Vantage_Features/New UI Scheme/SetUpL11_5vFlash_NewUI)
% 06/19/2020 - VTS-1641 bug fixes
% 10/10/2018 - updated SetUp file name and .mat name

clear all

P.numRays = 64;      % no. of Rays (1 for Flash transmit)
P.startDepth = 0;
P.endDepth = 80;    % Acquisition depth in wavelengths
P.txFocus = 40;     % initial value of P.txFocus

P.HV = 8;           % preset voltage
P.Elnum=64;         % num of el in transmit aperture, must be no more than Trans.numelements

% Specify system parameters.
Resource.Parameters.numTransmit = 64;       % number of transmit channels.
Resource.Parameters.numRcvChannels = 64;    % number of receive channels.
Resource.Parameters.speedOfSound = 6320;
Resource.Parameters.Connector = 1;          % allows muti-entries. 1 - upper hypertronics conn 2- lower hypertronics conn 3- lemos
Resource.Parameters.waitForProcessing = 1;  % DMA transfers wait for processing.
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.
Resource.System.UTA = '160-SI';

% Specify Trans structure array.
Trans.name = 'Olympus-5MHz-64El';
Trans.units = 'mm';                                         % required in Gen3 to prevent default to mm units
Trans.numelements = 64;
Trans.frequency = 5;                                % nominal frequency in MHz
Trans.Bandwidth = [5*(1-0.67/2), 5*(1+0.67/2)];
Trans.type = 0;                                     % Array geometry is linear (x values only).
Trans.id = -1;
Trans.connType = 10;                                 % Ipex connector
Trans.elevationApertureMm = 10;                      % active elevation aperture in mm
Trans.elevationFocusMm = 1000000;                   % nominal elevation focus depth from lens on face of transducer
Trans.elementWidth = 0.52;                          % element width in mm; assumes 50 micron kerf
Trans.spacingMm = 0.6;                             % Spacing between elements in mm.
Trans.ElementPos = zeros(Trans.numelements,4);
Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
if ~isfield(Trans,'ElementSens')                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
    Theta = (-pi/2:pi/100:pi/2);
    Theta(51) = 0.0000001;                          % set to almost zero to avoid divide by zero.
    eleWidthWl = Trans.elementWidth * Trans.frequency/Resource.Parameters.speedOfSound;
    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
end
Trans.lensCorrection = 0;
Trans.impedance = 50;
Trans.maxHighVoltage = 50;

% Now convert all units as required, based on Trans.units
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000); % conversion factor from mm to wavelengths
% regardless of units, always provide spacing in wavelengths
Trans.spacing = Trans.spacingMm * scaleToWvl;   % Spacing between elements in wavelengths

P.theta = -pi/4;
P.rayDelta = 2*(-P.theta)/(P.numRays-1);
P.aperture = P.Elnum*Trans.spacing;
P.radius = (P.aperture/2)/tan(-P.theta); % dist. to virt. apex
% Specify PData structure array.
PData(1).PDelta = [0.5, 0, 0.5];  % x, y, z pdeltas
PData(1).Size(1) = 10 + ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));
PData(1).Size(2) = 10 + ceil(2*(P.endDepth + P.radius)*sin(-P.theta)/PData(1).PDelta(1));
PData(1).Size(3) = 1;      % single image page
PData(1).Origin = [-(PData(1).Size(2)/2)*PData(1).PDelta(1),0,P.startDepth];
% - specify 128 Region structures.
PData(1).Region = repmat(struct('Shape',struct( ...
                    'Name','SectorFT',...
                    'Position',[0,0,-P.radius],...
                    'z',P.startDepth,...
                    'r',P.radius+P.endDepth,...
                    'angle',P.rayDelta,...
                    'steer',0)),1,P.numRays);
% - set position of regions to correspond to beam spacing.
for i = 1:P.numRays
    PData(1).Region(i).Shape.steer(1) = P.theta + (i-1)*P.rayDelta;
end
PData(1).Region = computeRegions(PData(1));

% Specify Media.  Use point targets in middle of PData.
% Set up Media points
% - Uncomment for speckle
% Media.MP = rand(40000,4);
% Media.MP(:,2) = 0;
% Media.MP(:,4) = 0.04*Media.MP(:,4) + 0.04;  % Random amplitude
% Media.MP(:,1) = 2*halfwidth*(Media.MP(:,1)-0.5);
% Media.MP(:,3) = P.endDepth*Media.MP(:,3);
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
% Media.attenuation = -0.5;
Media.function = 'movePoints';

% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 2048*P.numRays; % This is for the max range on range slider.
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 10;
Resource.InterBuffer(1).numFrames = 1;  % 1 frame defined but no intermediate buffer needed.
Resource.ImageBuffer(1).numFrames = 10;
Resource.DisplayWindow(1).Title = ['LineMode' num2str(P.numRays) 'Acqs'];
Resource.DisplayWindow(1).pdelta = PData(1).PDelta(1)/2;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight =ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).numFrames = 20;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow.Colormap = gray(256);

% Specify TW structure array.
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,.67,2,1];

ElOffset = floor((Trans.numelements - P.Elnum)/2); if ElOffset<0, error('reset P.Elnum to a value larger than Trans.numelements.'); end

% Specify TPC structure.
TPC(1).hv = P.HV;

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'focus', P.txFocus, ...
                   'Steer', [0.0,0.0], ...
                   'Apod', [zeros(1,ElOffset), ones(1,P.Elnum), zeros(1,Trans.numelements-P.Elnum-ElOffset)], ...
                   'Delay', zeros(1,min([Trans.numelements, Resource.Parameters.numTransmit]))), 1, P.numRays);
% - Set event specific TX attributes.
Angles = P.theta:P.rayDelta:(P.theta + (P.numRays-1)*P.rayDelta);
TXorgs = P.radius*tan(Angles);
for n = 1:P.numRays   % P.numRays transmit events
    TX(n).Origin = [TXorgs(n),0.0,0.0];
    TX(n).Steer = [Angles(n),0.0];
    TX(n).Delay = computeTXDelays(TX(n));
end
%clear Angles TXorgs

% Specify Receive structure arrays.
% - We need P.numRays Receives for each frame.
% -- Compute the maximum receive path length, using the law of cosines.
maxAcqLength = ceil(sqrt(P.aperture^2 + P.endDepth^2 - 2*P.aperture*P.endDepth*cos(P.theta-pi/2)) - P.startDepth);
wlsPer128 = 128/(2*4); % wavelengths in 128 samples for 4 samplesPerWave
Receive = repmat(struct('Apod', [zeros(1,ElOffset), ones(1,P.Elnum), zeros(1,Trans.numelements-P.Elnum-ElOffset)], ...
                        'startDepth', P.startDepth, ...
                        'endDepth', P.startDepth + wlsPer128*ceil(maxAcqLength/wlsPer128), ...
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
TGC.CntrlPts = [0,581,696,709,769,830,890,950];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Recon structure array.
Recon = struct('senscutoff', 0.5, ...
               'pdatanum', 1, ...
               'rcvBufFrame',-1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...
               'RINums', 1:P.numRays);

% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 'replaceIntensity', ...  % replace data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 0), 1, P.numRays);
% - Set specific ReconInfo attributes.
for i = 1:P.numRays
    ReconInfo(i).txnum = i;
    ReconInfo(i).rcvnum = i;
    ReconInfo(i).regionnum = i;
end

% Specify Process structure array.
pers = 1;
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
SeqControl(1).command = 'jump'; %  - Jump back to start.
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % set time between rays
SeqControl(2).argument = 500;
SeqControl(3).command = 'timeToNextAcq';  % set time between frames
SeqControl(3).argument = 50000 - (P.numRays-1)*SeqControl(2).argument;
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
    Event(n-1).seqControl = [3,nsc]; % Replace last event's seqControl for frame time and transferToHost.
      SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
      nsc = nsc+1;

    Event(n).info = 'recon and process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 1;
    Event(n).seqControl = 4;
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
    'SliderMinMaxVal',[64,150,P.endDepth]*wls2mm,...
    'SliderStep',[0.1,0.2],'ValueFormat','%3.0f',...
    'Callback',@RangeChangeCallback);

% - Transmit focus change
UI(3).Control = VsSliderControl('LocationCode','UserB4',...
    'Label',['TX Focus (',AxesUnit,')'],...
    'SliderMinMaxVal',[20,150,P.txFocus]*wls2mm,...
    'SliderStep',[0.1,0.2],'ValueFormat','%3.0f',...
    'Callback',@TxFocusCallback);

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 1;

% Save all the structures to a .mat file.
save('MatFiles/NDT_LineMode_64_SI');

%% Automatic VSX Execution:
% Uncomment the following line to automatically run VSX every time you run
% this SetUp script (note that if VSX finds the variable 'filename' in the
% Matlab workspace, it will load and run that file instead of prompting the
% user for the file to be used):

% filename = 'NDT_LineMode_64_SI'; VSX;


% **** Callback routines used by UIControls (UI) ****

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

    PData = evalin('base','PData');
    PData(1).Size(1) = 10 + ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));
    % - specify 128 Region structures.
    PData(1).Region = repmat(struct('Shape',struct( ...
                        'Name','SectorFT',...
                        'Position',[0,0,-P.radius],...
                        'z',P.startDepth,...
                        'r',P.radius+P.endDepth,...
                        'angle',P.rayDelta,...
                        'steer',0)),1,P.numRays);
    % - set position of regions to correspond to beam spacing.
    for i = 1:P.numRays
        PData(1).Region(i).Shape.steer(1) = P.theta + (i-1)*PData(1).Region(i).Shape.angle;
    end
    PData(1).Region = computeRegions(PData(1));
    assignin('base','PData',PData);
    evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);');
    Receive = evalin('base', 'Receive');
    maxAcqLength = ceil(sqrt(P.aperture^2 + P.endDepth^2 - 2*P.aperture*P.endDepth*cos(P.theta-pi/2)) - P.startDepth);
    wlsPer128 = evalin('base','wlsPer128');
    for i = 1:size(Receive,2)
        Receive(i).endDepth = P.startDepth + wlsPer128*ceil(maxAcqLength/wlsPer128);
    end
    assignin('base','Receive',Receive);
    evalin('base','TGC.rangeMax = P.endDepth;');
    evalin('base','TGC.Waveform = computeTGCWaveform(TGC);');
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'PData','InterBuffer','Receive','TGC','Recon','DisplayWindow','ImageBuffer'};
    assignin('base','Control', Control);
    assignin('base', 'action', 'displayChange');
end

function TxFocusCallback(hObject, ~, UIValue)
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
        if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
            P.txFocus = UIValue*scaleToWvl;
        end
    end
    assignin('base','P',P);

    % - Redefine event specific TX attributes for the new focus.
    TX = evalin('base', 'TX');
    for n = 1:size(TX,2)   % 128 transmit events
        % write new focus value to TX
        TX(n).focus = P.txFocus;
        TX(n).Delay = computeTXDelays(TX(n));
    end
    assignin('base','TX', TX);
    % Set Control command to update TX
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'TX'};
    assignin('base','Control', Control);
end
