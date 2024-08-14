% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use.
%
% File name: SetUpP5_64vN_FMC.m - Example of imaging using full matrix capture
%                         and total focusing reconstruction
%
% Description:
%   Sequence programming file for P5-64vN array, using total
%   focusing reconstruction. For each acquisition, one channel is active
%   for transmit and 64 channels are active for receive.
%
% Last update:
% 11/02/2020

clear all

P.startDepth = 0;   % Acquisition depth in wavelengths
P.endDepth = 60;    % This should preferrably be a multiple of 128 samples.
P.HV = 25;          % preset voltage
P.Elnum = 64;       % num of el in transmit aperture, must be no more than Trans.numelements

% Specify system parameters.
Resource.Parameters.numTransmit = 64;     % number of transmit channels.
Resource.Parameters.numRcvChannels = 64;  % number of receive channels.
Resource.Parameters.speedOfSound = 6320;
Resource.Parameters.Connector = 1;
Resource.Parameters.waitForProcessing = 1; % set to 1 for synchronous acquisition and processing
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.
Resource.System.UTA = '160-SI';

% Specify Trans structure array.
Trans.name = 'P5-64vN';
Trans.units = 'mm'; % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);

scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000); % conversion factor from mm to wavelengths

P.theta = -pi/4;
P.rayDelta = 2*(-P.theta);
P.aperture = P.Elnum*Trans.spacing; % P.aperture in wavelengths
P.radius = (P.aperture/2)/tan(-P.theta); % dist. to virt. apex

% Set up PData structure.
PData(1).PDelta = [0.5, 0, 0.5];
PData(1).Size(1) = 10 + ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));
PData(1).Size(2) = 10 + ceil(2*(P.endDepth + P.radius)*sin(-P.theta)/PData(1).PDelta(1));
PData(1).Size(3) = 1;
PData(1).Origin = [-(PData(1).Size(2)/2)*PData(1).PDelta(1), 0, P.startDepth];
PData(1).Region = struct(...
            'Shape',struct('Name','SectorFT', ...
            'Position',[0,0,-P.radius], ...
            'z',P.startDepth, ...
            'r',P.radius+P.endDepth, ...
            'angle',P.rayDelta, ...
            'steer',0));
PData(1).Region = computeRegions(PData(1));

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
% Media.attenuation = -0.5;
Media.function = 'movePoints';

% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 1024*P.Elnum;
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 10;     % 10 frames for RF data.
Resource.InterBuffer(1).numFrames = 1;  % one intermediate buffer needed.
Resource.ImageBuffer(1).numFrames = 10;
Resource.DisplayWindow(1).Title = ['FMC' num2str(P.Elnum) 'Acqs'];
Resource.DisplayWindow(1).pdelta = PData(1).PDelta(1)/2;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).numFrames = 20;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = gray(256);

% Specify Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,.67,2,1];

ElOffset = floor((Trans.numelements - P.Elnum)/2); if ElOffset<0, error('reset P.Elnum to a value smaller than Trans.numelements.'); end

% Specify P.Elnum TX structure arrays. Transmit on element n in the array for event n.
TX = repmat(struct('waveform', 1, ...
                   'Origin', zeros(1,3), ...
                   'focus', 0, ...
                   'Steer', [0.0,0.0], ...
                   'Apod', zeros(1,Trans.numelements), ...
                   'Delay', zeros(1,Trans.numelements)), 1, P.Elnum);

% - Set event specific TX attributes.
for n = 1:P.Elnum   % P.Elnum transmit events
    % Set transmit Origins to positions of elements.
    TX(n).Origin = scaleToWvl*Trans.ElementPos(n,1:3);
    % Set transmit Apodization so that only one element is active.
    TX(n).Apod(n+ElOffset) = 1.0;    % Only one active transmitter for each TX.
end

% Specify TPC structure.
TPC(1).hv = P.HV;

% Specify TGC Waveform structure.
TGC.CntrlPts = [478, 602, 703, 812, 906, 973, 985, 1022];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Receive structure arrays.
% - We need P.Elnum Receive structures for each frame.
maxAcqLength = ceil(sqrt(P.aperture^2 + P.endDepth^2 - 2*P.aperture*P.endDepth*cos(P.theta-pi/2)) - P.startDepth);
wlsPer128 = 128/(2*2); % wavelengths in 128 samples for 2 samplesPerWave
Receive = repmat(struct('Apod', [zeros(1,ElOffset), ones(1,P.Elnum), zeros(1,Trans.numelements-P.Elnum-ElOffset)], ...
                        'startDepth', P.startDepth, ...
                        'endDepth', P.startDepth + wlsPer128*ceil(maxAcqLength/wlsPer128), ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode','BS100BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, P.Elnum*Resource.RcvBuffer(1).numFrames);
% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    Receive(P.Elnum*(i-1)+1).callMediaFunc = 1;
    for j = 1:P.Elnum
        Receive(P.Elnum*(i-1)+j).framenum = i;
        Receive(P.Elnum*(i-1)+j).acqNum = j;
    end
end

% Specify Recon structure arrays.
% - We need one Recon structure.  Each frame will use
%   P.Elnum ReconInfo structures, since we are using P.Elnum
%   synthetic aperture acquisitions.
Recon = struct('senscutoff', 0, ...
               'pdatanum', 1, ...
               'rcvBufFrame', -1, ...     % use most recently transferred frame
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...  % auto-increment ImageBuffer each recon
               'RINums', 1:P.Elnum);

% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 'accumIQ', ...  % accumulate IQ data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 1), 1, P.Elnum);
% - Set specific ReconInfo attributes.
if P.Elnum>1
    ReconInfo(1).mode = 'replaceIQ';  % replace IQ data
    for j = 1:P.Elnum  % For each row in the column
        ReconInfo(j).txnum = j;
        ReconInfo(j).rcvnum = j;
    end
    ReconInfo(P.Elnum).mode = 'accumIQ_replaceIntensity'; % accum & detect
else
    ReconInfo(1).mode = 'replaceIntensity';
end

% Specify Process structure array.
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',5.0,...            % pgain is image processing gain
                         'reject',2,...      % reject level
                         'persistMethod','simple',...
                         'persistLevel',1.0,...
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
SeqControl(1).command = 'jump'; % jump back to start.
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(2).argument = 100;  % 100 usec
SeqControl(3).command = 'returnToMatlab';
SeqControl(4).command = 'timeToNextAcq';  % time between frames
SeqControl(4).argument = 50000;  % 50 ms
nsc = 5; % nsc is count of SeqControl objects

% Specify Event structure array.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    for j = 1:P.Elnum  % Acquire frames
        Event(n).info = 'Acquisition.';
        Event(n).tx = j;
        Event(n).rcv = P.Elnum*(i-1)+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;
    end   
    Event(n-1).seqControl = [4 nsc];  % append transferToHost to last acquisition
    SeqControl(nsc).command = 'transferToHost';    
    nsc = nsc+1;

    Event(n).info = 'recon and process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 1;
    Event(n).seqControl = 3;
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
MinMaxVal = [10,150,P.endDepth]; % default unit is wavelength
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
frameRateFactor = 1;

% Save all the structures to a .mat file.
save('MatFiles/P5_64vN_FMC');

%% Automatic VSX Execution:
% Uncomment the following line to automatically run VSX every time you run
% this SetUp script (note that if VSX finds the variable 'filename' in the
% Matlab workspace, it will load and run that file instead of prompting the
% user for the file to be used):

% filename = 'P5_64vN_FMC'; VSX;


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
    PData(1).Region = struct(...
                'Shape',struct('Name','SectorFT', ...
                'Position',[0,0,-P.radius], ...
                'z',P.startDepth, ...
                'r',P.radius+P.endDepth, ...
                'angle',P.rayDelta, ...
                'steer',0));
    PData(1).Region = computeRegions(PData(1));
    assignin('base','PData',PData);
    evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);');

    Receive = evalin('base', 'Receive');
    wlsPer128 = evalin('base','wlsPer128');
    maxAcqLength = ceil(sqrt(P.aperture^2 + P.endDepth^2 - 2*P.aperture*P.endDepth*cos(P.theta-pi/2)) - P.startDepth);
    for i = 1:size(Receive,2)
        Receive(i).endDepth = P.startDepth + wlsPer128*ceil(maxAcqLength/wlsPer128);
    end
    assignin('base','Receive',Receive);
    evalin('base','TGC.rangeMax = P.endDepth;');
    evalin('base','TGC.Waveform = computeTGCWaveform(TGC);');
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'PData','InterBuffer','ImageBuffer','DisplayWindow','Receive','TGC','Recon'};
    assignin('base','Control', Control);
    assignin('base', 'action', 'displayChange');
end
