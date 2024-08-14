% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use.
%
% File name: SetUp_NDT_FMC_UTA64LEMO_2El.m - Example of full matrix
% capature using the 64 LEMO UTA with 2 LEMO connectors activated
%
%
% Description:
%   Sequence programming file for full matrix capture using the 64 LEMO UTA
%   with 2 LEMOs activated. For each acquisition, one LEMO connector is activated
%   for transmit and both LEMOs are activated for receive. The activated LEMOs
%   are set by the Trans.ConnectorES array. The display shows the RF data
%   acquired when the ith element transmits and jth elemetn receive. The i
%   and j are set by "Trans El No" and "Recv El No" on GUI.
%
% Last update:
% 06/23/2020 - Update to SW 4.3 format for new user UIControls and External function definitions (VTS 1691).
%   More info:(?/Example_Scripts/Vantage_Features/New UI Scheme/SetUpL11_5vFlash_NewUI)

clear all

P.startDepth = 0;   % Acquisition depth in wavelengths
P.endDepth = 80;    % This should preferrably be a multiple of 128 samples.

P.HV = 25;          % preset voltage

% Specify system parameters.
Resource.Parameters.numTransmit = 64;     % number of transmit channels.
Resource.Parameters.numRcvChannels = 64;  % number of receive channels.
Resource.Parameters.speedOfSound = 6320;
Resource.Parameters.simulateMode = 0;
Resource.Parameters.Connector = 1;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Specify Trans structure array.
Trans.name = 'Custom-5MHz-2El';
Trans.units = 'mm';                                         % required in Gen3 to prevent default to mm units
Trans.numelements = 2;
        Trans.frequency = 5;                                % nominal frequency in MHz
        Trans.Bandwidth = [5*(1-0.67/2), 5*(1+0.67/2)];
        Trans.type = 0;                                     % Array geometry is linear (x values only).
        Trans.id = hex2dec('0000');
        Trans.connType = 11;                                % 64 LEMO UTA
        Trans.elevationApertureMm = 4;                      % active elevation aperture in mm
        Trans.elevationFocusMm = 1000000;                   % nominal elevation focus depth from lens on face of transducer
        Trans.elementWidth = 0.30;                          % element width in mm; assumes 50 micron kerf
        Trans.spacingMm = 0.35;                             % Spacing between elements in mm.
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
        Trans.ConnectorES = [16 26]'; % LEMO 16 and 26 are activated. LEMO 16 is seen as element 1. LEMO 26 is seen as element 2.

        % Now convert all units as required, based on Trans.units
        scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000); % conversion factor from mm to wavelengths
        % regardless of units, always provide spacing in wavelengths
        Trans.spacing = Trans.spacingMm * scaleToWvl;   % Spacing between elements in wavelengths


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
Resource.RcvBuffer(1).rowsPerFrame = 2560*Trans.numelements;
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 10;     % 10 frames for RF data.

% Specify Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,.67,2,1];

% Specify P.Elnum TX structure arrays. Transmit on element n in the array for event n.
TX = repmat(struct('waveform', 1, ...
                   'Origin', zeros(1,3), ...
                   'focus', 0, ...
                   'Steer', [0.0,0.0], ...
                   'Apod', zeros(1,Trans.numelements), ...
                   'Delay', zeros(1,Trans.numelements)), 1, Trans.numelements);

% - Set event specific TX attributes.
for n = 1:Trans.numelements
    % Set transmit Origins to positions of elements.
    TX(n).Origin = scaleToWvl*Trans.ElementPos(n,1:3);
    % Set transmit Apodization so that only one element is active.
    TX(n).Apod(n) = 1.0;    % Only one active transmitter for each TX.
end

% Specify TPC structure.
TPC(1).hv = P.HV;

% Specify TGC Waveform structure.
TGC.CntrlPts = [500, 500, 500, 500, 500, 500, 500, 500];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Receive structure arrays.
% - We need P.Elnum Receive structures for each frame.
Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', P.endDepth, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode','NS200BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, Trans.numelements*Resource.RcvBuffer(1).numFrames);
% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    Receive(Trans.numelements*(i-1)+1).callMediaFunc = 0;
    for j = 1:Trans.numelements
        Receive(Trans.numelements*(i-1)+j).framenum = i;
        Receive(Trans.numelements*(i-1)+j).acqNum = j;
    end
end

% Specify Process structure array.
Process(1).classname = 'External';
Process(1).method = 'plotChannelRF';
Process(1).Parameters = {'srcbuffer','receive',...  % name of buffer to process.
                         'srcbufnum',1,...
                         'srcframenum',-1,...
                         'dstbuffer','none'};

% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; % jump back to start.
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(2).argument = 100;  % 100 usec
SeqControl(3).command = 'returnToMatlab';
nsc = 4; % nsc is count of SeqControl objects

% Specify Event structure array.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    for j = 1:Trans.numelements                % Acquire frames
        Event(n).info = 'Acquisition.';
        Event(n).tx = j;
        Event(n).rcv = Trans.numelements*(i-1)+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;
    end
    Event(n-1).seqControl = nsc; % modify last acquisition Event's seqControl
      SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
      nsc = nsc+1;

    Event(n).info = 'Display RF line';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
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

% - transmit channel display slider
UI(1).Control = VsSliderControl('LocationCode','UserB2',...
    'Label','Trans El No.',...
    'SliderMinMaxVal',[1,Trans.numelements,1],...
    'SliderStep',[1/(Trans.numelements-1),1/1],...
    'ValueFormat','%3.0f',...
    'Callback',@(~,~,UIValue) assignin('base','myPlotTransChnl',round(UIValue)));

% - recieve channel dispaly slider
UI(2).Control = VsSliderControl('LocationCode','UserB1',...
    'Label','Recv El No.',...
    'SliderMinMaxVal',[1,Trans.numelements,Trans.numelements],...
    'SliderStep',[1/(Trans.numelements-1),1/1],'ValueFormat','%3.0f',...
	'Callback',@(~,~,UIValue) assignin('base','myPlotRecvChnl',round(UIValue)));

% External function definition.
EF(1).Function = vsv.seq.function.ExFunctionDef('plotChannelRF', @plotChannelRF);

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 1;

% Save all the structures to a .mat file, set filename and run VSX
save('MatFiles/NDT_FMC_UTA64LEMO_2El');

%% Automatic VSX Execution:
% Uncomment the following line to automatically run VSX every time you run
% this SetUp script (note that if VSX finds the variable 'filename' in the
% Matlab workspace, it will load and run that file instead of prompting the
% user for the file to be used):

% filename = 'NDT_FMC_UTA64LEMO_2El'; VSX;


% **** Callback routines used by External function definition (EF) ****

function plotChannelRF(RData)
    persistent myHandle
    Receive = evalin('base','Receive');
    Trans = evalin('base','Trans');
    % If myPlotRecvChnl exists, read it for the Recv channel to plot.
    if evalin('base','exist(''myPlotRecvChnl'',''var'')')
        RecvChannel = evalin('base','myPlotRecvChnl');
    else
        RecvChannel = Trans.numelements;  % Recv Channel no. to plot
    end
    % If myPlotTransChnl exists, read it for the Trans channel to plot.
    if evalin('base','exist(''myPlotTransChnl'',''var'')')
        TransChannel = evalin('base','myPlotTransChnl');
    else
        TransChannel = 1;  % Trans Channel no. to plot
    end

    % Create the figure if it does not exist.
    if isempty(myHandle)||~ishandle(myHandle)
        figure('name','Receive Signal','NumberTitle','off');
        myHandle = axes('XLim',[0,Receive(1).endSample],'YLim',[-2048 2048], ...
                        'NextPlot','replacechildren');
    end
    % Find Lemo connector numbers using Trans.Connector
    txLemo = Trans.Connector(TransChannel);
    rcvLemo = Trans.Connector(RecvChannel);
    % Plot the element's RF data.
    plot(myHandle,RData(Receive(TransChannel).startSample:Receive(TransChannel).endSample,Trans.Connector(RecvChannel)));
    title(myHandle, ['Transmit Lemo ', num2str(txLemo), '  Receive Lemo ', num2str(rcvLemo)]);
    drawnow limitrate
end
