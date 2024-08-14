% SetUpL11_4vFlash_2SysDataTransfer_SyncModule_TCPIP.m -
% Example of phase synchronizing four systems by using Sync Module and
% transferring data using TCPIP communication
%
% Description:
%
%   This script shows the use of the Synchronization Module hardware and
%   software configuration needed to share a single 250 MHz system clock
%   between multiple Vantage systems, and thus enable synchronous and
%   phase-aligned operation of multiple hardware sequencers.
%
%   Two sets of event sequence are included in the script, including
%   1) regular plane wave image, one frame is enough and
%   2) high frame rate (HFR) data acquition along with data transfer
%
% Testing: Tested with software release 3.2.2 on two Vantage systems.
%
% Last update:
% 05/06/2020 - Update to SW 4.3 format for new user UIControls and External function definitions (VTS 1691).
%   More info:(?/Example_Scripts/Vantage_Features/New UI Scheme/SetUpL11_5vFlash_NewUI)
% 07/19/2017 - Finalize the script
%
% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use

clear all
close all
delete(instrfind)

% Master or Slave?
master = 1;   % =1 Master   =0 Slave  Identify if this is the Master computer or the Slave

% The tcpip object will be different between Master and Slave PC
% Please increase the Java Heap Memory to the maximum to prevent an error
% message about "not enough memory to create the input buffer"
% The IP address should be defined correctly.
masterIP = '192.168.0.1';
slaveIP = '192.168.0.2';
portNum = 11335;

if master
    tcpipMtoS = tcpip(slaveIP,portNum,'NetworkRole','server','TimeOut',30,'InputBufferSize',2^30);
else
    tcpipStoM = tcpip(masterIP,portNum,'NetworkRole','client','TimeOut',30,'OutputBufferSize',2^30);
end

% Specify numFrame used in HFR data acquisition and combination
numFrameHFR = 1000;

% Specify P structure array.
P.startDepth = 5;   % Acquisition depth in wavelengths
P.endDepth = 192;   % This should preferrably be a multiple of 128 samples.
filename = ('L11-5vFlash_2SysDataTransfer_SyncModule'); % assigning the matfile name to the variable 'filename' permits avoiding the VSX query

Resource.Parameters.SystemLED = {'running', 'pausedSync', 'activeTxAndOrRx', 'pausedOnMultiSysSync'};
Resource.VDAS.dmaTimeout = 10*1000; % (ms) time software sequencer will wait for 'transferToHost'
% set this time long enough to permit launching the Slave sequence
% and then launching the Master sequence

% Define system parameters.
Resource.Parameters.numTransmit = 128;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;    % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 1;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Specify Trans structure array.
Trans.name = 'L11-5v';
Trans.units = 'mm'; % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);  % L11-5v transducer is 'known' transducer so we can use computeTrans.
Trans.maxHighVoltage = 50;  % set maximum high voltage limit for pulser supply.

% Specify PData structure array.
PData(1).PDelta = [Trans.spacing, 0, 0.5];
PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3)); % startDepth, endDepth and pdelta set PData(1).Size.
PData(1).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(3) = 1;      % single image page
PData(1).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of upper lft crnr.
% No PData.Region specified, so a default Region for the entire PData array will be created by computeRegions.

% Specify Media object. 'pt1.m' script defines array of point targets.
pt1;
Media.function = 'movePoints';

% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 4096;   % this size allows for maximum range
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 1;
Resource.RcvBuffer(2).datatype = 'int16';
Resource.RcvBuffer(2).rowsPerFrame = 4096*numFrameHFR;   % this size allows for maximum range
Resource.RcvBuffer(2).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(2).numFrames = 1;
Resource.InterBuffer(1).numFrames = 1;  % one intermediate buffer needed.
Resource.ImageBuffer(1).numFrames = 10;
Resource.DisplayWindow(1).Title = 'L11-5vFlash';
Resource.DisplayWindow(1).pdelta = 0.35;
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
TW.type = 'parametric';
TW.Parameters = [Trans.frequency, .67, 2, 1];   % A, B, C, D

% Specify TX structure array.
TX(1).waveform = 1;            % use 1st TW structure.
TX(1).Origin = [0.0,0.0,0.0];  % flash transmit origin at (0,0,0).
TX(1).focus = 0;
TX(1).Steer = [0.0,0.0];       % theta, alpha = 0.
TX(1).Apod = ones(1,Trans.numelements);
TX(1).Delay = computeTXDelays(TX(1));

% Specify TGC Waveform structure.
TGC.CntrlPts = [0,200,344,452,606,747,870,920];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Receive structure arrays -
%   endDepth - add additional acquisition depth to account for some channels
%              having longer path lengths.
%   InputFilter - The same coefficients are used for all channels. The
%              coefficients below give a broad bandwidth bandpass filter.
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
    'startDepth', P.startDepth, ...
    'endDepth', maxAcqLength,...
    'TGC', 1, ...
    'bufnum', 1, ...
    'framenum', 1, ...
    'acqNum', 1, ...
    'sampleMode', 'NS200BW', ...
    'mode', 0, ...
    'callMediaFunc', 1),1,Resource.RcvBuffer(1).numFrames+Resource.RcvBuffer(2).numFrames*numFrameHFR);

for i = 1:Resource.RcvBuffer(2).numFrames
    for j = 1:numFrameHFR
        rcvNum = 1 + numFrameHFR*(i-1)+j;
        Receive(rcvNum).bufnum = 2;
        Receive(rcvNum).framenum = i;
        Receive(rcvNum).acqNum = j;
    end
end

% Specify Recon structure arrays.
Recon = struct('senscutoff', 0.6, ...
    'pdatanum', 1, ...
    'rcvBufFrame', 1, ...     % only one frame is required for realtime B-mode imaging
    'IntBufDest', [1,1], ...
    'ImgBufDest', [1,-1], ...  % auto-increment ImageBuffer each recon
    'RINums', 1);

% Define ReconInfo structures.
ReconInfo = struct('mode', 'replaceIntensity', ...
    'txnum', 1, ...
    'rcvnum', 1, ...
    'regionnum', 1);

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
    'interpMethod','4pt',...  %method of interp. (1=4pt)
    'grainRemoval','none',...
    'processMethod','none',...
    'averageMethod','none',...
    'compressMethod','power',...
    'compressFactor',40,...
    'mappingMethod','full',...
    'display',1,...      % display image after processing
    'displayWindow',1};

% Process(2) uses EF2 for timing control
TimingControl = 2;
Process(TimingControl).classname = 'External';
Process(TimingControl).Parameters = {'srcbuffer','none'};
if master
    Process(TimingControl).method = 'waitForSlave';
else
    Process(TimingControl).method = 'confirmReady';
end

% Process(2) uses EF3 for UI control
UIControl = 3;
Process(UIControl).classname = 'External';
Process(UIControl).Parameters = {'srcbuffer','none'};
if master
    Process(UIControl).method = 'saveSettings';
else
    Process(UIControl).method = 'loadSettings';
end

% Process(4) uses EF4 for data transfer
DataTransfer = 4;
Process(DataTransfer).classname = 'External';
Process(DataTransfer).Parameters = {'srcbuffer','receive',...
    'srcbufnum', 2,...
    'srcframenum',0,... % process the all frames in a RcvBuffer
    'dstbuffer','none'};
if master
    Process(DataTransfer).method = 'DataTransferMaster';
else
    Process(DataTransfer).method = 'DataTransferSlave';
end

% Specify SeqControl structure arrays.
nsc = 1;
SeqControl(nsc).command = 'jump'; % jump back to event in argument set below
SeqControl(nsc).condition = 'exitAfterJump'; % only required if it doesn't jump back to event 1
JUMP = nsc;     nsc = nsc+1;      % set a mnemonic for the jump command index

SeqControl(nsc).command = 'timeToNextAcq';
SeqControl(nsc).argument = 50;  % HFR for data combination
TTNA = nsc;     nsc = nsc+1;

SeqControl(nsc).command = 'multiSysSync';
SeqControl(nsc).condition = 'normal';       % 'normal' uses internal trigger; other options are 'BNC_Rising' and 'BNC_Falling'
MSSYNC= nsc;   nsc = nsc+1;

SeqControl(nsc).command = 'noop';
SeqControl(nsc).argument = 500/.2;  % 0.5 ms
NOOP = nsc;     nsc = nsc+1;

SeqControl(nsc).command = 'sync'; % synchronize SW and HW sequencers
SeqControl(nsc).argument = 5e6;     % 5s
SYNC = nsc;    nsc = nsc+1;

n = 1;
SeqControl(JUMP).argument = n; % set jump back to top of loop

%% Master Image Sequence ===
if master

    % HW and SW should be synchronized while jump back to the first event
    Event(n).info = 'Wait for SLAVE';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = TimingControl;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = 'Sync';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = [SYNC,NOOP];
    n = n+1;

    Event(n).info = 'Full aperture acquisition';
    Event(n).tx = 1;
    Event(n).rcv = 1;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = MSSYNC;
    n = n+1;

    Event(n).info = 'DMA';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = nsc;
    SeqControl(nsc).command = 'transferToHost';
    nsc = nsc + 1;
    n = n+1;

    Event(n).info = 'Reconstruct & Process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 1;
    Event(n).seqControl = 0;
    n = n+1;

    % === Save settings for the slave
    Event(n).info = 'Save UI settings';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = UIControl;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = 'Sync, pause HW';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = SYNC;
    n = n+1;

    Event(n).info = 'Jump back';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = JUMP;
    n = n+1;

    %% Master Saving Data Sequence ===

    nStartTransferRF = n;

    Event(n).info = 'Wait for SLAVE';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = TimingControl;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = 'Sync';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = [SYNC,NOOP];
    n = n+1;

    % Acquire all frames defined in RcvBuffer
    for i = 1:Resource.RcvBuffer(2).numFrames
        for j = 1:numFrameHFR
            Event(n).info = 'Full aperture acquisition';
            Event(n).tx = 1;
            Event(n).rcv = 1+numFrameHFR*(i-1)+j;
            Event(n).recon = 0;
            Event(n).process = 0;
            if i< Resource.RcvBuffer(2).numFrames*numFrameHFR
                Event(n).seqControl = [TTNA,MSSYNC];  % TTNA is used to control the frame rate within one DMA transfer
            else
                Event(n).seqControl = MSSYNC;  % TTNA is not required for the last frame
            end
            n = n+1;
        end

        Event(n).info = 'DMA';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = nsc;
        SeqControl(nsc).command = 'transferToHost';
        nsc = nsc + 1;
        n = n+1;
    end

    Event(n).info = 'Sync';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = SYNC;
    n = n+1;

    % === Data Transfer
    Event(n).info = 'Transfer RF to Master';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = DataTransfer;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = 'Jump back';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = JUMP;

    %% Slave Image Sequence
else

    % HW and SW should be synchronized while jump back to the first event
    Event(n).info = 'Sync, pause HW after jump back';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = SYNC;
    n = n+1;

    % Tell Master to start the sequence
    Event(n).info = 'Trigger master';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = TimingControl;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = 'Full aperture acquisition';
    Event(n).tx = 1;
    Event(n).rcv = 1;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = MSSYNC;
    n = n+1;

    Event(n).info = 'DMA';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = nsc;
    SeqControl(nsc).command = 'transferToHost';
    nsc = nsc + 1;
    n = n+1;

    Event(n).info = 'Reconstruct & Process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 1;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = 'Sync, pause HW';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = SYNC;
    n = n+1;

    % === Load settings
    Event(n).info = 'Load UI settings';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = UIControl;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = 'Jump back';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = JUMP;
    n = n+1;

    %% SLAVE: Saving Data Sequence ===

    nStartTransferRF = n;

    % HW and SW should be synchronized while jump back to the first event
    Event(n).info = 'Sync, pause HW';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 0;      % reconstruction
    Event(n).process = 0;    % processing
    Event(n).seqControl = SYNC;
    n = n+1;

    % Tell Master to start the sequence
    Event(n).info = 'Trigger Master';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = TimingControl;
    Event(n).seqControl = 0;
    n = n+1;

    % Acquire all frames defined in RcvBuffer
    for i = 1:Resource.RcvBuffer(2).numFrames

        for j = 1:numFrameHFR
            Event(n).info = 'Full aperture acquisition';
            Event(n).tx = 1;
            Event(n).rcv = 1+numFrameHFR*(i-1)+j;
            Event(n).recon = 0;
            Event(n).process = 0;
            Event(n).seqControl = MSSYNC;
            n = n+1;
        end

        Event(n).info = 'DMA';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = nsc;
        SeqControl(nsc).command = 'transferToHost';
        nsc = nsc + 1;
        n = n+1;

    end

    Event(n).info = 'Sync';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = SYNC;
    n = n+1;

    Event(n).info = 'Transfer Data to Master';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = DataTransfer;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = 'Jump back';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = JUMP;

end

%% User specified UI Control Elements

import vsv.seq.uicontrol.VsSliderControl

% - Sensitivity Cutoff
UI(1).Control =  VsSliderControl('LocationCode','UserB7',...
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
    'ValueFormat','%3.0f', 'Callback', @RangeChangeCallback);

if master
    UI(3).Control = VsSliderControl('LocationCode','UserC1',...
        'Label','Transfer RF','Callback',@TransferRFcallback);
end

% External function definitions

import vsv.seq.function.ExFunctionDef

if master
    EF(1).Function = vsv.seq.function.ExFunctionDef('FreezeFunctionMaster', @FreezeFunctionMaster);
else
    EF(1).Function = vsv.seq.function.ExFunctionDef('FreezeFunctionSlave', @FreezeFunctionSlave);
end
if master
    EF(2).Function = vsv.seq.function.ExFunctionDef('waitForSlave', @waitForSlave);
else
    EF(2).Function = vsv.seq.function.ExFunctionDef('confirmReady', @confirmReady);
end
if master
    EF(3).Function = vsv.seq.function.ExFunctionDef('saveSettings', @saveSettings);
else
    EF(3).Function = vsv.seq.function.ExFunctionDef('loadSettings', @loadSettings);
end
if master
    EF(4).Function = vsv.seq.function.ExFunctionDef('DataTransferMaster', @DataTransferMaster);
else
    EF(4).Function = vsv.seq.function.ExFunctionDef('DataTransferSlave', @DataTransferSlave);
end

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = Resource.RcvBuffer(1).numFrames;

if master
    fopen(tcpipMtoS);
else
    fopen(tcpipStoM);
end

% Save all the structures to a .mat file.
save(['MatFiles/',filename]); if ~master, pause(1), end
VSX % setting 'filename' variable permits running automatically upon running the setup script


%% **** Callback routines used by UIControls (UI) ****

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

    evalin('base','PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));');
    evalin('base','PData(1).Region = computeRegions(PData(1));');
    evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);');
    Receive = evalin('base', 'Receive');
    maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
    for i = 1:size(Receive,2)
        Receive(i).endDepth = maxAcqLength;
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

function TransferRFcallback(~, ~, ~)
    tcpipMtoS = evalin('base','tcpipMtoS');
    fwrite(tcpipMtoS,1); % pass one command
    fprintf(tcpipMtoS,'TransferRF = 1;');
    % assign TransferRF to the workspace to make the Slave update the start
    % event earlier than the Master
    assignin('base','TransferRF',1);
end


%% **** Callback routines used by External function definition (EF) ****

function FreezeFunctionMaster()
    evalin('base','stop(freezeTimer);');
    evalin('base','delete(freezeTimer);');
    evalin('base','clear freezeTimer;');

    tcpipMtoS = evalin('base','tcpipMtoS');

    while isequal(get(findobj('String','Freeze'),'Value'),1)
        pause(0.25),
        fprintf(tcpipMtoS,'freeze');
    end
    fprintf(tcpipMtoS,'unfreeze');
end

function FreezeFunctionSlave()
    tcpipStoM = evalin('base','tcpipStoM');
    if strcmp(fscanf(tcpipStoM,'%s'),'unfreeze')
        h = findobj('String','Freeze');
        evalin('base','stop(freezeTimer);');
        evalin('base','delete(freezeTimer);');
        evalin('base','clear freezeTimer;');
        set(h,'Value',0);
        feval(get(h,'Callback'),h);
    end
end

function saveSettings()
    tcpipMtoS = evalin('base','tcpipMtoS');
    command = [];

    % other Uicontrols
    UItagName = {'UI','ProcessTool','filterTool'};
    for n = 1:3
        h = findobj('tag',UItagName{n});
        cH = get(h,'CurrentObject'); % currentHandle
        if ~isempty(cH)&&~isequal(cH,h)
            break
        end
    end

    if ~isempty(cH)&&~isequal(cH,h)
        switch get(cH,'Style')
            case {'slider','popupmenu','togglebutton'}
                pause(0.5)
                command{1} = ['nameTag = ''',get(cH,'tag'),''';'];
                command{2} = 'cH = findobj(''Tag'',nameTag);';
                command{3} = ['set(cH,''Value'',',num2str(get(cH,'Value')),');'];
            case 'edit'
                command{1} = ['nameTag = ''',get(cH,'tag'),''';'];
                command{2} = 'cH = findobj(''Tag'',nameTag);';
                command{3} = ['set(cH,''String'',''',get(cH,'String'),''');'];
            otherwise
        end
        set(h,'CurrentObject',h); % set it back to default after execution
    end

    % close UI?
    if evalin('base','vsExit')
        command{1} = 'feval(get(h,''CloseRequestFcn''),h);';
    end

    % pass command to slave
    fwrite(tcpipMtoS,length(command)) % pass command length
    for i = 1:length(command)
        fprintf(tcpipMtoS,command{i});
    end

    % is Transfer RF pressed?
    if evalin('base','exist(''TransferRF'',''var'');')
        if evalin('base','TransferRF')
            nStart = evalin('base','nStartTransferRF');
            Control = evalin('base','Control');
            Control.Command = 'set&Run';
            Control.Parameters = {'Parameters',1,'startEvent',nStart};
            evalin('base',['Resource.Parameters.startEvent =',num2str(nStart),';']);
            assignin('base','Control',Control);
            evalin('base','clear TransferRF;');
            return
        end
    end
end

function loadSettings()
    tcpipStoM = evalin('base','tcpipStoM');
    h = findobj('tag','UI'); % will be used in eval(command)

    % read number of commands
    nCommand = fread(tcpipStoM,1);

    for i = 1:nCommand
        command = fscanf(tcpipStoM);
        if ~isempty(command(1:end-1))
            eval(command)
        end
    end

    % is Transfer RF pressed?
    if exist('TransferRF','var')
        nStart = evalin('base','nStartTransferRF');
        Control = evalin('base','Control');
        Control.Command = 'set&Run';
        Control.Parameters = {'Parameters',1,'startEvent',nStart};
        evalin('base',['Resource.Parameters.startEvent =',num2str(nStart),';']);
        assignin('base','Control',Control);
        clear TransferRF
        return
    end

    % execute callback
    if exist('cH','var')
        CB = get(cH,'Callback');
        if ~iscell(CB)
            feval(CB,cH);
        else
            feval(CB{1},cH,[],CB{2});
        end
    end
end

function DataTransferMaster(RDataMaster)
    tcpipMtoS = evalin('base','tcpipMtoS');

    Receive = evalin('base','Receive');
    RcvDepth = Receive(end).endSample;

    dataSize = size(RDataMaster);
    if isequal(nnz(dataSize),2)
        dataSize(3) = 1;
    end

    MasterBuffer = RDataMaster(1:RcvDepth,1:dataSize(2),1:dataSize(3));
    segmentSize = numel(MasterBuffer);

    % receive data from the Slave
    tic
    disp('Wait for data transfer....');
    SlaveBuffer = fread(tcpipMtoS,segmentSize,'int16');
    disp(['Duration_Of_Transfer = ',num2str(toc,'%2.3f'),' s']);

    % combine it for image reconstruction
    tic
    disp('Wait for data combination....');
    RDataSlave = zeros(dataSize,'int16');
    RDataSlave(1:RcvDepth,:) = reshape(SlaveBuffer,RcvDepth,dataSize(2),dataSize(3));
    combinedRF = [RDataMaster RDataSlave];
    assignin('base','combinedRF',combinedRF);
    disp(['Duration_Of_DataCombination = ',num2str(toc,'%2.3f'),' s']);
    disp('Please refer to combinedRF in the workspace');

    Control = evalin('base','Control');
    Control.Command = 'set';
    Control.Parameters = {'Parameters',1,'startEvent',1};
    evalin('base','Resource.Parameters.startEvent = 1;');
    assignin('base','Control',Control);

    fprintf(tcpipMtoS,'finished'); % tell the slave that data combination is completed
end

function DataTransferSlave(RDataSlave)

    tcpipStoM = evalin('base','tcpipStoM');
    clearBuffer = fread(tcpipStoM,1); % clear one left byte from previous transfer

    Receive = evalin('base','Receive');
    RcvDepth = Receive(end).endSample;

    dataSize = size(RDataSlave);
    if isequal(nnz(dataSize),2)
        dataSize(3) = 1;
    end
    SlaveBuffer = RDataSlave(1:RcvDepth,1:dataSize(2),1:dataSize(3));
    fwrite(tcpipStoM,SlaveBuffer(:),'int16');

    status = fscanf(tcpipStoM,'%s'); % wait Master to combine data

    Control = evalin('base','Control');
    Control.Command = 'set';
    Control.Parameters = {'Parameters',1,'startEvent',1};
    evalin('base','Resource.Parameters.startEvent = 1;');
    assignin('base','Control',Control);
end

function waitForSlave()
    tcpipMtoS = evalin('base','tcpipMtoS');
    if ~strcmp(fscanf(tcpipMtoS,'%s'),'ready')
        error('TCPIP communication is not ready!');
    end
end

function confirmReady()
    tcpipStoM = evalin('base','tcpipStoM');
    fprintf(tcpipStoM,'ready');
end