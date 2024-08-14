% SetUpL11_4vFlash_2SysSync.m - Example of phase synchronizing two systems by sharing the 250 MHz clock and triggers
%                                           Based on SetUpL11_5vFlash imaging with single plane wave transmit.
% Description:
%   Sequence programming file for L11-5v Linear array, using a plane wave
%   transmit and single receive acquisition. All 128 transmit and receive
%   channels are active for each acquisition. Processing is asynchronous
%   with respect to acquisition.
%
%   The original script was modified to use in a 2-system shared clock configuration, in which an HDMI cable passes the 250 MHz clock
%   from the "master" system clock OUT to the "slave" system clock IN port, and BNC cables are used to trigger the sequences.
%   The master system triggers the acquisition for both systems, using a delayed Trigger OUT command.
%   See "Application Note_Phase-Accurate Synchronization of Two Vantage Systems" for further documentation.
%   Note the sections of code below identified by the comment: "% === add for 2-system sync ==="
%
% Testing: Tested with software release 2.10.4 on two Vantage 256 systems.
%
% Last update:
% 05/06/2020 - Update to SW 4.3 format for new user UIControls and External function definitions (VTS 1691).
%   More info:(?/Example_Scripts/Vantage_Features/New UI Scheme/SetUpL11_5vFlash_NewUI)
% 12/09/2019 - 4.2.0 verif. testing
% 12/15/2015 - Update to SW 3.0 format
% 7/31/2015 - Add synchronization logic. Tested with Software 2.10.4.
% 9/3/2015  - Refined Master/Slave timing; Tested with 2.11.6
%
% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use

clear all

% Specify P structure array.
P.startDepth = 5;   % Acquisition depth in wavelengths
P.endDepth = 192;   % This should preferrably be a multiple of 128 samples.

filename = ('L11-5vFlash_2SysSync'); % assigning the matfile name to the variable 'filename' permits avoiding the VSX query

% === add for 2-system sync ===
master = 1; % =1 Master   =0 Slave  Identify if this is the Master computer or the Slave
if master
master2slaveDelay = 0.020;  % (microsecs) measured delay between the master and slave start times.
                            % Must measure this delay for every installation because it will vary depending on cable lengths
                            % and other hardware details. Quantized in 4 ns increments.
else
    master2slaveDelay = 0;
end
Resource.VDAS.dmaTimeout = 10*1000; % (ms) time software sequencer will wait for 'transferToHost'
                                    % set this time long enough to permit launching the Slave sequence
                                    % and then launching the Master sequence
% ===

% Define system parameters.
Resource.Parameters.numTransmit = 128;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;    % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
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
Resource.RcvBuffer(1).numFrames = 100;       % 100 frames used for RF cineloop.
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
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).numFrames = 20;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = gray(256);

% Specify Transmit waveform structure.
TW.type = 'parametric';
TW.Parameters = [Trans.frequency, .67, 2, 1];

% Specify TX structure array.
TX(1).waveform = 1;            % use 1st TW structure.
TX(1).Origin = [0.0,0.0,0.0];  % flash transmit origin at (0,0,0).
TX(1).focus = 0;
TX(1).Steer = [0.0,0.0];       % theta, alpha = 0.
TX(1).Apod = ones(1,Trans.numelements);
TX(1).Delay = computeTXDelays(TX(1));
% === add for 2-system sync ===
extraDelay = master2slaveDelay * Trans.frequency;
TX(1).Delay = TX(1).Delay + extraDelay; % extraDelay must be zero for the slave system

TX(2) = TX(1);                 % dummy transmit structure that does not turn on transmitters but does allow a trigger out
TX(2).Apod = zeros(1,Trans.numelements);
% ===

% Specify TGC Waveform structure.
TGC.CntrlPts = [0,200,344,452,606,747,870,920];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Receive structure arrays -
%   endDepth - add additional acquisition depth to account for some channels
%              having longer path lengths.
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
                        'callMediaFunc', 1),1,Resource.RcvBuffer(1).numFrames);

% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    % -- Acquisition for full frame.
    Receive(i).framenum = i;
end

% Specify Recon structure arrays.
Recon = struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame', -1, ...     % use most recently transferred frame
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
SeqControl(1).command = 'jump'; % jump back to event in argument set below
JUMP = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between frames
SeqControl(2).argument = 20*1000;  % 20 msec for 50 fps
TTNA = 2;
SeqControl(3).command = 'returnToMatlab';
RTNMatlab = 3;
nsc = 4;    % nsc-1 is count of SeqControl objects

% === add for 2-system sync ===
SeqControl(nsc).command = 'triggerOut';
SeqControl(nsc).condition = 'syncSYNC_CLK';
SeqControl(nsc).argument = 10*250; % 250 counts per microsecond (up to about 4100 microsecs)
TrigOUT = nsc;
nsc = nsc+1;

SeqControl(nsc).command = 'triggerIn';
SeqControl(nsc).condition = 'Trigger_1_Rising';
SeqControl(nsc).argument = 0; % (integer: 250 ns counts ... 8-bit counter)  =0 waits forever
TrigIN1 = nsc;
nsc = nsc+1;
% ===

% SeqControl(nsc).command = 'sync'; % synchronize SW and HW sequencers (optional)
% SYNC = nsc;
% nsc = nsc+1;

n = 1; % n is count of Events

% === add for 2-system sync ===
Event(n).info = 'Initial delayed Trigger OUT';
Event(n).tx = 2;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = TrigOUT;
n=n+1;

SeqControl(1).argument = n; % set jump back to top of loop
% ===

% Acquire all frames defined in RcvBuffer
for i = 1:Resource.RcvBuffer(1).numFrames

    Event(n).info = 'Full aperture acquisition';
    Event(n).tx = 1;
    Event(n).rcv = i;
    Event(n).recon = 0;
    Event(n).process = 0;
    % === add for 2-system sync ===
    if master
        Event(n).seqControl = [TrigIN1, TTNA];
                                                % === triggerIN cancels any running TTNA timer (but allows setting TTNA for next transmit)
    else
        Event(n).seqControl = [TrigIN1];
    end
    % ===
    n = n+1;

    Event(n).info = 'DMA';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = [nsc,nsc+1];
       SeqControl(nsc).command = 'transferToHost';
        % === add for 2-system sync ===
       SeqControl(nsc+1).command = 'waitForTransferComplete';   % used to inform the user that acquisition hasn't occured,
                                                                % presumably because a trigger was not received within DMA timeout
       SeqControl(nsc+1).argument = nsc;
       nsc = nsc + 2;
       % ===
    n = n+1;

    Event(n).info = 'Reconstruct & Process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 1;
    if floor(i/5) == i/5     % Exit to Matlab every 5th frame
        Event(n).seqControl = RTNMatlab;
    else
        Event(n).seqControl = 0;
    end
    n = n+1;

    % === add for 2-system sync ===
    Event(n).info = 'delayed Trigger OUT';
    Event(n).tx = 2;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = TrigOUT;

    %                                             % if  SYNC (optional) times out, DMA + Recon + Process must be taking longer than expected
    %                                             % SYNC should only be used on the Master, and only if the frame rate is
    %                                             % to be controlled by software processes (e.g. external functions)
    n = n+1;
    % ===

end

Event(n).info = 'Jump back';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = JUMP;


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
    'ValueFormat','%3.0f','Callback',@RangeChangeCallback);

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 5;

% Save all the structures to a .mat file.
save(['MatFiles/',filename]);
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