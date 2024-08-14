% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use.
%
% --- PhotoAcoustic application example script ---
% Generate .mat Sequence Object file for FlashAngles B-mode and PhotoAcoustic (PA) imaging with L11-5v transducer, and run VSX automatically.
% Modified from the SetupL11-5vFlashDoppler script, retaining some extra flexibility that is not required for PA imaging;
% for example, PData, P and other structures can be identical for both modes but are separate here for added flexibility.
%
% The PA data is displayed using the split palette colormap overlay method, using the Power Doppler colormap.
% Note that the compression slider does not work properly when using split palettes.
% When running in simulation (with oneway=0), the media points are imaged as though the medium moves between B-mode and PA
% acquisitions, and thus are slightly offset from each other to better visualize each result. The B-mode or PA mode images
% can be set to dominate the display by adjusting gain and compression and color priority controls.
%   - For the 2D (B-mode) image, 'na' flat focus transmits at steered angles are used.
%   - For the PA acquisition, 'ne' T/R events with output trigger for an external source (laser) are performed,
%       and use receive-only reconstruction (when oneway=1).
%       * The sequence waits for an input trigger (e.g. from laser flash lamp), then pauses while the laser pumps (flash2Qdelay microsecs),
%           then fires the Q-switch using a trigger out that also begins receive acquisition. This approach results in very low
%           jitter between the laser pulse and the ultrasound.
%       * Multiple PA acquisitions are reconstructed and then coherently averaged by accumulation in the Inter Buffer.
%           Incoherent averaging is possible by summing in the Image buffer instead, and relaxes jitter requirements.
%           Though not implemented here, the RF data could be accumulated in hardware prior to transfer
%           using the Receive.mode attribute; doing so also requires precise synchronization with the laser.
%       * The receive-only reconstruction mode is not currently supported in simulation, so to test the PA code we simulate by
%           turning on the transmitters for a conventional T/R acquisition and reconstruction. Toggling between modes is easily
%           done using the 'oneway' parameter (=0 turn on transmitters, =1 turn off transmitters for real PA acquisition) below.
%       * Two TGC structures are defined, but only TGC(1) is set up to be controlled by the GUI sliders. To use the TGC slider controls
%           for the PA acquisitions, use TGC(1) for PA and TGC(2) for 2D.
%
% Last update:
% Aug 10 2015: replaced 'pause' command with 'triggerIn' with 0.5 sec timeout
% 12/15/2015 - Update to SW 3.0 format

clear all

%% === Set commonly modified parameters ========================================

P(1).startDepth = 2;          % Acquisition start depth in wavelengths
P(1).endDepth = 192;          % Acquisition end depth
P(2).startDepth = 0;          % Acquisition start depth in wavelengths
P(2).endDepth = 128;          % Acquisition end depth may be different for PA

% Set 2D parameters
na = 7;      % Set na = number of flash angles for 2D.
    if na>1
        dtheta2D = (36*pi/180)/(na-1); % set dtheta2D to range over +/- 18 degrees.
    else
        dtheta2D = 0; % for a single angle, the angle will be zero.
    end

% Set PA parameters
oneway = 0;     % (logical) oneway=1 turns off the transmitters by setting TX.Apod to zero for all transmitters
flash2Qdelay = 200; % microseconds between trigger input and start of acquisition (which outputs a trigger pulse at time=0)
ne = 1;         % ne = number of acquisitions in PA ensemble for coherent addition in I/Q buffer.
PA_Angle = 0;   % angle of transmit plane wave for use in testing the PA mode in simulation
PA_PRF = 100;   % PA PRF in Hz. To be set in accordance with laser rep rate, when not using the input trigger mode.
                % When using the input trigger mode, remove the TTNA (SeqControl(4)) from the PA events to avoid
                % "missed TTNA" messages.

    if oneway==1
        disp(' *** PhotoAcoustic mode: Using one-way receive-only reconstruction ***')
    else
        disp(' *** Ultrasound Transmit mode: Using conventional T/R reconstruction ***')
    end
% ===============================================================================

%% Specify system parameters(Vantage 128)
Resource.Parameters.numTransmit = 128;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;    % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

%% Specify Trans structure array.
Trans.name = 'L11-5v';
Trans.units = 'mm';             % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);    % computeTrans is used for known transducers.
nElem = Trans.numelements;
Trans.maxHighVoltage = 50;      % set a reasonable high voltage limit.

%% Specify PData structure arrays.
% - 2D PData structure
PData(1).PDelta(1) = 1.0;
PData(1).PDelta(3) = 0.5;
PData(1).Size(1,1) = ceil((P(1).endDepth-P(1).startDepth)/PData(1).PDelta(3)); % rows
PData(1).Size(1,2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1)); % cols
PData(1).Size(1,3) = 1;             % single image page
PData(1).Origin = [-Trans.spacing*63.5,0,P(1).startDepth]; % x,y,z of uppr lft crnr.
% - PA PData structure
PData(2).PDelta(1) = 1.0;
PData(2).PDelta(3) = 0.5;
PData(2).Size(1,1) = ceil((P(2).endDepth-P(2).startDepth)/PData(2).PDelta(3)); % PA window rows
PData(2).Size(1,2) = ceil((Trans.numelements*Trans.spacing)/PData(2).PDelta(1));  % PA window columns
PData(2).Size(1,3) = 1;             % single image page
PData(2).Origin = [-Trans.spacing*63.5,0,P(2).startDepth]; % x,y,z of upper lft crnr.

%% Specify Media object and point displacement function
pt1;
Media.function = 'movePoints';

%% Specify Resources.
% - RcvBuffer(1) is for both 2D and PA acquisitions.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 4096*(na + ne);
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 20;           % 20 frames allocated for RF acqusitions.
% InterBuffer(1) is for 2D reconstructions.
Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).numFrames = 1;          % one intermediate frame needed for 2D.
% InterBuffer(2) is for PA reconstructions.
Resource.InterBuffer(2).datatype = 'complex';
Resource.InterBuffer(2).numFrames = 1;          % one intermediate frame needed for PA.
% ImageBuffer(1) is for 2D image.
Resource.ImageBuffer(1).datatype = 'double';    % image buffer for 2D
Resource.ImageBuffer(1).numFrames = 20;
% ImageBuffer(2) is for PA image.
Resource.ImageBuffer(2).datatype = 'double';    % image buffer for PA
Resource.ImageBuffer(2).numFrames = Resource.ImageBuffer(1).numFrames;
% DisplayWindow is for 2D combined with PA
Resource.DisplayWindow(1).Title = mfilename;
Resource.DisplayWindow(1).pdelta = 0.4;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)]; % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).numFrames = 20;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = grayscaleCFImap;
Resource.DisplayWindow(1).splitPalette = 1;

% ------Specify structures used in Events------
%% Specify Transmit waveforms structure
% - 2D transmit waveform
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,0.67,2,1];
% - PA transmit waveform
TW(2).type = 'parametric';
TW(2).Parameters = [Trans.frequency,0.67,6,1];

%% Specify Transmit beams structure
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'focus', 0.0, ...
                   'Steer', [0.0,0.0], ...
                   'Apod', ones(1,Trans.numelements), ...
                   'Delay', zeros(1,Trans.numelements)), 1, na+1); % na TXs for 2D + 1 for PA
% - Set event specific TX attributes.
if fix(na/2) == na/2       % if na even
    startAngle = (-(fix(na/2) - 1) - 0.5)*dtheta2D;
else
    startAngle = -fix(na/2)*dtheta2D;
end
for n = 1:na   % na transmit events for 2D
    TX(n).Steer = [(startAngle+(n-1)*dtheta2D),0.0];
    TX(n).Delay = computeTXDelays(TX(n));
end
% -- only one TX struct needed for PA
TX(na+1).waveform = 2;
if oneway
    TX(na+1).Apod = zeros(1,Trans.numelements);     % THIS COMMAND TURNS OFF ALL TRANSMITTERS AND INVOKES THE RECEIVE-ONLY BEAMFORMER
else
    TX(na+1).Apod =  ones(1,Trans.numelements);     % This is the conventional T/R condition and invokes the default beamformer
end
TX(na+1).Steer = [PA_Angle,0.0];            % only relevant when transmitters are active
TX(na+1).Delay = computeTXDelays(TX(na+1)); % only relevant when transmitters are active

%% Specify TPC structures.
TPC(1).name = '2D';
TPC(1).maxHighVoltage = 50;

% This allows one to use different transmit profile for PA ... only relevant if transmitters are active
% --- currently TPC(2) is not used ---
TPC(2).name = 'PA';
TPC(2).maxHighVoltage = 35;

%% Analog front end gain settings.
RcvProfile(1).LnaGain = 18;     % 12, 18, or 24 dB  (18=default)
RcvProfile(1).condition = 'immediate';

RcvProfile(2).LnaGain = 24;
RcvProfile(2).condition = 'immediate';

%% Specify Receive structure arrays.
%   We need to acquire all the 2D and PA data within a single RcvBuffer frame.  This allows
%   the transfer-to-host DMA after each frame to transfer a large amount of data, improving throughput.
% - We need 2*na Receives for a 2D frame and ne Receives for a PA frame.
maxAcqLngth2D = sqrt(P(1).endDepth^2 + (Trans.numelements*Trans.spacing)^2) - P(1).startDepth;
maxAcqLngthPA =  sqrt(P(2).endDepth^2 + (Trans.numelements*Trans.spacing)^2) - P(2).startDepth;
wl4sPer128 = 128/(4*2);  % wavelengths in a 128 sample block for 4 smpls per wave round trip.
% wl2sPer128 = 128/(2*2);  % wavelengths in a 128 sample block for 2 smpls per wave round trip.
Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', P(1).startDepth, ...
                        'endDepth', P(1).startDepth + wl4sPer128*ceil(maxAcqLngth2D/wl4sPer128), ...
                        'TGC', 1, ...           % TGC(1) is tied to the GUI sliders
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, (na+ne)*Resource.RcvBuffer(1).numFrames);

% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    k = (na + ne)*(i-1); % k keeps track of Receive index increment per frame.
    % - Set attributes for each frame.
    Receive(k+1).callMediaFunc = 1; % move points before doing ensemble of different angle plane waves
    % acquisitions for 2D
    for j = 1:na
        Receive(j+k).framenum = i;
        Receive(j+k).acqNum = j;
    end
    % PA acquisitions
    for j = (na+1):(na+ne);
        Receive(j+k).framenum = i;
        Receive(j+k).acqNum = j;        % PA acqNums continue after 2D
        Receive(j+k).startDepth = P(2).startDepth;
        Receive(j+k).endDepth = P(2).startDepth + wl4sPer128*ceil(maxAcqLngthPA/wl4sPer128);
        Receive(j+k).TGC = 2;           % TGC(1) is tied to the GUI sliders
        if j==na+1, Receive(j+k).callMediaFunc = 1; end % move points between 2D and PA to see difference in simulation
    end
end

%% Specify TGC Waveform structures.
% - 2D TGC
TGC.CntrlPts = [0,141,275,404,510,603,702,782];
TGC(1).rangeMax = P(1).endDepth;
TGC(1).Waveform = computeTGCWaveform(TGC(1));
% - PA TGC
TGC(2).CntrlPts =[ 0 272 662 662 662 662 662 662]; % TO BE MODIFIED HERE AFTER COLLECTING PA DATA
TGC(2).rangeMax = P(2).endDepth;
TGC(2).Waveform = computeTGCWaveform(TGC(2));

%% Specify Recon structure arrays.
% - We need two Recon structures, one for 2D, one for PA. These will be referenced in the same
%   event, so that they will use the same (most recent) acquisition frame.
Recon = repmat(struct('senscutoff', 0.7, ...
               'pdatanum', 1, ...
               'rcvBufFrame', -1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...
               'RINums', zeros(1,1)), 1, 2);
% - Set Recon values for 2D frame.
Recon(1).RINums(1,1:na) = 1:na;  % na ReconInfos needed for na angles
k = na + 1;
% - Set Recon values for PA ensemble.
Recon(2).pdatanum = 2;
Recon(2).IntBufDest = [2,1];
Recon(2).ImgBufDest = [2,-1];
Recon(2).RINums(1,1:ne) = k:(k+ne-1);   % 'ne' ReconInfos needed for PA ensemble.

%% Define ReconInfo structures.
% - For 2D, we need na ReconInfo structures for na steering angles.
% - For PA, we need ne ReconInfo structures.
ReconInfo = repmat(struct('mode', 'accumIQ', ...    % 4=accumulate IQ data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 1), 1, na + ne);
% - Set specific ReconInfo attributes.
%   - ReconInfos for 2D frame.
ReconInfo(1).mode = 'replaceIQ';          % 3=replace IQ data (expect to use mode 5 on last acquisition)
for j = 1:na
    ReconInfo(j).txnum = j;
    ReconInfo(j).rcvnum = j;
end
if na>1
    ReconInfo(na).mode = 'accumIQ_replaceIntensity';     % 5=Reconstruct IQ data, add values to InterBuffer and compute magnitude, replacing data in ImageBuffer.
else
    ReconInfo(na).mode = 'replaceIntensity';     % 0=replace IQ data, detect, and replace Intensity data in ImageBuffer. (single acquisition)
end

%  - ReconInfos for PA ensemble.
k = na;
for j = 1:ne
    if j==1, ReconInfo(k+j).mode = 'replaceIQ'; end
    ReconInfo(k+j).txnum = na + 1;
    ReconInfo(k+j).rcvnum = na + j;
end
if ne>1
    ReconInfo(na+ne).mode = 'accumIQ_replaceIntensity'; % 5=accum and detect
else
    ReconInfo(na+ne).mode = 'replaceIntensity'; % 0=replace IQ data, detect, and replace Intensity data;  1=Add the new reconstructed intensity data to the data in the ImageBuffer
end

%% Specify Process structure arrays.
cpt = 22;       % define here so we can use in UIControl below
cpers = 80;     % define here so we can use in UIControl below

Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',1.0,...            % pgain is image processing gain
                         'reject',10,...      % reject level
                         'persistMethod','simple',...
                         'persistLevel',50,...
                         'interpMethod','4pt',...
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',40,...
                         'mappingMethod','lowerHalf',...
                         'display',0,...      % don't display image after processing
                         'displayWindow',1};

Process(2).classname = 'Image';
Process(2).method = 'imageDisplay';
Process(2).Parameters = {'imgbufnum',2,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',2,...    % number of PData structure to use
                         'pgain',1.0,...            % pgain is image processing gain
                         'reject',10,...      % reject level
                         'persistMethod','dynamic',...
                         'persistLevel',cpers,...
                         'interpMethod','4pt',...
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',40,...
                         'mappingMethod','upperHalf',...
                         'threshold',cpt,...
                         'display',1,...      % display image after processing
                         'displayWindow',1};


%% Specify SeqControl structure arrays.
% -- Time between 2D flash angle acquisitions
SeqControl(1).command = 'timeToNextAcq';
SeqControl(1).argument = 200;
% -- Change to Profile 2 (PA)
SeqControl(2).command = 'setTPCProfile';
SeqControl(2).condition = 'next';
SeqControl(2).argument = 2;
% -- Time between 2D acquisition and PA ensemble. Set to allow time for profile change.
SeqControl(3).command = 'timeToNextAcq';
SeqControl(3).argument = 7000; % time in usec
% -- PRF for PA ensemble
SeqControl(4).command = 'timeToNextAcq';
SeqControl(4).argument = round(1/(PA_PRF*1e-06)); % (10 msecs for PA_PRF=100 Hz)
% -- Change to Profile 1 (2D)
SeqControl(5).command = 'setTPCProfile';
SeqControl(5).condition = 'next';
SeqControl(5).argument = 1;
% -- Time between PA and next 2D acquisition. Set to allow time for profile change.
SeqControl(6).command = 'timeToNextAcq';
SeqControl(6).argument = 7000; % time in usec
% -- Jump back to start.
SeqControl(7).command = 'jump';
SeqControl(7).argument = 1;
% set receive profile
SeqControl(8).command = 'setRcvProfile';
SeqControl(8).argument = 1;
SeqControl(9).command = 'setRcvProfile';
SeqControl(9).argument = 2;
% output trigger
SeqControl(10).command = 'triggerOut';
% input trigger
SeqControl(11).command = 'triggerIn';
SeqControl(11).condition = 'Trigger_1_Rising'; % Trigger input 1, enable with rising edge
SeqControl(11).argument = 2; % 500 msec timeout delay
    % (Timeout range is 1:255 in 250 msec steps; 0 means timeout disabled)
% noop delay between trigger in and start of acquisition
SeqControl(12).command = 'noop';
SeqControl(12).argument = fix(flash2Qdelay)*5; % noop counts are in 0.2 microsec increments
% sync command
SeqControl(13).command = 'sync';

% - The remainder of the SeqControl structures are defined dynamically in the sequence events.
%   The variable nsc keeps track of the next SeqControl structure index.
nsc = 14;  % next SeqControl number

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 3;

%% Specify Event structure arrays.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    % Acquire 2D frame
    for j = 1:na
        Event(n).info = 'Acquire 2D flash angle';
        Event(n).tx = j;
        Event(n).rcv = (na+ne)*(i-1)+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 1;
        n = n+1;
    end
    Event(n-1).seqControl = [3,9];   % replace last 2D acquisition Event's seqControl (longer TTNA and new RCV profile)

    % Acquire PA ensemble.
    for j = (na+1):(na+ne)
        % Wait for input trigger from flash lamp firing
        Event(n).info = 'Wait for Trigger IN';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 11;
        n = n+1;

        % Pause for optical buildup
        Event(n).info = 'noop and sync';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = [12,13];
        n = n+1;

        % send trigger output at start of every PA acquisition to fire Q-switch
        Event(n).info = 'Acquire PA event';
        Event(n).tx = na+1;
        Event(n).rcv = (na+ne)*(i-1)+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 10;

        n = n+1;
    end
    Event(n-1).seqControl = [10,6,8]; % replace last PA acquisition Event's seqControl with longer TTNA and RCV profile change

    Event(n).info = 'Transfer Data';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = nsc;
    n = n+1;
    SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer (each needs a different value of nsc)
      nsc = nsc+1;

    Event(n).info = 'recons and 2D process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = [1,2];
    Event(n).process = 1;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = 'PA image display';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 2;
    Event(n).seqControl = 0;
    if floor(i/frameRateFactor) == i/frameRateFactor     % Exit to Matlab only every 3rd frame to prevent slowdown
        Event(n).seqControl = nsc;
        SeqControl(nsc).command = 'returnToMatlab';
        nsc = nsc+1;
    end
    n = n+1;
end

Event(n).info = 'Jump back';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 7;


%% User specified UI Control Elements
% - Sensitivity Cutoff
UI(1).Control =  {'UserB7','Style','VsSlider','Label','Sens. Cutoff',...
                  'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                  'SliderStep',[0.025,0.1],'ValueFormat','%1.3f'};
UI(1).Callback = text2cell('%-UI#1Callback');

% - Color Priority Threshold Slider
UI(2).Control = {'UserB2','Style','VsSlider','Label','Color Priority','SliderMinMaxVal',[0,255,cpt],...
                 'SliderStep',[1/255,0.1],'ValueFormat','%3.0f'};
UI(2).Callback = text2cell('%-UI#2Callback');

% - Color Persistence Slider
UI(3).Control = {'UserB1','Style','VsSlider','Label','Color Persistence','SliderMinMaxVal',[0,100,cpers],...
                 'SliderStep',[1/100,0.1],'ValueFormat','%3.0f'};
UI(3).Callback = text2cell('%-UI#3Callback');


%% Save all the structures to a .mat file, and run VSX automatically
filename = ('L11-5vFlashPA');   % define variable 'filename' to permit VSX to skip user query for matfile
save (['MatFiles/',filename])                 % save the structures to a matfile
% VSX                             % invoke VSX automatically when running this Setup script

return


%% **** Callback routines to be encoded by text2cell function. ****
% ---------------------------------------------------------------------
%-UI#1Callback - Sensitivity cutoff change
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
%-UI#1Callback


%-UI#2Callback - Color Threshold change
    % Set the value in the Process structure for use in cineloop playback.
    Process = evalin('base','Process');
    for k = 1:2:length(Process(2).Parameters)
        if strcmp(Process(2).Parameters{k},'threshold'), Process(2).Parameters{k+1} = UIValue; end
    end
    assignin('base','Process',Process);
    % Set Control.Command to set Image.threshold.
    Control = evalin('base','Control');
    Control.Command = 'set&Run';
    Control.Parameters = {'Process',2,'threshold',UIValue};
    assignin('base','Control', Control);
%-UI#2Callback

%-UI#3Callback - Color Persistence change
    % Set the value in the Process structure for use in cineloop playback.
    Process = evalin('base','Process');
    for k = 1:2:length(Process(2).Parameters)
        if strcmp(Process(2).Parameters{k},'persistLevel'), Process(2).Parameters{k+1} = UIValue; end
    end
    assignin('base','Process',Process);
    % Set Control.Command to set Image.persistLevel.
    Control = evalin('base','Control');
    Control.Command = 'set&Run';
    Control.Parameters = {'Process',2,'persistLevel',UIValue};
    assignin('base','Control', Control);
%-UI#3Callback
