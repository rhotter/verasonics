% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use.
%
% File name SetUpL7_4ShearWaveImaging.m:
% Generate .mat Sequence Object file for L7-4 Linear flash transmit with Shearwave visualization.
% An example for Changing startevent on conditional coding
%
% All 128 transmit channels are used, with receive on the 64 central elements.
% Update:
% Peter K finished the pushtest using Ron's IQ processing function
% YT, 08/29/2014: Focus adjusment and movie save is added
% YT, 12/02/2014: SWI on/off is added for customer to have regular imaging
% YT, 02/10/2015: add dragable ROI ability to move ROI around
% YT, 02/18/2015: modify the first set of events to regular flash imaging
% YT, 11/23/2015: make it compatible with 3.0
% April 2019 VTS-1142 IQBuffer replaced with separate I and Q Buffers
% YT, 09/20/2019: make one script working for multiple configurations
% 05/06/2020 - Update to SW 4.3 format for new UIControls and function definitions (VTS 1691). 
%   More info:(.../Example_Scripts/Vantage_Features/New UI Scheme/SetUpL11_5vFlash_NewUI)

clear all

%% System parameters
filename = ('L7-4ShearWave');

HIFU.port = 'COM4'; % not required for a non-HIFU system

na          = 50;      % Set na = number of detect acquisitions.
SWIFrames   = 10;
BmodeFrames = 20;

powermax   = 250;      % scaling of the display function
pushCycle  = 1000;
maxVoltage = 65;

% Define ROI for IQ processing, can be changed in GUI
SWIFocusX = 0;
SWIFocusZ = 50;
SWIROI.width = 80;
SWIROI.height = 50;  % [width and height in wavelength]

%% Define system parameters.
Resource.Parameters.numTransmit = 128;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;    % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

SysConfig = hwConfigCheck(1);

V64LE = SysConfig.LEsys;

% Specify Trans structure array.
Trans.name = 'L7-4';
Trans.units = 'mm';
Trans = computeTrans(Trans);  % L7-4 transducer is 'known' transducer so we can use computeTrans.
Trans.maxHighVoltage = maxVoltage;  % set maximum high voltage limit for pulser supply.

%% Imaging parameters
P.startDepth = 5;   % Acquisition depth in wavelengths
P.endDepth = 160;   % This should preferrably be a multiple of 128 samples.

% Specify PData(1) structure array for Bmode Imaging
PData(1).PDelta = [Trans.spacing,0,0.5];
PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3)); % startDepth, endDepth and pdelta set PData(1).Size.
PData(1).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(3) = 1;      % single image page
PData(1).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of upper lft crnr.

% Specify PData(2) structure array for Shearwave visulization
PData(2) = PData(1);
PData(2).PDelta = [0.5,0,0.25];
PData(2).Size(1) = ceil((P.endDepth-P.startDepth)/PData(2).PDelta(3)); % startDepth, endDepth and pdelta set PData(1).Size.
PData(2).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(2).PDelta(1));
PData(2).Region.Shape = struct(...
    'Name','Rectangle',...
    'Position',[SWIFocusX,0,SWIFocusZ-SWIROI.height/2],...
    'width', SWIROI.width,...
    'height', SWIROI.height);
PData(2).Region = computeRegions(PData(2));

% PData(2).Size(1) = ceil(SWIROI.height/PData(2).PDelta(3)); % startDepth, endDepth and pdelta set PData(2).Size.
% PData(2).Size(2) = ceil(SWIROI.width/PData(2).PDelta(1));
% PData(2).Size(3) = 1;      % single image page
% PData(2).Origin = [(SWIFocusX-SWIROI.width)/2,0,SWIFocusZ-SWIROI.height/2]; % x,y,z of upper lft crnr in wavelength

% Specify Media object. 'pt1.m' script defines array of point targets.
pt1;
Media.function = 'movePoints';

%% Specify Resources.
% RcvBuffer for all raw data
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 4096;
Resource.RcvBuffer(1).colsPerFrame = 128;
Resource.RcvBuffer(1).numFrames = BmodeFrames;

% RcvBuffer for all raw data
Resource.RcvBuffer(2).datatype = 'int16';
Resource.RcvBuffer(2).rowsPerFrame = na*4096;
Resource.RcvBuffer(2).colsPerFrame = 128;
Resource.RcvBuffer(2).numFrames = SWIFrames;

% InterBuffer for Bmode (not required)
Resource.InterBuffer(1).numFrames = 1;

% InterBuffer for SWI visualizaion, colsPerFrame must be larger enough for
% PData(2) update
Resource.InterBuffer(2).numFrames = 1;
Resource.InterBuffer(2).pagesPerFrame = na;

% ImageBuffer for reference Bmode image
Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).numFrames = BmodeFrames;

Resource.DisplayWindow(1).Title = filename;
Resource.DisplayWindow(1).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...
    DwWidth, DwHeight];  % lower left corner position
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)]; % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Colormap = gray(256);
Resource.DisplayWindow(1).numFrames = BmodeFrames;
Resource.DisplayWindow(1).Type = 'Matlab';

%% Transmit parameters
% Specify Transmit waveform structure.
% - detect waveform
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,0.67,2,1];

% - Push waveform.
TW(2).type = 'parametric';
TW(2).Parameters = [Trans.frequency,1,pushCycle*2,1];

% Configure for use with HIFU system if required
if SysConfig.p5req == 2
    Resource.HIFU.externalHifuPwr = 1;
    Resource.HIFU.extPwrComPortID = HIFU.port;
    Resource.HIFU.psType = 'QPX600DP'; % set to 'QPX600DP' or 'XG40-38' to match supply being used
end

% Set TPC profile 5 high voltage limit.
TPC(5).maxHighVoltage = Trans.maxHighVoltage;

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
    'Origin', [0.0,0.0,0.0], ...
    'focus', 0.0, ...
    'Steer', [0.0,0.0], ...
    'Apod', ones(1,Trans.numelements), ...
    'Delay', zeros(1,Trans.numelements)), 1, 2);

% - Set event specific TX attributes for push.
TX(2).waveform = 2;
TX(2).focus = SWIFocusZ;       % wavelength, can be changed in the GUI
TX(2).focusX = SWIFocusX;      % can be changed in the GUI
TX(2).pushElements = 32;       % can be changed in the GUI

% Apod = 1 for pushing elements
TX(2).Apod = zeros(1,Trans.numelements);
TX(2).Apod(64-TX(2).pushElements/2 +1 : 64+TX(2).pushElements/2) = 1;
TX(2).Delay = computeTXDelays(TX(2));

%% Receive parameters
% Specify Receive structure arrays.
% -- Receive will be changed after SWI on.
% -- Compute the maximum receive path length, using the law of cosines.
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
wlsPer128 = 128/(4*2); % wavelengths in 128 samples for 4 samplesPerWave

% Total (na+1)*RcvBuffer.numFrames  (1 regular flash and na detections)
if V64LE
    RcvApod = [zeros(1,32),ones(1,64),zeros(1,32)];
else
    RcvApod = ones(1,Trans.numelements);
end
Receive = repmat(struct('Apod', RcvApod, ...
    'startDepth', P.startDepth, ...
    'endDepth', P.startDepth + wlsPer128*ceil(maxAcqLength/wlsPer128), ...
    'TGC', 1, ...
    'bufnum', 1, ...
    'framenum', 1, ...
    'acqNum', 1, ...
    'sampleMode', 'NS200BW', ...
    'mode', 0, ...
    'callMediaFunc', 0), 1, BmodeFrames+na*SWIFrames);

% % Receive(1) to Receive(20) are regular flash imaging
for i = 1:BmodeFrames
    % -- Acquisition for full frame.
    Receive(i).callMediaFunc = 1;  % make media move per frame
    Receive(i).framenum = i;
end

d = i; % detection starts from i+1, here is 21;

% - Set event specific Receive attributes for each frame.
for i = 1:SWIFrames
    for j = 1:na  % na acquisitions per frame
        Receive(na*(i-1)+j+d).callMediaFunc = 1;  % make media move per frame
        Receive(na*(i-1)+j+d).bufnum = 2;
        Receive(na*(i-1)+j+d).framenum = i;
        Receive(na*(i-1)+j+d).acqNum = j;
    end
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [500,590,650,710,770,800,850,950];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

%% Reconstruction parameters
% Specify Recon structure arrays.
% - We need a Recon structure for the 2D image which will be used for each frame.
% Recon(1) is used in regular flash imaging
% Recon(2) is used in Bmode imaging after SWI on
% - 10 IQData frames will be stored after Recon(3) to Recon(12)

Recon = repmat(struct(...
    'senscutoff', 0.6, ...
    'pdatanum', 1, ...
    'rcvBufFrame',-1, ...
    'IntBufDest', [2,1], ...
    'ImgBufDest', [0,0], ...
    'RINums',1),1,3);

Recon(1).IntBufDest = [1,1];
Recon(1).ImgBufDest = [1,-1];
Recon(1).RINums = 1;

Recon(2).IntBufDest = [0,0];
Recon(2).ImgBufDest = [1,-1];
Recon(2).RINums = 2;

Recon(3).pdatanum = 2;
Recon(3).IntBufDest = [2,1];
Recon(3).ImgBufDest = [0,0];
Recon(3).RINums = 3:na+2;

% Define ReconInfo structures.
% - ReconInfo for 2D frame.
ReconInfo(1) = struct('mode','replaceIntensity', ...  % intensity output.
    'txnum',1, ...
    'rcvnum',1, ...
    'regionnum',1);

ReconInfo(2) = struct('mode','replaceIntensity', ...  % intensity output.
    'txnum',1, ...
    'rcvnum',d+1, ...                % the first Acq will be shown in the display window
    'regionnum',1);

k = 2; % k keeps track of index of last ReconInfo defined
% We need na ReconInfo structures for IQ reconstructions.

ReconInfo((k+1):(k+na)) = repmat(struct('mode', 'replaceIQ', ... % IQ output
    'txnum', 1, ...
    'rcvnum', 1, ...
    'regionnum', 1), 1, na);

% - Set specific ReconInfo attributes.
for j = 1:na  % For each row in the column
    ReconInfo(k+j).txnum = 1;
    ReconInfo(k+j).rcvnum = j+d;
    ReconInfo(k+j).pagenum = j;
end


%% Process parameters
% Specify Process structure array. (1) is used for B-mode imaging
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

% EF1 is external function for UI control
Process(2).classname = 'External';
Process(2).method = 'UIControl';
Process(2).Parameters = {'srcbuffer','none'};

% EF2 is external function for shearwave visualization
Process(3).classname = 'External';
Process(3).method = 'processIQ';
Process(3).Parameters = {'srcbuffer','inter',... % name of buffer to process.
    'srcbufnum',2,...
    'srcframenum',1,...
    'dstbuffer','none'};

%% SeqControl and Events for shearwave generation

% - Change to Profile 1 (low power)
SeqControl(1).command = 'setTPCProfile';
SeqControl(1).condition = 'immediate';
SeqControl(1).argument = 1;
% - Change to Profile 5 (high power)
SeqControl(2).command = 'setTPCProfile';
SeqControl(2).condition = 'immediate';
SeqControl(2).argument = 5;
% - Noop to allow time for charging external cap.
SeqControl(3).command = 'noop';
SeqControl(3).argument = 500000; % wait 100 msec.

% - time between regular flash imaging
SeqControl(4).command = 'timeToNextAcq';
SeqControl(4).argument = 10000;               % 10 ms
% - time between push and detect acquisitions
SeqControl(5).command = 'timeToNextAcq';
SeqControl(5).argument = 500;               % 500usec
afterpush = 5;
% - time between detect acquisitions
SeqControl(6).command = 'timeToNextAcq';
SeqControl(6).argument = 100;               % 100usec
PRF=6;
% - time between frames
SeqControl(7).command = 'timeToNextEB';    % set time between extended bursts
SeqControl(7).argument = 200000;            % 200000usec = 200msec (~ 5 fps)
TTNEB=7;

% - Return to Matlab
SeqControl(8).command = 'returnToMatlab';
% - Trigger out
SeqControl(9).command = 'triggerOut';

% - Jump back to start will be defined in the event due to the conditional
% event coding

nsc = length(SeqControl)+1;

% Specify Event structure arrays.
n = 1;

Event(n).info = 'ext func for UI control';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 2;
Event(n).seqControl = 0;
n = n+1;

%% Regular flash imaging starts fron event(nStartFlash)
nStartFlash = n;

% Switch to TPC profile 1 (low power) for flash imaging
Event(n).info = 'Switch to profile 1.';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 1;
n = n+1;

nStartAcqFlash = n;

for i = 1:BmodeFrames
    Event(n).info = 'trigger out';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = 9;
    n = n+1;

    Event(n).info = 'Full aperture.';
    Event(n).tx = 1;
    Event(n).rcv = i;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = [4,nsc];
    SeqControl(nsc).command = 'transferToHost';
    nsc = nsc + 1;
    n = n+1;

    Event(n).info = 'recon and process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 1;
    Event(n).seqControl = 8;
    n = n+1;
end

Event(n).info = 'Jump back';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = nsc;
SeqControl(nsc).command = 'jump';
SeqControl(nsc).argument = nStartAcqFlash;

n = n+1;
nsc = nsc+1;

%% Shear Wave Imaging starts from event(nStartPush)
nStartPush = n;

% Switch to TPC profile 5 (high power) and allow time for charging ext. cap.
Event(n).info = 'Switch to profile 5.';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 2;
n = n+1;

Event(n).info = 'noop for charging ext. cap.';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 3;
n = n+1;

nStartAcqPush = n;

for i = 1:SWIFrames
    Event(n).info = 'trigger out';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = 9;
    n = n+1;

    % Push transmit
    Event(n).info = 'Push transmit';
    Event(n).tx = 2;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = [afterpush, TTNEB];
    n = n+1;

    for j = 1:na                      % Acquire frame
        Event(n).info = 'Acquire data';
        Event(n).tx = 1;
        Event(n).rcv = na*(i-1)+j+d;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = PRF;
        n = n+1;
    end

    Event(n-1).seqControl = 0; % do not want a TTNA here, since the next transmit is a EB

    Event(n).info = 'transfer data to Host';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = nsc;
    SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
    nsc = nsc+1;
    n = n+1;

    Event(n).info = 'recon and process for SWI';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = [2,3];
    Event(n).process = 1;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = 'ext func to process IQ';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 3;
    Event(n).seqControl = 8;
    n = n+1;

end

Event(n).info = 'Jump back';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = nsc;
SeqControl(nsc).command = 'jump';
SeqControl(nsc).argument = nStartAcqPush;


%% User specified UI Control Elements
import vsv.seq.uicontrol.VsSliderControl; 
import vsv.seq.uicontrol.VsButtonControl

% Define UIPos, which contains the default GUI positions - three columns of 10 controls. The x,y
%    locations increment up columns, with each column being a separate page. The origin
%    specified by UIPos is the lower left corner of a virtual box that encloses the control.
UIPos = zeros(10,2,3);
UIPos(:,1,1) = 0.0625;
UIPos(:,1,2) = 0.375;
UIPos(:,1,3) = 0.6875;
UIPos(:,2,1) = 0.0:0.1:0.9;
UIPos(:,2,2) = 0.0:0.1:0.9;
UIPos(:,2,3) = 0.0:0.1:0.9;

% Define slider group offsets and sizes. All units are normalized.
SG = struct('TO',[0.0,0.0975],...   % title offset
    'TS',[0.25,0.025],...   % title size
    'TF',0.8,...            % title font size
    'SO',[0.0,0.06],...     % slider offset
    'SS',[0.25,0.031],...   % slider size
    'EO',[0.075,0.031],...   % edit box offset
    'ES',[0.11,0.031]);     % edit box size

% - Sensitivity Cutoff
UI(1).Control = VsSliderControl('LocationCode','UserB7','Label','Sens. Cutoff',...
                  'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                  'SliderStep',[0.025,0.1],'ValueFormat','%1.3f', ...
                  'Callback', @SensCutoffCallback );  


% - Range Change
MinMaxVal = [64,300,P.endDepth]; % default unit is wavelength
AxesUnit = 'wls';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
        AxesUnit = 'mm';
        MinMaxVal = MinMaxVal * (Resource.Parameters.speedOfSound/1000/Trans.frequency);
    end
end
UI(2).Control = VsSliderControl('LocationCode','UserA1','Label',['Range (',AxesUnit,')'],...
                 'SliderMinMaxVal',MinMaxVal,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f', ...
                 'Callback', @RangeChangeCallback);  


% - Focus Adjustment, off: move Focus = 1, on: move ROI = 2
UI(3).Control = vsv.seq.uicontrol.VsButtonGroupControl('LocationCode','UserB5','Title','move Focus / ROI',...
                 'PossibleCases',   {'move Focus','move ROI'},...
                 'Callback', @MoveFocusMoveROICallback);  



% - Push Elements Adjustment
TX(2).oldElements = TX(2).pushElements;
UI(4).Control = VsSliderControl('LocationCode','UserB4','Label','Push Elements',...
                  'SliderMinMaxVal',[20,100,TX(2).pushElements],...
                  'SliderStep',[1/80,5/80],'ValueFormat','%3.0f',...
                  'Callback', @NumPushElemCallback );  


% - IQ Caxis Adjustment
UI(5).Control = VsSliderControl('LocationCode','UserC3','Label','IQ Caxis',...
                  'SliderMinMaxVal',[100,600,powermax],...
                  'SliderStep',[10/500,50/500],'ValueFormat','%3.0f');  
UI(5).Callback = @(src,evt, UIValue)assignin('base', 'powermax', UIValue);


% - IQ loop fps adjustment, only works at freeze and replay status
replay = 'off';
loopfps = 15;
UI(6).Control = VsSliderControl('LocationCode','UserC2','Label','IQ loop fps',...
                  'SliderMinMaxVal',[1,31,loopfps],...
                  'SliderStep',[1/30,5/30],'ValueFormat','%3.0f',...
                  'Callback', @IQLoopFpsCallback );  


% The follow UIs are not using User##
Pos = UIPos(2,:,3);

% text for IQ replay, no callback
UI(7).Control = {'Style','text',...
    'String','IQ replay',...
    'Units','normalized',...
    'Position',[Pos+SG.TO,SG.TS],...
    'FontUnits','normalized',...
    'FontSize',0.8,...
    'FontWeight','bold'};

% Replay button
UI(8).Control = {'Style','pushbutton',...
    'String','Replay',...
    'Units','normalized',...
    'Position',[Pos+[0 0.05],0.12,0.04],...
    'FontUnits','normalized',...
    'FontSize',0.5,...
    'Callback', @replayIQ};

% Stop button, only visible after clicking "replay"
UI(9).Control = {'Style','pushbutton',...
    'String','Stop',...
    'Visible','off',...
    'Units','normalized',...
    'Position',[Pos+[0 0.05],0.12,0.04],...
    'FontUnits','normalized',...
    'FontSize',0.5};
UI(9).Callback=@(src,evt, UIValue)assignin('base', 'replay', 'off');


% Save shavewave imaging.avi with desired fps in the PC
UI(10).Control = {'Style','pushbutton',...
    'String','Save',...
    'Units','normalized',...
    'Position',[Pos+[0.13 0.05],0.12,0.04],...
    'FontUnits','normalized',...
    'FontSize',0.5,...
    'Callback',@saveIQ};

% ROI adjustment
UI(11).Control = VsSliderControl('LocationCode','UserB2','Label','ROI Width',...
                  'SliderMinMaxVal',[20,100,SWIROI.width],...
                  'SliderStep',[2/40,10/40],'ValueFormat','%3.0f',...
                  'Callback', @ROIWidthCallback );  


UI(12).Control = VsSliderControl('LocationCode','UserB1','Label','ROI Height',...
                  'SliderMinMaxVal',[20,120,SWIROI.height],...
                  'SliderStep',[2/40,10/40],'ValueFormat','%3.0f',...
                  'Callback', @ROIHeigthCallback );  


% SWI switch button
UI(13).Control  = VsButtonControl( 'LocationCode',   'UserB6', 'Label','SWI off');

UI(14).Control  = VsButtonControl( 'LocationCode',   'UserB6', 'Label','SWI on');

UI(15).Control = VsSliderControl('LocationCode','UserB3','Label','Push cycles',...
                  'SliderMinMaxVal',[50,2000,pushCycle],...
                  'SliderStep',[10/1500,50/1500],'ValueFormat','%3.0f',...
                  'Callback', @NumPushCyclesCallback );  


EF(1).Function = vsv.seq.function.ExFunctionDef('UIControl', @UIControl);
EF(2).Function = vsv.seq.function.ExFunctionDef('processIQ', @processIQ);

% ROI position
ROIpos = [TX(2).focusX-SWIROI.width/2,TX(2).focus-SWIROI.height/2,SWIROI.width,SWIROI.height];

% ElementLocation is used for focus correction
Location.Element = (0.5 + PData(1).Origin(1)):Trans.spacing/2:(PData(1).Origin(1)+DwWidth*Resource.DisplayWindow.pdelta - 0.5);
Location.PosForOddPush = Location.Element(1:2:end); % center of each element
Location.PosForEvenPush = Location.Element(2:2:end); % 127 spaces among 128 elements

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 1;

% Save all the structures to a .mat file.
save(['MatFiles/',filename]);

% Uncomment the following line to automatically run VSX every time you run
% this SetUp script (note that if VSX finds the variable 'filename' in the
% Matlab workspace, it will load and run that file instead of prompting the
% user for the file to be used):

% VSX

return

%% **** Callback routines used by UIControls (UI) ****

function SensCutoffCallback(~,~,UIValue)
%Sensitivity cutoff change
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
%Range change
    simMode = evalin('base','Resource.Parameters.simulateMode');
    % No range change if in simulate mode 2.
    if simMode == 2
        set(hObject,'Value',evalin('base','P.endDepth'));
        return
    end
    P = evalin('base','P');
    P.endDepth = UIValue;
    assignin('base','P',P);
    evalin('base','PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));');
    evalin('base','PData(1).Region = computeRegions(PData(1));');
    evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);');
    Receive = evalin('base', 'Receive');
    Trans = evalin('base', 'Trans');
    maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
    for i = 1:size(Receive,2)
        Receive(i).P.endDepth = maxAcqLength;
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

function MoveFocusMoveROICallback(~,~,UIState)
%move Focus or move ROI
    Resource = evalin('base','Resource');
    ROIHandle = evalin('base','ROIHandle');
    fh = Resource.DisplayWindow(1).figureHandle;

    if UIState == 1
        assignin('base','moveTarget','Focus');
        set(ROIHandle,'ButtonDownFcn',[]);
        set(fh,'WindowButtonDownFcn',@moveFocus);
        set(fh,'WindowButtonMotionFcn',[]);
    else
        assignin('base','moveTarget','SWIROI');
        set(fh,'WindowButtonDownFcn',[]);
        set(fh,'WindowButtonMotionFcn',{@moveMouse,Resource.DisplayWindow(1).imageHandle});
        set(ROIHandle,'ButtonDownFcn',@moveROI);
    end
end

function NumPushElemCallback(~,~,UIValue)
%PushElements adjustment
    TX = evalin('base','TX');
    UI = evalin('base','UI');
    freeze = evalin('base','freeze');

    assignin('base','moveTarget','Focus');
    pushElements = round(UIValue); %(get(UI(4).handle(2),'Value'));

    % If focusAdj is off, or at freeze, or Shearwave fig is closed, change value back to the old one
    if freeze == 1
        set(UI(4).handle(2),'Value',TX(2).pushElements);
        set(UI(4).handle(3),'String',TX(2).pushElements);
    else% on, focussCorrection and TXupdate
        TX(2).oldElements = TX(2).pushElements;
        TX(2).pushElements = pushElements;
        assignin('base','TX',TX);

        msg = ['Change push elements to ',num2str(pushElements),' elements...'];

        focusCorrection();
        TXupdate(msg,'TX');
    end
end

function IQLoopFpsCallback(~,~,UIValue)
%fps adjustment for replay shearwave imaging
    UI = evalin('base','UI');
    loopfps = evalin('base','loopfps');
    freeze = evalin('base','freeze');

    if freeze
        assignin('base','loopfps', UIValue);
    else
        set(UI(6).handle(2),'Value',loopfps);
        set(UI(6).handle(3),'String',loopfps);
    end
end

function replayIQ(varargin)
    freeze = evalin('base','freeze');
    % figClose = evalin('base','figClose');

    if freeze
        replay = 'on';
        assignin('base','replay','on');

        % need copyBuffers to access IQData
        Control.Command = 'copyBuffers';
        runAcq(Control); % NOTE:  If runAcq() has an error, it reports it then exits MATLAB.

        UI = evalin('base','UI');
        TX = evalin('base','TX');
        ROIpos = evalin('base','ROIpos');
        MovieData = evalin('base','MovieData');
        SWIfigHandle = evalin('base','SWIfigHandle');

        set(UI(8).handle,'Visible','off');
        set(UI(9).handle,'Visible','on');

        disp('replay shearwave imaging....');

        x = round(TX(2).focusX);
        z = TX(2).focus;

        IQBuffer = complex(IData{2},QData{2});

        Depth = size(IQBuffer,1);
        Width = size(IQBuffer,2);

        na = evalin('base','na');

        IQMovie(na-1) = struct('cdata',[],'colormap',[]);
        SWIimageHandle = evalin('base','SWIimageHandle');

        evalin('base','clear IQMovie');

        while strcmp(replay,'on')

            for i = 1:na-1
                tic
                replay = evalin('base','replay');
                freeze = evalin('base','freeze');
                loopfps = evalin('base','loopfps');
                powermax = evalin('base', 'powermax');

                % if GUI is closed, return
                if evalin('base','vsExit == 1')
                    return
                end

                % unfreeze will stop replay
                if ~freeze || strcmp(replay,'off') == 1
                    set(UI(8).handle,'Visible','on');
                    set(UI(9).handle,'Visible','off');
                    assignin('base','replay','off');
                    return
                end

                set(SWIimageHandle,'CData',squeeze(MovieData(:,:,i)));
                drawnow
                caxis(get(SWIfigHandle,'currentAxes'),[50,powermax]);

                frameData = getframe(SWIfigHandle);
                % padding to a multiple of two for MPEG-4 format
                [row,col,~]=size(frameData.cdata);
                if ~isequal(mod(row,2),0)
                    frameData.cdata(end,:,:) = [];
                end
                if ~isequal(mod(col,2),0)
                    frameData.cdata(:,end,:) = [];
                end
                IQMovie(i) = frameData;
                pause(1/loopfps-toc);
            end

            assignin('base','IQMovie',IQMovie);
        end
    end

    end

function saveIQ(varargin)
    UI = evalin('base','UI');
    freeze = evalin('base','freeze');

    if evalin('base','exist(''IQMovie'',''var'')')
        if freeze
            assignin('base','replay','off');
            set(UI(8).handle,'Visible','on');
            set(UI(9).handle,'Visible','off');

            IQMovie = evalin('base','IQMovie');
            loopfps = evalin('base','loopfps');
            [fn,pn,filterindex] = uiputfile('*.mp4','Save Shearwave movie as');
            if ~isequal(fn,0) % fn will be zero if user hits cancel
                fn = strrep(fullfile(pn,fn), '''', '''''');

                vidObj = VideoWriter(fn,'MPEG-4');
                vidObj.Quality = 100;
                vidObj.FrameRate = loopfps;
                open(vidObj);
                writeVideo(vidObj,IQMovie);
                close(vidObj);

                fprintf('The shearwave movie has been saved at %s \n',fn);
            else
                disp('The shearwave movie is not saved.');
            end
        end
    else
        msgbox('replay is not finished!');
    end
end

function ROIWidthCallback(~,~,UIValue)
%ROI Width change
    PData = evalin('base','PData');
    SWIROI = evalin('base','SWIROI');
    Location = evalin('base','Location');
    ROIHandle = evalin('base','ROIHandle');
    SWIfigHandle = evalin('base','SWIfigHandle');

    ROIpos = get(ROIHandle,'Position');

    newX = ROIpos(1) - (UIValue-ROIpos(3))/2;
    if newX < Location.Element(1)-0.5, newX = Location.Element(1)-0.5;
    elseif newX > Location.Element(end)+0.5-UIValue, newX = Location.Element(end)+0.5-UIValue; end

    ROIpos(1) = newX;
    ROIpos(3) = UIValue;
    SWIROI.width = UIValue;
    assignin('base','ROIpos',ROIpos);
    assignin('base','SWIROI',SWIROI);

    PData(2).Region(1).Shape.width = SWIROI.width;
    PData(2).Region(1) = computeRegions(PData(2));

    assignin('base','PData',PData);

    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'PData','Recon'};
    assignin('base','Control', Control);
    assignin('base','newIQ',1);
    assignin('base','moveTarget','SWIROI');

    pos = get(SWIfigHandle,'Position');
    set(SWIfigHandle,'Position',[pos(1:2),[SWIROI.width+3,SWIROI.height]*7]);

end

function ROIHeigthCallback(~,~,UIValue)
%ROI Height change
    P = evalin('base','P');
    PData = evalin('base','PData');
    SWIROI = evalin('base','SWIROI');
    ROIHandle = evalin('base','ROIHandle');
    SWIfigHandle = evalin('base','SWIfigHandle');

    ROIpos = get(ROIHandle,'Position');

    newY = ROIpos(2) - (UIValue-ROIpos(4))/2;
    if newY < P.startDepth, newY = P.startDepth;
    elseif newY > P.endDepth - UIValue, newY = P.endDepth - UIValue; end

    ROIpos(2) = newY;
    ROIpos(4) = UIValue;
    SWIROI.height = UIValue;
    assignin('base','ROIpos',ROIpos);
    assignin('base','SWIROI',SWIROI);

    PData(2).Region(1).Shape.Position(3) = newY;
    PData(2).Region(1).Shape.height = SWIROI.height;
    PData(2).Region(1) = computeRegions(PData(2));

    assignin('base','PData',PData);

    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'PData','Recon'};
    assignin('base','Control', Control);
    assignin('base','newIQ',1);
    assignin('base','moveTarget','SWIROI');

    pos = get(SWIfigHandle,'Position');
    set(SWIfigHandle,'Position',[pos(1:2),[SWIROI.width+3,SWIROI.height]*7]);
end

function NumPushCyclesCallback(~,~,UIValue)
%pushCycle change
    TW = evalin('base', 'TW');
    TX = evalin('base', 'TX');
    TW(2).oldCycles = evalin('base','pushCycle');

    pushCycle = round(UIValue);

    % TW(2) is used for push
    if strcmp(TW(2).type,'parametric') || strcmp(TW(2).type,'envelop')
        TW(2).Parameters(3) = pushCycle*2;
        evalin('base',['pushCycle = ',num2str(pushCycle),';']);
    else
        set(UI(15).handle(2),'Value',TW.oldCycles);
        set(UI(15).handle(3),'String',num2str(TW.oldCycles));
    end

    [~,~,~,~,TW(2)] = computeTWWaveform(TW(2));
    TX(2).Bdur = TW(2).Bdur;

    msg = ['Change push cycles to ',num2str(pushCycle,'%4.0f'),' cycles...'];

    assignin('base','TW',TW);
    assignin('base','TX',TX);

    TXupdate(msg,'TW');

end


%% **** Callback routines used by External function definition (EF) ****

function UIControl(varargin)

    UI = evalin('base','UI');
    f = evalin('base','f');

    for i = 3:15
        set(UI(i).handle,'Visible','off');
    end

    hv2Handle(1) = findobj('String','High Voltage P5');
    hv2Handle(2) = findobj('tag','hv2Sldr');
    hv2Handle(3) = findobj('tag','hv2Value');
    hv2Handle(4) = findobj('tag','hv2Actual');
    set(hv2Handle,'Visible','off');

    % modify font of buttongroup
    % f = findobj('tag','UI');
    bkgrnd = get(f,'Color');
    set(UI(3).handle(1),'FontSize',0.2,'BackgroundColor',bkgrnd);
    if ispc
        pos = [0.58 0.08]; % for Windows
    else
        pos = [0.52 0.05]; % for Mac
    end

    set(UI(3).handle(2),'FontSize',0.9,'Position',[0.05,pos(1),0.9,0.34]);
    set(UI(3).handle(3),'FontSize',0.9,'Position',[0.05,pos(2),0.9,0.34]);

    set(UI(13).handle,'Callback',@SWIoff);
    set(UI(14).handle,'Callback',@SWIon);
    set(UI(14).handle,'Visible','on');

end

function processIQ(IBuffer, QBuffer)
    IQBuffer = complex(IBuffer, QBuffer); % VTS-1142 separate I and Q buffers
    %processIQFunction: Computes power estimates from IQData
    %		Im = I(k) * Q(k+1) - I(k+1) * Q(k)
    %		Re = I(k) * I(k+1) + Q(k) * Q(k+1)
    %		Power = sqrt(Im*Im + Re*Re);

    persistent recHandle markHandle myHandle recTrans

    SWIROI = evalin('base','SWIROI');
    ROIpos = evalin('base','ROIpos');
    SWIfigHandle = evalin('base','SWIfigHandle');
    bmodeFigHandle = evalin('base','Resource.DisplayWindow(1).figureHandle');

    na = evalin('base','na');
    TX = evalin('base','TX');
    PData = evalin('base','PData');
    newIQ = evalin('base','newIQ');
    Trans = evalin('base','Trans');
    powermax = evalin('base', 'powermax');
    moveTarget = evalin('base','moveTarget');

    ROILA = (PData(2).Region.PixelsLA+1);

    Depth = floor(SWIROI.height/PData(2).PDelta(3));%size(IQBuffer,1);%PData(2).Size(1);
    Width = length(ROILA)/Depth;

    if ~isequal(mod(Width,1),0)
        Width = floor(SWIROI.width/PData(2).PDelta(1));%size(IQBuffer,2);%PData(2).Size(2);
        Depth = length(ROILA)/Width;
    end

    MovieData = zeros(Depth,Width,na-1);

    x = TX(2).focusX;
    z = TX(2).focus;

    % Focus mark, SWIROI,  and Transducer rect on bmode figure
    if ishandle(bmodeFigHandle)
        if  isempty(recHandle) || ~ishandle(recHandle)
            figure(bmodeFigHandle), hold on,
            recHandle = rectangle('Position',[x-SWIROI.width/2,z-SWIROI.height/2,SWIROI.width,SWIROI.height],'EdgeColor','w','LineWidth',2);
            recTrans  = rectangle('Position',[x-TX(2).pushElements*Trans.spacing/2,4,TX(2).pushElements*Trans.spacing,4],'EdgeColor','r','FaceColor','r');
            markHandle = plot(x,z,'xr','MarkerFaceColor','r','MarkerSize',8,'Linewidth',2);hold off;
            assignin('base','TransHandle',recTrans);
            assignin('base','ROIHandle',recHandle);
            assignin('base','markHandle',markHandle);
        else
            switch moveTarget
                case 'Focus'
                    set(markHandle,'XData',x);
                    set(markHandle,'YData',z);
                    set(recTrans,'Position',[x-TX(2).pushElements*Trans.spacing/2,4,TX(2).pushElements*Trans.spacing,4]);
                    assignin('base','moveTarget','jump');
                case 'SWIROI'
                    set(recHandle,'Position',ROIpos);
                    assignin('base','moveTarget','jump');
                case 'jump'
            end
        end
    end

    % kernel2D is used for 2D filter to smooth SWI
    kernel2D = ...
    [   0.0073    0.0208    0.0294    0.0208    0.0073;
        0.0208    0.0589    0.0833    0.0589    0.0208;
        0.0294    0.0833    0.1179    0.0833    0.0294;
        0.0208    0.0589    0.0833    0.0589    0.0208;
        0.0073    0.0208    0.0294    0.0208    0.0073;];

    if newIQ == 1  % Need to replot focus IQ data after changing ROI
        assignin('base','newIQ',0);

        axisChannel = linspace(ROIpos(1),ROIpos(1)+ROIpos(3),Width);
        axisDepth   = linspace(ROIpos(2),ROIpos(2)+ROIpos(4),Depth);

        IQ3 = IQBuffer(:,:,1,3); IQ4 = IQBuffer(:,:,1,4);

        ImMean = (imag(IQ3(ROILA)) + imag(IQ4(ROILA)))/2;
        ReMean = (real(IQ3(ROILA)) + real(IQ4(ROILA)))/2;
        Im = (imag(IQ3(ROILA))-ImMean) .* (real(IQ4(ROILA))-ReMean) - ...
            (imag(IQ4(ROILA))-ImMean) .* (real(IQ3(ROILA))-ReMean);
        Re = (imag(IQ3(ROILA))-ImMean) .* (imag(IQ4(ROILA))-ImMean) + ...
            (real(IQ3(ROILA))-ReMean) .* (real(IQ4(ROILA))-ReMean);
        Power = reshape((Im .* Im + Re .* Re).^0.125,Depth,Width);
        buffer = filter2(Power,kernel2D,'full');
        bufferS = size(buffer);
        ind1 = (bufferS(1)-Depth)/2;
        ind2 = (bufferS(2)-Width)/2;
        Power = rot90(buffer(ind1+1:end-ind1,ind2+1:end-ind2),2);
        figure(SWIfigHandle), if ishandle(myHandle), delete(myHandle); end
        myHandle = imagesc(axisChannel,axisDepth,Power);
        title('Shear Wave Visualization with L7-4','FontSize',14,'Fontweight','bold');
        xlabel('Wavelength','FontSize',12,'Fontweight','bold');
        ylabel('Wavelength','FontSize',12,'Fontweight','bold');
        set(gca,'FontSize',12,'Fontweight','bold');
        axis tight equal, colormap('gray')
        MovieData(:,:,1) = Power;
        assignin('base','SWIimageHandle',myHandle);
    end

    caxis(get(SWIfigHandle,'currentAxes'),[50,powermax]);
    % The size of IQData here is [nRows, nCols, nFrames, nPages]
    for i = 2:na-1 % for all combinations of 2 pages
        IQ1 = IQBuffer(:,:,1,i); IQ2 = IQBuffer(:,:,1,i+1);
        ImMean = (imag(IQ1(ROILA)) + imag(IQ2(ROILA)))/2;
        ReMean = (real(IQ1(ROILA)) + real(IQ2(ROILA)))/2;
        Im = (imag(IQ1(ROILA))-ImMean) .* (real(IQ2(ROILA))-ReMean) - ...
            (imag(IQ2(ROILA))-ImMean) .* (real(IQ1(ROILA))-ReMean);
        Re = (imag(IQ1(ROILA))-ImMean) .* (imag(IQ2(ROILA))-ImMean) + ...
            (real(IQ1(ROILA))-ReMean) .* (real(IQ2(ROILA))-ReMean);
        Power = reshape((Im .* Im + Re .* Re).^0.125,Depth,Width);
        buffer = filter2(Power,kernel2D,'full');
        bufferS = size(buffer);
        ind1 = (bufferS(1)-Depth)/2;
        ind2 = (bufferS(2)-Width)/2;
        Power = rot90(buffer(ind1+1:end-ind1,ind2+1:end-ind2),2);
        set(myHandle,'CData',Power);
        MovieData(:,:,i) = Power;
        drawnow
    end

    assignin('base','MovieData',MovieData);


end

%% Other external functions used in callback functions

function moveFocus(varargin)

    freeze = evalin('base','freeze');
    vsExit = evalin('base','vsExit');

    if freeze == 0 && vsExit == 0  % no response if freeze or exit

        UI = evalin('base','UI');
        TX = evalin('base','TX');
        Trans = evalin('base','Trans');
        SWIROI = evalin('base','SWIROI');
        Location = evalin('base','Location');

        Location.oldX = TX(2).focusX;
        Location.oldZ = TX(2).focus;

        bmodeFigHandle = evalin('base','Resource.DisplayWindow(1).figureHandle');
        bmodeAxes = get(bmodeFigHandle,'currentAxes');
        currentPos = get(bmodeAxes,'CurrentPoint');
        YLim = get(bmodeAxes,'YLim');

        Location.newX = currentPos(1);
        Location.newZ = currentPos(3);

        if Location.newZ > YLim(2)-5, Location.newZ = YLim(2)-5;
        elseif Location.newZ < 10, Location.newZ = 10;
        end

        assignin('base','Location',Location);
        assignin('base','moveTarget','Focus');
        focusCorrection();

        Location = evalin('base','Location');

        msg = ['Move focus to (',num2str(Location.newX,'%3.2f'),',',num2str(Location.newZ,'%3.2f'),')....'];
        TXupdate(msg,'TX');

    end
end

function moveROI(varargin)

    freeze = evalin('base','freeze');
    vsExit = evalin('base','vsExit');

    if freeze == 0 && vsExit == 0  % no response if freeze or exit

        ax = axis;
        xrange = ax(2)-ax(1); % range of x-axis
        yrange = ax(4)-ax(3); % range of y-axis

        gcfUnits = get(gcf,'Units');
        pltUnits = get(gca,'Units');
        set(gcf,'Units','pixels')
        set(gca,'Units','normalized');

        theRect = get(gcbo,'Position'); % position of clicked rectangle (figure units)
        theFig = get(gcf,'Position'); % position of figure on desktop (pixel units)
        thePlot = get(gca,'Position'); % position of plot within figure (percentages)

        pt1 = get(gcf,'CurrentPoint'); % button down detected (pixel units)
        pt2 = get(gca,'CurrentPoint'); % button down detected (figure units)

        set(gcf,'Units',gcfUnits);
        set(gca,'Units',pltUnits);

        xoffset = pt2(1,1)-theRect(1);
        if isequal( get(gca,'Ydir'), 'reverse' )
            yoffset = theRect(2)-pt2(1,2)+theRect(4);
        else
            yoffset = pt2(1,2)-theRect(2);
        end

        x = pt1(1) - xoffset/xrange * theFig(3)*thePlot(3); % calc. x location in pixels
        y = pt1(2) - yoffset/yrange * theFig(4)*thePlot(4); % calc. y location in pixels
        wt = theRect(3)/xrange * theFig(3) * thePlot(3); % calc. width in pixels
        ht = theRect(4)/yrange * theFig(4) * thePlot(4); % calc. height in pixels

        beforeRect = [x y wt ht]; % starting rectange (pixels)
        afterRect = dragrect(beforeRect); % move around
        pixDiff = afterRect(1:2)-beforeRect(1:2); % find the change in x/y pixel units
        figDiff = pixDiff .* ... % convert change to figure units
            [xrange/(theFig(3)*thePlot(3)) yrange/(theFig(4)*thePlot(4))];

        if isequal( get(gca,'Ydir'), 'reverse' )
            newY = pt2(1,2) - (figDiff(2)-yoffset)-theRect(4);
        else
            newY = pt2(1,2) + figDiff(2) - yoffset;
        end

        newX = pt2(1,1) + figDiff(1) - xoffset;

        P = evalin('base','P');
        Location = evalin('base','Location');

        if newX < Location.Element(1)-0.5, newX = Location.Element(1)-0.5;
        elseif newX > Location.Element(end)+0.5-theRect(3), newX = Location.Element(end)+0.5-theRect(3); end

        if newY < P.startDepth, newY = P.startDepth;
        elseif newY > P.endDepth - theRect(4), newY = P.endDepth - theRect(4); end

        set(gcbo,'Position',[newX newY theRect(3) theRect(4)]) % move rectangle to the new location
        assignin('base','ROIpos',[newX newY theRect(3) theRect(4)]);

        % Trans and PData are required for computeRegions
        PData = evalin('base','PData');
        SWIROI = evalin('base','SWIROI');

        PData(2).Region.Shape.Position = [newX+SWIROI.width/2,0,newY];
        PData(2).Region = computeRegions(PData(2));

        assignin('base','PData',PData);
        assignin('base','newIQ',1);

        Control = evalin('base','Control');
        Control.Command = 'update&Run';
        Control.Parameters = {'PData','Recon'};
        assignin('base','Control', Control);

    end
end

function TXupdate(finalMsg,adjustCase)

    TX = evalin('base','TX');
    Trans = evalin('base','Trans');
    Location = evalin('base','Location');

    % Location.newX won't be assigned if new focus is not determined
    if ~isfield(Location,'pushStartEle')
        Ind = find(TX(2).Apod == 1);
        Location.pushStartEle = Ind(1);
    end

    oldApod = TX(2).Apod;
    TX(2).Apod = zeros(1,Trans.numelements);
    TX(2).Apod(Location.pushStartEle:Location.pushStartEle+TX(2).pushElements-1) = 1;
    fprintf([finalMsg, '\n']);

    switch adjustCase
        case 'TX'
            TX(2).focusX = Location.newX;
            TX(2).Origin = Location.newX;
            TX(2).focus = Location.newZ;
            TX(2).oldElements = TX(2).pushElements;
            TX(2).Delay = computeTXDelays(TX(2));
            assignin('base','TX', TX);
            Control.Parameters = {'TX'};
        case 'TW'
            evalin('base','TW(2).oldCycles = pushCycle;');
            Control.Parameters = {'TW'};
    end

    Control.Command = 'update&Run';
    assignin('base','Control', Control);

end

function closeIQfig(varargin)

    if evalin('base','vsExit')
        delete(gcf)
    else
        @SWIoff;
    end

end

function SWIoff(varargin)

    % if the movie is replay, no response
    replay = evalin('base','replay');

    if strcmp(replay,'on')
        msgbox('please stop reply first');

    elseif strcmp(replay,'off')

        % make related UI controls unvisible
        UI = evalin('base','UI');
        for i = 3:15
            set(UI(i).handle,'Visible','off');
        end
        set(UI(2).handle,'Visible','on');
        set(UI(14).handle,'Visible','on');

        % Set Focus Adjustment back to off
        set(UI(3).handle(2),'Value',1);
        set(UI(3).handle(3),'Value',0);

        % make SWI figure, blue region, and red mark unvisible
        set(evalin('base','SWIfigHandle'),'Visible','off');

        if evalin('base','exist(''ROIHandle'',''var'')')
            if evalin('base','ishandle(ROIHandle)')
                set(evalin('base','ROIHandle'),'Visible','off');
                set(evalin('base','markHandle'),'Visible','off');
                set(evalin('base','TransHandle'),'Visible','off');
            else
                evalin('base','clear ROIHandle markHandle TransHandle');
            end
        end

        % make P1 slider visible and P5 slider unvisible
        hv1Handle(1) = findobj('String','High Voltage P1');
        hv1Handle(2) = findobj('tag','hv1Sldr');
        hv1Handle(3) = findobj('tag','hv1Value');
        set(hv1Handle,'Visible','on');

        hv2Handle(1) = findobj('String','High Voltage P5');
        hv2Handle(2) = findobj('tag','hv2Sldr');
        hv2Handle(3) = findobj('tag','hv2Value');
        hv2Handle(4) = findobj('tag','hv2Actual');
        set(hv2Handle,'Visible','off');

        % set hv2 back to 1.6v
        hv2 = 1.6;
        set(hv2Handle(2),'Value',hv2);
        feval(hv2Handle(2).Callback{1},hv2Handle(2),[],hv2Handle(3),hv1Handle(2),hv1Handle(3));

        % change the start event
        nStart = evalin('base','nStartFlash');
        Control = evalin('base','Control');
        if isempty(Control(1).Command), n=1; else n=length(Control)+1; end
        Control(n).Command = 'set&Run';
        Control(n).Parameters = {'Parameters',1,'startEvent',nStart};
        evalin('base',['Resource.Parameters.startEvent =',num2str(nStart),';']);
        assignin('base','Control',Control);

        % no WindowButtonDownFcn
        Resource = evalin('base','Resource');
        set(Resource.DisplayWindow(1).figureHandle,'WindowButtonDownFcn',[]);

    end

end

function SWIon(varargin)

    Resource = evalin('base','Resource');

    % Handle for Shearwave figure
    if ~evalin('base','exist(''SWIfigHandle'',''var'')')

        SWIROI = evalin('base','SWIROI');

        SWIfigHandle = figure('Name','ShearWaveVisulization',...
            'NumberTitle','off','Visible','on',...
            'Position',[Resource.DisplayWindow(1).Position(1)+Resource.DisplayWindow(1).Position(3), ... % left edge
            Resource.DisplayWindow(1).Position(2), ... % bottom
            [SWIROI.width+3,SWIROI.height]*7], ...            % width, height
            'CloseRequestFcn',@closeIQfig);
        assignin('base','SWIfigHandle',SWIfigHandle);
    else
        set(evalin('base','SWIfigHandle'),'Visible','on');
    end

    if evalin('base','exist(''ROIHandle'',''var'')')
        if evalin('base','ishandle(ROIHandle)')
            % make ROI, recTrans and mark visible
            set(evalin('base','ROIHandle'),'Visible','on');
            set(evalin('base','markHandle'),'Visible','on');
            set(evalin('base','TransHandle'),'Visible','on');
        end
    end

    % newIQ is used for check whether new SWI axes is required
    evalin('base','newIQ = 1;');

    % make related UI controls visible
    UI = evalin('base','UI');
    for i = 3:15
        set(UI(i).handle,'Visible','on','Interruptible','off');
    end
    set(UI(2).handle,'Visible','off');
    set(UI(9).handle,'Visible','off');
    set(UI(14).handle,'Visible','off');
    set(UI(8).handle,'Interruptible','on');

    % set moveTarget and WindowButtonDownFcn
    assignin('base','moveTarget','Focus');
    set(Resource.DisplayWindow(1).figureHandle,'WindowButtonDownFcn',@moveFocus);
    set(Resource.DisplayWindow(1).figureHandle,'WindowButtonMotionFcn',[]);


    % make P1 slider unvisible and P5 slider visible
    hv1Handle(1) = findobj('String','High Voltage P1');
    hv1Handle(2) = findobj('tag','hv1Sldr');
    hv1Handle(3) = findobj('tag','hv1Value');
    set(hv1Handle,'Visible','off');

    hv2Handle(1) = findobj('String','High Voltage P5');
    hv2Handle(2) = findobj('tag','hv2Sldr');
    hv2Handle(3) = findobj('tag','hv2Value');
    hv2Handle(4) = findobj('tag','hv2Actual');
    set(hv2Handle,'Visible','on');

    focusCorrection();

    % change the start event
    nStart = evalin('base','nStartPush');
    Control = evalin('base','Control');
    if isempty(Control(1).Command), n=1; else n=length(Control)+1; end
    Control(n).Command = 'set&Run';
    Control(n).Parameters = {'Parameters',1,'startEvent',nStart};
    evalin('base',['Resource.Parameters.startEvent =',num2str(nStart),';']);
    assignin('base','Control',Control);

end

function focusCorrection(varargin)

    TX = evalin('base','TX');
    Trans = evalin('base','Trans');
    Location = evalin('base','Location');

    pushElements = TX(2).pushElements;

    if ~isfield(Location,'newX')
        Location.oldX = TX(2).focusX;
        Location.oldZ = TX(2).focus;
        Location.newX = TX(2).focusX;
        Location.newZ = TX(2).focus;
        Ind = find(TX(2).Apod == 1);
        Location.pushStartEle = Ind(1);
    end

    x = Location.newX;
    z = Location.newZ;

    % Location.PosForEvenPush indicates the x-coor of the right edge of each element
    % LOcation.PosForOddPush indicates the x-coor of the center of each element
    if mod(pushElements,2) == 0
        focusLimit = [Location.PosForEvenPush(pushElements/2),Location.PosForEvenPush(Trans.numelements-pushElements/2)];
        if x < focusLimit(1), x = focusLimit(1); elseif x > focusLimit(2), x = focusLimit(2); end
        [value,eleInd] = min(abs(x - Location.PosForEvenPush));
        Location.newX = Location.PosForEvenPush(eleInd);
        Location.pushStartEle = eleInd - pushElements/2 + 1;
    else
        focusLimit = [Location.PosForOddPush((pushElements+1)/2),Location.PosForOddPush(Trans.numelements-(pushElements-1)/2)];
        if x < focusLimit(1), x = focusLimit(1);
        elseif x > focusLimit(2), x = focusLimit(2); end
        [value,eleInd] = min(abs(x - Location.PosForOddPush));
        Location.newX = Location.PosForOddPush(eleInd);
        Location.pushStartEle = eleInd - (pushElements-1)/2;
    end

    Location.pushStartCoor = Location.PosForEvenPush(Location.pushStartEle);
    assignin('base','Location',Location);

end

function moveMouse(varargin)

    h = varargin{3};
    p = get(h.Parent,'CurrentPoint');
    x = p(1,1);
    y = p(1,2);

    ROIHandle = evalin('base','ROIHandle');
    ROIpos = get(ROIHandle,'Position');
    x1 = ROIpos(1);
    x2 = ROIpos(1)+ROIpos(3);
    y1 = ROIpos(2);
    y2 = ROIpos(2)+ROIpos(4);
    d = 1;

    if (x1<x && x<x2 && y1-d<y && y<y1+d) || ...
       (x1<x && x<x2 && y2-d<y && y<y2+d) || ...
       (y1<y && y<y2 && x1-d<x && x<x1+d) || ...
       (y1<y && y<y2 && x2-d<x && x<x2+d)
        set(ROIHandle.Parent.Parent,'pointer','hand');
    else
        set(ROIHandle.Parent.Parent,'pointer','arrow');
    end

end



