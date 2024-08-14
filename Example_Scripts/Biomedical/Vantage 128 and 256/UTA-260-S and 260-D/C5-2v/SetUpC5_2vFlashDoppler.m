% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: SetUpC5_2vFlashDoppler.m - Example of curved array doppler
%                                       imaging with flash transmit
% Description:
%   Sequence programming file for C5-2v curved array, using flash transmits
%   for 2D (B-mode) and doppler imaging.
%   - for the 2D (B-mode) imaging, na steered angles are used.
%   - for the doppler ensemble, ne pulses are transmitted at a steering
%     angle of dopAngle radians.
%   - all 128 transmit and receive channels are active for each acquisition.
%   - processing is asynchronous with respect to acquisition.
%
% Last update:
% 09/11/2020 - Update to SW 4.3 format for new user UIControls and External function definitions (VTS-1691).
%   More info:(?/Example_Scripts/Vantage_Features/New UI Scheme/SetUpL11_5vFlash_NewUI)

clear all

% P(1) is used for Bmode, P(2) is used for Doppler
P(1).startDepth = 2;   % Acquisition depth in wavelengths
P(1).endDepth = 192;   % This should preferrably be a multiple of 128 samples.
P(2).startDepth = 2;   % Acquisition depth in wavelengths
P(2).endDepth = 192;   % This should preferrably be a multiple of 128 samples.

na = 7;      % Set na = number of angles.
if (na > 1)
    dtheta2D = (24*pi/180)/(na-1);
    startAngle = -24*pi/180/2;
else
    dtheta2D=0;
    startAngle=0;
end % set dtheta2D to range over +/- 12 degrees.

% Set Doppler parameters
ne = 14;     % Set ne = number of acquisitions in Doppler ensemble.
dopPRF = 1.0e+03; % Doppler PRF in Hz.
pwrThres = 0.5;

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
Trans.name = 'C5-2v';
Trans.units = 'wavelengths'; % required in Gen3 to prevent default to mm units
Trans.frequency = 2.976;  % center frequency set lower for Doppler.
% note nominal center frequence from computeTrans is 4.464 MHz
Trans = computeTrans(Trans);  % C5-2v transducer is 'known' transducer so we can use computeTrans.
Trans.maxHighVoltage = 50;  % set maximum high voltage limit for pulser supply.
radius = Trans.radius;
scanangle = Trans.numelements*Trans.spacing/radius;

% Specify PData structure array.
% - 2D PData structure
PData(1).PDelta = [1.0, 0, 0.5];  % x, y and z pdeltas
sizeRows = 10 + ceil((P(1).endDepth + radius - (radius * cos(scanangle/2)))/PData(1).PDelta(3));
sizeCols = 10 + ceil(2*(P(1).endDepth + radius)*sin(scanangle/2)/PData(1).PDelta(1));
PData(1).Size = [sizeRows,sizeCols,1];     % size, origin and pdelta set region of interest.
PData(1).Origin(1,1) = (P(1).endDepth+radius)*sin(-scanangle/2) - 5;
PData(1).Origin(1,2) = 0;
PData(1).Origin(1,3) = ceil(radius * cos(scanangle/2)) - radius - 5;
PData(1).Region = struct(...
    'Shape',struct('Name','Sector',...
                   'Position',[0,0,-radius],...
                   'r1',radius+P(1).startDepth,...
                   'r2',radius+P(1).endDepth,...
                   'angle',scanangle));
PData(1).Region = computeRegions(PData(1));

% - Doppler PData structure
PData(2).PDelta = [1.0, 0, 1.0];  % x, y and z pdeltas
sizeRows = 10 + ceil((P(1).endDepth + radius - (radius * cos(-(scanangle/4))))/PData(2).PDelta(3));
sizeCols = 10 + ceil(2*(P(1).endDepth + radius)*sin(scanangle/4)/PData(2).PDelta(1));
PData(2).Size = [sizeRows,sizeCols,1];     % size, origin and pdelta set region of interest.
PData(2).Origin(1,1) = (P(1).endDepth+radius)*sin(-(scanangle/4)) - 5;
PData(2).Origin(1,2) = 0;
PData(2).Origin(1,3) = ceil(radius * cos(scanangle/4)) - radius - 5;
PData(2).Region = struct(...
    'Shape',struct('Name','Sector',...
                   'Position',[0,0,-radius],...
                   'r1',radius+P(1).startDepth,...
                   'r2',radius+P(1).endDepth,...
                   'angle',scanangle/2));
PData(2).Region = computeRegions(PData(2));

%  Media points for curved array.
% - Uncomment for speckle
% Media.numPoints = (20000);
% Media.MP = rand(Media.numPoints,4);
% Media.MP(:,2) = 0;
% Media.MP(:,4) = 0.01 + 0.04*Media.MP(:,4);  % Random low amplitude
% RandR = P.endDepth *Media.MP(:,1)+radius;
% RandTheta = scanangle*(Media.MP(:,3)-0.5);
% Media.MP(:,1) = RandR.*sin(RandTheta);
% Media.MP(:,3) = RandR.*cos(RandTheta)-radius;
% - Define points
%Media.MP(1,:) = [0,0,70,1.0];
Media.MP(1,:) = [0,0,10,1.0];
Media.MP(2,:) = [(radius+10)*sin(-0.2608),0,(radius+10)*cos(-0.2608)-radius,1.0];
Media.MP(3,:) = [(radius+10)*sin(0.2608),0,(radius+10)*cos(0.2608)-radius,1.0];
Media.MP(4,:) = [(radius+10)*sin(-0.5267),0,(radius+10)*cos(-0.5267)-radius,1.0];
Media.MP(5,:) = [(radius+10)*sin(0.5267),0,(radius+10)*cos(0.5267)-radius,1.0];
Media.MP(6,:) = [0,0,40,1.0];
Media.MP(7,:) = [0,0,70,1.0];
Media.MP(8,:) = [(radius+70)*sin(-0.2608),0,(radius+70)*cos(-0.2608)-radius,1.0];
Media.MP(9,:) = [(radius+70)*sin(0.2608),0,(radius+70)*cos(0.2608)-radius,1.0];
Media.MP(10,:) = [(radius+70)*sin(-0.5267),0,(radius+70)*cos(-0.5267)-radius,1.0];
Media.MP(11,:) = [(radius+70)*sin(0.5267),0,(radius+70)*cos(0.5267)-radius,1.0];
Media.MP(12,:) = [0,0,100,1.0];
Media.MP(13,:) = [0,0,130,1.0];
Media.MP(14,:) = [(radius+130)*sin(-0.2608),0,(radius+130)*cos(-0.2608)-radius,1.0];
Media.MP(15,:) = [(radius+130)*sin(0.2608),0,(radius+130)*cos(0.2608)-radius,1.0];
Media.MP(16,:) = [(radius+130)*sin(-0.5267),0,(radius+130)*cos(-0.5267)-radius,1.0];
Media.MP(17,:) = [(radius+130)*sin(0.5267),0,(radius+130)*cos(0.5267)-radius,1.0];
Media.MP(18,:) = [0,0,160,1.0];
Media.MP(19,:) = [0,0,190,1.0];
Media.MP(20,:) = [(radius+190)*sin(-0.2608),0,(radius+190)*cos(-0.2608)-radius,1.0];
Media.MP(21,:) = [(radius+190)*sin(0.2608),0,(radius+190)*cos(0.2608)-radius,1.0];
Media.MP(22,:) = [(radius+190)*sin(-0.5267),0,(radius+190)*cos(-0.5267)-radius,1.0];
Media.MP(23,:) = [(radius+190)*sin(0.5267),0,(radius+190)*cos(0.5267)-radius,1.0];
Media.function = 'movePoints';

% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 2048*(na + ne);
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 20;     % Two frames needed for double buffer.
% InterBuffer(1) is for 2D reconstructions.
Resource.InterBuffer(1).numFrames = 1;          % one intermediate frame needed for 2D.
% InterBuffer(2) is for Doppler reconstructions.
Resource.InterBuffer(2).numFrames = 1;          % one intermediate frame needed for Doppler.
Resource.InterBuffer(2).pagesPerFrame = ne;     % ne pages per ensemble
% ImageBuffer(1) is for 2D image.
Resource.ImageBuffer(1).numFrames = 20;
% ImageBuffer(2) is for Doppler image.
Resource.ImageBuffer(2).numFrames = 20;
% DisplayWindow is for 2D combined with Doppler
Resource.DisplayWindow(1).Title = 'C5-2vFlashDoppler';
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

% Specify Transmit waveform structure.
% - 2D transmit waveform
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,0.67,2,1];
% - Doppler transmit waveform, transmit frequency should be equivalent to
% supported demodFrequency
TW(2).type = 'parametric';
TW(2).Parameters = [3.7,0.67,8,1];

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,radius], ...
                   'focus', 0.0, ...
                   'Steer', [0.0,0.0], ...
                   'Apod', ones(1,Trans.numelements), ...
                   'Delay', zeros(1,Trans.numelements)), 1, na+1);
% - Set event specific TX attributes.
for n = 1:na   % na transmit events
    TX(n).Steer = [(startAngle+(n-1)*dtheta2D),0.0];
    TX(n).Delay = computeTXDelays(TX(n));
end
% -- only one TX struct needed for Doppler
TX(na+1).waveform = 2;
TX(na+1).Steer = [0,0.0];
TX(na+1).Delay = computeTXDelays(TX(na+1));

% Specify TPC structures.
TPC(1).name = '2D';
TPC(1).maxHighVoltage = 50;
TPC(2).name = 'Doppler';
TPC(2).maxHighVoltage = 35;

RcvProfile(1).LnaGain = 18; % Profile used for imaging
RcvProfile(2).LnaGain = 24; % Profile used for Doppler
RcvProfile(2).LnaZinSel = 31; % Force high-Z state for best Doppler sensitivity

% Specify TGC Waveform structures.
% - 2D TGC
TGC.CntrlPts = [0,141,275,404,510,603,702,782];
TGC(1).rangeMax = P(1).endDepth;
TGC(1).Waveform = computeTGCWaveform(TGC(1));
% - Doppler TGC
TGC(2).CntrlPts =[ 0 272 662 662 662 662 662 662];
TGC(2).rangeMax = P(2).endDepth;
TGC(2).Waveform = computeTGCWaveform(TGC(2));

% Specify Receive structure arrays.
%   We need to acquire all the 2D and Doppler data within a single RcvBuffer frame.  This allows
%   the transfer-to-host DMA after each frame to transfer a large amount of data, improving throughput.
% - We need na Receives for a 2D frame and ne Receives for a Doppler frame.
maxAcqLength2D = sqrt((P(1).endDepth+radius)^2 + radius^2 - ...
                     2*(P(1).endDepth+radius)*radius*cos(scanangle)) - P(1).startDepth;
maxAcqLengthDop = sqrt((P(2).endDepth+radius)^2 + radius^2 - ...
                     2*(P(2).endDepth+radius)*radius*cos(scanangle/2)) - P(2).startDepth;
wl4sPer128 = 128/(4*2);  % wavelengths in a 128 sample block for 4 smpls per wave round trip.
wl2sPer128 = 128/(2*2);  % wavelengths in a 128 sample block for 2 smpls per wave round trip.
Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', P(1).startDepth, ...
                        'endDepth', P(1).startDepth + wl4sPer128*ceil(maxAcqLength2D/wl4sPer128), ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ... % 200% Bandwidth for 2D
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, (na+ne)*Resource.RcvBuffer(1).numFrames);
%   'InputFilter', [-0.0141,0.0420,-0.0825,0.1267,-0.1612,0.1783], ...  % bandpass filter for Doppler (~ 20% BW)
% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    k = (na + ne)*(i-1); % k keeps track of Receive index increment per frame.
    % - Set attributes for each frame.
    Receive(k+1).callMediaFunc = 1;
    for j = 1:na  % acquisitions for 2D).
        % 2D acquisition
        Receive(j+k).framenum = i;
        Receive(j+k).acqNum = j;
    end
    for j = (na+1):(na+ne)
        % Doppler acquisition
        Receive(j+k).startDepth = P(2).startDepth;
        Receive(j+k).endDepth = P(2).startDepth + wl2sPer128*ceil(maxAcqLengthDop/wl2sPer128);
        Receive(j+k).sampleMode = 'BS100BW';
        Receive(j+k).demodFrequency = TW(2).Parameters(1);
        Receive(j+k).TGC = 2;
        Receive(j+k).framenum = i;
        Receive(j+k).acqNum = j;        % Doppler acqNums continue after 2D
    end
end

% Specify Recon structure arrays.
% - We need two Recon structures, one for 2D, one for Doppler. These will be referenced in the same
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
% - Set Recon values for Doppler ensemble.
Recon(2).pdatanum = 2;
Recon(2).IntBufDest = [2,1];
Recon(2).ImgBufDest = [0,0];
Recon(2).RINums(1,1:ne) = k:(k+ne-1);   % ne ReconInfos needed for Doppler ensemble.

% Define ReconInfo structures.
% - For 2D, we need na ReconInfo structures for na steering angles.
% - For Doppler, we need ne ReconInfo structures.
ReconInfo = repmat(struct('mode', 'accumIQ', ...    % accumulate IQ data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'pagenum',1, ...
                   'regionnum', 1), 1, na + ne);
% - Set specific ReconInfo attributes.
%   - ReconInfos for 2D frame.
if na>1
    ReconInfo(1).mode = 'replaceIQ'; % replace IQ data
    for j = 1:na
        ReconInfo(j).txnum = j;
        ReconInfo(j).rcvnum = j;
    end
    ReconInfo(na).mode = 'accumIQ_replaceIntensity'; % accum and detect
else
    ReconInfo(1).mode = 'replaceIntensity';
end

%  - ReconInfos for Doppler ensemble.
k = na;
for j = 1:ne
    ReconInfo(k+j).mode = 'replaceIQ';
    ReconInfo(k+j).txnum = na + 1;
    ReconInfo(k+j).rcvnum = na + j;
    ReconInfo(k+j).pagenum = j;
end

% Specify Process structure arrays.
cpt = 28;  % define here so we can use in UIControl below
persf = 80;
persp = 90;
DopState = 'freq';

Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',1.0,...            % pgain is image processing gain
                         'reject',2,...      % reject level
                         'persistMethod','simple',...
                         'persistLevel',20,...
                         'interpMethod','4pt',...
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',40,...
                         'mappingMethod','lowerHalf',...
                         'display',0,...      % don't display image after processing
                         'displayWindow',1};

Process(2).classname = 'Doppler';                   % process structure for 1st Doppler ensemble
Process(2).method = 'computeCFIFreqEst';
Process(2).Parameters = {'IntBufSrc',[2,1],...          % number of Interbuffer to process.
                         'SrcPages',[3,ne-2],...        % start frame number in source buffer
                         'ImgBufDest',[2,-1],...
                         'pdatanum',2,...           % number of PData structure
                         'prf',dopPRF,...           % Doppler PRF in Hz
                         'wallFilter','regression',...  % 1 -> quadratic regression
                         'pwrThreshold',pwrThres,...
                         'maxPower',50,...
                         'postFilter',1};

Process(3).classname = 'Image';
Process(3).method = 'imageDisplay';
Process(3).Parameters = {'imgbufnum',2,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'srcData','signedColor',... % type of data to display.
                         'pdatanum',2,...    % number of PData structure to use
                         'pgain',1.0,...            % pgain is image processing gain
                         'reject',2,...      % reject level
                         'persistMethod','dynamic',...
                         'persistLevel',persf,...
                         'interpMethod','4pt',...
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',40,...
                         'mappingMethod','upperHalf',...
                         'threshold',cpt,...
                         'display',1,...      % don't display image after processing
                         'displayWindow',1};

% Specify SeqControl structure arrays.
% -- Time between 2D flash angle acquisitions
SeqControl(1).command = 'timeToNextAcq';
SeqControl(1).argument = 250;
% -- Change to Profile 2 (Doppler)
SeqControl(2).command = 'setTPCProfile';
SeqControl(2).condition = 'next';
SeqControl(2).argument = 2;
% -- Time between 2D acquisition and Doppler ensemble. Set to allow time for profile change.
SeqControl(3).command = 'timeToNextAcq';
SeqControl(3).argument = 7000; % time in usec
% -- PRF for Doppler ensemble
SeqControl(4).command = 'timeToNextAcq';
SeqControl(4).argument = round(1/(dopPRF*1e-06)); % (for 3KHz dopPRF & 14 ensemble = 4.7 msecs)
% -- Change to Profile 1 (2D)
SeqControl(5).command = 'setTPCProfile';
SeqControl(5).condition = 'next';
SeqControl(5).argument = 1;
% -- Time between Doppler and next 2D acquisition. Set to allow time for profile change.
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
% - The remainder of the SeqControl structures are defined dynamically in the sequence events.
%   The variable nsc keeps track of the next SeqControl structure index.
nsc = 10;  % next SeqControl number

% Specify Event structure arrays.
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
        if j == 1
            Event(n).seqControl = [1,8];
        end
        n = n+1;
    end
    Event(n-1).seqControl = [2,3];   % replace last 2D acquisition Event's seqControl
    % Acquire Doppler ensemble.
    for j = (na+1):(na+ne)
        Event(n).info = 'Acquire Doppler ensemble';
        Event(n).tx = na+1;
        Event(n).rcv = (na+ne)*(i-1)+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 4;
        if j == na+1
            Event(n).seqControl = [4,9];
        end
        n = n+1;
    end
    Event(n-1).seqControl = [5,6,nsc]; % replace last Doppler acquisition Event's seqControl
      SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
      nsc = nsc+1;

    Event(n).info = 'recons and 2D process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = [1,2];
    Event(n).process = 1;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = 'Doppler processing';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 2;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = 'Doppler image display';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 3;
    Event(n).seqControl = 0;
    if floor(i/3) == i/3     % Exit to Matlab every 3rd frame
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


% User specified UI Control Elements

import vsv.seq.uicontrol.VsSliderControl
import vsv.seq.uicontrol.VsButtonGroupControl

% - Sensitivity Cutoff
UI(1).Control = VsSliderControl('LocationCode','UserB7',...
    'Label','Sens. Cutoff',...
    'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
    'SliderStep',[0.025,0.1],'ValueFormat','%1.3f',...
    'Callback',@SensCutoffCallback);

% - Doppler Mode Button Group
UI(2).Control = VsButtonGroupControl('LocationCode','UserB4',...
    'Title','Doppler Mode','NumButtons',2,...
    'Labels',{'Velocity','Power'},...
    'Callback',@DopplerModeChange);

% - Doppler Power Threshold Slider
UI(3).Control = VsSliderControl('LocationCode','UserB3',...
    'Label','DopPwrThres','SliderMinMaxVal',[0.0,1.0,pwrThres],...
    'SliderStep',[0.02,0.1],'ValueFormat','%3.2f',...
    'Callback',@DopplerPowerChange);

% - Color Priority Threshold Slider
UI(4).Control = VsSliderControl('LocationCode','UserB2',...
    'Label','Color Priority','SliderMinMaxVal',[0,255,cpt],...
    'SliderStep',[1/255,0.1],'ValueFormat','%3.0f',...
    'Callback',@ColorPriorityChange);

% - Color Persistence Slider
UI(5).Control = VsSliderControl('LocationCode','UserB1',...
    'Label','Color Persistence','SliderMinMaxVal',[0,100,persf],...
    'SliderStep',[1/100,0.1],'ValueFormat','%3.0f',...
    'Callback',@ColorPersistenceChange);

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 3;

% Save all the structures to a .mat file.
save('MatFiles/C5-2vFlashDoppler');

%% Automatic VSX Execution:
% Uncomment the following line to automatically run VSX every time you run
% this SetUp script (note that if VSX finds the variable 'filename' in the
% Matlab workspace, it will load and run that file instead of prompting the
% user for the file to be used):

filename = 'C5-2vFlashDoppler';  VSX;


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

function DopplerModeChange(~, ~, UIState)
    Control = repmat(struct('Command','set&Run','Parameters',[]),1,4);
    Process = evalin('base','Process');
    Resource = evalin('base','Resource');

    switch UIState
       case 1  % Velocity mode
          newMap = grayscaleCFImap;
          newMap(1:128,:) = Resource.DisplayWindow(1).Colormap(1:128,:);
          Resource.DisplayWindow(1).Colormap = newMap;
          assignin('base','persp',get(findobj('Tag','UserB1Slider'),'Value'));
          persf = evalin('base','persf'); persValue = persf(1);
          Control(1).Parameters = {'Process',2,'method','computeCFIFreqEst'};
          Control(2).Parameters = {'Process',3,'srcData','signedColor','persistMethod','dynamic','persistLevel',persValue};
          Control(3).Parameters = {'DisplayWindow',1,'colormap',newMap};
          Control(4).Parameters = {'ImageBuffer',1,'lastFrame',0};
          set(findobj('tag','UserB1Edit'),'String',num2str(persValue,'%3.0f'));
          set(findobj('tag','UserB1Slider'),'Value',persValue);
          assignin('base','DopState','freq');
          % Set modified Process attributes in base Matlab environment.
          Process(2).method = 'computeCFIFreqEst';
          for k = 1:2:length(Process(3).Parameters)
              if strcmp(Process(3).Parameters{k},'srcData'), Process(3).Parameters{k+1} = 'signedColor';
              elseif strcmp(Process(3).Parameters{k},'persistMethod'), Process(3).Parameters{k+1} = 'dynamic';
              elseif strcmp(Process(3).Parameters{k},'persistLevel'), Process(3).Parameters{k+1} = persf;
              end
          end
       case 2  % Power mode
          newMap = grayscaleCPAmap;
          newMap(1:128,:) = Resource.DisplayWindow(1).Colormap(1:128,:);
          Resource.DisplayWindow(1).Colormap = newMap;
          for k = 1:2:length(Process(3).Parameters)
              if strcmp(Process(3).Parameters{k},'persistLevel'), persf = Process(3).Parameters{k+1}; end
          end
          assignin('base','persf',persf);
          persValue = evalin('base','persp');
          Control(1).Parameters = {'Process',2,'method','computeCFIPowerEst'};
          Control(2).Parameters = {'Process',3,'srcData','unsignedColor','persistMethod','simple','persistLevel',persValue};
          Control(3).Parameters = {'DisplayWindow',1,'colormap',newMap};
          Control(4).Parameters = {'ImageBuffer',1,'lastFrame',0};
          set(findobj('tag','UserB1Edit'),'String',num2str(persValue,'%3.0f'));
          set(findobj('tag','UserB1Slider'),'Value',persValue);
          assignin('base','DopState','power');
          Process(2).method = 'computeCFIPowerEst';
          for k = 1:2:length(Process(3).Parameters)
              if strcmp(Process(3).Parameters{k},'srcData'), Process(3).Parameters{k+1} = 'unsignedColor';
              elseif strcmp(Process(3).Parameters{k},'persistMethod'), Process(3).Parameters{k+1} = 'simple';
              elseif strcmp(Process(3).Parameters{k},'persistLevel'), Process(3).Parameters{k+1} = persValue;
              end
          end
    end

    assignin('base','newMap',newMap);
    evalin('base','Resource.DisplayWindow(1).Colormap = newMap;');
    assignin('base','Process',Process);
    assignin('base','Control', Control);

    % If PTool window is open, adjust all uicontrols
    hPTool = findobj('tag','ProcessTool');
    if ishandle(hPTool)
        posPTool = get(hPTool,'position');
        PTool;
        set(findobj('tag','ProcessTool'),'position',posPTool);
    end

    % If ColorMapTool is open, close it.
    hCMTool = findobj('tag','ColorMapTool');
    if ishandle(hCMTool)
        delete(hCMTool);
        set(findobj('tag','toolsMenu'),'Value',1); % set tools selection back to none
    end
end

function DopplerPowerChange(~, ~, UIValue)
    Process = evalin('base','Process');
    for k = 1:2:length(Process(2).Parameters)
        if strcmp(Process(2).Parameters{k},'pwrThreshold'), Process(2).Parameters{k+1} = UIValue; end
    end
    assignin('base','Process',Process);
    % Set Control.Command to set Doppler threshold.
    Control = evalin('base','Control');
    Control.Command = 'set&Run';
    Control.Parameters = {'Process',2,'pwrThreshold',UIValue};
    assignin('base','Control', Control);
end

function ColorPriorityChange(~, ~, UIValue)
    % Set the value in the Process structure for use in cineloop playback.
    Process = evalin('base','Process');
    for k = 1:2:length(Process(3).Parameters)
        if strcmp(Process(3).Parameters{k},'threshold')
            Process(3).Parameters{k+1} = UIValue;
        end
    end
    assignin('base','Process',Process);
    % Set Control.Command to set Image.threshold.
    Control = evalin('base','Control');
    Control.Command = 'set&Run';
    Control.Parameters = {'Process',3,'threshold',UIValue};
    assignin('base','Control', Control);
end

function ColorPersistenceChange(~, ~, UIValue)
    % Set the value in the Process structure for use in cineloop playback.
    Process = evalin('base','Process');
    for k = 1:2:length(Process(3).Parameters)
        if strcmp(Process(3).Parameters{k},'persistLevel'), Process(3).Parameters{k+1} = UIValue; end
    end
    assignin('base','Process',Process);
    % Set Control.Command to set Image.persistLevel.
    Control = evalin('base','Control');
    Control.Command = 'set&Run';
    Control.Parameters = {'Process',3,'persistLevel',UIValue};
    assignin('base','Control', Control);

    % If PTool window is open, adjust persistLevel1 in Process(3)
    hPTool = findobj('tag','ProcessTool');
    if ishandle(hPTool)
        if isequal(get(findobj('tag','processNum'),'Value'),3)
            set(findobj('tag','persistSlider1'),'Value',UIValue);
            set(findobj('tag','persistValue1'),'String',num2str(UIValue));
        end
    end
end