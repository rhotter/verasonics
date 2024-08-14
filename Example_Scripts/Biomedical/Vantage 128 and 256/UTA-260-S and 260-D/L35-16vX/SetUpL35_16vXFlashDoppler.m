% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use.
%
% File name SetUpL35-16vXFlashDoppler.m:
%    Generate .mat Sequence Object file for L35-16vX Linear array, with 2D using
%    multiple plane wave steering angles with interleaved 4X sampling. Doppler
%    uses a single ensemble of plane waves at a steered angle with 67% BW
%    sampling.
%
% Description: This script uses the "interleaved sampling" mechanism for 2D RF
%    acquisition.  The L35-16vX is operated at an offset center frequency for
%    Recon processing of 31.25 MHz, with 4X sampling so the required RF data
%    acquisition sample rate is 125 MHz.  Since the HW system cannot acquire
%    data at sample rates above 62.5 MHz, we acquire two sets of receive data
%    for each line, at an A/D sample rate of 62.5 MHz.  These pairs of
%    acquisition data lines are then interleaved in a 'Pre' routine for the
%    reconstruction processing to produce a single composite line sampled at
%    125 MHz.  The data acquired in the second acquistion of each pair
%    must have its data shifted by half a sample period, so it can be
%    interleaved with the first line.  This is accomplished by adding an
%    offset to the transmit delays for the first acquisition equal to that
%    half-sample interval; this in effect shifts each receive sample's data
%    earlier by the desired one half sample period.
%        For Doppler acquisitions the RF data are sampled at 4/3 of the demod
%    frequency of 26.7858, or 35.7143 MHz, which is 250/7. This provides 67%
%    bandwith, supporting frequencies from 17.813 to 35.714 MHz. To transmit
%    at an effective center frequency of 26.7858, an arbitrary waveform is used,
%    which is loaded from the file 'AwDop25.mat'.
%
% Last update: 10/30/17 Working with releases >= 3.3

clear all

% Specify parameters for settings and presets.
% - 2D settings
P(1).startDepth = 5;
P(1).endDepth = 192;
P(1).numAngles = 8;
% -- set angle increment, dtheta, to space lines over +/- 18 degrees.
if (P(1).numAngles > 1)
    dtheta = (36*pi/180)/(P(1).numAngles-1);
    startAngle = -36*pi/180/2;
else dtheta = 0;
    startAngle=0;
end
% - Doppler settings
P(2).startDepth = 5;
P(2).endDepth = 128;
ne = 16;     % Set ne = number of acquisitions in Doppler ensemble.
dopAngle = 12 * pi/180;
dopPRF = 3.0e+03; % Doppler PRF in Hz.
pwrThres = 0.5;

% Specify system parameters.
Resource.Parameters.connector = 1;
Resource.Parameters.numTransmit = 128;  % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;  % number of receive channels.
Resource.Parameters.speedOfSound = 1540;
Resource.Parameters.simulateMode = 0; % 0 means no simulation, if hardware is present.
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

Resource.Parameters.fakeScanhead = 1;

% Specify Trans structure array.
Trans.name = 'L35-16vX';
Trans.units = 'mm';
Trans = computeTrans(Trans);
mm2wl = 1000*Trans.frequency/Resource.Parameters.speedOfSound; % conversion factor from mm to wavelengths

% Set up receive profile 1 for 2D and 2 for Doppler
RcvProfile(1).LnaZinSel = 31; % put LNA in "high-z" input state for best sensitivity for 2D
RcvProfile(1).AntiAliasCutoff = 35; % 2D acquisition is at 25 MHz

RcvProfile(2).LnaZinSel = 0; % enable the highpass filter for 4/3 sampling Doppler acquisition
RcvProfile(2).AntiAliasCutoff = 30; % Doppler acquisition is at 26.8 MHz with narrow bandwidth

% Specify PData structure array.
% - specify 2D PData
PData(1).PDelta = [Trans.spacing/2, 0, 0.5]; % [pdeltaX, pdeltaY, pdeltaZ]
PData(1).Size(1) = ceil((P(1).endDepth-P(1).startDepth)/PData(1).PDelta(3)); % range, origin and pdelta set PData(1).Size.
PData(1).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(3) = 1;      % single image page
PData(1).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of uppr lft crnr.
% - No PData(1).Region specified, so a default Region for the entire PData(1) array will be created by computeRegions.
% - Doppler PData structure
PData(2).PDelta = [Trans.spacing, 0, 0.5];
PData(2).Size(1) = ceil((P(2).endDepth-P(2).startDepth)/PData(2).PDelta(3)); % flow window rows
PData(2).Size(2) = ceil((128*Trans.spacing)/PData(2).PDelta(1));  % flow window columns
PData(2).Size(3) = 1;             % single image page
PData(2).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P(2).startDepth]; % x,y,z of upper lft crnr.
PData(2).Region = struct('Shape',struct(...
    'Name','Parallelogram',...
    'Position',[0,0,P(2).startDepth],...
    'width', (Trans.numelements*Trans.spacing),...
    'height', P(2).endDepth-P(2).startDepth,...
    'angle',dopAngle));

% Specify Media object.
pt1;
Media.function = 'movePoints';

% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = P(1).numAngles*6144+ne*2048; % this size allows for maximum range
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 20;    % frames stored in RcvBuffer.
% InterBuffer(1) is for 2D reconstructions.
Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).numFrames = 1;   % one intermediate buffer needed.
% InterBuffer(2) is for Doppler reconstructions.
Resource.InterBuffer(2).numFrames = 1;          % one intermediate frame needed for Doppler.
Resource.InterBuffer(2).pagesPerFrame = ne;     % ne pages per ensemble
% ImageBuffer(1) is for 2D image.
Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).numFrames = 10;
% ImageBuffer(2) is for Doppler image.
Resource.ImageBuffer(2).numFrames = 10;
Resource.DisplayWindow(1).Title = 'L35-16vX Flash Doppler, Interleaved 2D Sampling';
Resource.DisplayWindow(1).pdelta = 0.3;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).numFrames = 20;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = grayscaleCFImap;

% Specify Transmit waveform structures.
% - Doppler transmit waveform
load('AwDop26_79.mat','TW');
% - 2D transmit waveform
TW(2).type = 'parametric';
TW(2).Parameters = [25,1,3,1];   % A, B, C, D

TW(3).type = 'pulseCode';
TW(3).PulseCode =  [ 1   3   2  -3   1; ...
                     1   4   0  -5   1; ...
                     0   5   0  -4   1; ...
                     0   5   0  -5   1; ...
                     0   4   0  -5   1; ...
                     0   5   0  -4   1; ...
                     0   5   0  -5   1; ...
                     0   4   1  -3   1; ...
                     2   3   1   0   1 ];


% Specify TX structure array.  We need two TX definitions for each interleave 2D acquisition and one for Doppler.
TX = repmat(struct('waveform', 2, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'focus', 0.0, ...
                   'Steer', [0.0,0.0], ...
                   'Apod', ones(1,Trans.numelements), ...
                   'Delay', zeros(1,Trans.numelements)), 1, 2*P(1).numAngles+1);
% - Set event specific TX attributes.
%   Define two transmit events at each of the desired 2D steering angles.
%   The two transmits will be identical with the exception of the first
%   having an extra delay of 0.25 wavelengths. TX(1) added delay will be
%   1/4 of the wls, based on the demodulation frequency.
%   TX.delay will be converted to time based on the Trans.frequency
for n = 1:2:2*P(1).numAngles
    TX(n).Steer = [(startAngle+((n-1)/2)*dtheta),0.0];
    TX(n+1).Steer = TX(n).Steer;
    TX(n+1).Delay = computeTXDelays(TX(n+1));
    TX(n).Delay = TX(n+1).Delay + 0.25*(Trans.frequency/31.25);
end
% -- only one TX struct needed for Doppler
k = 2*P(1).numAngles;
TX(k+1).waveform = 1;
TX(k+1).Steer = [dopAngle,0.0];
TX(k+1).Delay = computeTXDelays(TX(k+1));

% Specify Receive structure arrays.
% - We need 2*P.numAngles Receives for every 2D acquisition frame, since two
% acquisitions are needed to produce one interleaved acquisition line of data.
% For the interleaved acquisition approach with 2X interleave, we must take
% into account the actual sample rate at which the digital filters in the
% hardware will be operating: twice the center frequency set by
% Trans.frequency, not the typical 4X factor.  This means that the Nyquist
% limit for the filters will be at Trans.frequency; the higher half of
% the transducer signal frequency spectrum will be folded over the lower
% half due to aliasing.  (The subsequent interleave combination of the two
% acquisition events will unfold this aliasing).  Therefore the filter
% actually used for the Input Filter must be defined as a high-pass
% filter;  the net effect after interleave will be a symmetric bandpass
% filter centered at 31.25 MHz.
% - A Highpass coefficient array has been defined below, representing a
% fractional bandwidth of 100%, relative to Fc at 31.25 MHz.  The coefficient
% array listed below was developed using the matlab fdatool, with a Hamming
% window.
% - Note that for the L28 we would ideally use a bandpass filter centered
% near the transducer's nominal center frequency of 28 MHz, not the 31.25
% MHz forced on the CGD filters.
HighPassCoef100 = [-0.0000    0.0014   -0.0000   -0.0024   -0.0000    0.0046   -0.0000...
                   -0.0081   -0.0000    0.0136   -0.0000   -0.0217   -0.0000    0.0341...
                   -0.0000   -0.0551   -0.0000    0.1009   -0.0000   -0.3169    0.5006];
% For Doppler, we use a 35% BW bandpass filter center at the aliased center frequency for 4/3 sampling.
BandPassCoef = [-0.00021 +0.00000 +0.00058 +0.00000 -0.00027 +0.00000 -0.00250 ...
                +0.00000 +0.01065 +0.00000 -0.02682 +0.00000 +0.05188 +0.00000 ...
                -0.08319 +0.00000 +0.11462 +0.00000 -0.13812 +0.00000 +0.14679];

maxAcqLength = ceil(sqrt(P(1).endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
maxAcqLngthDop = ceil(sqrt(P(2).endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', P(1).startDepth, ...
                        'endDepth', maxAcqLength, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BWI', ...
                        'demodFrequency', 31.25, ...
                        'InputFilter', HighPassCoef100, ...
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, (2*P(1).numAngles+ne)*Resource.RcvBuffer(1).numFrames);
% - Set event specific Receive attributes for each frame.
for i = 1:Resource.RcvBuffer(1).numFrames
    k = (2*P(1).numAngles+ne)*(i-1);
    Receive(k+1).callMediaFunc = 1;
    for j = 1:2:2*P(1).numAngles
        % pairs of acquisition receives to buffer 1 for interleaved acquisition at 1/2 sample rate:
        % First of each pair
        Receive(k+j).framenum = i;
        Receive(k+j).acqNum = j;
        % second of each pair
        Receive(k+j+1).framenum = i;
        Receive(k+j+1).acqNum = j+1;
    end
    for j = (2*P(1).numAngles+1):(2*P(1).numAngles+ne);
        % Doppler acquisition
        Receive(j+k).startDepth = P(2).startDepth;
        Receive(j+k).endDepth = maxAcqLngthDop;
        Receive(j+k).TGC = 2;
        Receive(j+k).framenum = i;
        Receive(j+k).acqNum = j;        % Doppler acqNums continue after 2D
        Receive(j+k).sampleMode = 'BS67BW';
        Receive(j+k).demodFrequency = 26.7858;
        Receive(j+k).InputFilter = BandPassCoef;
    end
end

% Specify TGC Waveform structures.
% - 2D TGC
TGC.CntrlPts = [0,141,275,404,510,603,702,782];
TGC(1).rangeMax = P(1).endDepth;
TGC(1).Waveform = computeTGCWaveform(TGC(1));
% - Doppler TGC
TGC(2).CntrlPts =[ 0 272 662 662 662 662 662 662];
TGC(2).rangeMax = P(2).endDepth;
TGC(2).Waveform = computeTGCWaveform(TGC(2));

% Specify Recon structure arrays.
% - We need two Recon structures, one for 2D, one for Doppler. These will be referenced in the same
%   event, so that they will use the same (most recent) acquisition frame.
na = P(1).numAngles;
Recon = repmat(struct('senscutoff', 0.7, ...
               'pdatanum', 1, ...
               'rcvBufFrame', -1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...
               'RINums', 1:na), 1, 2);
% - Set Recon values for 2D frame.
k = na + 1;
% - Set Recon values for Doppler ensemble.
Recon(2).pdatanum = 2;
Recon(2).IntBufDest = [2,1];
Recon(2).ImgBufDest = [0,0];
Recon(2).RINums = (k:(k+ne-1));   % ne ReconInfos needed for Doppler ensemble.

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
        ReconInfo(j).txnum = 2*j-1;
        ReconInfo(j).rcvnum = 2*j-1;
    end
    ReconInfo(na).mode = 'accumIQ_replaceIntensity'; % accum and detect
else
    ReconInfo(1).mode = 'replaceIntensity';
end

%  - ReconInfos for Doppler ensemble.
k = 2*na;
for j = 1:ne
    ReconInfo(na+j).mode = 'replaceIQ';
    ReconInfo(na+j).txnum = k+1;
    ReconInfo(na+j).rcvnum = k + j;
    ReconInfo(na+j).pagenum = j;
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
                         'interpMethod','4pt',...  %method of interp. (1=4pt)
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',40,...
                         'mappingMode','lowerHalf',...
                         'display',0,...      % don't display image after processing
                         'displayWindow',1};

Process(2).classname = 'Doppler';                   % process structure for 1st Doppler ensemble
Process(2).method = 'computeCFIFreqEst';
Process(2).Parameters = {'IntBufSrc',[2,1],...          % number of (Inter?)buffer to process.
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
                         'interpMethod','4pt',...  %method of interp. (1=4pt)
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',40,...
                         'mappingMode','upperHalf',...
                         'threshold',cpt,...
                         'display',1,...      % don't display image after processing
                         'displayWindow',1};

% Specify SeqControl structure arrays.
% -- Time between 2D flash angle acquisitions
SeqControl(1).command = 'timeToNextAcq';
SeqControl(1).argument = 100;  % time in usec
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
% -- set receive profile 1 for 2D
SeqControl(8).command = 'setRcvProfile';
SeqControl(8).argument = 1;
% -- set receive profile 2 for Doppler
SeqControl(9).command = 'setRcvProfile';
SeqControl(9).argument = 2;
% - The remainder of the SeqControl structures are defined dynamically in the sequence events.
%   The variable nsc keeps track of the next SeqControl structure index.
nsc = 10;  % next SeqControl number

% Specify Event structure arrays.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    % Acquire 2D frame
    k = (2*P(1).numAngles+ne)*(i-1);
    for j = 1:2:2*P(1).numAngles  % Acquire frame
        Event(n).info = 'Acquire angle, 1st half interleave samples';
        Event(n).tx = j;   % use steered TX structure.
        Event(n).rcv = k+j;
        Event(n).recon = 0;      % no reconstruction.
        Event(n).process = 0;    % no processing
        Event(n).seqControl = 1;
        n = n+1;

        Event(n).info = 'Acquire angle, 2nd half interleave samples';
        Event(n).tx = j+1;   % use TX structure with interleave offset delay added.
        Event(n).rcv = k+j+1;
        Event(n).recon = 0;      % no reconstruction.
        Event(n).process = 0;    % no processing
        Event(n).seqControl = 1;
        n = n+1;
    end
    Event(n-1).seqControl = [2,3];   % replace last 2D acquisition Event's seqControl

    % switch to the Doppler receive profile
    Event(n).info = 'Set Doppler rcv profile';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = 9;
    n = n+1;

    % Acquire Doppler ensemble.
    for j = (2*P(1).numAngles+1):(2*P(1).numAngles+ne)
        Event(n).info = 'Acquire Doppler ensemble';
        Event(n).tx = 2*P(1).numAngles+1;      % use next TX structure after 2D.
        Event(n).rcv = k+j;
        Event(n).recon = 0;      % no reconstruction.
        Event(n).process = 0;    % no processing
        Event(n).seqControl = 4; % TTNA for Doppler
        n = n+1;
    end
    Event(n-1).seqControl = [5,6,nsc]; % replace last Doppler acquisition Event's seqControl
      SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
      nsc = nsc+1;

    Event(n).info = 'recons and 2D process';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = [1,2];  % reconstruction for both 2D and Doppler
    Event(n).process = 1;    % process 2D
    Event(n).seqControl = 8;  % select the 2D receive profile for next frame
    n = n+1;

    Event(n).info = 'Doppler processing';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 0;      % no reconstruction
    Event(n).process = 2;    % IQ to estimate processing
    Event(n).seqControl = 0; % no seqCntrl
    n = n+1;

    Event(n).info = 'Doppler image display';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 0;      % no reconstruction
    Event(n).process = 3;    % display processing
    Event(n).seqControl = 0; % no seqCntrl
    if floor(i/3) == i/3     % Exit to Matlab every 3rd frame
        Event(n).seqControl = nsc;
        SeqControl(nsc).command = 'returnToMatlab';
        nsc = nsc+1;
    end
    n = n+1;
end

Event(n).info = 'Jump back';
Event(n).tx = 0;        % no TX
Event(n).rcv = 0;       % no Rcv
Event(n).recon = 0;     % no Recon
Event(n).process = 0;
Event(n).seqControl = 7;


% User specified UI Control Elements
% - Sensitivity Cutoff
UI(1).Control =  {'UserB7','Style','VsSlider','Label','Sens. Cutoff',...
                  'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                  'SliderStep',[0.025,0.1],'ValueFormat','%1.3f'};
UI(1).Callback = text2cell('%SensCutoffCallback');

% - Doppler Mode Button Group
UI(2).Control = {'UserB4','Style','VsButtonGroup','Title','Doppler Mode','NumButtons',2,'Labels',{'Velocity','Power'}};
UI(2).Callback = text2cell('%-UI#2Callback');

% - Doppler Power Threshold Slider
UI(3).Control = {'UserB3','Style','VsSlider','Label','DopPwrThres','SliderMinMaxVal',[0.0,1.0,pwrThres],...
                 'SliderStep',[0.02,0.1],'ValueFormat','%3.2f'};
UI(3).Callback = text2cell('%-UI#3Callback');

% - Color Priority Threshold Slider
UI(4).Control = {'UserB2','Style','VsSlider','Label','Color Priority','SliderMinMaxVal',[0,255,cpt],...
                 'SliderStep',[1/255,0.1],'ValueFormat','%3.0f'};
UI(4).Callback = text2cell('%-UI#4Callback');

% - Color Persistence Slider
UI(5).Control = {'UserB1','Style','VsSlider','Label','Color Persistence','SliderMinMaxVal',[0,100,persf],...
                 'SliderStep',[1/100,0.1],'ValueFormat','%3.0f'};
UI(5).Callback = text2cell('%-UI#5Callback');

% - Doppler Waveform Slider
UI(6).Control = {'UserC3','Style','VsSlider','Label','Doppler Waveform','SliderMinMaxVal',[1,2,1],...
                 'SliderStep',[1, 1],'ValueFormat','%3.0f'};
UI(6).Callback = text2cell('%-UI#6Callback');

% - Doppler Zin Slider
UI(7).Control = {'UserC2','Style','VsSlider','Label','Doppler LNA Zin','SliderMinMaxVal',[0,31,0],...
                 'SliderStep',[1/31, 1/31],'ValueFormat','%3.0f'};
UI(7).Callback = text2cell('%-UI#7Callback');

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 3;

% Save all the structures to a .mat file.
save('MatFiles/L35-16vXFlashDoppler');

%% Automatic VSX Execution:
% Uncomment the following line to automatically run VSX every time you run
% this SetUp script (note that if VSX finds the variable 'filename' in the
% Matlab workspace, it will load and run that file instead of prompting the
% user for the file to be used):

% filename = 'L35-16vXFlashDoppler';  VSX;

return


% **** Callback routines to be encoded by text2cell function. ****
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

%-UI#2Callback - Doppler mode change
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
if ishandle(hPTool),
    posPTool = get(hPTool,'position');
    PTool;
    set(findobj('tag','ProcessTool'),'position',posPTool);
end

% If ColorMapTool is open, close it.
hCMTool = findobj('tag','ColorMapTool');
if ishandle(hCMTool),
    delete(hCMTool);
    set(findobj('tag','toolsMenu'),'Value',1); % set tools selection back to none
end

return
%-UI#2Callback


%-UI#3Callback - Doppler Power change
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
%-UI#3Callback

%-UI#4Callback - Color Priority change
% Set the value in the Process structure for use in cineloop playback.
Process = evalin('base','Process');
for k = 1:2:length(Process(3).Parameters)
    if strcmp(Process(3).Parameters{k},'threshold'), Process(3).Parameters{k+1} = UIValue; end
end
assignin('base','Process',Process);
% Set Control.Command to set Image.threshold.
Control = evalin('base','Control');
Control.Command = 'set&Run';
Control.Parameters = {'Process',3,'threshold',UIValue};
assignin('base','Control', Control);
%-UI#4Callback

%-UI#5Callback - Color Persistence change
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
%-UI#5Callback

%-UI#6Callback - Doppler waveform select
TX = evalin('base','TX');
wvsel = 2*round(UIValue) - 1; % map to select either 1 or 3
TX(17).waveform = wvsel;
assignin('base','TX',TX);
% Set Control.Command to update TX
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'TX'};
assignin('base','Control', Control);
%-UI#6Callback

%-UI#7Callback - Doppler rcv profile Zin select
RcvProfile = evalin('base','RcvProfile');
RcvProfile(2).LnaZinSel = round(UIValue);
RcvProfile = computeRcvProfile( RcvProfile );
assignin('base','RcvProfile',RcvProfile);
% Set Control.Command to update RcvProfile.
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'RcvProfile'};
assignin('base','Control', Control);
%-UI#7Callback
