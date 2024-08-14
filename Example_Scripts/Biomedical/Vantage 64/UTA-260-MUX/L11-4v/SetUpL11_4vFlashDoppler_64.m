% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: SetUpL11_4vFlashDoppler_64.m -  Example of 3-1 synthetic aperture doppler
%                                     imaging woth flash transmits
% Description:
%   Sequence programming file for L11-4v Linear array, using 3-1
%   synthetic aperture plane wave transmits and receive acquisitions on
%   64 channels system. 64 transmit channels and 42 or 44 receive channels
%   are active and positioned as follows (each char represents 2 elements)
%   for each of the 3 synthetic apertures over na steering angles.
%
%   Element Nos.                                                               1
%                               3   4           6         8    9               2
%               1               3   2           5         6    6               8
%   Aperture 1: |               |   |           |         |    |               |
%               tttttttttttttttttttttttttttttttt--------------------------------
%               rrrrrrrrrrrrrrrrrrrrr-------------------------------------------
%               |               |    |          |         |    |               |
%   Aperture 2: |               |    |          |         |    |               |
%               ----------------tttttttttttttttttttttttttttttttt----------------
%               ---------------------rrrrrrrrrrrrrrrrrrrrrr---------------------
%               |               |    |          |         |    |               |
%   Aperture 3: |               |    |          |         |    |               |
%               --------------------------------tttttttttttttttttttttttttttttttt
%               -------------------------------------------rrrrrrrrrrrrrrrrrrrrr
%               |               |    |          |         |    |               |
%
%   The receive data from each of these apertures are stored under different
%   acqNums in the Receive buffer. The reconstruction sums the IQ data from the
%   3 apertures and multiple angles and computes intensity values to produce
%   the full frame.
%
%   For Doppler, we use a 64 element aperture starting at P(2).firstElem (e.g. 20)
%
%   Element Nos.                                                                1
%                         2                     6         8                     2
%               1         0                     5         3                     8
%   Aperture 1: |         |                     |         |                     |
%               ----------ttttttttttttttttttttttttttttttttt----------------------
%               ----------rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr----------------------
%
% Notes:
%   - for the Doppler ensemble, ne pulses are transmitted at a steering
%     angle of dopAngle radians.
%   - this version uses TPC profiles to provide a seperate high voltage
%     level for Doppler, and RcvProfiles to set different AFE preamp
%     settings for 2D and Doppler.
%   - processing is asynchronous with respect to acquisition.
%
% Last update:
%   08/30/2019 - Updated for 4.1.1 release, VTS-1401 bug fix.

clear all

% Set 2D parameters - P(1) is used for 2D, P(2) is used for Doppler
P(1).startDepth = 5;
P(1).endDepth = 192;   % Acquisition depth in wavelengths
na = 7;      % Set na = number of flash angles for 2D.
if (na > 1)
    dtheta2D = (30*pi/180)/(na-1);
    startAngle = -15*pi/180;
else
    dtheta2D = 0;
    startAngle=0;
end % set dtheta2D to range over +/- 15 degrees.

% Set Doppler parameters in P(2)
ne = 14;     % Set ne = number of acquisitions in Doppler ensemble.
P(2).startDepth = 5;
P(2).endDepth = 160;   % Acquisition depth in wavelengths
P(2).firstElem = 20; % first element in 64 element Doppler aperture
P(2).dopAngle = 12 * pi/180;
P(2).dopPRF = 3.0e+03; % Doppler PRF in Hz.
P(2).pwrThres = 0.5;

% Specify system parameters.
Resource.Parameters.numTransmit = 64;          % number of transmit channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.
Resource.System.SoftwareVersion = [3 5 0]; % Minimum software release for this script.
Resource.System.UTA = '260-MUX'; % This script requires the 260-MUX UTA.

% Specify Trans structure array.
Trans.name = 'L11-4v';
Trans.units = 'wavelengths'; % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);  % computeTrans is used for known transducers.
Trans.maxHighVoltage = 50;    % set a reasonable high voltage limit.
Trans = computeUTAMux64(Trans); % add the HVMux structure for UTA module

% Specify PData structure. 2D uses full PData array while Doppler uses parallelogram region.
PData(1).PDelta = [Trans.spacing, 0, 0.5];
PData(1).Size(1) = ceil((P(1).endDepth-P(1).startDepth)/PData(1).PDelta(3)); % rows
PData(1).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1)); % cols
PData(1).Size(3) = 1;             % single image page
PData(1).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P(1).startDepth]; % x,y,z of uppr lft crnr.
% - 2D Region definition
PData(1).Region = struct('Shape',struct('Name','PData'));
PData(1).Region = computeRegions(PData(1));
% Specify Doppler PData structure. The PData shape is the same as for 2D, but PDelta is different.
PData(2).PDelta = [Trans.spacing, 0, 1.0];
PData(2).Size(1) = ceil((P(1).endDepth-P(1).startDepth)/PData(2).PDelta(3)); % rows
PData(2).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(2).PDelta(1)); % cols
PData(2).Size(3) = 1;             % single image page
PData(2).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P(1).startDepth]; % x,y,z of uppr lft crnr.
% - Doppler PData region definition
PData(2).Region = struct('Shape',struct('Name','Parallelogram', ...
                               'Position', [PData(2).Origin(1)+(P(2).firstElem+31)*Trans.spacing,0,P(2).startDepth], ...
                               'width', 64*Trans.spacing, ...
                               'height', P(2).endDepth - P(2).startDepth, ...
                               'angle', P(2).dopAngle));
PData(2).Region = computeRegions(PData(2));

% Specify Media object.
pt1;
Media.attenuation = 0; %-0.5;
Media.function = 'movePoints';

% Specify Resources.
% - RcvBuffer(1) is for both 2D and Doppler acquisitions.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 4096*(na*3 + ne);
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numTransmit;
Resource.RcvBuffer(1).numFrames = 20;           % 20 frames allocated for RF acqusitions.
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
Resource.DisplayWindow(1).Title = 'L11-4vFlashDoppler_64';
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).pdelta = 0.4;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)]; % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).numFrames = 20;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = grayscaleCFImap;
Resource.DisplayWindow(1).splitPalette = 1;

% ------Specify structures used in Events------
% Specify Transmit waveform structures.
% - 2D transmit waveform
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,0.67,2,1];
% - Doppler transmit waveform, transmit frequency should be equivalent to
% supported demodFrequency
TW(2).type = 'parametric';
TW(2).Parameters = [6.25,0.67,6,1];

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'focus', 0.0, ...
                   'Steer', [0.0,0.0], ...
                   'Apod', ones(1,Resource.Parameters.numTransmit), ...
                   'Delay', zeros(1,Resource.Parameters.numTransmit)), 1, (3*na)+1); % na TXs for 2D + 1 for Doppler

% - Set event specific TX attributes.
angle = startAngle;
for j = 1:3:3*na
    TX(j).aperture = 1;% Use the tx aperture that starts at element 1.
    TX(j).Steer = [angle,0.0];
    TX(j).Delay = computeTXDelays(TX(j));
    TX(j+1).aperture = 33;% Use the tx aperture that starts at element 33.
    TX(j+1).Steer = [angle,0.0];
    TX(j+1).Delay = computeTXDelays(TX(j+1));
    TX(j+2).aperture = 65; % Use the tx aperture that starts at element 65.
    TX(j+2).Steer = [angle,0.0];
    TX(j+2).Delay = computeTXDelays(TX(j+2));
    angle = angle + dtheta2D;
end

% -- only one TX struct needed for Doppler
k = 3*na + 1;
TX(k).aperture = 20;
TX(k).waveform = 2;
TX(k).Steer = [P(2).dopAngle,0.0];
TX(k).Delay = computeTXDelays(TX(k));

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
TGC.CntrlPts = [0,345,530,711,816,858,855,884];
TGC(1).rangeMax = P(1).endDepth;
TGC(1).Waveform = computeTGCWaveform(TGC(1));
% - Doppler TGC
TGC(2).CntrlPts =[ 0 310 488 587 587 587 587 538];
TGC(2).rangeMax = P(2).endDepth;
TGC(2).Waveform = computeTGCWaveform(TGC(2));

% Specify Receive structure arrays.
%   We need to acquire all the 2D and Doppler data within a single RcvBuffer frame.  This allows
%   the transfer-to-host DMA after each frame to transfer a large amount of data, improving throughput.
% - We need 3*na Receives for a 2D frame and ne Receives for a Doppler frame.
maxAcqLngth2D = sqrt(P(1).endDepth^2 + (Trans.numelements*Trans.spacing)^2);
a = 64*Trans.spacing;
b = P(2).endDepth/cos(P(2).dopAngle);
maxAcqLngthDop =  sqrt(a^2 + b^2 - 2*a*b*cos(P(2).dopAngle+(pi/2)));
Receive = repmat(struct('Apod', zeros(1,Resource.Parameters.numTransmit), ...
                        'startDepth', P(1).startDepth, ...
                        'endDepth', maxAcqLngth2D, ...
                        'aperture', 1, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'InputFilter', [], ...
                        'sampleMode', 'NS200BW', ... % 200% Bandwidth for 2D
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, (3*na+ne)*Resource.RcvBuffer(1).numFrames);
% Define bandpass filter for Doppler - 6.25 center frequency, 30% bandwidth.
BPF30 = [-0.00061 +0.00000 +0.00275 +0.00000 -0.00665 +0.00000 +0.01047 ...
         +0.00000 -0.00912 +0.00000 -0.00443 +0.00000 +0.03537 +0.00000 ...
         -0.08228 +0.00000 +0.13477 +0.00000 -0.17654 +0.00000 +0.19257];

% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    k = (3*na + ne)*(i-1); % k keeps track of Receive index increment per frame.
    % - Set attributes for each frame.
    Receive(k+1).callMediaFunc = 1;
    for j = 1:3:3*na  % acquisitions for 2D (3 acqs for each angle).
        % -- 2D acquisition, 1st synthetic aperture acquisition for aperture 1.
        Receive(j+k).Apod(1:42) = 1.0;
        Receive(j+k).aperture = 1;
        Receive(j+k).framenum = i;
        Receive(j+k).acqNum = j;      % three acquisitions per frame
        % -- 2D acquisition, 2nd synthetic aperture acquisition for aperture 33.
        Receive(j+k+1).Apod(11:54) = 1.0;
        Receive(j+k+1).aperture = 33;
        Receive(j+k+1).framenum = i;
        Receive(j+k+1).acqNum = j+1;  % three acquisitions per frame
        % -- 2D acquisition, 3rd synthetic aperture acquisition for aperture 65.
        Receive(j+k+2).Apod(23:64) = 1.0;
        Receive(j+k+2).aperture = 65;
        Receive(j+k+2).framenum = i;
        Receive(j+k+2).acqNum = j+2;  % three acquisitions per frame
    end
    for j = (3*na+1):(3*na+ne)
        % Doppler acquisition
        Receive(j+k).Apod(:) = 1.0;
        Receive(j+k).aperture = 20;  % use aperture starting at P(2).firstElem
        Receive(j+k).startDepth = P(2).startDepth;
        Receive(j+k).endDepth = maxAcqLngthDop;
        Receive(j+k).InputFilter = BPF30;
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
Recon(1).RINums = 1:3*na;  % 3*na ReconInfos needed for na angles and 3 apertures.
k = 3*na + 1;
% - Set Recon values for Doppler ensemble.
Recon(2).pdatanum = 2;
Recon(2).IntBufDest = [2,1];
Recon(2).ImgBufDest = [0,0];
Recon(2).RINums = k:(k+ne-1);   % ne ReconInfos needed for Doppler ensemble.

% Define ReconInfo structures.
% - For 2D, we need 3*na ReconInfo structures for na steering angles and 3 apertures.
% - For Doppler, we need ne ReconInfo structures.
ReconInfo = repmat(struct('mode', 'accumIQ', ...    % accumulate IQ data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'pagenum',1, ...
                   'scaleFactor', 1.0, ...
                   'regionnum', 1), 1, 3*na + ne);
% - Set specific ReconInfo attributes.
%   - ReconInfos for 2D frame.
if na>1
    ReconInfo(1).mode = 'replaceIQ'; % replace IQ data
    for j = 1:3:3*na
        ReconInfo(j).txnum = j;
        ReconInfo(j).rcvnum = j;
        ReconInfo(j+1).txnum = j+1;
        ReconInfo(j+1).rcvnum = j+1;
        ReconInfo(j+2).txnum = j+2;
        ReconInfo(j+2).rcvnum = j+2;
    end
    ReconInfo(3*na).mode = 'accumIQ_replaceIntensity'; % accum and detect
else
    ReconInfo(1).mode = 'replaceIntensity'; % accum and detect
end
%  - ReconInfos for Doppler ensemble.
k = 3*na;
for j = 1:ne
%   ReconInfo(k+j).mode = 'replaceIQ_normalize';
    ReconInfo(k+j).mode = 'replaceIQ';
    ReconInfo(k+j).txnum = k + 1;
    ReconInfo(k+j).rcvnum = k + j;
    ReconInfo(k+j).pagenum = j;
    ReconInfo(k+j).scaleFactor = 2.0;
end

% Specify Process structure arrays.
cpt = 50;  % define here so we can use in UIControl below
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
                         'compressFactor',45,...
                         'mappingMethod','lowerHalf',...
                         'display',0,...      % don't display image after processing
                         'displayWindow',1};
%
% Process(2).classname = 'Doppler';                   % process structure for 1st Doppler ensemble
% Process(2).method = 'computeCFIFreqEst';
% Process(2).Parameters = {'IntBufSrc',[2,1],...      % number of Interbuffer to process.
%                          'SrcPages',[3,ne-2],...    % start frame number in source buffer
%                          'ImgBufDest',[2,-1],...
%                          'pdatanum',2,...           % number of PData structure
%                          'prf',P(2).dopPRF,...      % Doppler PRF in Hz
%                          'wallFilter','regression',...  % 1 -> quadratic regression
%                          'pwrThreshold',P(2).pwrThres,...
%                          'maxPower',50,...
%                          'estAvg',1,...
%                          'postFilter',1};

Process(2).classname = 'Doppler';                   % process structure for 1st Doppler ensemble
Process(2).method = 'computeCFIFreqEst';
Process(2).Parameters = {'IntBufSrc',[2,1],...      % number of Interbuffer to process.
                         'SrcPages',[3,ne-2],...    % start frame number in source buffer
                         'ImgBufDest',[2,-1],...
                         'pdatanum',2,...           % number of PData structure
                         'prf',P(2).dopPRF,...      % Doppler PRF in Hz
                         'wallFilter','regression',...  % 1 -> quadratic regression
                         'pwrThreshold',P(2).pwrThres,...
                         'maxPower',50};

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
                         'compressFactor',45,...
                         'mappingMethod','upperHalf',...
                         'threshold',cpt,...
                         'display',1,...      % display image after processing
                         'displayWindow',1};

% Specify SeqControl structure arrays.
% -- Time between 2D flash angle acquisitions
SeqControl(1).command = 'timeToNextAcq';
SeqControl(1).argument = 200;
% -- Change to Profile 2 (Doppler)
SeqControl(2).command = 'setTPCProfile';
SeqControl(2).condition = 'next';
SeqControl(2).argument = 2;
% -- Time between 2D acquisition and Doppler ensemble. Set to allow time for profile change.
SeqControl(3).command = 'timeToNextAcq';
SeqControl(3).argument = 7000; % time in usec
% -- PRF for Doppler ensemble
SeqControl(4).command = 'timeToNextAcq';
SeqControl(4).argument = round(1/(P(2).dopPRF*1e-06)); % (for 3KHz dopPRF & 14 ensemble = 4.7 msecs)
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
    k = (3*na+ne);
    % Acquire 2D frame
    for j = 1:3:3*na
        Event(n).info = '1st aperture.';
        Event(n).tx = j;
        Event(n).rcv = k*(i-1)+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 1;
        if j == 1
            Event(n).seqControl = [1,8];
        end
        n = n+1;

        Event(n).info = '2nd aperture.';
        Event(n).tx = j+1;
        Event(n).rcv = k*(i-1)+j+1;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 1;
        n = n+1;

        Event(n).info = '3rd aperture.';
        Event(n).tx = j+2;
        Event(n).rcv = k*(i-1)+j+2;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 1;
        n=n+1;
    end

    Event(n-1).seqControl = [2,3];   % replace last 2D acquisition Event's seqControl

    % Acquire Doppler ensemble.
    for j = (3*na+1):(3*na+ne)   % Acquire Doppler ensemble.
        Event(n).info = 'Doppler ensemble';
        Event(n).tx = (3*na)+1;
        Event(n).rcv = (3*na+ne)*(i-1)+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 4;
        if j == 3*na+1
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
% - Sensitivity Cutoff
UI(1).Control =  {'UserB7','Style','VsSlider','Label','Sens. Cutoff',...
                  'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                  'SliderStep',[0.025,0.1],'ValueFormat','%1.3f'};
UI(1).Callback = text2cell('%SensCutoffCallback');

% - Doppler Mode Button Group
UI(2).Control = {'UserB4','Style','VsButtonGroup','Title','Doppler Mode','NumButtons',2,'Labels',{'Velocity','Power'}};
UI(2).Callback = text2cell('%-UI#2Callback');

% - Doppler Power Threshold Slider
UI(3).Control = {'UserB3','Style','VsSlider','Label','DopPwrThres','SliderMinMaxVal',[0.0,1.0,P(2).pwrThres],...
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

% - MinPower factor
UI(6).Control = {'UserC4','Style','VsSlider','Label','minPower','SliderMinMaxVal',[0,1000,500],...
                 'SliderStep',[1/100,0.1],'ValueFormat','%4.0f'};
UI(6).Callback = text2cell('%-UI#6Callback');

% - MaxPower factor
UI(7).Control = {'UserC3','Style','VsSlider','Label','maxPower','SliderMinMaxVal',[20,100,50],...
                 'SliderStep',[1/100,0.1],'ValueFormat','%3.0f'};
UI(7).Callback = text2cell('%-UI#7Callback');

% - Power compression factor
UI(8).Control = {'UserC2','Style','VsSlider','Label','pwrCompression','SliderMinMaxVal',[0.1,1.0,0.5],...
                 'SliderStep',[1/100,0.1],'ValueFormat','%1.2f'};
UI(8).Callback = text2cell('%-UI#8Callback');

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 3;

% Save all the structures to a .mat file.
save('MatFiles/L11-4vFlashDoppler_64');
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

% If PTool window is open, adjust persistLevel1 in Process(3)
hPTool = findobj('tag','ProcessTool');
if ishandle(hPTool)
    hPNum = findobj('tag','processNum');
    if isequal(get(findobj('tag','processNum'),'Value'),3)
        set(findobj('tag','persistSlider1'),'Value',UIValue);
        set(findobj('tag','persistValue1'),'String',num2str(UIValue));
    end
end
return
%-UI#5Callback

%-UI#6Callback
% MinPower
% Set the value in the Doppler Process structure.
Process = evalin('base','Process');
flag = 0;
k = length(Process(2).Parameters);
for j = 1:2:k
    if strcmp(Process(2).Parameters{j},'minPower')
        Process(2).Parameters{j+1} = UIValue;
        flag = 1;
        break;
    end
end
% In case attribute was not defined.
if flag == 0
    Process(2).Parameters{k+1} = 'minPower';
    Process(2).Parameters{k+2}= UIValue;
end
assignin('base','Process',Process);
% Set Control.Command to set Doppler attribute.
Control = evalin('base','Control');
Control.Command = 'set&Run';
Control.Parameters = {'Process',2,'minPower',UIValue};
assignin('base','Control', Control);
%-UI#6Callback

%-UI#7Callback
% Maxpower
% Set the value in the Doppler Process structure.
Process = evalin('base','Process');
flag = 0;
k = length(Process(2).Parameters);
for j = 1:2:k
    if strcmp(Process(2).Parameters{j},'maxPower')
        Process(2).Parameters{j+1} = UIValue;
        flag = 1;
        break;
    end
end
% In case attribute was not defined.
if flag == 0
    Process(2).Parameters{k+1} = 'maxPower';
    Process(2).Parameters{k+2}= UIValue;
end
assignin('base','Process',Process);
% Set Control.Command to set Doppler attribute.
Control = evalin('base','Control');
Control.Command = 'set&Run';
Control.Parameters = {'Process',2,'maxPower',UIValue};
assignin('base','Control', Control);
%-UI#7Callback

%-UI#8Callback
% Power Compression
% Set the pwrCompression value in the Doppler Process structure.
Process = evalin('base','Process');
flag = 0;
k = length(Process(2).Parameters);
for j = 1:2:k
    if strcmp(Process(2).Parameters{j},'pwrCompression')
        Process(2).Parameters{j+1} = UIValue;
        flag = 1;
        break;
    end
end
% In case attribute was not defined.
if flag == 0
    Process(2).Parameters{k+1} = 'pwrCompression';
    Process(2).Parameters{k+2}= UIValue;
end
assignin('base','Process',Process);
% Set Control.Command to set Doppler attribute.
Control = evalin('base','Control');
Control.Command = 'set&Run';
Control.Parameters = {'Process',2,'pwrCompression',UIValue};
assignin('base','Control', Control);
%-UI#8Callback
