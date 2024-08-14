% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: SetUpL22_14vFLashDoppler_32LE.m - Example of plane wave Doppler imaging, using FlashAngles for B-mode image
%
% Description:
%   Sequence programming file for L22-14v Linear array, using plane wave
%   transmits for 2D (B-mode) and Doppler imaging.
%   - for the 2D (B-mode) imaging, na steered angles are used.
%   - for the Doppler ensemble, ne pulses are transmitted at a steering angle of dopAngle radians.
%   - this version uses TPC profiles to provide a seperate high voltage level for Doppler.
%   - processing is asynchronous with respect to acquisition.
%   - This script uses 4X sampling with A/D sample rate of 62.5 MHz for a
%     15.625 MHz processing center frequency.  Transmit is at 17.8 MHz and
%     receive bandpass filter has been shifted to 18 MHz center frequency,
%     13.9 MHz bandwidth (-3 dB) to support the 12 MHz bandwidth of the
%     L22-14v (18 MHz center frequency, ~67% bandwidth).
%
%   For Imaging:
%   Sequence programming file for L22-4v Linear array, using 4-1
%   synthetic aperture plane wave transmits and receive acquisitions on
%   32 channels system. 64 transmit channels and 32 receive channels
%   are active and positioned as follows (each char represents 2 elements)
%   for each of the 3 synthetic apertures.
%
%   Element Nos.                                                               1
%                               3   4           6         8    9               2
%               1               3   2           5         6    6               8
%   Aperture 1: |               |   |           |         |    |               |
%               tttttttttttttttttttttttttttttttt--------------------------------
%               rrrrrrrrrrrrrrrrr-----------------------------------------------
%   Aperture 2: |               |   |           |         |    |               |
%               |               |   |           |         |    |               |
%               tttttttttttttttttttttttttttttttt--------------------------------
%               ----------------rrrrrrrrrrrrrrrr--------------------------------
%               |               |    |          |         |    |               |
%   Aperture 3: |               |    |          |         |    |               |
%               --------------------------------tttttttttttttttttttttttttttttttt
%               --------------------------------rrrrrrrrrrrrrrr-----------------
%               |               |    |          |         |    |               |
%   Aperture 3: |               |    |          |         |    |               |
%               --------------------------------tttttttttttttttttttttttttttttttt
%               -----------------------------------------------rrrrrrrrrrrrrrrrr
%               |               |    |          |         |    |               |
%
%   The receive data from each of these apertures are stored under
%   different acqNums in the Receive buffer. The reconstruction sums the
%   IQ data from the 4 apertures and multiple angles and computes intensity values to produce
%   the full frame. Processing is asynchronous with respect to acquisition.
%
%   For Doppler
%
%   Element Nos.                                                               1
%                               3   4           6         8    9               2
%               1               3   2           5         6    6               8
%   Aperture 1: |               |   |           |         |    |               |
%               ----------------tttttttttttttttttttttttttttttttt----------------
%               ----------------------rrrrrrrrrrrrrrrrrrrrr----------------------
%
% Last update:
% 09/03/2018 - script improvement for 32LE

clear all
% P(1) is used for Bmode, P(2) is used for Doppler
P(1).startDepth = 0;   % Acquisition depth in wavelengths
P(1).endDepth = 150;   % This should preferrably be a multiple of 128 samples.
P(2).startDepth = 0;   % Acquisition depth in wavelengths
P(2).endDepth = 128;   % This should preferrably be a multiple of 128 samples.

% Set 2D parameters
na = 7;      % Set na = number of flash angles for 2D.
if (na > 1)
    dtheta2D = (30*pi/180)/(na-1);
    startAngle = -30*pi/180/2;
else
    dtheta2D = 0;
    startAngle=0;
end % set dtheta2D to range over +/- 15 degrees.

% Set Doppler parameters
ne = 14;     % Set ne = number of acquisitions in Doppler ensemble.
dopAngle = 12 * pi/180;
dopPRF = 2.0e+03; % Doppler PRF in Hz.
pwrThres = 0.3;

% Specify system parameters.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.
Resource.System.SoftwareVersion = [3 5 0]; % Minimum software release for this script.
Resource.System.UTA = '260-MUX'; % This script requires the 260-MUX UTA.
Resource.System.transmitChannels = 64;
Resource.System.receiveChannels = 32;

% Specify Trans structure array.
Trans.name = 'L22-14v';
Trans.units = 'wavelengths';
Trans = computeTrans(Trans);
Trans = computeUTAMux64(Trans);
Trans.maxHighVoltage = 25; % mfr data sheet lists 30 Volt limit

% Specify PData structure arrays.
% - 2D PData structure
PData(1).PDelta = [0.4, 0, 0.35];
PData(1).Size(1) = ceil((P(1).endDepth-P(1).startDepth)/PData(1).PDelta(3)); % rows
PData(1).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1)); % cols
PData(1).Size(3) = 1;             % single image page
PData(1).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P(1).startDepth]; % x,y,z of uppr lft crnr.
% - Doppler PData structure
PData(2).PDelta = [0.6, 0, 0.6];
PData(2).Size(1) = ceil((P(2).endDepth-P(2).startDepth)/PData(2).PDelta(3)); % flow window rows
PData(2).Size(2) = ceil((64*Trans.spacing)/PData(2).PDelta(1));  % flow window columns
PData(2).Size(3) = 1;             % single image page
PData(2).Origin = [-Trans.spacing*(64-1)/2,0,P(2).startDepth]; % x,y,z of upper lft crnr.

% Specify Media object.
pt1;
Media.attenuation = -0.5;
Media.function = 'movePoints';

% Specify Resources.
% - RcvBuffer(1) is for both 2D and Doppler acquisitions.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 2048*(3*na + ne);
Resource.RcvBuffer(1).colsPerFrame = Resource.System.receiveChannels;
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
Resource.DisplayWindow(1).Title = 'L22-14v Flash Doppler, 4X sampling at 62.5 MHz';
Resource.DisplayWindow(1).pdelta = 0.3;
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
% Specify Transmit waveform structures.
% - 2D transmit waveform
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,0.67,2,1];
% - Doppler transmit waveform, transmit frequency should be equivalent to
% supported demodFrequency
TW(2).type = 'parametric';
TW(2).Parameters = [15.625,0.67,6,1];

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'focus', 0.0, ...
                   'Steer', [0.0,0.0], ...
                   'Apod', ones(1,Resource.System.transmitChannels), ...
                   'Delay', zeros(1,Resource.System.transmitChannels)), 1, (2*na)+1); % na TXs for 2D + 1 for Doppler
% - Set event specific TX attributes.
angle = startAngle;
for j = 1:2:2*na
    TX(j).aperture = 1;% Use the tx aperture that starts at element 1.
    TX(j).Steer = [angle,0.0];
    TX(j).Delay = computeTXDelays(TX(j));
    TX(j+1).aperture = 65;% Use the tx aperture that starts at element 65.
    TX(j+1).Steer = [angle,0.0];
    TX(j+1).Delay = computeTXDelays(TX(j+1));
    angle = angle + dtheta2D;
end
% -- only one TX struct needed for Doppler
TX((2*na)+1).aperture = 33;
TX((2*na)+1).waveform = 2;
TX((2*na)+1).Steer = [dopAngle,0.0];
TX((2*na)+1).Delay = computeTXDelays(TX((2*na)+1));

% Specify TPC structures.
TPC(1).name = '2D';
TPC(1).maxHighVoltage = 25;
TPC(2).name = 'Doppler';
TPC(2).maxHighVoltage = 25;

RcvProfile(1).LnaGain = 18; % Profile used for imaging
RcvProfile(2).LnaGain = 24; % Profile used for Doppler
RcvProfile(2).LnaZinSel = 31; % Force high-Z state for best Doppler sensitivity

% Specify TGC Waveform structures.
% - Doppler TGC
TGC(1).CntrlPts = [200,511,716,920,1023,1023,1023,1023]; %[400,625,750,825,825,875,875,875];
TGC(1).rangeMax = P(2).endDepth;
TGC(1).Waveform = computeTGCWaveform(TGC(1));
% - 2D TGC
TGC(2).CntrlPts = [200,511,716,920,1023,1023,1023,1023]; %[400,500,625,700,750,800,850,950];
TGC(2).rangeMax = P(1).endDepth;
TGC(2).Waveform = computeTGCWaveform(TGC(2));

% Specify Receive structure arrays.
%   We need to acquire all the 2D and Doppler data within a single RcvBuffer frame.  This allows
%   the transfer-to-host DMA after each frame to transfer a large amount of data, improving throughput.
% - We need na Receives for a 2D frame and ne Receives for a Doppler frame.

% sampling center frequency is 15.625, but we want the bandpass filter
% centered on the actual transducer center frequency of 18 MHz with 67%
% bandwidth, or 12 to 24 MHz.  Coefficients below were set using
% "G3_BPFdevelopment" with normalized cf=1.15 (18 MHz), bw=0.85,
% xsn wdth=0.41 resulting in -3 dB 0.71 to 1.6 (11.1 to 25 MHz), and
% -20 dB 0.57 to 1.74 (8.9 to 27.2 MHz)
%
BPF1 = [ -0.00009 -0.00128 +0.00104 +0.00085 +0.00159 +0.00244 -0.00955 ...
         +0.00079 -0.00476 +0.01108 +0.02103 -0.01892 +0.00281 -0.05206 ...
         +0.01358 +0.06165 +0.00735 +0.09698 -0.27612 -0.10144 +0.48608 ];

% For Doppler, use narrow bandwidth coefficients (50% BW) centered at
% 15.625; this is a copy of update function's default coef array for 1
% samples per wave
BPFDop = [ -0.00162 +0.00000 +0.00568 +0.00000 -0.01065 +0.00000 +0.01349 ...
           +0.00000 -0.00858 +0.00000 -0.00955 +0.00000 +0.04312 +0.00000 ...
           -0.08841 +0.00000 +0.13550 +0.00000 -0.17130 +0.00000 +0.18463 ];

maxAcqLngth2D = sqrt(P(1).endDepth^2 + (Trans.numelements*Trans.spacing)^2) - P(1).startDepth;
maxAcqLngthDop =  sqrt(P(2).endDepth^2 + (96*Trans.spacing)^2) - P(2).startDepth;
wl4sPer128 = 128/(4*2);  % wavelengths in a 128 sample block for 4 smpls per wave round trip.
wl2sPer128 = 128/(2*2);  % wavelengths in a 128 sample block for 2 smpls per wave round trip.
Receive = repmat(struct('Apod', zeros(1,Resource.System.transmitChannels), ...
                        'startDepth', P(1).startDepth, ...
                        'endDepth', P(1).startDepth + wl4sPer128*ceil(maxAcqLngth2D/wl4sPer128), ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ... % 200% Bandwidth for 2D
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, (4*na+ne)*Resource.RcvBuffer(1).numFrames);

% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    k = (4*na + ne)*(i-1); % k keeps track of Receive index increment per frame.
    % - Set attributes for each frame.
    Receive(k+1).callMediaFunc = 1;
    for j = 1:4:4*na  % acquisitions for 2D
         % -- 2D acquisition, 1st synthetic aperture acquisition for aperture 1.
        Receive(j+k).Apod(1:32) = 1.0;
        Receive(j+k).aperture = 1;
        Receive(j+k).framenum = i;
        Receive(j+k).acqNum = j;      % 4 acquisitions per frame
        % -- 2D acquisition, 2nd synthetic aperture acquisition for aperture 33.
        Receive(j+k+1).Apod(33:64) = 1.0;
        Receive(j+k+1).aperture = 1;
        Receive(j+k+1).framenum = i;
        Receive(j+k+1).acqNum = j+1;  % 4 acquisitions per frame
        % -- 2D acquisition, 3rd synthetic aperture acquisition for aperture 65.
        Receive(j+k+2).Apod(1:32) = 1.0;
        Receive(j+k+2).aperture = 65;
        Receive(j+k+2).framenum = i;
        Receive(j+k+2).acqNum = j+2;  % 4 acquisitions per frame
         % -- 2D acquisition, 4rd synthetic aperture acquisition for aperture 65.
        Receive(j+k+3).Apod(33:64) = 1.0;
        Receive(j+k+3).aperture = 65;
        Receive(j+k+3).framenum = i;
        Receive(j+k+3).acqNum = j+3;  % 4 acquisitions per frame

        Receive(j+k).InputFilter = BPF1; % 18 MHz for 2D
        Receive(j+k+1).InputFilter = BPF1;
        Receive(j+k+2).InputFilter = BPF1;
        Receive(j+k+3).InputFilter = BPF1;
    end
    for j = (4*na+1):(4*na+ne)
        % Doppler acquisition
        Receive(j+k).Apod(:) = [zeros(1,16),ones(1,32),zeros(1,16)];
        Receive(j+k).aperture = 33;
        Receive(j+k).startDepth = P(2).startDepth;
        Receive(j+k).endDepth = P(2).startDepth + wl2sPer128*ceil(maxAcqLngthDop/wl2sPer128);
        Receive(j+k).sampleMode = 'BS100BW';
        Receive(j+k).demodFrequency = TW(2).Parameters(1);
        Receive(j+k).TGC = 2;
        Receive(j+k).framenum = i;
        Receive(j+k).acqNum = j;        % Doppler acqNums continue after 2D
        Receive(j+k).InputFilter = BPFDop; % as defined above for Doppler
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
Recon(1).RINums(1,1:4*na) = (1:4*na);  % na ReconInfos needed for na angles
k = 4*na + 1;
% - Set Recon values for Doppler ensemble.
Recon(2).pdatanum = 2;
Recon(2).IntBufDest = [2,1];
Recon(2).ImgBufDest = [0,0];
Recon(2).RINums(1,1:ne) = (k:(k+ne-1));   % ne ReconInfos needed for Doppler ensemble.

% Define ReconInfo structures.
% - For 2D, we need 3*na ReconInfo structures for na steering angles.
% - For Doppler, we need ne ReconInfo structures.
ReconInfo = repmat(struct('mode', 'accumIQ', ...    % accumulate IQ data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'pagenum',1, ...
                   'regionnum', 1), 1, 4*na + ne);
% - Set specific ReconInfo attributes.
%   - ReconInfos for 2D frame.
if na>1
    ReconInfo(1).mode = 'replaceIQ'; % replace IQ data
    for j = 1:4:4*na
        ReconInfo(j).txnum = (j+1)/2;
        ReconInfo(j).rcvnum = j;
        ReconInfo(j+1).txnum = (j+1)/2;
        ReconInfo(j+1).rcvnum = j+1;
        ReconInfo(j+2).txnum = (j+1)/2+1;
        ReconInfo(j+2).rcvnum = j+2;
        ReconInfo(j+3).txnum = (j+1)/2+1;
        ReconInfo(j+3).rcvnum = j+3;
    end
    ReconInfo(4*na).mode = 'accumIQ_replaceIntensity'; % accum and detect
else
    ReconInfo(1).mode = 'replaceIntensity';
end

%  - ReconInfos for Doppler ensemble.
k = 4*na;
for j = 1:ne
    ReconInfo(k+j).mode = 'replaceIQ';
    ReconInfo(k+j).txnum = k/2+ 1;
    ReconInfo(k+j).rcvnum = k + j;
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
                         'pgain',3.0,...            % pgain is image processing gain
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
    k = 4*na+ne;
    for j = 1:4:4*na
        Event(n).info = '1st aperture.';
        Event(n).tx = (j+1)/2;
        Event(n).rcv = k*(i-1)+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 1;
        if j == 1
            Event(n).seqControl = [1,8];
        end
        n = n+1;

        Event(n).info = '2nd aperture.';
        Event(n).tx = (j+1)/2;
        Event(n).rcv = k*(i-1)+j+1;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 1;
        n = n+1;

        Event(n).info = '3rd aperture.';
        Event(n).tx = (j+1)/2+1;
        Event(n).rcv = k*(i-1)+j+2;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 1;
        n=n+1;

        Event(n).info = '4rd aperture.';
        Event(n).tx = (j+1)/2+1;
        Event(n).rcv = k*(i-1)+j+3;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 1;
        n=n+1;
    end

    Event(n-1).seqControl = [2,3];   % replace last 2D acquisition Event's seqControl

    % Acquire Doppler ensemble.
    for j = (4*na+1):(4*na+ne)
        Event(n).info = 'Acquire Doppler ensemble';
        Event(n).tx = 2*na+1;
        Event(n).rcv = (4*na+ne)*(i-1)+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 4;
        if j == 4*na+1
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

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 4;

% Save all the structures to a .mat file.
save('MatFiles/L22-14vFlashDoppler_32LE');
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
