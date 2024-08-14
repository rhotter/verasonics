% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: SetUpP4_2gHWideBeamDoppler.m - Example of WideBeam doppler
% imaging used in Clinical Application
%
% Description:
%   Sequence programming file for P4-2gH Phased array, using WideBeam
%   transmits for 2D (B-mode) and Doppler imaging.
%   - this version uses TPC profiles to provide a seperate high voltage
%     level for Doppler.
%   - processing is asynchronous with respect to acquisition.
%   - The image depth, angle, and acquisition PRF can be adjusted in script
%   - An outline on the display shows the area undergoing Doppler processing
%
% Last update:
%   11/01/2020 - Adapted from SetUpP4_2vWideBeamDoppler.m script.

clear all

% Specify system parameters.
Resource.Parameters.numTransmit = 128;          % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;        % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Specify Trans structure array.
Trans.name = 'P4-2gH';
Trans.units = 'wavelengths'; % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);

mm2wl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);

% Set 2D parameters
P.txFocusMm = -600;
P.startDepthMm = 0;
P.endDepthMm = 150;
P.maxDepthMm = 200;
P.numRays = 48; % no. of rays in frame
P.txFocus = P.txFocusMm*mm2wl;
P.startDepth = P.startDepthMm*mm2wl;  % startDepth in wavelength
P.endDepth = P.endDepthMm*mm2wl;
P.maxDepth = P.maxDepthMm*mm2wl;
P.theta = -(pi/360)*75;  % 75 degree sector, theta is angle to left edge
P.rayDelta = 2*(-P.theta)/(P.numRays-1);
aperture = 96*Trans.spacing; % aperture based on 96 elements
dapex = (aperture/2)/tan(-P.theta); % dist. to virt. apex

% Set Doppler parameters
P.dopDispEle = 64; % middle 64 elements width will be used for Doppler display
P.dopNumTx = 96;  % no. of elements in TX aperture
P.dopNumRays = 3; % no. of rays in frame
P.dopStartDepth = P.startDepth;
P.dopEndDepth = P.endDepth;
P.dopPRIs = 14;
P.dopPRF = 3.0e+03; % Doppler PRF in Hz.
P.regionSteerAngle = 0;
P.pwrThres = 0.4;
P.persfreq = 50;
P.perspwr = 70;
P.cpl = 70;  % level at which 2D overwrites Doppler
demodFreq = 2.5;
m = P.dopNumRays;

% Specify PData structure arrays.
% - 2D PData structure
beamWidth = 0.13; % approximate beam width in radians
Angle = (P.theta-beamWidth):(2*(-P.theta+beamWidth)/(P.numRays-1)):(-P.theta+beamWidth);
PData(1).PDelta = [1.0, 0, 0.5];
PData(1).Size(1) = 10 + ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));
PData(1).Size(2) = 10 + ceil(2*(P.endDepth + dapex)*sin(-P.theta)/PData(1).PDelta(1));
PData(1).Size(3) = 1;
PData(1).Origin = [-(PData(1).Size(2)/2)*PData(1).PDelta(1),0,P.startDepth];
PData(1).Region = repmat(struct(...
            'Shape',struct('Name','SectorFT', ...
            'Position',[0,0,-dapex], ...
            'z',P.startDepth, ...
            'r',dapex+P.endDepth, ...
            'angle',P.rayDelta*10, ...
            'steer',0,...
            'andWithPrev', 1)),1,P.numRays+1);

PData(1).Region(1).Shape.angle = -P.theta*2;
PData(1).Region(1).Shape.andWithPrev = 0;

for n = 2:P.numRays+1
    PData(1).Region(n).Shape.steer = Angle(n-1);%P.theta + (n-1)*P.rayDelta;
end
PData(1).Region = computeRegions(PData(1));

% - Doppler PData structure for TXPD calculation
zVal = 25; rVal = 75;
P.dopTheta = -(pi/360)*30;  % 30 degree sector
P.dopRadius = (aperture/3)/tan(-P.dopTheta); % dist. to virt. apex
P.dopRayDelta = 2*(-P.dopTheta)/(P.dopNumRays);
dopZ = zVal*mm2wl;
dopR = rVal*mm2wl;
dopBeamSteer = -P.dopTheta/m*(-2:2:2);

PData(2) = PData(1); PData(2).Region = [];
PData(2).PDelta(1) = Trans.spacing/2;
PData(2).Size(2) = 10 + ceil(2*(P.endDepth + P.dopRadius)*sin(-P.dopTheta)/PData(2).PDelta(1));
PData(2).Origin = [-(PData(2).Size(2)/2)*PData(2).PDelta(1),0,P.startDepth];
PData(2).Region = repmat(struct(...
            'Shape',struct('Name','SectorFT', ...
            'Position',[0,0,-P.dopRadius], ...
            'z',dopZ, ...
            'r',P.dopRadius+dopR, ...
            'angle',P.dopRayDelta, ...
            'steer',0)),1,P.dopNumRays);

for n = 1:m
    PData(2).Region(n).Shape.steer = dopBeamSteer(n);
end

PData(2).Region = computeRegions(PData(2));

% - PData(3) is used to outline the Doppler area
PData(3) = PData(2); PData(3).Region = [];
PData(3).Region = struct(...
            'Shape',struct('Name','SectorFT', ...
            'Position',[0,0,-P.dopRadius], ...
            'z',dopZ, ...
            'r',P.dopRadius+dopR, ...
            'angle',-2*P.dopTheta, ...
            'steer',0));

PData(3).Region = computeRegions(PData(3));

% Specify Media object.
pt1;
Media.attenuation = -0.5; % not required for doppler script
Media.function = 'movePoints';

% Specify Resources.
% - RcvBuffer(1) is for both 2D and Doppler acquisitions.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 4096*(P.numRays + m*P.dopPRIs);
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 20;           % 20 frames allocated for RF acqusitions.
% InterBuffer(1) is for 2D reconstructions.
Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).numFrames = 1;          % one intermediate frame needed for 2D.
% InterBuffer(2) is for Doppler reconstructions.
Resource.InterBuffer(2).datatype = 'complex';
Resource.InterBuffer(2).numFrames = 1;          % one intermediate frame needed for Doppler.
Resource.InterBuffer(2).pagesPerFrame = P.dopPRIs;     % P.dopPRIs pages per ensemble
% ImageBuffer(1) is for 2D image.
Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).numFrames = 20;
% ImageBuffer(2) is for Doppler image.
Resource.ImageBuffer(2).datatype = 'double';    % image buffer for Doppler
Resource.ImageBuffer(2).numFrames = 20;
% DisplayWindow is for 2D combined with Doppler
Resource.DisplayWindow(1).Title = 'P4-2gHWideBeamDoppler';
Resource.DisplayWindow(1).pdelta = 0.3;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)]; % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).numFrames = 50;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).Colormap = grayscaleCFImap;
Resource.DisplayWindow(1).splitPalette = 1;

% ------Specify structures used in Events------
% Specify Transmit waveform structures.
% - 2D transmit waveform
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,0.67,1,1];
% - Doppler transmit waveform, transmit frequency should be equivalent to
% supported demodFrequency
TW(2).type = 'parametric';
TW(2).Parameters = [demodFreq,0.67,12,1];

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'focus', P.txFocus, ...
                   'Steer', [0.0,0.0], ...
                   'Apod', ones(1,Trans.numelements), ...
                   'Delay', zeros(1,Trans.numelements),...
                   'TXPD', [], ...
                   'peakCutOff', 0.5, ...
                   'peakBLMax',15), 1, P.numRays+m); % P.numRays TXs for 2D + m for Doppler

% -- P.numRays TX structs needed for 2D
beamWidth = 0.15; % approximate beam width in radians
Angle = (P.theta-beamWidth):(2*(-P.theta+beamWidth)/(P.numRays-1)):(-P.theta+beamWidth);
for n = 1:P.numRays
    TX(n).Origin = [dapex*tan(Angle(n)), 0.0, 0.0];
    TX(n).Steer = [Angle(n),0.0];
end

% -- P.dopNumRays TX structs needed for Doppler
winNum = 192;
W = hann(winNum)'; % 128 elements for hann window
ApodRange = (-P.dopNumTx/2+1:P.dopNumTx/2);

% dopBeamWidth = 1.8; % approximate beam width in radians
% dopBeamSteer = (P.dopTheta-dopBeamWidth):(2*(-P.dopTheta+dopBeamWidth)/(P.dopNumRays-1)):(-P.dopTheta+dopBeamWidth);
dopTxSteer = dopBeamSteer; %-P.dopTheta/m*(-3.5:3.5:3.5);
TxDopOrgX = P.dopRadius*tan(dopTxSteer);
for n = 1:m
    TxInd = P.numRays+n;
    TX(TxInd).Origin(1) = TxDopOrgX(n);
    TX(TxInd).Apod = W(winNum/2+ApodRange);
    TX(TxInd).waveform = 2;
    TX(TxInd).Steer = [dopTxSteer(n),0.0];
end

% TX delay and TXPD calculation
h = waitbar(0,'Program TX parameters, please wait!');
for i = 1:P.numRays
    TX(i).Delay = computeTXDelays(TX(i));
    TX(i).TXPD = computeTXPD(TX(i),PData(1));
    waitbar(i/(P.numRays+m))
end
for i = P.numRays+1:P.numRays+m
    TX(i).Delay = computeTXDelays(TX(i));
    TX(i).TXPD = computeTXPD(TX(i),PData(2));
    waitbar(i/(P.numRays+m))
end
close(h)

% Specify TPC structures.
TPC(1).name = '2D';
TPC(1).maxHighVoltage = Trans.maxHighVoltage;
TPC(2).name = 'Doppler';
TPC(2).maxHighVoltage = Trans.maxHighVoltage;

RcvProfile(1).LnaGain = 18; % Profile used for imaging
RcvProfile(2).LnaGain = 24; % Profile used for Doppler
RcvProfile(2).LnaZinSel = 31; % Force high-Z state for best Doppler sensitivity

% Specify TGC Waveform structures.
% - 2D TGC
TGC(1).CntrlPts = [0,243,428,557,663,756,804,884];
TGC(1).rangeMax = P.endDepth;
TGC(1).Waveform = computeTGCWaveform(TGC(1));
% - Doppler TGC
TGC(2).CntrlPts = [0,272,662,662,662,662,662,662];
TGC(2).rangeMax = P.dopEndDepth;
TGC(2).Waveform = computeTGCWaveform(TGC(2));

% Specify Receive structure arrays.
%   We need to acquire all the 2D and Doppler data within a single RcvBuffer frame.  This allows
%   the transfer-to-host DMA after each frame to transfer a large amount of data, improving throughput.
% - We need P.numRays Receives for a 2D frame and P.dopPRIs Receives for a Doppler frame.
maxAcqLngth2D = ceil(sqrt(aperture^2 + P.endDepth^2 - 2*aperture*P.endDepth*cos(P.theta-pi/2)) - P.startDepth);
maxAcqLngthDop = maxAcqLngth2D;
samplesPerWaveDop = 2*TW(2).Parameters(1)/Trans.frequency;

% 30% BW BPF
BPF = ...
[+0.00150 +0.00000 -0.00348 +0.00000 +0.00116 +0.00000 +0.00873 ...
 +0.00000 -0.02283 -0.00003 +0.02713 +0.00000 -0.00284 +0.00000 ...
 -0.05923 -0.00003 +0.14578 +0.00003 -0.22281 +0.00000 +0.25385];

wl4sPer128 = 128/(4*2);  % wavelengths in a 128 sample block for 4 smpls per wave round trip.
wl2sPer128 = 128/(2*samplesPerWaveDop);  % wavelengths in a 128 sample block for 2 smpls per wave round trip.
Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', wl4sPer128*ceil(maxAcqLngth2D/wl4sPer128)+P.startDepth, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ... % 200% Bandwidth for 2D
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, (P.numRays+m*P.dopPRIs)*Resource.RcvBuffer(1).numFrames);

% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    k = (P.numRays + m*P.dopPRIs)*(i-1); % k keeps track of Receive index increment per frame.
    % - Set attributes for each frame.
    Receive(k+1).callMediaFunc = 1;
    for j = 1:P.numRays  % acquisitions for 2D
        Receive(j+k).framenum = i;
        Receive(j+k).acqNum = j;
    end
    for j = (P.numRays+1):(P.numRays+P.dopPRIs*m)
        % Doppler acquisition
        Receive(j+k).startDepth = P.dopStartDepth;
        Receive(j+k).endDepth = wl2sPer128*ceil(maxAcqLngthDop/wl2sPer128)+P.dopStartDepth;
        Receive(j+k).sampleMode = 'BS100BW';
        Receive(j+k).demodFrequency = TW(2).Parameters(1);
        Receive(j+k).InputFilter = BPF;
        Receive(j+k).TGC = 2;
        Receive(j+k).framenum = i;
        Receive(j+k).acqNum = j;        % Doppler acqNums continue after 2D
    end
end

% Specify Recon structure arrays.
% - We need two Recon structures, one for 2D, one for Doppler. These will be referenced in the same
%   event, so that they will use the same (most recent) acquisition frame.
Recon = repmat(struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame', -1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...
               'RINums', zeros(1,1)), 1, 2);
% - Set Recon values for 2D frame.
Recon(1).RINums(1,1:P.numRays) = (1:P.numRays);  % P.numRays ReconInfos needed for P.numRays angles
k = P.numRays + 1;
% - Set Recon values for Doppler ensemble.
Recon(2).pdatanum = 2;
Recon(2).IntBufDest = [2,1];
Recon(2).ImgBufDest = [0,0];
Recon(2).RINums(1,1:m*P.dopPRIs) = (k:(k+m*P.dopPRIs-1));   % P.dopPRIs ReconInfos needed for Doppler ensemble.

% Define ReconInfo structures.
% - For 2D, we need P.numRays ReconInfo structures for P.numRays steering angles.
% - For Doppler, we need m*P.dopPRIs ReconInfo structures.
ReconInfo = repmat(struct('mode', 'accumIQ', ...    % accumulate IQ data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'pagenum',1, ...
                   'scaleFactor', 0.2, ...
                   'regionnum', 1,...
                   'threadSync',1), 1, P.numRays + m*P.dopPRIs);
% - Set specific ReconInfo attributes.
%   - ReconInfos for 2D frame.
ReconInfo(1).Pre = 'clearInterBuf';
for j = 1:P.numRays
    ReconInfo(j).txnum = j;
    ReconInfo(j).rcvnum = j;
    ReconInfo(j).regionnum = j+1;
end
ReconInfo(P.numRays).Post = 'IQ2IntensityImageBuf';

%  - ReconInfos for Doppler ensemble.
ReconInfo(P.numRays+1).Pre = 'clearInterBuf';
for n = 1:m
    for j = 1:P.dopPRIs
        ReconInfo(P.numRays+j+(n-1)*P.dopPRIs).mode = 'replaceIQ';
        ReconInfo(P.numRays+j+(n-1)*P.dopPRIs).txnum = P.numRays + n;
        ReconInfo(P.numRays+j+(n-1)*P.dopPRIs).rcvnum = P.numRays+j+(n-1)*P.dopPRIs;
        ReconInfo(P.numRays+j+(n-1)*P.dopPRIs).regionnum = n;
        ReconInfo(P.numRays+j+(n-1)*P.dopPRIs).pagenum = j;
        ReconInfo(P.numRays+j+(n-1)*P.dopPRIs).scaleFactor = 1;
    end
end

% Specify Process structure arrays.
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
                         'SrcPages',[3,P.dopPRIs-2],...        % start frame number in source buffer
                         'ImgBufDest',[2,-1],...
                         'pdatanum',2,...           % number of PData structure
                         'prf',P.dopPRF,...           % Doppler PRF in Hz
                         'wallFilter','regression',...  % 1 -> quadratic regression
                         'pwrThreshold',P.pwrThres,...
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
                         'persistLevel',P.persfreq,...
                         'interpMethod','4pt',...
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',40,...
                         'mappingMethod','upperHalf',...
                         'threshold',P.cpl,...
                         'display',1,...      % display image after processing
                         'displayWindow',1};

% EF1 is external function for ROI plot
Process(4).classname = 'External';
Process(4).method = 'ROIplot';
Process(4).Parameters = {'srcbuffer','none'};

% Specify SeqControl structure arrays.
% - set receive profile for 2D
SeqControl(1).command = 'setRcvProfile';
SeqControl(1).argument = 1;
% - time between wide beams - frame time = 48*120usec
SeqControl(2).command = 'timeToNextAcq';
SeqControl(2).argument = 230;
% -- Time between 2D acquisition and Doppler ensemble. Set to allow time for profile change.
SeqControl(3).command = 'timeToNextAcq';
SeqControl(3).argument = 2000; % time in usec
% -- Change to TPC Profile 2 (Doppler)
SeqControl(4).command = 'setTPCProfile';
SeqControl(4).condition = 'next';
SeqControl(4).argument = 2;
% - set Receive profile for Doppler
SeqControl(5).command = 'setRcvProfile';
SeqControl(5).argument = 2;
% - time between Doppler ensumbles (dopPRF)
SeqControl(6).command = 'timeToNextAcq';
SeqControl(6).argument = round(1/(P.dopPRF*1e-06));
%---------------check Doppler PRF--------------------
tof = ceil(2*Receive(P.numRays+1).endDepth/demodFreq);
if SeqControl(6).argument < tof
    SeqControl(6).argument = tof;
    P.dopPRF = round(1/(tof*1e-06));
    SeqControl(6).argument = round(1/(P.dopPRF*1e-06));
    fprintf(['"timeToNextAcq" is adjusted to ' num2str(tof) '\n']);
    fprintf(['"P.dopPRF" is adjusted to ' num2str(P.dopPRF) '\n']);
    for k = 1:2:length(Process(2).Parameters)
        if strcmp(Process(2).Parameters{k},'prf'), Process(2).Parameters{k+1} = P.dopPRF; end
    end
end
%---------------check Doppler PRF--------------------
% -- Change to Profile 1 (2D)
SeqControl(7).command = 'setTPCProfile';
SeqControl(7).condition = 'next';
SeqControl(7).argument = 1;
% - time between frames (15000 usec = 15msec)
SeqControl(8).command = 'timeToNextAcq';
SeqControl(8).argument = 2000;
% - return to Matlab
SeqControl(9).command = 'returnToMatlab';
% - Jump back to start.
SeqControl(10).command = 'jump';
SeqControl(10).argument = 2;
nsc = 11; % nsc is count of SeqControl objects

% Specify Event structure arrays.
n = 1;

Event(n).info = 'ROI plot';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 4;
Event(n).seqControl = 0;
n = n+1;

for i = 1:Resource.RcvBuffer(1).numFrames
    rcvInd = (P.numRays + m*P.dopPRIs)*(i-1); % rcvInd keeps track of Receive index increment per frame.

    Event(n).info = 'Set RcvProfile for 2D';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = 1;
    n = n+1;

    % Acquire 2D frame
    for j = 1:P.numRays
        Event(n).info = 'Acquire 2D Wide Beam';
        Event(n).tx = j;
        Event(n).rcv = rcvInd+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;
    end
    Event(n-1).seqControl = [3,4];   % replace last 2D acquisition Event's seqControl

    Event(n).info = 'Set RcvProfile for Doppler';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = 5;
    n = n+1;

    % Acquire Doppler ensemble.
    for j = 1:m
        for k = 1:P.dopPRIs
            Event(n).info = 'Acquire Doppler ensemble';
            Event(n).tx = P.numRays+j;
            Event(n).rcv = rcvInd+P.numRays+k+(j-1)*P.dopPRIs;
            Event(n).recon = 0;
            Event(n).process = 0;
            Event(n).seqControl = 6;
            n = n+1;
        end
    end
    Event(n-1).seqControl = [7,8,nsc]; % replace last Doppler acquisition Event's seqControl
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
        Event(n).seqControl = 9;
    end
    n = n+1;

end

Event(n).info = 'Jump back';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 10;


% User specified UI Control Elements
import vsv.seq.uicontrol.VsSliderControl; 
import vsv.seq.uicontrol.VsButtonGroupControl; 

% - Sensitivity Cutoff
UI(1).Control = VsSliderControl('LocationCode','UserB7','Label','Sens. Cutoff',...
                  'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                  'SliderStep',[0.025,0.1],'ValueFormat','%1.3f', ...
                  'Callback', @SensCutoffCallback );

% - Doppler Mode Button Group
UI(2).Control = VsButtonGroupControl('LocationCode','UserB4','Title','Doppler Mode',...
                  'PossibleCases', {'Velocity','Power'},'Callback',@DopplerModeCallback);

% - Doppler Power Threshold Slider
UI(3).Control = VsSliderControl('LocationCode','UserB3','Label','DopPwrThres',...
                  'SliderMinMaxVal',[0.0,1.0,P.pwrThres],...
                  'SliderStep',[0.02,0.1],'ValueFormat','%3.2f', ...
                  'Callback', @DopPowerThreshold);

% - Color Priority Threshold Slider
UI(4).Control = VsSliderControl('LocationCode','UserB2','Label','Color Priority',...
                  'SliderMinMaxVal',[0,255,P.cpl],...
                  'SliderStep',[1/255,0.1],'ValueFormat','%3.0f', ...
                  'Callback', @ColorPriorityLevel);
            
% - Color Persistence Slider
UI(5).Control = VsSliderControl('LocationCode','UserB1','Label','Color Persistence',...
                  'SliderMinMaxVal',[0,100,P.persfreq],...
                  'SliderStep',[1/100,0.1],'ValueFormat','%3.0f', ...
                  'Callback', @ColorPersistenceLevel);

% - Replot ROI for doppler
UI(6).Control = VsButtonGroupControl('LocationCode','UserC4','Title','Doppler display',...
                  'PossibleCases', {'on','off'},'Callback',@ReplotROI);

% - PRF adjustment
PRFmin = 250; PRFmax = 4250; stepDiff = PRFmax-PRFmin;
UI(7).Control = VsSliderControl('LocationCode','UserC3','Label','PRF (Hz)',...
                  'SliderMinMaxVal',[PRFmin,PRFmax,P.dopPRF],...
                  'SliderStep',[50/stepDiff,500/stepDiff],'ValueFormat','%3.0f', ...
                  'Callback', @PRFChange);

% - Range Change
MinMaxMm = [20,P.maxDepthMm]; % min max in mm
if isfield(Resource.DisplayWindow(1),'AxesUnits') && strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
    AxesUnit = 'mm';
    MinMaxVal = [MinMaxMm, P.dopEndDepth/mm2wl];
else
    AxesUnit = 'wls';
    MinMaxVal = [MinMaxMm*mm2wl, P.dopEndDepth];
end
stepBase = MinMaxVal(2)-MinMaxVal(1);
UI(8).Control = VsSliderControl('LocationCode','UserA1','Label',['Range (',AxesUnit,')'],...
                  'SliderMinMaxVal',MinMaxVal,...
                  'SliderStep',[5/stepBase,10/stepBase],'ValueFormat','%3.1f', ...
                  'Callback', @RangeChangeCallback); 

% - Doppler Region z
MinMax = [P.startDepthMm,P.maxDepthMm];
if strcmp(AxesUnit,'wls')
    MinMax = MinMax*mm2wl;
    zVal = zVal*mm2wl;
    rVal = rVal*mm2wl;
end
stepBase = MinMax(2)-MinMax(1);
UI(9).Control = VsSliderControl('LocationCode','UserC2','Label','Dop Region z',...
                  'SliderMinMaxVal',[MinMax,zVal],...
                  'SliderStep',[5/stepBase,10/stepBase],'ValueFormat','%3.0f', ...
                  'Callback', @DopRegionZ);
          
% - Doppler Region r
UI(10).Control = VsSliderControl('LocationCode','UserC1','Label','Dop Region r',...
                  'SliderMinMaxVal',[MinMax,rVal],...
                  'SliderStep',[5/stepBase,10/stepBase],'ValueFormat','%3.0f', ...
                  'Callback', @DopRegionR);


% - External function for ROIplot
EF(1).Function = vsv.seq.function.ExFunctionDef('ROIplot', @ROIplot);

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 3;

% Save all the structures to a .mat file.
save('MatFiles/P4-2gHWideBeamDoppler');

% Automatic VSX Execution:
% Uncomment the following line to automatically run VSX every time you run
% this SetUp script (note that if VSX finds the variable 'filename' in the
% Matlab workspace, it will load and run that file instead of prompting the
% user for the file to be used):

% filename = 'P4-2gHWideBeamDoppler';  VSX;

% **** Callback routines used by UIControls (UI)  ****
function SensCutoffCallback(~,~,UIValue)
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
    
end

function DopplerModeCallback(~,~,UIState)
%DopplerModeCallback - Doppler mode change

    P = evalin('base','P');
    Control = repmat(struct('Command','set&Run','Parameters',[]),1,4);
    Process = evalin('base','Process');
    Resource = evalin('base','Resource');

    switch UIState
       case 1  % Velocity mode
          newMap = grayscaleCFImap;
          newMap(1:128,:) = Resource.DisplayWindow(1).Colormap(1:128,:);
          P.perspwr = get(findobj('Tag','UserB1Slider'),'Value');
          Control(1).Parameters = {'Process',2,'method','computeCFIFreqEst'};
          Control(2).Parameters = {'Process',3,'srcData','signedColor','persistMethod','dynamic','persistLevel',P.persfreq};
          Control(3).Parameters = {'DisplayWindow',1,'colormap',newMap};
          Control(4).Parameters = {'ImageBuffer',1,'lastFrame',0};
          set(findobj('tag','UserB1Edit'),'String',num2str(P.persfreq,'%3.0f'));
          set(findobj('tag','UserB1Slider'),'Value',P.persfreq);
          assignin('base','DopState','freq');
          % Set modified Process attributes in base Matlab environment.
          Process(2).method = 'computeCFIFreqEst';
          for k = 1:2:length(Process(3).Parameters)
              if strcmp(Process(3).Parameters{k},'srcData'), Process(3).Parameters{k+1} = 'signedColor';
              elseif strcmp(Process(3).Parameters{k},'persistMethod'), Process(3).Parameters{k+1} = 'dynamic';
              elseif strcmp(Process(3).Parameters{k},'persistLevel'), Process(3).Parameters{k+1} = P.persfreq;
              end
          end
       case 2  % Power mode
          newMap = grayscaleCPAmap;
          newMap(1:128,:) = Resource.DisplayWindow(1).Colormap(1:128,:);
          for k = 1:2:length(Process(3).Parameters)
              if strcmp(Process(3).Parameters{k},'persistLevel'), P.persfreq = Process(3).Parameters{k+1}; end
          end
          Control(1).Parameters = {'Process',2,'method','computeCFIPowerEst'};
          Control(2).Parameters = {'Process',3,'srcData','unsignedColor','persistMethod','simple','persistLevel',P.perspwr};
          Control(3).Parameters = {'DisplayWindow',1,'colormap',newMap};
          Control(4).Parameters = {'ImageBuffer',1,'lastFrame',0};
          set(findobj('tag','UserB1Edit'),'String',num2str(P.perspwr,'%3.0f'));
          set(findobj('tag','UserB1Slider'),'Value',P.perspwr);
          assignin('base','DopState','power');
          Process(2).method = 'computeCFIPowerEst';
          for k = 1:2:length(Process(3).Parameters)
              if strcmp(Process(3).Parameters{k},'srcData'), Process(3).Parameters{k+1} = 'unsignedColor';
              elseif strcmp(Process(3).Parameters{k},'persistMethod'), Process(3).Parameters{k+1} = 'simple';
              elseif strcmp(Process(3).Parameters{k},'persistLevel'), Process(3).Parameters{k+1} = P.perspwr;
              end
          end
    end
    assignin('base','P',P);
    assignin('base','Process',Process);
    assignin('base','Control', Control);

    % If PTool window is open, adjust all uicontrols
    hPTool = findobj('tag','ProcessTool');
    if ishandle(hPTool)
        posPTool = get(hPTool,'position');
        PTool;
        set(findobj('tag','ProcessTool'),'position',posPTool);
    end
    
end

function DopPowerThreshold(~,~,UIValue)
%DopPowerThreshold - Doppler Power change

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

function ColorPriorityLevel(~,~,UIValue)
%ColorPriorityLevel - Color Priority change

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
    
end

function ColorPersistenceLevel(~,~,UIValue)
%ColorPersistenceLevel - Color Persistence change

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

function ReplotROI(~,~,UIState)
%ReplotROI

    Process = evalin('base','Process');
    VsType = evalin('base','Resource.DisplayWindow(1).Type');

    if UIState == 1
        evalin('base','drawRegionOutline(hROI,''on'');');
        for k = 1:2:length(Process(1).Parameters)
            if strcmp(Process(1).Parameters{k},'display'), Process(1).Parameters{k+1} = 0; end
        end
        for k = 1:2:length(Process(3).Parameters)
            if strcmp(Process(3).Parameters{k},'display'), Process(3).Parameters{k+1} = 1; end
        end
        % Set Control.Command to set Image.persistLevel.
        Control = evalin('base','Control');
        Control(1).Command = 'set&Run';
        Control(2).Command = 'set&Run';
        Control(1).Parameters = {'Process',1,'display',0};
        Control(2).Parameters = {'Process',3,'display',1};
    else
        evalin('base','drawRegionOutline(hROI,''off'');');
        for k = 1:2:length(Process(1).Parameters)
            if strcmp(Process(1).Parameters{k},'display'), Process(1).Parameters{k+1} = 1; end
        end
        for k = 1:2:length(Process(3).Parameters)
            if strcmp(Process(3).Parameters{k},'display'), Process(3).Parameters{k+1} = 0; end
        end
        % Set Control.Command to set Image.persistLevel.
        Control = evalin('base','Control');
        Control(1).Command = 'set&Run';
        Control(2).Command = 'set&Run';
        Control(1).Parameters = {'Process',1,'display',1};
        Control(2).Parameters = {'Process',3,'display',0};
    end
    assignin('base','Control', Control);
    
end

function PRFChange(~,~,UIValue)
%PRFChange - PRF

    P = evalin('base','P');
    P.dopPRF = UIValue;
    Trans = evalin('base','Trans');
    Process = evalin('base','Process');
    demodFreq = evalin('base','demodFreq');
    SeqControl = evalin('base','SeqControl');
    m = P.dopNumRays;
    SeqControl(6).argument = round(1/(P.dopPRF*1e-06));

    %---------------check Doppler PRF--------------------
    currentDepth = evalin('base','Receive(P.numRays+1).endDepth');
    tof = ceil(2*currentDepth/demodFreq);
    if SeqControl(6).argument < tof
        SeqControl(6).argument = tof;
        P.dopPRF = round(1/(tof*1e-06));
        SeqControl(6).argument = round(1/(P.dopPRF*1e-06));
        fprintf(['"timeToNextAcq" is adjusted to ' num2str(tof) '\n']);
        fprintf(['"dopPRF" is adjusted to ' num2str(P.dopPRF) '\n']);
        UI = evalin('base','UI');
        set(UI(7).handle(2),'Value',P.dopPRF);
        set(UI(7).handle(3),'String',num2str(P.dopPRF));
    end
    %---------------check Doppler PRF--------------------

    for k = 1:2:length(Process(2).Parameters)
        if strcmp(Process(2).Parameters{k},'prf'), Process(2).Parameters{k+1} = P.dopPRF; end
    end

    assignin('base','P',P);
    assignin('base','Process',Process);
    assignin('base','SeqControl',SeqControl);

    Control(1).Command = 'set&Run';
    Control(1).Parameters = {'Process',2,'prf',P.dopPRF};
    Control(2).Command = 'update&Run';
    Control(2).Parameters = {'SeqControl'};
    assignin('base','Control',Control);

    
end

function RangeChangeCallback(hObject,~,UIValue)
%RangeChangeCallback - range change

    simMode = evalin('base','Resource.Parameters.simulateMode');
    % No range change if in simulate mode 2.
    if simMode == 2
        set(hObject,'Value',evalin('base','P.endDepth'));
        return
    end

    % P.endDepth
    Trans = evalin('base','Trans');
    Resource = evalin('base','Resource');
    mm2wl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);
    P = evalin('base','P');
    P.endDepth = UIValue;
    if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
        if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
            P.endDepth = UIValue*mm2wl;
        end
    end
    P.dopEndDepth = P.endDepth;

    % PData
    PData = evalin('base','PData');
    m = P.dopNumRays;
    
    aperture = 96*Trans.spacing; 
    dapex = (aperture/2)/tan(-P.theta);

    PData(1).Size(1) = 10 + ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));
    PData(1).Size(2) = 10 + ceil(2*(P.endDepth + dapex)*sin(-P.theta)/PData(1).PDelta(1));
    PData(1).Origin = [-(PData(1).Size(2)/2)*PData(1).PDelta(1),0,P.startDepth];
    % - Doppler PData structure for TXPD calculation
    PData(2).Size(1) = 10 + ceil((P.endDepth-P.startDepth)/PData(2).PDelta(3));
    PData(2).Size(2) = 10 + ceil(2*(P.endDepth + P.dopRadius)*sin(-P.dopTheta)/PData(2).PDelta(1));
    PData(2).Origin = [-(PData(2).Size(2)/2)*PData(2).PDelta(1),0,P.startDepth];

    for n = 1:P.numRays
        PData(1).Region(n).Shape.r = dapex + P.endDepth;
    end

    PData(3).Size = PData(2).Size;
    PData(3).Origin = PData(2).Origin;

    % check Doppler Region
    dopRegionR = P.dopRadius + P.dopEndDepth;
    if PData(3).Region.Shape.r > dopRegionR
        UI = evalin('base','UI');
        for n = 1:m
            PData(2).Region(n).Shape.r = dopRegionR;
        end
        PData(3).Region.Shape.r = dopRegionR;
        set(UI(10).handle(2),'Value',UIValue);
        set(UI(10).handle(3),'String',num2str(UIValue,'%3.0f'));

        if (PData(3).Region.Shape.z - PData(3).Region.Shape.Position(3)) > PData(3).Region.Shape.r - 5
            PData(3).Region.Shape.z = PData(3).Region.Shape.r - 5 + PData(3).Region.Shape.Position(3);
            for n = 1:m
                PData(2).Region(n).Shape.z = PData(3).Region.Shape.z;
            end
            z = PData(3).Region.Shape.z;
            if strcmp(evalin('base','AxesUnit'),'mm')
                z = z/mm2wl;
            end
            set(UI(9).handle(2),'Value',z);
            set(UI(9).handle(3),'String',num2str(z,'%3.0f'));
        end
    end

    PData(1).Region = computeRegions(PData(1));
    PData(2).Region = computeRegions(PData(2));
    PData(3).Region = computeRegions(PData(3));

    assignin('base','PData',PData);
    evalin('base','Resource.DisplayWindow(1).Position(3) = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);');
    evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);');
    evalin('base','Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];');

    % Update TXPD data of TX structures.
    TX = evalin('base','TX');
    h = waitbar(0,'Program TX parameters, please wait!');

    for i = 1:P.numRays
        TX(i).TXPD = computeTXPD(TX(i),PData(1));
        waitbar(i/(P.numRays+m))
    end
    for i = P.numRays+1:P.numRays+m
        TX(i).TXPD = computeTXPD(TX(i),PData(2));
        waitbar(i/(P.numRays+m))
    end
    close(h)
    assignin('base','TX',TX);

    % Receive
    Receive = evalin('base', 'Receive');
    dopTWfreq = evalin('base','TW(2).Parameters(1)');

    maxAcqLngth2D = ceil(sqrt(aperture^2 + P.endDepth^2 - 2*aperture*P.endDepth*cos(P.theta-pi/2)) - P.startDepth);
    maxAcqLngthDop = maxAcqLngth2D;
    samplesPerWaveDop = 2*dopTWfreq/Trans.frequency;

    wl4sPer128 = 128/(4*2);  % wavelengths in a 128 sample block for 4 smpls per wave round trip.
    wl2sPer128 = 128/(2*samplesPerWaveDop);  % wavelengths in a 128 sample block for 2 smpls per wave round trip.
    for i = 1:Resource.RcvBuffer(1).numFrames
        k = (P.numRays + m*P.dopPRIs)*(i-1); % k keeps track of Receive index increment per frame.
        for j = 1:P.numRays  % acquisitions for 2D
            Receive(j+k).endDepth = P.startDepth + wl4sPer128*ceil(maxAcqLngth2D/wl4sPer128);
        end
        for j = (P.numRays+1):(P.numRays+P.dopPRIs*m) % Doppler acquisition
            Receive(j+k).endDepth = P.dopStartDepth + wl2sPer128*ceil(maxAcqLngthDop/wl2sPer128);
        end
    end
    assignin('base','Receive',Receive);

    % PRF
    %---------------check Doppler PRF--------------------
    Process = evalin('base','Process');
    SeqControl = evalin('base','SeqControl');
    SeqControl(6).argument = round(1/(P.dopPRF*1e-06));
    currentDepth = evalin('base','Receive(P.numRays+1).endDepth');
    tof = ceil(2*currentDepth/dopTWfreq);
    if SeqControl(6).argument < tof
        SeqControl(6).argument = tof;
        P.dopPRF = round(1/(tof*1e-06));
        SeqControl(6).argument = round(1/(P.dopPRF*1e-06));
        fprintf(['"timeToNextAcq" is adjusted to ' num2str(tof) '\n']);
        fprintf(['"dopPRF" is adjusted to ' num2str(P.dopPRF) '\n']);
        UI = evalin('base','UI');
        set(UI(7).handle(2),'Value',P.dopPRF);
        set(UI(7).handle(3),'String',num2str(P.dopPRF));
        for k = 1:2:length(Process(2).Parameters)
            if strcmp(Process(2).Parameters{k},'prf'), Process(2).Parameters{k+1} = P.dopPRF; end
        end
        assignin('base','Process',Process);
    end
    %---------------check Doppler PRF--------------------

    assignin('base','P',P);
    assignin('base','SeqControl',SeqControl);

    % TGC
    TGC = evalin('base','TGC');
    TGC(1).rangeMax = P.endDepth;
    TGC(1).Waveform = computeTGCWaveform(TGC(1));
    TGC(2).rangeMax = P.dopEndDepth;
    TGC(2).Waveform = computeTGCWaveform(TGC(2));
    assignin('base','TGC',TGC);

    Control(1).Command = 'set&Run';
    Control(1).Parameters = {'Process',2,'prf',P.dopPRF};
    Control(2).Command = 'update&Run';
    Control(2).Parameters = {'PData','InterBuffer','ImageBuffer','DisplayWindow','SeqControl','TX','Receive','TGC','Recon'};
    assignin('base','Control', Control);
    assignin('base', 'action', 'displayChange');
end

function DopRegionZ(~,~,UIValue)
%DopRegionZ

    P = evalin('base','P');
    UI = evalin('base','UI');
    PData = evalin('base','PData');
    mm2wl = evalin('base','mm2wl');

    if strcmp(evalin('base','AxesUnit'),'mm')
        z = UIValue*mm2wl;
    else
        z = UIValue;
    end

    if z > PData(3).Region.Shape.r-P.dopRadius-10
        z = PData(3).Region.Shape.r-P.dopRadius-10;
    end

    for n = 1: P.dopNumRays
        PData(2).Region(n).Shape.z = z;
    end
    PData(3).Region.Shape.z = z;

    PData(2).Region = computeRegions(PData(2));
    PData(3).Region = computeRegions(PData(3));

    if strcmp(evalin('base','AxesUnit'),'mm')
        z = z/mm2wl;
    end
    set(UI(9).handle(2),'Value',z);
    set(UI(9).handle(3),'String',num2str(z,'%3.0f'));

    assignin('base','PData',PData);
    Control.Command = 'update&Run';
    Control.Parameters = {'PData','InterBuffer','ImageBuffer','Recon'};
    assignin('base','Control', Control);
end

function DopRegionR(~,~,UIValue)
%DopRegionR

    P = evalin('base','P');
    UI = evalin('base','UI');
    PData = evalin('base','PData');
    mm2wl = evalin('base','mm2wl');

    if strcmp(evalin('base','AxesUnit'),'mm')
        r = UIValue*mm2wl;
        currentDepth = get(UI(8).handle(2),'Value')*mm2wl;
    else
        r = UIValue;
        currentDepth = get(UI(8).handle(2),'Value');
    end

    if r < PData(3).Region.Shape.z+5
        r = PData(3).Region.Shape.z+5;
    elseif r > currentDepth
        r = currentDepth;
    end

    for n = 1: P.dopNumRays
        PData(2).Region(n).Shape.r = r+P.dopRadius;
    end
    PData(3).Region.Shape.r = r+P.dopRadius;

    PData(2).Region = computeRegions(PData(2));
    PData(3).Region = computeRegions(PData(3));

    if strcmp(evalin('base','AxesUnit'),'mm')
        r = r/mm2wl;
    end

    set(UI(10).handle(2),'Value',r);
    set(UI(10).handle(3),'String',num2str(r,'%3.0f'));

    assignin('base','PData',PData);
    Control.Command = 'update&Run';
    Control.Parameters = {'PData','InterBuffer','ImageBuffer','Recon'};
    assignin('base','Control',Control);
end

% **** Callback routines used by External function definition (EF) ****

function ROIplot(varargin)

    persistent drawROI

    % change the txt for voltage slider
    hv1txt = findobj('tag','hv1txt');
    hv2txt = findobj('tag','hv2txt');
    set(hv1txt,'String','Bmode Voltage');
    set(hv2txt,'String','Doppler Voltage');

    % drawRegionOutline(WinNum,PDataNum,RegionNum) creates an outline from
    % PData(PDataNum).Region(RegionNum) on displayWindow(WinNum) with default color - white.
    if isempty(drawROI)
        drawROI = 1;
        evalin('base','hROI = drawRegionOutline(1,3,1);')
    else
        evalin('base','drawRegionOutline(hROI,1,3,1);')
    end
end
%-EF#1
