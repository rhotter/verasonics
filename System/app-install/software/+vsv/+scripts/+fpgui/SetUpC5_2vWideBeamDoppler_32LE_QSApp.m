% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: SetUpC5_2vWideBeamDoppler_32LE.m - Example of WideBeam doppler
% imaging used in Clinical Application
%
% Description:
%   Sequence programming file for C5-2v curved array, using WideBeam
%   transmits for 2D (B-mode) and doppler imaging.
%   - this version uses TPC profiles to provide a seperate high voltage
%     level for Doppler.
%   - numTx channels are used for transmit, with numTx/2 transmitters on each side of the center element
%     (where possible).
%   - processing is asynchronous with respect to acquisition.
%   - The image depth, angle, and acquisition PRF can be adjused in Advanced script
%   - An outline on the display shows the area undergoing Doppler processing
%
% Last update:
% 08/06/2019 - modified for software release 4.1.
% 05/06/2020 - Update to SW 4.3 format for new UIControls and function definitions (VTS 1691).
%   More info:(.../Example_Scripts/Vantage_Features/New UI Scheme/SetUpL11_5vFlash_NewUI)

%clear[ 	]+all

demodFreq = 2.9762;

P.txFocusMm = 800;
P.startDepthMm = 2;
P.endDepthMm = 120;
P.maxDepthMm = 200;

% Specify system parameters.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.
Resource.System.SoftwareVersion = [4 1 0]; % Minimum software release for this script.
Resource.System.UTA = '260-MUX'; % This script requires the 260-MUX UTA.
Resource.Parameters.numTransmit = 64;

% Specify Trans structure array.
Trans.name = 'C5-2v';
Trans.units = 'wavelengths'; % required in Gen3 to prevent default to mm units
Trans = computeTrans(Trans);  % C5-2v transducer is 'known' transducer so we can use computeTrans.
Trans.maxHighVoltage = 50;  % set maximum high voltage limit for pulser supply.
radius = Trans.radius;
Trans = computeUTAMux64(Trans); % Add HVMux field for use with UTA 260-Mux

% Convert mm to wavelength
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);
% Set 2D parameters
P.numTx = 48;   % no. of elements in TX aperture. Should no larger than 63
P.numRays = 64; % no. of rays in frame
P.txFocus = P.txFocusMm*scaleToWvl;
P.startDepth = P.startDepthMm*scaleToWvl;  % startDepth in wavelength
P.endDepth = P.endDepthMm*scaleToWvl;
P.maxDepth = P.maxDepthMm*scaleToWvl;

% Set Doppler parameters
P.dopNumTx = 64;  % no. of elements in TX aperture (might be adjusted later)
P.dopNumRays = 4; % no. of rays in frame
P.dopStartDepth = P.startDepth;
P.dopEndDepth = P.endDepth;
P.dopPRIs = 14;
P.dopPRF = 3.0e+03; % Doppler PRF in Hz.
P.regionSteerAngle = 0;
P.pwrThres = 0.3;
P.persfreq = 70;
P.perspwr = 80;
P.cpl = 70;  % level at which 2D overwrites Doppler

m = P.dopNumRays;

% Specify PData structure array.
% - 2D PData structure
scanangle = Trans.numelements*Trans.spacing/radius;
dtheta = scanangle/P.numRays; % angle between rays
theta = -(scanangle/2) + 0.5*dtheta; % angle to left edge from centerline
Angle = theta:dtheta:(-theta);

PData(1).PDelta = [1.0, 0, 0.5];  % x, y and z pdeltas
sizeRows = 10 + ceil((P.endDepth + radius - (radius * cos(scanangle/2)))/PData(1).PDelta(3));
sizeCols = 10 + ceil(2*(P.endDepth + radius)*sin(scanangle/2)/PData(1).PDelta(1));
PData(1).Size = [sizeRows,sizeCols,1];     % size, origin and pdelta set region of interest.
PData(1).Origin(1,1) = (P.endDepth+radius)*sin(-scanangle/2) - 5;
PData(1).Origin(1,2) = 0;
PData(1).Origin(1,3) = ceil(radius * cos(scanangle/2)) - radius - 5;

% Define PData Regions for numRays scanlines
PData(1).Region = repmat(struct(...
        'Shape',struct('Name','Sector',...
        'Position',[0,0,-radius],...
        'r1',radius+P.startDepth,...
        'r2',radius+P.endDepth,...
        'angle',dtheta*7,...
        'steer',0,...
        'andWithPrev', 0)),1,P.numRays+1);

PData(1).Region(1).Shape.angle = scanangle;
% AND with Region 1 to have sharp edge in Bmode
for n = 1:5
    PData(1).Region(n+1).Shape.andWithPrev = 1;
    PData(1).Region(end+1-n).Shape.andWithPrev = 1;
end

% Change the steer angle
for n = 2:P.numRays+1
    PData(1).Region(n).Shape.steer = Angle(n-1);
end

PData(1).Region = computeRegions(PData(1));

% - Doppler PData structure
r1Val = 40;  r2Val = 100;
dopScanangle = pi/180*30;%scanangle/2;
dopDtheta = dopScanangle/P.dopNumRays; % angle between rays
dopTheta = -(dopScanangle/2) + 0.5*dopDtheta; % angle to left edge from centerline
dopBeamSteer = dopTheta:dopDtheta:(-dopTheta);
dopR1 = r1Val*scaleToWvl;
dopR2 = r2Val*scaleToWvl;

PData(2) = PData(1); PData(2).Region = [];
PData(2).PDelta(3) = 1.0;
PData(2).Size(1) = 10 + ceil((P.dopEndDepth + radius - (radius * cos(scanangle/2)))/PData(2).PDelta(3));
PData(2).Size(2) = 10 + ceil(2*(P.dopEndDepth + radius)*sin(scanangle/2)/PData(2).PDelta(1));
PData(2).Region = repmat(struct(...
    'Shape',struct('Name','Sector',...
    'Position',[0,0,-radius],...
    'r1',radius+dopR1,...
    'r2',radius+dopR2,...
    'angle',dopDtheta,...
    'steer',dopBeamSteer(1))),1,m);
for n = 1:m
    PData(2).Region(n).Shape.steer = dopBeamSteer(n);
end

PData(2).Region = computeRegions(PData(2));

% - PData(3) is used to outline the doppler area
PData(3) = PData(2); PData(3).Region = [];
PData(3).Region.Shape = struct(...
    'Name','Sector',...
    'Position',[0,0,-radius],...
    'r1',radius+dopR1,...
    'r2',radius+dopR2,...
    'angle',dopScanangle);
PData(3).Region = computeRegions(PData(3));

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

maxBufLength = ceil(sqrt((P.maxDepth+radius)^2 + radius^2 - ...
    2*(P.maxDepth+radius)*radius*cos(scanangle)));
maxBufSizePerAcq = 128*ceil(maxBufLength*8/128); % wavelengths in a 128 sample block for 4 smpls per wave round trip.
% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = maxBufSizePerAcq*(P.numRays + m*P.dopPRIs);
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numTransmit;
Resource.RcvBuffer(1).numFrames = numRcvFrames;     % Two frames needed for double buffer.
% InterBuffer(1) is for 2D reconstructions.
Resource.InterBuffer(1).numFrames = 1;          % one intermediate frame needed for 2D.
% InterBuffer(2) is for Doppler reconstructions.
Resource.InterBuffer(2).numFrames = 1;          % one intermediate frame needed for Doppler.
Resource.InterBuffer(2).pagesPerFrame = P.dopPRIs;     % ne pages per ensemble
% ImageBuffer(1) is for 2D image.
Resource.ImageBuffer(1).numFrames = 20;
% ImageBuffer(2) is for Doppler image.
Resource.ImageBuffer(2).numFrames = 20;
% DisplayWindow is for 2D combined with Doppler
Resource.DisplayWindow(1).Title = 'C5-2vWideBeamDoppler_32LE';
Resource.DisplayWindow(1).pdelta = 0.5;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
    DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)]; % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).numFrames = numImageFrames;
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
TW(2).Parameters = [demodFreq,0.67,12,1];

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
    'Origin', [0.0,0.0,0.0], ...
    'focus', P.txFocus, ...
    'Steer', [0.0,0.0], ...
    'Apod', zeros(1,Trans.numelements), ...
    'Delay', zeros(1,Trans.numelements), ...
    'TXPD', [], ...
    'peakCutOff', 0.1, ...
    'peakBLMax', 15.0), 1, P.numRays+m);

% - Set event specific TX attributes.
%    numTx/2 is the number of elements to include on each side of the
%    center element, for the specified focus and sensitivity cutoff.
L = P.numTx+1;
W1 = blackman(2*L);
W = W1(round(L/2):round(L/2)+L);
for n = 1:P.numRays   % numRays transmit events
    % Set transmit Origins.
    TX(n).Origin = [radius*sin(Angle(n)), 0.0, radius*cos(Angle(n))-radius];
    ce = round(1+(Trans.numelements-1)*(Angle(n) - theta)/(-2*theta));
    % compute active elements for 64 element aperture around center element, ce
    ActiveElements = zeros(1,Trans.numelements);
    lft = ce - 31;
    if lft < 1, lft = 1; end
    rt = ce + 32;
    if rt > Trans.numelements, rt = Trans.numelements; end
    ActiveElements(lft:rt) = 1;
    % Compute MUX aperture for all active elements.
    TX(n).aperture = computeMuxAperture(ActiveElements, Trans);
    % Set transmit Apodization so that a maximum of numTx + 1 transmitters are active.
    w1 = 1;
    w2 = L;
    lft = round(ce - P.numTx/2);
    if lft < 1, w1 = -lft+2; lft = 1; end
    rt = round(ce + P.numTx/2);
    if rt > Trans.numelements, w2 = L - (rt - Trans.numelements); rt = Trans.numelements; end
    TX(n).Apod(lft:rt) = W(w1:w2);
end

% -- P.dopNumRays TX structs needed for Doppler
% TX.Steer is not required for curved array
winNum = 192;
W = hann(winNum);

[~,ind] = min(abs(radius*sin(dopBeamSteer(1)-abs(P.regionSteerAngle))-Trans.ElementPos(:,1)));
while ind <= P.dopNumTx/2
    P.dopNumTx = P.dopNumTx-2;
end
fprintf('The aperture size of each doppler beam: %g elements \n', P.dopNumTx);

ApodRange = (-P.dopNumTx/2+1:P.dopNumTx/2);
for n = 1:P.dopNumRays
    TxInd = P.numRays+n;
    TX(TxInd).Origin = [radius*sin(dopBeamSteer(n)+P.regionSteerAngle),...
        0.0, radius*cos(dopBeamSteer(n)+P.regionSteerAngle)-radius];
    [~,ind] = min(abs(TX(TxInd).Origin(1) - Trans.ElementPos(:,1)));
    P.dopTxOrgChnl(n) = ind;
    TX(TxInd).Apod(P.dopTxOrgChnl(n)+ApodRange) = W(winNum/2+ApodRange);
    TX(TxInd).aperture = computeMuxAperture(TX(TxInd).Apod, Trans);
    TX(TxInd).waveform = 2;
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
TPC(2).maxHighVoltage = 30; % VTS-1007 reduce Doppler TX to safe level for mux chips

RcvProfile(1).LnaGain = 18; % Profile used for imaging
RcvProfile(2).LnaGain = 24; % Profile used for Doppler
RcvProfile(2).LnaZinSel = 31; % Force high-Z state for best Doppler sensitivity

% Specify TGC Waveform structures.
% - 2D TGC
TGC(1).CntrlPts = [0,267,522,668,830,921,1023,1023];
TGC(1).rangeMax = P.endDepth;
TGC(1).Waveform = computeTGCWaveform(TGC(1));
% - Doppler TGC
TGC(2).CntrlPts = [0,272,662,662,662,662,662,662];
TGC(2).rangeMax = P.dopEndDepth;
TGC(2).Waveform = computeTGCWaveform(TGC(2));

% Specify Receive structure arrays.
%   We need to acquire all the 2D and Doppler data within a single RcvBuffer frame.  This allows
%   the transfer-to-host DMA after each frame to transfer a large amount of data, improving throughput.
% - We need na Receives for a 2D frame and ne Receives for a Doppler frame.
maxAcqLength2D = sqrt((P.endDepth+radius)^2 + radius^2 - ...
    2*(P.endDepth+radius)*radius*cos(scanangle)) - P.startDepth;
maxAcqLengthDop = sqrt((P.dopEndDepth+radius)^2 + radius^2 - ...
    2*(P.dopEndDepth+radius)*radius*cos(dopScanangle)) - P.dopStartDepth;
samplesPerWaveDop = 2*demodFreq/Trans.frequency;

% 30% BW BPF
BPF = ...
    [+0.00021 +0.00000 -0.00131 +0.00000 +0.00443 +0.00000 -0.01099 ...
     +0.00000 +0.02228 +0.00000 -0.03876 +0.00000 +0.05957 +0.00000 ...
     -0.08209 +0.00000 +0.10257 +0.00000 -0.11697 +0.00000 +0.12213];

wl4sPer128 = 128/(4*2);  % wavelengths in a 128 sample block for 4 smpls per wave round trip.
wl2sPer128 = 128/(2*samplesPerWaveDop);  % wavelengths in a 128 sample block for 2 smpls per wave round trip.
Receive = repmat(struct('Apod', zeros(1,Trans.numelements), ...
    'startDepth', P.startDepth, ...
    'endDepth', P.startDepth + wl4sPer128*ceil(maxAcqLength2D/wl4sPer128), ...
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
        Receive(j+k).aperture = TX(j).aperture;
        ce = round(1+(Trans.numelements-1)*(Angle(j) - theta)/(-2*theta));
        lft = ce - 15;
        if lft < 1, lft = 1; end
        rt = ce + 16;
        if rt > Trans.numelements, rt = Trans.numelements; end
        Receive(j+k).Apod(lft:rt) = 1;
    end
    for j = 1:m
        for n = 1:P.dopPRIs
        % Doppler acquisition
            ind = (j-1)*P.dopPRIs + n + P.numRays;
            Receive(ind+k).aperture = TX(j+P.numRays).aperture;
            lft = P.dopTxOrgChnl(j) - 15;
            if lft < 1, lft = 1; end
            rt = P.dopTxOrgChnl(j) + 16;
            if rt > Trans.numelements, rt = Trans.numelements; end
            Receive(ind+k).Apod(lft:rt) = 1;
            Receive(ind+k).startDepth = P.dopStartDepth;
            Receive(ind+k).endDepth = wl2sPer128*ceil(maxAcqLengthDop/wl2sPer128)+P.dopStartDepth;
            Receive(ind+k).sampleMode = 'BS100BW';
            Receive(ind+k).demodFrequency = demodFreq;
            Receive(ind+k).InputFilter = BPF;
            Receive(ind+k).TGC = 2;
            Receive(ind+k).framenum = i;
            Receive(ind+k).acqNum = ind;        % Doppler acqNums continue after 2D
        end
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
        ReconInfo(P.numRays+n+(j-1)*m).mode = 'replaceIQ';
        ReconInfo(P.numRays+n+(j-1)*m).txnum = P.numRays + n;
        ReconInfo(P.numRays+n+(j-1)*m).rcvnum = P.numRays+j+(n-1)*P.dopPRIs;
        ReconInfo(P.numRays+n+(j-1)*m).regionnum = n;
        ReconInfo(P.numRays+n+(j-1)*m).pagenum = j;
        ReconInfo(P.numRays+n+(j-1)*m).scaleFactor = 1;
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
SeqControl(2).argument = 350;
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

%% User specified UI Control Elements
import vsv.seq.uicontrol.VsSliderControl
import vsv.seq.uicontrol.VsButtonGroupControl

% - Sensitivity Cutoff
UI(1).Control = VsSliderControl('LocationCode','UserB7','Label','Sens. Cutoff',...
                   'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                   'SliderStep',[0.025,0.1],'ValueFormat','%1.3f',...
                   'Callback', @SensCutoffCallback);


% - Doppler Mode Button Group
UI(2).Control = VsButtonGroupControl('LocationCode','UserB4','Title','Doppler Mode',...
                 'PossibleCases',   {'Velocity','Power'},...
                 'Callback', @DopplerModeCallback);


% - Doppler Power Threshold Slider
UI(3).Control = VsSliderControl('LocationCode','UserB3','Label','DopPwrThres',...
                   'SliderMinMaxVal',[0.0,1.0,P.pwrThres],...
                   'SliderStep',[0.02,0.1],'ValueFormat','%3.2f',...
                   'Callback', @DopPowerThresholdCallback);

% - Color Priority Threshold Slider
UI(4).Control = VsSliderControl('LocationCode','UserB2','Label','Color Priority',...
                   'SliderMinMaxVal',[0,255,P.cpl],...
                   'SliderStep',[1/255,0.1],'ValueFormat','%3.0f',...
                   'Callback', @ColorPriorityLevelCallback);

% - Color Persistence Slider
UI(5).Control = VsSliderControl('LocationCode','UserB1','Label','Color Persistence',...
                   'SliderMinMaxVal',[0,100,P.persfreq],...
                   'SliderStep',[1/100,0.1],'ValueFormat','%3.0f',...
                   'Callback', @ColorPersistenceLevelCallback);


% - Replot ROI for doppler
UI(6).Control = VsButtonGroupControl('LocationCode','UserC4','Title','Doppler Display',...
                 'PossibleCases',   {'on','off'},...
                 'Callback', @ReplotROICallback);

% - PRF adjustment
PRFmin = 250; PRFmax = 4250; stepDiff = PRFmax-PRFmin;
UI(7).Control = VsSliderControl('LocationCode','UserC3','Label','PRF (Hz)',...
                   'SliderMinMaxVal',[PRFmin,PRFmax,P.dopPRF],...
                   'SliderStep',[50/stepDiff,500/stepDiff],'ValueFormat','%3.0f',...
                   'Callback', @PRFChangeCallback);

% - Range Change
MinMaxMm = [20,P.maxDepthMm]; % min max in mm
if isfield(Resource.DisplayWindow(1),'AxesUnits') && strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
    AxesUnit = 'mm';
    MinMaxVal = [MinMaxMm, P.dopEndDepth/scaleToWvl];
else
    AxesUnit = 'wls';
    MinMaxVal = [MinMaxMm*scaleToWvl, P.dopEndDepth];
end
stepBase = MinMaxVal(2)-MinMaxVal(1);
UI(8).Control = VsSliderControl('LocationCode','UserA1','Label',['Range (',AxesUnit,')'],...
                   'SliderMinMaxVal',MinMaxVal,...
                   'SliderStep',[5/stepBase,10/stepBase],'ValueFormat','%3.0f',...
                   'Callback', @RangeChangeCallback);


% - Doppler Region Radius r1
MinMax = [P.startDepthMm,P.maxDepthMm];
if strcmp(AxesUnit,'wls')
    MinMax = MinMax*scaleToWvl;
    r1Val = r1Val*scaleToWvl;
    r2Val = r2Val*scaleToWvl;
end
stepBase = MinMax(2)-MinMax(1);
UI(9).Control = VsSliderControl('LocationCode','UserC2','Label','Dop Region r1',...
                   'SliderMinMaxVal',[MinMax,r1Val],...
                   'SliderStep',[5/stepBase,10/stepBase],'ValueFormat','%3.0f',...
                   'Callback', @DopRegionR1Callback);


% - Doppler Region Radius r2
UI(10).Control = VsSliderControl('LocationCode','UserC1','Label','Dop Region r2',...
                   'SliderMinMaxVal',[MinMax,r2Val],...
                   'SliderStep',[5/stepBase,10/stepBase],'ValueFormat','%3.0f',...
                   'Callback', @DopRegionR2Callback);


% - External function for ROIplot
EF(1).Function = vsv.seq.function.ExFunctionDef('ROIplot',@ROIplot);

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 3;

% Save all the structures to a .mat file.
filename = fullfile( vsv.file.getVSXDir, 'MatFiles/SetUpC5_2vWideBeamDoppler_32LE_QSApp.mat');

save(filename);
return

%% **** Callback routines used by UIControls (UI) ****
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

function DopPowerThresholdCallback(~,~,UIValue)
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

function ColorPriorityLevelCallback(~,~,UIValue)
%ColorPriorityLevel - Color Priority change
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
end

function ColorPersistenceLevelCallback(~,~,UIValue)
%ColorPersistenceLevel - Color Persistence change
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
end

function ReplotROICallback(~,~,UIState)
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

function PRFChangeCallback(~,~,UIValue)
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
    scaleToWvl = evalin('base','scaleToWvl');
    P = evalin('base','P');
    P.endDepth = UIValue;
    if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
        if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
            P.endDepth = UIValue*scaleToWvl;
        end
    end
    P.dopEndDepth = P.endDepth;
    P.txFocus = 4*P.dopEndDepth;

    % PData
    PData = evalin('base','PData');
    radius = evalin('base','radius');
    scanangle = evalin('base','scanangle');
    m = P.dopNumRays;

    sizeRows = 10 + ceil((P.endDepth + radius - (radius * cos(scanangle/2)))/PData(1).PDelta(3));
    sizeCols = 10 + ceil(2*(P.endDepth + radius)*sin(scanangle/2)/PData(1).PDelta(1));
    PData(1).Size = [sizeRows,sizeCols,1];     % size, origin and pdelta set region of interest.
    PData(1).Origin(1,1) = (P.endDepth+radius)*sin(-scanangle/2) - 5;

    for n = 1:P.numRays+1
        PData(1).Region(n).Shape.r2 = radius + P.endDepth;
    end

    PData(2).Size(1) = 10 + ceil((P.dopEndDepth + radius - (radius * cos(scanangle/2)))/PData(2).PDelta(3));
    PData(2).Size(2) = 10 + ceil(2*(P.dopEndDepth + radius)*sin(scanangle/2)/PData(2).PDelta(1));
    PData(2).Origin = PData(1).Origin;

    PData(3).Size = PData(2).Size;
    PData(3).Origin = PData(2).Origin;

    % check Doppler Region
    if PData(3).Region.Shape.r2 > radius + P.endDepth
        UI = evalin('base','UI');
        for n = 1:m
            PData(2).Region(n).Shape.r2 = radius + P.endDepth;
        end
        PData(3).Region.Shape.r2 = radius + P.endDepth;
        set(UI(10).handle(2),'Value',UIValue);
        set(UI(10).handle(3),'String',num2str(UIValue,'%3.0f'));

        if PData(3).Region.Shape.r1 > PData(3).Region.Shape.r2 - 5
            PData(3).Region.Shape.r1 = PData(3).Region.Shape.r2 - 5;
            for n = 1:m
                PData(2).Region(n).Shape.r1 = PData(3).Region.Shape.r1;
            end
            r1 = PData(3).Region.Shape.r1 - radius;
            if strcmp(evalin('base','AxesUnit'),'mm')
                r1 = r1/scaleToWvl;
            end
            set(UI(9).handle(2),'Value',r1);
            set(UI(9).handle(3),'String',num2str(r1,'%3.0f'));
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
        TX(i).focus = P.txFocus;
        TX(i).Delay = computeTXDelays(TX(i));
        TX(i).TXPD = computeTXPD(TX(i),PData(1));
        waitbar(i/(P.numRays+m))
    end
    for i = P.numRays+1:P.numRays+m
        TX(i).focus = P.txFocus;
        TX(i).Delay = computeTXDelays(TX(i));
        TX(i).TXPD = computeTXPD(TX(i),PData(2));
        waitbar(i/(P.numRays+m))
    end
    close(h)
    assignin('base','TX',TX);

    % Receive
    Receive = evalin('base', 'Receive');
    demodFreq = evalin('base','demodFreq');
    dopScanangle = evalin('base','dopScanangle');

    maxAcqLength2D = sqrt((P.endDepth+radius)^2 + radius^2 - ...
        2*(P.endDepth+radius)*radius*cos(scanangle)) - P.startDepth;
    maxAcqLengthDop = sqrt((P.dopEndDepth+radius)^2 + radius^2 - ...
        2*(P.dopEndDepth+radius)*radius*cos(dopScanangle)) - P.dopStartDepth;
    samplesPerWaveDop = 2*demodFreq/Trans.frequency;
    wl4sPer128 = 128/(4*2);  % wavelengths in a 128 sample block for 4 smpls per wave round trip.
    wl2sPer128 = 128/(2*samplesPerWaveDop);  % wavelengths in a 128 sample block for 2 smpls per wave round trip.
    for i = 1:Resource.RcvBuffer(1).numFrames
        k = (P.numRays + m*P.dopPRIs)*(i-1); % k keeps track of Receive index increment per frame.
        for j = 1:P.numRays  % acquisitions for 2D
            Receive(j+k).endDepth = P.startDepth + wl4sPer128*ceil(maxAcqLength2D/wl4sPer128);
        end
        for j = (P.numRays+1):(P.numRays+P.dopPRIs*m) % Doppler acquisition
            Receive(j+k).endDepth = P.dopStartDepth + wl2sPer128*ceil(maxAcqLengthDop/wl2sPer128);
        end
    end

    assignin('base','Receive',Receive);

    % PRF
    %---------------check Doppler PRF--------------------
    Process = evalin('base','Process');
    SeqControl = evalin('base','SeqControl');
    SeqControl(6).argument = round(1/(P.dopPRF*1e-06));
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

function DopRegionR1Callback(~,~,UIValue)
    m = evalin('base','P.dopNumRays');
    PData = evalin('base','PData');
    radius = evalin('base','radius');
    scaleToWvl = evalin('base','scaleToWvl');

    if strcmp(evalin('base','AxesUnit'),'mm')
        r1 = UIValue*scaleToWvl;
    else
        r1 = UIValue;
    end

    if r1 > PData(3).Region.Shape.r2-radius-5
        r1 = PData(3).Region.Shape.r2-radius-5;
    end

    for n = 1:m
        PData(2).Region(n).Shape.r1 = r1+radius;
    end
    PData(3).Region.Shape.r1 = r1+radius;

    PData(2).Region = computeRegions(PData(2));
    PData(3).Region = computeRegions(PData(3));

    if strcmp(evalin('base','AxesUnit'),'mm')
        r1 = r1/scaleToWvl;
    end
    UI = evalin('base','UI');
    set(UI(9).handle(2),'Value',r1);
    set(UI(9).handle(3),'String',num2str(r1,'%3.0f'));

    assignin('base','PData',PData);
    Control.Command = 'update&Run';
    Control.Parameters = {'PData','InterBuffer','ImageBuffer','Recon'};
    assignin('base','Control', Control);

end

function DopRegionR2Callback(~,~,UIValue)
    m = evalin('base','P.dopNumRays');
    UI = evalin('base','UI');
    PData = evalin('base','PData');
    radius = evalin('base','radius');
    scaleToWvl = evalin('base','scaleToWvl');

    if strcmp(evalin('base','AxesUnit'),'mm')
        r2 = UIValue*scaleToWvl;
        currentDepth = get(UI(8).handle(2),'Value')*scaleToWvl;
    else
        r2 = UIValue;
        currentDepth = get(UI(8).handle(2),'Value');
    end

    if r2 < PData(3).Region.Shape.r1-radius+5
        r2 = PData(3).Region.Shape.r1-radius+5;
    elseif r2 > currentDepth
        r2 = currentDepth;
    end

    for n = 1:m
        PData(2).Region(n).Shape.r2 = r2+radius;
    end
    PData(3).Region.Shape.r2 = r2+radius;

    PData(2).Region = computeRegions(PData(2));
    PData(3).Region = computeRegions(PData(3));

    if strcmp(evalin('base','AxesUnit'),'mm')
        r2 = r2/scaleToWvl;
    end

    set(UI(10).handle(2),'Value',r2);
    set(UI(10).handle(3),'String',num2str(r2,'%3.0f'));

    assignin('base','PData',PData);
    Control.Command = 'update&Run';
    Control.Parameters = {'PData','InterBuffer','ImageBuffer','Recon'};
    assignin('base','Control',Control);

end

%% **** Callback routines used by External function definition (EF) ****
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

