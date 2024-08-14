% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: SetUpL11_5gHWideBeamDoppler_64LE.m - Example of WideBeam 2D amd Doppler
%   imaging.
%
% Description:
%   Sequence programming file for L11-5gH Linear array, using WideBeam transmits
%   for 2D (B-mode) and Doppler imaging. For 2D imaging, P.numRays overlapping wide
%   beams are used, with a typical beam width, P.numTx, of 32 elements and a focal
%   point, P.txFocus, far below the depth of interest. For Doppler imaging, four
%   overlapping wide beams are used to cover the Doppler region of interest. The
%   central portion of the Doppler wide beams are processed to provide the Doppler
%   flow image.  The acquisition sequence acquires all the 2D wide beams, followed
%   by the Doppler wide beams.
%   - This version uses TPC profiles to provide a seperate high voltage
%     level for 2D and Doppler.
%   - Acquisition of Doppler ensembles is interleaved - e.g. for the 4 Doppler beams,
%     each beam is acquired at 1x, 2x or 4x the Doppler PRF, repeated for the length
%     of the ensemble.
%   - processing is asynchronous with respect to acquisition.
%   - The image depth, Doppler angle, and Doppler PRF can be adjusted.
%   - An outline on the display shows the area undergoing Doppler processing
%
%   WideBeam Imaging with 64 element apertures:
%   In the following aperture examples, each space represents 2 elements.
%   Wide beam 1 (P.numTx=32):
%   tttttttttttttttt------------------------------------------------
%   rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr--------------------------------
%   Wide beam 24 (P.numTx=32, P.numRays=48):
%   ------------------------tttttttttttttttt------------------------
%   ----------------rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr----------------
%   Wide beam 48 (P.numTx=32, P.numRays=48):
%   ------------------------------------------------tttttttttttttttt
%   --------------------------------rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
%
%
% Last update:
% 09/17/21 Create for new L11-5gH transducer

clear all

% Adjustable 2D parameters:
P.startDepth = 2; % startDepth in wavelengths
P.endDepth = 200; % endDepth in wavelengths - range 100(~20mm) - 300(~60mm).
P.numRays = 48;   % no. of 2D rays in frame
P.numTx = 32;     % no. of elements used for 2D wide beams.
P.txFocus = 4*P.endDepth; % widebeam focal point (below depth of image)
% Adjustable Doppler parameters:
P.startDepthDop = 2;
P.endDepthDop = 0.75*P.endDepth; % initial endDepth of Doppler color box.
P.angleDop = 12*pi/180;
P.PRIsDop = 16;   % No. of Doppler acquisitions in ensemble
P.PRFDop = 3e+03; % Doppler PRF in Hz - range 500 - 4000
P.DopColorBox = [-10.0,P.startDepthDop,P.endDepthDop,...
                 (P.endDepthDop-P.startDepthDop),P.angleDop]; % [x,z,width,height,angle]
% Fixed Doppler parameters.
numRaysDop = 4;
demodFreqDop = 5.208; % Doppler demod frequency

% Specify system parameters.
Resource.Parameters.numTransmit = 128;
Resource.Parameters.numRcvChannels = 64;
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.
Resource.System.transmitChannels = 128;
Resource.System.receiveChannels = 64;

% Specify Trans structure array.
Trans.name = 'L11-5gH';
Trans.units = 'wavelengths'; % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);  % computeTrans is used for known transducers.
mm2wl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);

% Specify PData structure arrays.
% - 2D PData structure
PData(1).PDelta = [Trans.spacing, 0, 0.5];
PData(1).Size(1,1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3)); % range, origin and pdelta set PData(1).Size.
PData(1).Size(1,2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(1,3) = 1;      % single image page
PData(1).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of uppr lft crnr.
% Define P.numRays rectangular regions centered on TX beam origins and 2/3 of wide beam width.
TxOrgX = (-63.5*Trans.spacing):(127*Trans.spacing/(P.numRays-1)):(63.5*Trans.spacing);
PData(1).Region = repmat(struct('Shape',struct('Name','Rectangle',...
                                               'Position',[0,0,P.startDepth],...
                                               'width',Trans.spacing*(P.numTx-round(P.numTx/3)),...
                                               'height',P.endDepth-P.startDepth)),1,P.numRays);
% - set position of regions to correspond to transmit origins and beam spacing.
for n = 1:P.numRays, PData(1).Region(n).Shape.Position(1) = TxOrgX(n); end
PData(1).Region = computeRegions(PData(1));

% - Doppler PData structure
PData(2) = PData(1);
PData(2).PDelta = [Trans.spacing, 0, Trans.spacing];
PData(2).Size(1) = ceil((P.endDepth-P.startDepth)/PData(2).PDelta(3)); % range, origin and pdelta set PData(1).Size.
PData(2).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(2).PDelta(1)); % cols
% - Define the x position of the centers of the Doppler beams.
dopRegWidth = P.DopColorBox(3)/numRaysDop;  % (color box width)/(no. of beams) in wavelengths
TxDopOrgX = (P.DopColorBox(1)-dopRegWidth*(numRaysDop/2-0.5)):dopRegWidth:(P.DopColorBox(1)+dopRegWidth*(numRaysDop/2-0.5));
numTxDop = round(1.7*dopRegWidth/Trans.spacing); % define Doppler wide beam width 1.7 times Doppler region width
PData(2).Region = repmat(struct('Shape',struct(...
    'Name','Parallelogram',...
    'Position',[0,0,0],... %will be changed later
    'width', dopRegWidth,...
    'height', P.DopColorBox(4),...
    'angle',P.DopColorBox(5))),1,numRaysDop+1); % define numRaysDop regions and one extra for full window
for n = 1:numRaysDop
    PData(2).Region(n).Shape.Position(1) = TxDopOrgX(n);
end
% - PData(2).Region(numRaysDop+1) is used to outline the Doppler area
PData(2).Region(numRaysDop+1).Shape.Position = [P.DopColorBox(1),0,P.DopColorBox(2)];
PData(2).Region(numRaysDop+1).Shape.width = P.DopColorBox(3);
PData(2).Region = computeRegions(PData(2));

% Specify Media object.
pt1;
Media.attenuation = -0.0; % not required for doppler script
Media.function = 'movePoints';

% Specify Resources.
% RcvBuffer(1) is for both 2D and Doppler acquisitions.
% - Compute rowsPerFrame from max 2D depth and no. of widebeams + max Dop. depth and
%   no. of ensemble PRIs * no. of ensembles(4).
maxAcqLngth2D = ceil(sqrt(300^2 + ((Trans.numelements-1)*Trans.spacing)^2));
% - Calculate maxAcqLngthdop as the diagonal of the color box from the law of cosines.
%      c^2 = a^2 + b^2 - 2*a*b*cos(angle)
a = P.DopColorBox(3); % width of top of color box
b = P.endDepthDop/cos(P.angleDop); % length of right side of color box
angle = P.angleDop + pi/2;  % angle between top and rt. side (P.dopAngle positive).
maxAcqLngthDop =  sqrt(a^2 + b^2 - 2*a*b*cos(angle))-P.startDepthDop;
% - Adjust maxAcqLngthDop, since hardware needs multiple of 128 samples.
nSmpls = 2*maxAcqLngthDop * numRaysDop/2; % no. of acquisition samples for BS100BW
if abs(round(nSmpls/128) - nSmpls/128) < .01 % round to nearest 128 sample boundary for hardware.
    numRcvSamples = 128*round(nSmpls/128);
else
    numRcvSamples = 128*ceil(nSmpls/128); % round up to next 128 sample boundary.
end
if numRcvSamples == 0, numRcvSamples = 128; end
maxAcqLngthDop = P.startDepthDop + (numRcvSamples/2)/2; % /2 for wvlngths for BS100BW acquisitions
max2DSmpls = maxAcqLngth2D*2*4*P.numRays; % maxAcqLngth*2(round trip)*4(4 smpls/wvlngth)*P.numRays
maxDopSmpls = maxAcqLngthDop*2*2*P.PRIsDop*numRaysDop; % maxAcqLngth*2(round trip)*2(2 smpls/wvlngth)*...
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = max2DSmpls + maxDopSmpls;
Resource.RcvBuffer(1).colsPerFrame = Resource.System.receiveChannels;
Resource.RcvBuffer(1).numFrames = 20;           % 20 frames allocated for RF acqusitions.
% InterBuffer(1) is for 2D reconstructions.
Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).numFrames = 1;          % one intermediate frame needed for 2D.
% InterBuffer(2) is for Doppler reconstructions.
Resource.InterBuffer(2).datatype = 'complex';
Resource.InterBuffer(2).numFrames = 1;          % one intermediate frame needed for Doppler.
Resource.InterBuffer(2).pagesPerFrame = P.PRIsDop; % P.PRIsDop pages per ensemble
% ImageBuffer(1) is for 2D image.
Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).numFrames = 10;
% ImageBuffer(2) is for Doppler image.
Resource.ImageBuffer(2).datatype = 'double';    % image buffer for Doppler
Resource.ImageBuffer(2).numFrames = 10;
% DisplayWindow is for 2D combined with Doppler
Resource.DisplayWindow(1).Title = 'L11-5gHWideBeamDoppler';
Resource.DisplayWindow(1).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)]; % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).numFrames = 50;
Resource.DisplayWindow(1).AxesUnits = 'mm';
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
TW(2).Parameters = [demodFreqDop,0.67,6,1];

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'focus', P.txFocus, ...
                   'Steer', [0.0,0.0], ...
                   'Apod', zeros(1,Resource.System.transmitChannels), ...
                   'Delay', zeros(1,Resource.System.transmitChannels),...
                   'TXPD', [], ...
                   'peakCutOff', 0.5, ...
                   'peakBLMax',15), 1, P.numRays+numRaysDop); % P.numRays TXs for 2D + numRaysDop for Doppler

% - Set event specific TX attributes.
% -- P.numRays TX structs needed for 2D
if mod(P.numTx,2), P.numTx=P.numTx-1; end % need even P.numTx for ApodFunc
ApodFunc = kaiser(P.numTx,1);
k = P.numTx/2;
Ce = fix(1:(127/(P.numRays-1)):128); % calculate center element numbers.
for n = 1:P.numRays
    TX(n).Origin(1) = TxOrgX(n);
    % Set P.numTX channel transmit apodization centered around center element, Ce(n), if possible
    for j = 1:128
        if (j>(Ce(n)-k))&&(j<=(Ce(n)+k))
            TX(n).Apod(j) = ApodFunc(j-(Ce(n)-k));
        end
    end
    % Compute transmit delays
    TX(n).Delay = computeTXDelays(TX(n));
    % Compute transmit pixel data
    TX(n).TXPD = computeTXPD(TX(n),PData(1));
end

% -- numRaysDop TX structs needed for Doppler
k = P.numRays;
if mod(numTxDop,2), numTxDop=numTxDop-1; end % need even P.numTxDop for ApodFunc
kk = numTxDop/2;
Ced = zeros(1,numRaysDop);
for n = 1:numRaysDop
    TX(n+k).waveform = 2;
    TX(n+k).Origin = [TxDopOrgX(n),0,0];
    % Compute closest element number to transmit origin
    [Dummy,Ced(n)] = min(abs(Trans.ElementPos(:,1)-TX(n+k).Origin(1))); % ce is closest ele to cntr of aper.
    % Compute TX(n).Apod.
    for j = 1:128
        if (j>(Ced(n)-kk))&&(j<=(Ced(n)+kk))
            TX(n+k).Apod(j) = 1;
        end
    end
    TX(n+k).Steer = [P.angleDop,0.0];
    TX(n+k).Delay = computeTXDelays(TX(n+k));
    TX(n+k).TXPD = computeTXPD(TX(n+k),PData(2));
end

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
TGC(1).CntrlPts = [0,241,428,608,714,756,804,833];
TGC(1).rangeMax = P.endDepth;
TGC(1).Waveform = computeTGCWaveform(TGC(1));
% - Doppler TGC
TGC(2).CntrlPts = [0,425,662,713,713,713,662,610];
TGC(2).rangeMax = P.endDepthDop;
TGC(2).Waveform = computeTGCWaveform(TGC(2));

% Specify Receive structure arrays.
%   We need to acquire all the 2D and Doppler data within a single RcvBuffer frame.  This allows
%   the transfer-to-host DMA after each frame to transfer a large amount of data, improving throughput.
% - We need P.numRays Receives for a 2D frame and P.PRIsDop Receives for a Doppler frame.
%
% Define 30% BW InputFilter for Doppler
BPF = ...
[+0.00021 +0.00000 -0.00131 +0.00000 +0.00443 +0.00000 -0.01099 ...
 +0.00000 +0.02228 +0.00000 -0.03876 +0.00000 +0.05957 +0.00000 ...
 -0.08209 +0.00000 +0.10257 +0.00000 -0.11697 +0.00000 +0.12213];

maxAcqLngth2D = sqrt(P.endDepth^2 + (Trans.numelements*Trans.spacing)^2)-P.startDepth;
Receive = repmat(struct('Apod', zeros(1,Resource.Parameters.numTransmit), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLngth2D, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...  % 200% Bandwidth for 2D
                        'demodFrequency', Trans.frequency, ... % Trans.frequency used for 2D
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, (P.numRays+numRaysDop*P.PRIsDop)*Resource.RcvBuffer(1).numFrames);

% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    k = (P.numRays + numRaysDop*P.PRIsDop)*(i-1); % k keeps track of Receive index increment per frame.
    % - Set attributes for each frame.
    Receive(k+1).callMediaFunc = 1;
    for j = 1:P.numRays  % acquisitions for 2D
        % Set 64 channel aperture centered around center element, Ce(j), if possible
        lft = Ce(j) - 31;
        if lft < 1, lft = 1; end
        if lft > 65, lft = 65; end
        rt = Ce(j) + 31;
        if rt < 64, rt = 64; end
        if rt > Trans.numelements, rt = Trans.numelements; end
        Receive(j+k).Apod(lft:rt) = 1;
        Receive(j+k).framenum = i;
        Receive(j+k).acqNum = j;
    end
    % - Doppler acquisitions. Specify as series of single ensemble acquisitions, but
    %   acquisition may interleave 2 or more ensembles.
    for j = 1:numRaysDop
        ind = P.PRIsDop*(j-1) + P.numRays + k;
        % Set 64 channel aperture centered around center element, Ced(j), if possible
        lft = Ced(j) - 31;
        if lft < 1, lft = 1; end
        if lft > 65, lft = 65; end
        rt = Ced(j) + 31;
        if rt < 64, rt = 64; end
        if rt > Trans.numelements, rt = Trans.numelements; end
        for n = 1:P.PRIsDop  % acquire Doppler PRIs for each Doppler beam
            Receive(ind+n).Apod(lft:rt) = 1;
            Receive(ind+n).startDepth = P.startDepthDop;
            Receive(ind+n).endDepth = maxAcqLngthDop;
            Receive(ind+n).sampleMode = 'BS100BW';
            Receive(ind+n).demodFrequency = demodFreqDop;
            Receive(ind+n).InputFilter = BPF;
            Receive(ind+n).TGC = 2;
            Receive(ind+n).framenum = i;
            Receive(ind+n).acqNum = P.numRays+P.PRIsDop*(j-1)+n;  % Doppler acqNums continue after 2D
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
               'RINums', []), 1, 2);
% - Set Recon values for 2D frame.
Recon(1).RINums = 1:P.numRays;  % P.numRays ReconInfos needed for 2D
% - Set Recon values for Doppler ensemble.
Recon(2).pdatanum = 2;
Recon(2).IntBufDest = [2,1];
Recon(2).ImgBufDest = [0,0];
Recon(2).RINums = (P.numRays+1):(P.numRays + numRaysDop*P.PRIsDop); % numRaysDop*P.PRIsDop ReconInfos needed for Doppler ensembles.

% Define ReconInfo structures.
% - For 2D, we need P.numRays ReconInfo structures for P.numRays steering angles.
% - For Doppler, we need P.PRIsDop ReconInfo structures.
ReconInfo = repmat(struct('mode', 'accumIQ', ...    % accumulate IQ data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'pagenum', 1, ...
                   'Pre', [], ...
                   'Post', [], ...
                   'scaleFactor', 0.2, ...
                   'regionnum', 1,...
                   'threadSync',1), 1, P.numRays + numRaysDop*P.PRIsDop);
% - Set specific ReconInfo attributes.
%   - ReconInfos for 2D frame.
ReconInfo(1).Pre = 'clearInterBuf';
for j = 1:P.numRays
    ReconInfo(j).txnum = j;
    ReconInfo(j).rcvnum = j;
    ReconInfo(j).regionnum = j;
end
ReconInfo(P.numRays).Post = 'IQ2IntensityImageBuf';

%  - ReconInfos for Doppler ensemble.
ReconInfo(P.numRays+1).Pre = 'clearInterBuf';
for n = 1:numRaysDop
    k = P.numRays+(n-1)*P.PRIsDop;
    for j = 1:P.PRIsDop
        ReconInfo(k+j).mode = 'replaceIQ_normalize';
        ReconInfo(k+j).txnum = P.numRays + n;
        ReconInfo(k+j).rcvnum = k+j;
        ReconInfo(k+j).regionnum = n;
        ReconInfo(k+j).pagenum = j;
        ReconInfo(k+j).scaleFactor = 1;
    end
end

% Specify Process structure arrays.
DopState = 'freq';
pwrThres = 0.2;
persFreq = 70;
persPwr = 90;
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
Process(2).Parameters = {'IntBufSrc',[2,1],...          % number of Inter buffer to process.
                         'SrcPages',[3,P.PRIsDop-2],...        % start frame number in source buffer
                         'ImgBufDest',[2,-1],...
                         'pdatanum',2,...           % number of PData structure
                         'prf',P.PRFDop,...           % Doppler PRF in Hz
                         'wallFilter','regression',...
                         'pwrThreshold',pwrThres,...
                         'maxPower',20,...
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
                         'persistLevel',persFreq,...
                         'interpMethod','4pt',...
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',40,...
                         'mappingMethod','upperHalf',...
                         'threshold',20,...
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
% - time between wide beams
SeqControl(2).command = 'timeToNextAcq';
SeqControl(2).argument = 200;
% -- Time between 2D acquisition and Doppler ensemble. Set to allow time for profile change.
SeqControl(3).command = 'timeToNextAcq';
SeqControl(3).argument = 5000; % time in usec
% -- Change to TPC Profile 2 (Doppler)
SeqControl(4).command = 'setTPCProfile';
SeqControl(4).condition = 'next';
SeqControl(4).argument = 2;
% - set Receive profile for Doppler
SeqControl(5).command = 'setRcvProfile';
SeqControl(5).argument = 2;
% - Time between interleaved Doppler acquisitions
%   (PRI1beam1,PRI1beam2,...,PRI1beamN,PRI2beam1,PRI2Beam2,...
% - Use Doppler PRF to determine interleave factor (4, 2 or 1)
tof = ceil(2*Receive(P.numRays+1).endDepth/numRaysDop); % round trip time of flight for Doppler (in usecs)
pri = round(1e+06/P.PRFDop);  % Doppler PRI in usecs
if 4*tof < pri
    numIntLv = 4;
elseif 2*tof < pri
    numIntLv = 2;
else
    numIntLv = 1;
end
SeqControl(6).command = 'timeToNextAcq';
SeqControl(6).argument = round(1/(P.PRFDop*numIntLv) * 1e+06);
% -- Change to Profile 1 (2D)
SeqControl(7).command = 'setTPCProfile';
SeqControl(7).condition = 'next';
SeqControl(7).argument = 1;
% - time between frames
SeqControl(8).command = 'timeToNextAcq';
SeqControl(8).argument = 2000;
% - return to Matlab
SeqControl(9).command = 'returnToMatlab';
% - Jump back to start.
SeqControl(10).command = 'jump';
SeqControl(10).argument = 2;
nsc = 11; % nsc is count of SeqControl objects
startDopEvents = zeros(1,Resource.RcvBuffer(1).numFrames); % keeps track of events for Doppler acquisition.

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
    rcvInd = (P.numRays + numRaysDop*P.PRIsDop)*(i-1); % rcvInd keeps track of Receive index increment per frame.

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
    startDopEvents(i) = n; % need start of Doppler acqs for PRF/numIntLv change.
    % Acquire Doppler ensembles.
    if numIntLv == 4
        for j = 1:P.PRIsDop
            for k = 1:numRaysDop
                Event(n).info = 'Acquire four interleaved Doppler ensemble';
                Event(n).tx = P.numRays+k;
                Event(n).rcv = rcvInd+P.numRays+(k-1)*P.PRIsDop+j;
                Event(n).recon = 0;
                Event(n).process = 0;
                Event(n).seqControl = 6;
                n = n+1;
            end
        end
    elseif numIntLv == 2
        % Aquire ensembles from Doppler beams 1 & 2
        for j = 1:P.PRIsDop
            for k = 1:2
                Event(n).info = 'Acquire two interleaved Doppler ensemble';
                Event(n).tx = P.numRays+k;
                Event(n).rcv = rcvInd+P.numRays+(k-1)*P.PRIsDop+j;
                Event(n).recon = 0;
                Event(n).process = 0;
                Event(n).seqControl = 6;
                n = n+1;
            end
        end
        % Aquire ensembles from Doppler beams 3 & 4
        for j = 1:P.PRIsDop
            for k = 3:4
                Event(n).info = 'Acquire two interleaved Doppler ensemble';
                Event(n).tx = P.numRays+k;
                Event(n).rcv = rcvInd+P.numRays+(k-1)*P.PRIsDop+j;
                Event(n).recon = 0;
                Event(n).process = 0;
                Event(n).seqControl = 6;
                n = n+1;
            end
        end
    else % numIntLv = 1
        % Acquire non-interleaved ensembles from each beam.
        for k = 1:numRaysDop
            for j = 1:P.PRIsDop
                Event(n).info = 'Acquire non-interleaved Doppler ensemble';
                Event(n).tx = P.numRays+k;
                Event(n).rcv = rcvInd+P.numRays+(k-1)*P.PRIsDop+j;
                Event(n).recon = 0;
                Event(n).process = 0;
                Event(n).seqControl = 6;
                n = n+1;
            end
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
                   'SliderMinMaxVal',[0.0,1.0,pwrThres],...
                   'SliderStep',[0.02,0.1],'ValueFormat','%3.2f',...
                   'Callback', @DopPowerThresholdCallback);


% - Color Priority Threshold Slider
UI(4).Control = VsSliderControl('LocationCode','UserB2','Label','Color Priority',...
                   'SliderMinMaxVal',[0,255,20],...
                   'SliderStep',[1/255,0.1],'ValueFormat','%3.0f',...
                   'Callback', @ColorPriorityLevelCallback);


% - Color Persistence Slider
UI(5).Control = VsSliderControl('LocationCode','UserB1','Label','Color Persistence',...
                   'SliderMinMaxVal',[0,100,persFreq],...
                   'SliderStep',[1/100,0.1],'ValueFormat','%3.0f',...
                   'Callback', @ColorPersistenceLevelCallback);


% - Replot ROI for doppler
UI(6).Control = VsButtonGroupControl('LocationCode','UserC4','Title','Doppler Display',...
                 'PossibleCases',   {'on','off'},...
                 'Callback', @ReplotROICallback);


% - PRF adjustment
PRFmin = 500; PRFmax = 4000; stepDiff = PRFmax-PRFmin;
UI(7).Control = VsSliderControl('LocationCode','UserC3','Label','PRF (Hz)',...
                   'SliderMinMaxVal',[PRFmin,PRFmax,P.PRFDop],...
                   'SliderStep',[50/stepDiff,500/stepDiff],'ValueFormat','%3.0f',...
                   'Callback', @PRFChangeCallback);


% - Steer Angle adjustment
UI(8).Control = VsSliderControl('LocationCode','UserC2','Label','Doppler Angle',...
                   'SliderMinMaxVal',[0,20,round(P.angleDop*180/pi)],...
                   'SliderStep',[1/40,5/40],'ValueFormat','%3.0f',...
                   'Callback', @SteerAngleCallback);


% - Doppler depth adjustment
if isfield(Resource.DisplayWindow(1),'AxesUnits') && strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
    MinMaxVal = [10,40,P.endDepthDop/mm2wl];
else
    MinMaxVal = [10*mm2wl,40*mm2wl,P.endDepthDop];
end
UI(9).Control = VsSliderControl('LocationCode','UserC1','Label','Doppler Depth',...
                   'SliderMinMaxVal',MinMaxVal,...
                   'SliderStep',[1/20,1/10],'ValueFormat','%3.1f',...
                   'Callback', @DopplerDepthCallback);


% - Range Change
MinMaxMm = [10,60]; % min max in mm
mm2wl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);
if isfield(Resource.DisplayWindow(1),'AxesUnits') && strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
    AxesUnit = 'mm';
    MinMaxVal = [MinMaxMm, P.endDepth/mm2wl];
else
    AxesUnit = 'wls';
    MinMaxVal = [MinMaxMm*mm2wl, P.endDepth];
end
stepBase = MinMaxVal(2)-MinMaxVal(1);
UI(10).Control = VsSliderControl('LocationCode','UserA1','Label',['Range (',AxesUnit,')'],...
                   'SliderMinMaxVal',MinMaxVal,...
                   'SliderStep',[5/stepBase,10/stepBase],'ValueFormat','%3.0f',...
                   'Callback', @RangeChangeCallback);


% - External function for ROIplot
EF(1).Function = vsv.seq.function.ExFunctionDef('ROIplot',@ROIplot);

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 3;

% Save all the structures to a .mat file.
filename = 'MatFiles/L11-5gHWideBeamDoppler_64LE';
save(filename);
return

%% **** Callback routines used by UIControls (UI) ****
function SensCutoffCallback(~,~,UIValue)
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
    P = evalin('base','P');
    persFreq = evalin('base','persFreq');
    persPwr = evalin('base','persPwr');
    Control = repmat(struct('Command','set&Run','Parameters',[]),1,4);
    Process = evalin('base','Process');
    Resource = evalin('base','Resource');

    switch UIState
        case 1  % Velocity mode
            newMap = grayscaleCFImap;
            newMap(1:128,:) = Resource.DisplayWindow(1).Colormap(1:128,:);
            persPwr = get(findobj('Tag','UserB1Slider'),'Value');
            Control(1).Parameters = {'Process',2,'method','computeCFIFreqEst'};
            Control(2).Parameters = {'Process',3,'srcData','signedColor','persistMethod','dynamic','persistLevel',persFreq};
            Control(3).Parameters = {'DisplayWindow',1,'colormap',newMap};
            Control(4).Parameters = {'ImageBuffer',1,'lastFrame',0};
            set(findobj('tag','UserB1Edit'),'String',num2str(persFreq,'%3.0f'));
            set(findobj('tag','UserB1Slider'),'Value',persFreq);
            assignin('base','DopState','freq');
            % Set modified Process attributes in base Matlab environment.
            Process(2).method = 'computeCFIFreqEst';
            for k = 1:2:length(Process(3).Parameters)
                if strcmp(Process(3).Parameters{k},'srcData'), Process(3).Parameters{k+1} = 'signedColor';
                elseif strcmp(Process(3).Parameters{k},'persistMethod'), Process(3).Parameters{k+1} = 'dynamic';
                elseif strcmp(Process(3).Parameters{k},'persistLevel'), Process(3).Parameters{k+1} = persFreq;
                end
            end
        case 2  % Power mode
            newMap = grayscaleCPAmap;
            newMap(1:128,:) = Resource.DisplayWindow(1).Colormap(1:128,:);
            for k = 1:2:length(Process(3).Parameters)
                if strcmp(Process(3).Parameters{k},'persistLevel'), persFreq = Process(3).Parameters{k+1}; end
            end
            Control(1).Parameters = {'Process',2,'method','computeCFIPowerEst'};
            Control(2).Parameters = {'Process',3,'srcData','unsignedColor','persistMethod','simple','persistLevel',persPwr};
            Control(3).Parameters = {'DisplayWindow',1,'colormap',newMap};
            Control(4).Parameters = {'ImageBuffer',1,'lastFrame',0};
            set(findobj('tag','UserB1Edit'),'String',num2str(persPwr,'%3.0f'));
            set(findobj('tag','UserB1Slider'),'Value',persPwr);
            assignin('base','DopState','power');
            Process(2).method = 'computeCFIPowerEst';
            for k = 1:2:length(Process(3).Parameters)
                if strcmp(Process(3).Parameters{k},'srcData'), Process(3).Parameters{k+1} = 'unsignedColor';
                elseif strcmp(Process(3).Parameters{k},'persistMethod'), Process(3).Parameters{k+1} = 'simple';
                elseif strcmp(Process(3).Parameters{k},'persistLevel'), Process(3).Parameters{k+1} = persPwr;
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
    if UIState == 1
        evalin('base','drawRegionOutline(hROI,''on'');');
        for k = 1:2:length(Process(1).Parameters)
            if strcmp(Process(1).Parameters{k},'display'), Process(1).Parameters{k+1} = 0; end
        end
        for k = 1:2:length(Process(3).Parameters)
            if strcmp(Process(3).Parameters{k},'display'), Process(3).Parameters{k+1} = 1; end
        end
        % Set Control.Command.
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
        % Set Control.Command.
        Control = evalin('base','Control');
        Control(1).Command = 'set&Run';
        Control(2).Command = 'set&Run';
        Control(1).Parameters = {'Process',1,'display',1};
        Control(2).Parameters = {'Process',3,'display',0};
    end
    assignin('base','Control', Control);
end

function PRFChangeCallback(~,~,UIValue)
    P = evalin('base','P');
    P.PRFDop = UIValue;
    Resource = evalin('base','Resource');
    Receive = evalin('base','Receive');
    Process = evalin('base','Process');
    SeqControl = evalin('base','SeqControl');
    Event = evalin('base','Event');
    numRaysDop = evalin('base','numRaysDop');
    demodFreqDop = evalin('base','demodFreqDop');
    % - Use Doppler PRF to determine interleave factor (4, 2 or 1)
    tof = ceil(2*Receive(P.numRays+1).endDepth/demodFreqDop); % round trip time of flight for Doppler (in usecs)
    pri = round(1e+06/P.PRFDop);  % Doppler PRI in usecs
    if 4*tof < pri
        numIntLv = 4;
    elseif 2*tof < pri
        numIntLv = 2;
    else
        numIntLv = 1;
    end
    SeqControl(6).command = 'timeToNextAcq';
    SeqControl(6).argument = round(1/(P.PRFDop*numIntLv) * 1e+06);
    % If numIntLv changed, modify Event structures for Doppler acquisition.
    if numIntLv ~= evalin('base','numIntLv')
        startDopEvents = evalin('base','startDopEvents');
        nsc = 11;
        for i = 1:Resource.RcvBuffer(1).numFrames
            rcvInd = (P.numRays + numRaysDop*P.PRIsDop)*(i-1);
            n = startDopEvents(i);
            if numIntLv == 4
                for j = 1:P.PRIsDop
                    for k = 1:numRaysDop
                        Event(n).info = 'Acquire four interleaved Doppler ensemble';
                        Event(n).tx = P.numRays+k;
                        Event(n).rcv = rcvInd+P.numRays+(k-1)*P.PRIsDop+j;
                        Event(n).recon = 0;
                        Event(n).process = 0;
                        Event(n).seqControl = 6;
                        n = n+1;
                    end
                end
            elseif numIntLv == 2
                % Aquire ensembles from Doppler beams 1 & 2
                for j = 1:P.PRIsDop
                    for k = 1:2
                        Event(n).info = 'Acquire two interleaved Doppler ensemble';
                        Event(n).tx = P.numRays+k;
                        Event(n).rcv = rcvInd+P.numRays+(k-1)*P.PRIsDop+j;
                        Event(n).recon = 0;
                        Event(n).process = 0;
                        Event(n).seqControl = 6;
                        n = n+1;
                    end
                end
                % Aquire ensembles from Doppler beams 3 & 4
                for j = 1:P.PRIsDop
                    for k = 3:4
                        Event(n).info = 'Acquire two interleaved Doppler ensemble';
                        Event(n).tx = P.numRays+k;
                        Event(n).rcv = rcvInd+P.numRays+(k-1)*P.PRIsDop+j;
                        Event(n).recon = 0;
                        Event(n).process = 0;
                        Event(n).seqControl = 6;
                        n = n+1;
                    end
                end
            else % numIntLv = 1
                % Acquire non-interleaved ensembles from each beam.
                for k = 1:numRaysDop
                    for j = 1:P.PRIsDop
                        Event(n).info = 'Acquire non-interleaved Doppler ensemble';
                        Event(n).tx = P.numRays+k;
                        Event(n).rcv = rcvInd+P.numRays+(k-1)*P.PRIsDop+j;
                        Event(n).recon = 0;
                        Event(n).process = 0;
                        Event(n).seqControl = 6;
                        n = n+1;
                    end
                end
            end
            Event(n-1).seqControl = [7,8,nsc]; % replace last Doppler acquisition Event's seqControl
            SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
            nsc = nsc+1;
        end
        assignin('base','Event',Event);
    end
    % Set Doppler PRF in Doppler Process structure in Matlab workspace for reference.
    for k = 1:2:length(Process(2).Parameters)
        if strcmp(Process(2).Parameters{k},'prf'), Process(2).Parameters{k+1} = P.PRFDop; end
    end

    assignin('base','P',P);
    assignin('base','Process',Process);
    assignin('base','SeqControl',SeqControl);

    Control(1).Command = 'set&Run';
    Control(1).Parameters = {'Process',2,'prf',P.PRFDop};
    Control(2).Command = 'update&Run';
    if numIntLv ~= evalin('base','numIntLv')
        Control(2).Parameters = {'SeqControl','Event'};
    else
        Control(2).Parameters = {'SeqControl'};
    end
    assignin('base','numIntLv',numIntLv);
    assignin('base','Control',Control);
end

function SteerAngleCallback(~,~,UIValue)
%SteerAngle - Doppler Angle will change TX beam and PData

    P = evalin('base','P');
    TX = evalin('base','TX');
    PData = evalin('base','PData');

    P.angleDop = UIValue * pi/180;
    assignin('base','P',P);
    numRaysDop = evalin('base','numRaysDop');
    for n = 1:numRaysDop+1
        PData(2).Region(n).Shape.angle = P.angleDop;
    end
    PData(2).Region = computeRegions(PData(2));
    assignin('base','PData',PData);

    for n = 1:numRaysDop
        TX(P.numRays+n).Steer = [P.angleDop,0.0];
        TX(P.numRays+n).Delay = computeTXDelays(TX(P.numRays+n));
        TX(P.numRays+n).TXPD = computeTXPD(TX(P.numRays+n),PData(2));
    end
    assignin('base','TX',TX);

    Control.Command = 'update&Run';
    Control.Parameters = {'PData','InterBuffer','ImageBuffer','TX','Recon'};
    assignin('base','Control',Control);
end

function DopplerDepthCallback(hObject,~,UIValue)
    Trans = evalin('base','Trans');
    Resource = evalin('base','Resource');
    mm2wl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);
    P = evalin('base','P');
    % Determine if axis units in mm, instead of default units in wavelengths.
    if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
        if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
            axesUnits = 'mm';
            endDepthDop = UIValue*mm2wl;
        else
            axesUnits = 'wl';
            endDepthDop = UIValue;
        end
    end
    % No depth change if in simulate mode 2.
    if Resource.Parameters.simulateMode == 2
        if strcmp(axesUnits,'mm')
            set(hObject,'Value',P.endDepthDop/mm2wl);
        else
            set(hObject,'Value',P.endDepthDop);
        end
        return
    end
    if endDepthDop > P.endDepth
        endDepthDop = P.endDepth;
        if strcmp(axesUnits,'mm')
            set(hObject,'Value',endDepthDop/mm2wl);
        else
            set(hObject,'Value',endDepthDop);
        end
    end
    P.endDepthDop = endDepthDop;
    % PData
    PData = evalin('base','PData');
    numRaysDop = evalin('base','numRaysDop');
    for n = 1:numRaysDop+1, PData(2).Region(n).Shape.height = P.endDepthDop - P.startDepthDop; end
    PData(2).Region = computeRegions(PData(2));
    assignin('base','PData',PData);
    % Update Doppler Receive structures.
    Receive = evalin('base', 'Receive');
    % - Calculate maxAcqLngthdop as the diagonal of the color box from the law of cosines.
    %      c^2 = a^2 + b^2 - 2*a*b*cos(angle)
    a = PData(2).Region(numRaysDop+1).Shape.width;
    b = P.endDepthDop/cos(P.angleDop);
    angle = P.angleDop + pi/2;  % for P.dopAngle positive.
    maxAcqLngthDop =  sqrt(a^2 + b^2 - 2*a*b*cos(angle))-P.startDepthDop;
    % - Adjust maxAcqLngthDop, since hardware needs multiple of 128 samples.
    nSmpls = 2*maxAcqLngthDop * numRaysDop/2; % no. of acquisition samples for BS100BW
    if abs(round(nSmpls/128) - nSmpls/128) < .01 % round to nearest 128 sample boundary for hardware.
        numRcvSamples = 128*round(nSmpls/128);
    else
        numRcvSamples = 128*ceil(nSmpls/128); % round up to next 128 sample boundary.
    end
    if numRcvSamples == 0, numRcvSamples = 128; end
    maxAcqLngthDop = P.startDepthDop + (numRcvSamples/2)/2; % /2 for wvlngths for BS100BW acquisitions
    for i = 1:Resource.RcvBuffer(1).numFrames
        k = (P.numRays + numRaysDop*P.PRIsDop)*(i-1); % k keeps track of Receive index increment per frame.
        for j = (P.numRays+1):(P.numRays+P.PRIsDop*numRaysDop) % Doppler acquisition
            Receive(j+k).endDepth = maxAcqLngthDop;
        end
    end
    assignin('base','Receive',Receive);
    %---------------check if Doppler numIntLv needs modifying--------------------
    Process = evalin('base','Process');
    SeqControl = evalin('base','SeqControl');
    Event = evalin('base','Event');
    demodFreqDop = evalin('base','demodFreqDop');
    % - Use Doppler PRF to determine interleave factor (4, 2 or 1)
    tof = ceil(2*Receive(P.numRays+1).endDepth/demodFreqDop); % round trip time of flight for Doppler (in usecs)
    pri = round(1e+06/P.PRFDop);  % Doppler PRI in usecs
    if 4*tof < pri
        numIntLv = 4;
    elseif 2*tof < pri
        numIntLv = 2;
    else
        numIntLv = 1;
    end
    SeqControl(6).command = 'timeToNextAcq';
    SeqControl(6).argument = round(1/(P.PRFDop*numIntLv) * 1e+06);
    % If numIntLv changed, modify Event structures for Doppler acquisition.
    if numIntLv ~= evalin('base','numIntLv')
        startDopEvents = evalin('base','startDopEvents');
        nsc = 11;
        for i = 1:Resource.RcvBuffer(1).numFrames
            rcvInd = (P.numRays + numRaysDop*P.PRIsDop)*(i-1);
            n = startDopEvents(i);
            if numIntLv == 4
                for j = 1:P.PRIsDop
                    for k = 1:numRaysDop
                        Event(n).info = 'Acquire four interleaved Doppler ensemble';
                        Event(n).tx = P.numRays+k;
                        Event(n).rcv = rcvInd+P.numRays+(k-1)*P.PRIsDop+j;
                        Event(n).recon = 0;
                        Event(n).process = 0;
                        Event(n).seqControl = 6;
                        n = n+1;
                    end
                end
            elseif numIntLv == 2
                % Aquire ensembles from Doppler beams 1 & 2
                for j = 1:P.PRIsDop
                    for k = 1:2
                        Event(n).info = 'Acquire two interleaved Doppler ensemble';
                        Event(n).tx = P.numRays+k;
                        Event(n).rcv = rcvInd+P.numRays+(k-1)*P.PRIsDop+j;
                        Event(n).recon = 0;
                        Event(n).process = 0;
                        Event(n).seqControl = 6;
                        n = n+1;
                    end
                end
                % Aquire ensembles from Doppler beams 3 & 4
                for j = 1:P.PRIsDop
                    for k = 3:4
                        Event(n).info = 'Acquire two interleaved Doppler ensemble';
                        Event(n).tx = P.numRays+k;
                        Event(n).rcv = rcvInd+P.numRays+(k-1)*P.PRIsDop+j;
                        Event(n).recon = 0;
                        Event(n).process = 0;
                        Event(n).seqControl = 6;
                        n = n+1;
                    end
                end
            else % numIntLv = 1
                % Acquire non-interleaved ensembles from each beam.
                for k = 1:numRaysDop
                    for j = 1:P.PRIsDop
                        Event(n).info = 'Acquire non-interleaved Doppler ensemble';
                        Event(n).tx = P.numRays+k;
                        Event(n).rcv = rcvInd+P.numRays+(k-1)*P.PRIsDop+j;
                        Event(n).recon = 0;
                        Event(n).process = 0;
                        Event(n).seqControl = 6;
                        n = n+1;
                    end
                end
            end
            Event(n-1).seqControl = [7,8,nsc]; % replace last Doppler acquisition Event's seqControl
            SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
            nsc = nsc+1;
        end
        assignin('base','Event',Event);
    end
    % Set Doppler PRF in Doppler Process structure in Matlab workspace for reference.
    for k = 1:2:length(Process(2).Parameters)
        if strcmp(Process(2).Parameters{k},'prf'), Process(2).Parameters{k+1} = P.PRFDop; end
    end
    assignin('base','P',P);
    assignin('base','Process',Process);
    assignin('base','SeqControl',SeqControl);
    % TGC
    TGC = evalin('base','TGC');
    TGC(2).rangeMax = P.endDepthDop;
    TGC(2).Waveform = computeTGCWaveform(TGC(2));
    assignin('base','TGC',TGC);

    Control(1).Command = 'set&Run';
    Control(1).Parameters = {'Process',2,'prf',P.PRFDop};
    Control(2).Command = 'update&Run';
    if numIntLv ~= evalin('base','numIntLv')
        Control(2).Parameters = {'PData','InterBuffer','ImageBuffer','Receive','TGC','Recon','SeqControl','Event'};
        assignin('base','numIntLv',numIntLv);
    else
        Control(2).Parameters = {'PData','InterBuffer','ImageBuffer','Receive','TGC','Recon','SeqControl'};
    end
    assignin('base','Control', Control);
end

function RangeChangeCallback(hObject,~,UIValue)
%RangeChangeCallback - range change
% Only the 2D depth is changed, unless the 2D depth becomes less than the Doppler
% depth, in which case both are changed.

    Trans = evalin('base','Trans');
    Resource = evalin('base','Resource');
    mm2wl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);
    P = evalin('base','P');
    % Determine if axis units in mm, instead of default units in wavelengths.
    if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
        if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
            axesUnits = 'mm';
            endDepth = UIValue*mm2wl;
        else
            axesUnits = 'wl';
            endDepth = UIValue;
        end
    end
    % No depth change if in simulate mode 2.
    if Resource.Parameters.simulateMode == 2
        if strcmp(axesUnits,'mm')
            set(hObject,'Value',P.endDepth/mm2wl);
        else
            set(hObject,'Value',P.endDepth);
        end
        return
    end
    P.endDepth = endDepth;
    % If P.endDepthDop > P.endDepth, modify P.endDepthDop to equal P.endDepth
    modDopDepth = 0;
    if P.endDepthDop > P.endDepth
        P.endDepthDop = P.endDepth;
        h1 = findobj('Tag','UserC1Slider');
        h2 = findobj('Tag','UserC1Edit');
        if strcmp(axesUnits,'mm')
            set(h1,'Value',P.endDepthDop/mm2wl);
            set(h2,'String',num2str(P.endDepthDop/mm2wl,get(h2,'UserData')));
        else
            set(h1,'Value',P.endDepthDop);
            set(h2,'String',num2str(P.endDepthDop,get(h2,'UserData')));
        end
        modDopDepth = 1;
    end
    P.txFocus = 4*P.endDepth;
    numRaysDop = evalin('base','numRaysDop');

    % PData
    PData = evalin('base','PData');
    PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3)); % rows
    for n = 1:P.numRays, PData(1).Region(n).Shape.height = P.endDepth - P.startDepth; end
    PData(1).Region = computeRegions(PData(1));

    PData(2).Size(1) = PData(1).Size(1);
    for n = 1:numRaysDop+1, PData(2).Region(n).Shape.height = P.endDepthDop - P.startDepthDop; end
    PData(2).Region = computeRegions(PData(2));

    assignin('base','PData',PData);
    evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);');

    % Update TXPD data of TX structures.
    TX = evalin('base','TX');
    h = waitbar(0,'Program TX parameters, please wait!');
    for i = 1:P.numRays+numRaysDop
        TX(i).focus = P.txFocus;
        TX(i).Delay = computeTXDelays(TX(i));
        if i <= P.numRays
            TX(i).TXPD = computeTXPD(TX(i),PData(1));
        else
            TX(i).TXPD = computeTXPD(TX(i),PData(2));
        end
        waitbar(i/(P.numRays+numRaysDop))
    end
    close(h)
    assignin('base','TX',TX);
    % Receive
    Receive = evalin('base', 'Receive');
    maxAcqLngth2D = sqrt(P.endDepth^2 + (Trans.numelements*Trans.spacing)^2)-P.startDepth;
    if modDopDepth==1
        % - Calculate maxAcqLngthdop as the diagonal of the color box from the law of cosines.
        %      c^2 = a^2 + b^2 - 2*a*b*cos(angle)
        a = PData(2).Region(numRaysDop+1).Shape.width;
        b = P.endDepthDop/cos(P.angleDop);
        angle = P.angleDop + pi/2;  % for P.dopAngle positive.
        maxAcqLngthDop =  sqrt(a^2 + b^2 - 2*a*b*cos(angle))-P.startDepthDop;
    end
    for i = 1:Resource.RcvBuffer(1).numFrames
        k = (P.numRays + numRaysDop*P.PRIsDop)*(i-1); % k keeps track of Receive index increment per frame.
        for j = 1:P.numRays  % acquisitions for 2D
            Receive(j+k).endDepth = maxAcqLngth2D;
        end
        if modDopDepth==1
            for j = (P.numRays+1):(P.numRays+P.PRIsDop*numRaysDop) % Doppler acquisition
                Receive(j+k).endDepth = maxAcqLngthDop;
            end
        end
    end
    assignin('base','Receive',Receive);
    if modDopDepth==1
        %------check if numIntLv needs modifying, due to decreased Doppler depth-------------
        Process = evalin('base','Process');
        SeqControl = evalin('base','SeqControl');
        Event = evalin('base','Event');
        demodFreqDop = evalin('base','demodFreqDop');
        % - Use Doppler PRF to determine interleave factor (4, 2 or 1)
        tof = ceil(2*Receive(P.numRays+1).endDepth/demodFreqDop); % round trip time of flight for Doppler (in usecs)
        pri = round(1e+06/P.PRFDop);  % Doppler PRI in usecs
        if 4*tof < pri
            numIntLv = 4;
        elseif 2*tof < pri
            numIntLv = 2;
        else
            numIntLv = 1;
        end
        SeqControl(6).command = 'timeToNextAcq';
        SeqControl(6).argument = round(1/(P.PRFDop*numIntLv) * 1e+06);
        % If numIntLv changed, modify Event structures for Doppler acquisition.
        if numIntLv ~= evalin('base','numIntLv')
            startDopEvents = evalin('base','startDopEvents');
            nsc = 11;
            for i = 1:Resource.RcvBuffer(1).numFrames
                rcvInd = (P.numRays + numRaysDop*P.PRIsDop)*(i-1);
                n = startDopEvents(i);
                if numIntLv == 4
                    for j = 1:P.PRIsDop
                        for k = 1:numRaysDop
                            Event(n).info = 'Acquire four interleaved Doppler ensemble';
                            Event(n).tx = P.numRays+k;
                            Event(n).rcv = rcvInd+P.numRays+(k-1)*P.PRIsDop+j;
                            Event(n).recon = 0;
                            Event(n).process = 0;
                            Event(n).seqControl = 6;
                            n = n+1;
                        end
                    end
                elseif numIntLv == 2
                    % Aquire ensembles from Doppler beams 1 & 2
                    for j = 1:P.PRIsDop
                        for k = 1:2
                            Event(n).info = 'Acquire two interleaved Doppler ensemble';
                            Event(n).tx = P.numRays+k;
                            Event(n).rcv = rcvInd+P.numRays+(k-1)*P.PRIsDop+j;
                            Event(n).recon = 0;
                            Event(n).process = 0;
                            Event(n).seqControl = 6;
                            n = n+1;
                        end
                    end
                    % Aquire ensembles from Doppler beams 3 & 4
                    for j = 1:P.PRIsDop
                        for k = 3:4
                            Event(n).info = 'Acquire two interleaved Doppler ensemble';
                            Event(n).tx = P.numRays+k;
                            Event(n).rcv = rcvInd+P.numRays+(k-1)*P.PRIsDop+j;
                            Event(n).recon = 0;
                            Event(n).process = 0;
                            Event(n).seqControl = 6;
                            n = n+1;
                        end
                    end
                else % numIntLv = 1
                    % Acquire non-interleaved ensembles from each beam.
                    for k = 1:numRaysDop
                        for j = 1:P.PRIsDop
                            Event(n).info = 'Acquire non-interleaved Doppler ensemble';
                            Event(n).tx = P.numRays+k;
                            Event(n).rcv = rcvInd+P.numRays+(k-1)*P.PRIsDop+j;
                            Event(n).recon = 0;
                            Event(n).process = 0;
                            Event(n).seqControl = 6;
                            n = n+1;
                        end
                    end
                end
                Event(n-1).seqControl = [7,8,nsc]; % replace last Doppler acquisition Event's seqControl
                SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
                nsc = nsc+1;
            end
            assignin('base','Event',Event);
        end
        % Set Doppler PRF in Doppler Process structure in Matlab workspace for reference.
        for k = 1:2:length(Process(2).Parameters)
            if strcmp(Process(2).Parameters{k},'prf'), Process(2).Parameters{k+1} = P.PRFDop; end
        end
    end
    assignin('base','P',P);
    % TGC
    TGC = evalin('base','TGC');
    TGC(1).rangeMax = P.endDepth;
    TGC(1).Waveform = computeTGCWaveform(TGC(1));
    assignin('base','TGC',TGC);
    % Set Control structure
    Control(1).Command = 'update&Run';
    if (modDopDepth==1)&&(numIntLv ~= evalin('base','numIntLv')) % if interleave factor changed
        Control(1).Parameters = {'PData','InterBuffer','ImageBuffer','DisplayWindow','TX','Receive','TGC','Recon','SeqControl','Event'};
        assignin('base','numIntLv',numIntLv);
    else
        Control(1).Parameters = {'PData','InterBuffer','ImageBuffer','DisplayWindow','TX','Receive','TGC','Recon'};
    end
    assignin('base','Control', Control);
    assignin('base', 'action', 'displayChange');
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
        evalin('base','hROI = drawRegionOutline(1,2,numRaysDop+1);')
    else
        evalin('base','drawRegionOutline(hROI,1,2,numRaysDop+1);')
    end

end

