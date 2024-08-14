% Notice: 
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility 
%   for its use.
%
% File name SetUpL39_21gDWideBeamDoppler.m:
%    Generate .mat Sequence Object file for L39-21gD Linear array 
%    128 transmit and 128 receive channels are used.
%
% Description: This script uses the "interleaved sampling" mechanism for 2D RF
%    data acquisition.  The L39-21gD is operated at a nominal center frequency
%    for Recon processing of 31.25 MHz, with 4X sampling so the required RF data
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
%       Doppler acquisition uses a custom pulse code that transmits at 26.785MHz,
%    which is 4/3 of the A/D sample rate of 35.7143.
%
% Last update:
%    April 2021 VTS-2152 computeTrans entry for L39-21gD
%    03-11-2021 script for new Daxsonics L39-21gD
%    08-03-2018 Adopted color box programming from L38-22v widebeam Doppler script.
%    04-27-2018 Modified custom pulse code. Added new Doppler attributes for
%       optimization.

%clear[ 	]+all

% Set 2D parameters
P.numTx = 40;   % no. of elements in TX aperture.
P.numRays = 64; % no. of rays in frame
P.startDepth = 0;   % Acquisition depth in wavelengths
P.endDepth = 192;  % This should preferrably be a multiple of 128 samples (4 samples/wv).
P.txFocus = P.endDepth+10; % for 2D

% Set Doppler parameters
P.dopStartDepth = P.startDepth;
P.dopEndDepth = P.endDepth;
P.dopNumRays = 4; % no. of rays in frame
P.dopRegPos = [-10.0, 5]; % (x,z)
P.dopRegWidth = 100;  % L30 aperture is 114 wavelengths
P.dopAngle = 5*pi/180; %12*pi/180;
P.dopNumTx = 40;  % no. of elements in TX aperture
P.dopFocus = -3*P.endDepth; % negative focus for slightly diverging beam
P.dopPRIs = 18;
P.dopPRF = 1.5e+03; % Doppler PRF in Hz.
P.pwrThres = 0.1;
P.persfreq = 70;
P.perspwr = 90;
P.cpl = 40;  % level at which 2D overwrites Doppler

m = P.dopNumRays;
xd = P.dopRegPos(1);
zd = P.dopRegPos(2);
widthd = P.dopRegWidth;
heightd = P.dopEndDepth - P.dopStartDepth;
angle = P.dopAngle;

% Specify system parameters.
Resource.Parameters.connector = 1;
Resource.Parameters.numTransmit = 128;  % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;  % number of receive channels.
Resource.Parameters.speedOfSound = 1540;
Resource.Parameters.simulateMode = 0; % 0 means no simulation, if hardware is present.
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Specify Trans structure array.
Trans.name = 'L39-21gD';
Trans.units = 'mm';
Trans = computeTrans(Trans);

% Specify PData structure arrays.
% - 2D PData structure
PData(1).PDelta = [Trans.spacing, 0, 0.5];
PData(1).Size(1,1) = ceil((P.endDepth-5)/PData(1).PDelta(3)); % range, origin and pdelta set PData(1).Size.
PData(1).Size(1,2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(1,3) = 1;      % single image page
PData(1).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,5]; % x,y,z of uppr lft crnr.
% Define P.numRays rectangular regions centered on TX beam origins.
txdx = 127*Trans.spacing/(P.numRays-1);
TxOrgX = (-63.5*Trans.spacing):(127*Trans.spacing/(P.numRays-1)):(63.5*Trans.spacing);
rayDelta = Trans.numelements*Trans.spacing/P.numRays;
% - specify P.numRays default Region structures.
PData(1).Region = repmat(struct('Shape',struct('Name','Trapezoid',...
                                               'Position',[0,0,0],...
                                               'top',P.numTx*Trans.spacing,...
                                               'bottom',2*txdx,...
                                               'height',P.endDepth-5)),1,P.numRays);
% - set position and width of regions to correspond to transmit origins and beam spacing.
for n = 1:P.numRays, PData(1).Region(n).Shape.Position(1) = TxOrgX(n); end
PData(1).Region = computeRegions(PData(1));

% - Doppler PData structure
PData(2) = PData(1);
PData(2).PDelta = [Trans.spacing, 0, Trans.spacing];
PData(2).Size(1) = ceil((P.endDepth-5)/PData(2).PDelta(3)); % rows
PData(2).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(2).PDelta(1)); % cols
TxDopOrgX = (xd-widthd*(m-1)/(2*m)):(widthd/m):(xd+widthd*(m-1)/(2*m));
PData(2).Region = repmat(struct('Shape',struct(...
    'Name','Parallelogram',...
    'Position',[xd,0,zd],... %will be changed later
    'width', widthd/m,...
    'height', heightd,...
    'angle',angle)),1,m); % define m regions and one extra for full window
for n = 1:m    
    PData(2).Region(n).Shape.Position(1) = TxDopOrgX(n);
end
PData(2).Region = computeRegions(PData(2));

% - PData(3) is used to outline the doppler area
PData(3) = PData(2); PData(3).Region = [];
PData(3).Region.Shape = struct(...
    'Name','Parallelogram',...
    'Position',[TxDopOrgX(1)-PData(2).Region(1).Shape.width/2+44*Trans.spacing,0,P.dopStartDepth],...
    'width', 88*Trans.spacing,...
    'height', P.dopEndDepth-P.dopStartDepth,...
    'angle',P.dopAngle);
PData(3).Region = computeRegions(PData(3));

% Specify Media object.
pt1;
Media.attenuation = 0.0;
Media.function = 'movePoints';

% Specify Resources.
% - RcvBuffer(1) is for both 2D and Doppler acquisitions.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 2*400*8*(P.numRays + m*P.dopPRIs); 
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = numRcvFrames;           % 20 frames allocated for RF acqusitions.
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
Resource.DisplayWindow(1).Title = 'L39-21gD WideBeamDoppler';
Resource.DisplayWindow(1).pdelta = 0.4;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Type = 'Matlab';
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)]; % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).numFrames = numImageFrames;
Resource.DisplayWindow(1).AxesUnits = 'wavelengths';
Resource.DisplayWindow(1).Colormap = grayscaleCFImap;
Resource.DisplayWindow(1).splitPalette = 1;

% Specify Transmit waveform structures.  
% - 2D transmit waveform
TW(1).type = 'parametric';
TW(1).Parameters = [31.25,0.67,2,1];   % A, B, C, D
% - Doppler transmit waveform, transmit frequency should be equivalent to
% supported demodFrequency
TW(2).type = 'pulseCode'; % Pulse code for ~26.786 MHz (3/4)*(250/7)
TW(2).PulseCode =  [ 0   5   1  -4   1; ...
                     0   4   0  -5   1; ... 
                     0   5   0  -4   1; ... 
                     1   4   0  -5   1; ... 
                     0   5   0  -4   1; ... 
                     0   4   0  -5   1; ... 
                     0   5   1  -4   1; ... 
                     0   4   0  -5   1; ... 
                     0   5   0  -4   1; ... 
                     1   4   0  -5   1; ... 
                     0   5   0  -4   1];

scaleToWvl = 1;
if strcmp(Trans.units, 'mm')
    scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);
end 

% Specify TX structure array.
%   2*P.numRays TXs are required for 2D, and P.dopNumRays (m) are needed for Doppler.
%   The Doppler beams TXPD data are computed with PData(2).
%   The TXs are computed as follows:
%   WB(1)a,WB(2)a,WB(3)a,...,WB(P.numRays)a,WB(1)b,WB(2)b,WB(3)b,...,WB(P.numRays)b,
%      DB(1),DB(2),DB(3),...,DB(P.dopNumRays)
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'Apod', zeros(1,Resource.Parameters.numTransmit), ...
                   'focus', P.txFocus, ...
                   'Steer', [0.0,0.0], ...
                   'Delay', zeros(1,Resource.Parameters.numTransmit),...
                   'TXPD', [], ...
                   'peakCutOff', 0.3, ...
                   'peakBLMax', 27.0), 1, 2*P.numRays+m);

h = waitbar(0,'Program TX parameters, please wait!');                            
% - Set event specific TX attributes.
for n = 1:P.numRays  % specify P.numRays transmit events
    TX(n).Origin(1) = TxOrgX(n);
    % Compute closest element number to transmit origin
    [Dummy,ce] = min(abs(scaleToWvl*Trans.ElementPos(:,1)-TX(n).Origin(1))); % ce is closest ele to cntr of aper.
    % Compute TX.Apod.
    lft = round(ce - P.numTx/2);
    if lft < 1, lft = 1; end
    rt = round(ce + P.numTx/2);
    if rt > 128, rt = 128; end
    TX(n).Apod(lft:rt) = 1;
    % Apply apodization function.
    [RIndices,CIndices,V] = find(TX(n).Apod);
    V = kaiser(size(V,2),1);
    TX(n).Apod(CIndices) = V;
    % Define 2nd TXs for interleave.
    TX(n+P.numRays) = TX(n);
    % Compute transmit delays for normal and interleave
    TX(n+P.numRays).Delay = computeTXDelays(TX(n+P.numRays));
    % Add the delay offset representing one half of A/D sample period (1/62.5MHz = 16ns) to first acq's delay.
    TX(n).Delay = TX(n+P.numRays).Delay + 0.5*Trans.frequency/62.5; % 0.24 wls of Trans.frequency=31.25 equivalent to 8 nsec. 
    % Compute transmit pixel data only for the first TX in the interleave pair.
    TX(n).TXPD = computeTXPD(TX(n),PData(1));
    waitbar(n/P.numRays);            
end
close(h)
% -- P.dopNumRays TX structs needed for Doppler
k = 2*P.numRays;
% TX(1) to TX(2*P.numRays) are used for 2D
% TX(2*P.numRays+1) to TX(2*P.numRays+2*m) are used for Doppler
for n = 1:m
    TX(n+k).waveform = 2;
    TX(n+k).Origin = [TxDopOrgX(n),0.0,0.0];
    TX(n+k).focus = P.dopFocus;
    % Compute closest element number to transmit origin
    [Dummy,ce] = min(abs(scaleToWvl*Trans.ElementPos(:,1)-TX(n+k).Origin(1))); % ce is closest ele to cntr of aper.
    % Compute TX.Apod.
    lft = round(ce - P.dopNumTx/2);
    if lft < 1, lft = 1; end
    rt = round(ce + P.dopNumTx/2);
    if rt > 128, rt = 128; end
    TX(n+k).Apod(lft:rt) = 1;
%     % Apply apodization function.
%     [RIndices,CIndices,V] = find(TX(n+k).Apod);
%     V = kaiser(size(V,2),1);
%     TX(n+k).Apod(CIndices) = V;
    TX(n+k).Steer = [P.dopAngle,0.0];
    TX(n+k).Delay = computeTXDelays(TX(n+k));
    TX(n+k).TXPD = computeTXPD(TX(n+k),PData(2));
end

% Specify TPC structures.
TPC(1).name = '2D';
TPC(1).maxHighVoltage = 40;
TPC(2).name = 'Doppler';
TPC(2).maxHighVoltage = 40;

RcvProfile(1).LnaGain = 18; % Profile used for imaging
RcvProfile(2).LnaGain = 24; % Profile used for Doppler
RcvProfile(2).PgaGain = 30;
RcvProfile(2).LnaZinSel = 15; % enable high pass filter

% Specify TGC Waveform structure.
% - 2D TGC
TGC(1).CntrlPts = [0, 271, 498, 617, 767, 903, 1000, 1023];
TGC(1).rangeMax = P.endDepth;
TGC(1).Waveform = computeTGCWaveform(TGC(1));
% - Doppler TGC
TGC(2).CntrlPts =[ 150 570 713 713 662 610 610 610];
TGC(2).rangeMax = P.dopEndDepth;
TGC(2).Waveform = computeTGCWaveform(TGC(2));

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
% - Note that for the L35-16v we would ideally use a bandpass filter centered
% near the transducer's nominal center frequency of 28 MHz, not the 31.25
% MHz forced on the CGD filters.
maxAcqLngth2D = sqrt(P.endDepth^2 + (Trans.numelements*Trans.spacing)^2);
maxAcqLngthDop =  sqrt((P.dopEndDepth/cos(P.dopAngle))^2 + (Trans.numelements*Trans.spacing)^2);

% interleave needs a high pass filter instead of a bandpass one 
HPF = [-0.00082 -0.00101 +0.00299 -0.00150 -0.00366 +0.00656 -0.00110 ...
       -0.00922 +0.01138 +0.00208 -0.01913 +0.01678 +0.01099 -0.03647 ...
       +0.02179 +0.03351 -0.07281 +0.02533 +0.12036 -0.28638 +0.35950];
% For Doppler, we use a 35% BW bandpass filter center at the aliased center frequency for 4/3 sampling.
BandPassCoef1 = [-0.00021 +0.00000 +0.00058 +0.00000 -0.00027 +0.00000 -0.00250 ...
                 +0.00000 +0.01065 +0.00000 -0.02682 +0.00000 +0.05188 +0.00000 ...
                 -0.08319 +0.00000 +0.11462 +0.00000 -0.13812 +0.00000 +0.14679];
% 20% BW filter
BandPassCoef2 = [-0.0001 -0.0000 -0.0012 +0.0002 +0.0045 -0.0006 -0.0115 ...
                 +0.0011 +0.0235 -0.0017 -0.0405 +0.0023 +0.0610 -0.0025 ...
                 -0.0825 +0.0023 +0.1013 -0.0016 -0.1143 +0.0006 +0.1189];

Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLngth2D, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BWI', ...
                        'demodFrequency', 31.25, ...
                        'InputFilter', HPF, ... 
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, (2*P.numRays+m*P.dopPRIs)*Resource.RcvBuffer(1).numFrames);
% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    k = (2*P.numRays + m*P.dopPRIs)*(i-1); % k keeps track of Receive index increment per frame.
    % - Set attributes for each frame.
    Receive(k+1).callMediaFunc = 1;
    for j = 1:2:2*P.numRays  % acquisitions for 2D
        % first Receive of interleaved samples pair.
        Receive(k+j).framenum = i;
        Receive(k+j).acqNum = j;
        % 2nd Receive of interleaved samples pair
        Receive(k+j+1).framenum = i;
        Receive(k+j+1).acqNum = j+1;        
    end
    n = 2*P.numRays;
    % Doppler acquisitions
    for j = 1:P.dopPRIs*m
        % Doppler acquisition      
        Receive(k+n+j).startDepth = P.dopStartDepth;
        Receive(k+n+j).endDepth = maxAcqLngthDop;
        Receive(k+n+j).TGC = 2;
        Receive(k+n+j).sampleMode = 'BS67BW';
        Receive(k+n+j).demodFrequency = 26.7858;
        Receive(k+n+j).InputFilter = BandPassCoef2;
        Receive(k+n+j).framenum = i;
        Receive(k+n+j).acqNum = n+j;        % Doppler acqNums continue after 2D
    end
end                    

% Specify Recon structure arrays.
% - We need two Recon structures, one for 2D, one for Doppler. These will be referenced in the same
%   event, so that they will use the same (most recent) acquisition frame.
Recon = repmat(struct('senscutoff', 0.5, ...
               'pdatanum', 1, ...
               'rcvBufFrame', -1, ...     % use most recently transferred frame
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...  % auto-increment ImageBuffer each recon
               'RINums', zeros(1,1)), 1, 2);
% - Set ReconInfo for 2D frame.
Recon(1).RINums = (1:P.numRays);  % P.numRays ReconInfos needed for 2D
% - Set ReconInfo for Doppler ensemble.
Recon(2).pdatanum = 2;
Recon(2).IntBufDest = [2,1];
Recon(2).ImgBufDest = [0,0];
Recon(2).RINums = ((P.numRays+1):(P.numRays+m*P.dopPRIs));

% Define ReconInfo structures.
% - For 2D, we need P.numRays ReconInfo structures for P.numRays.
% - For Doppler, we need m*P.dopPRIs ReconInfo structures.
ReconInfo = repmat(struct('mode', 'accumIQ', ...    % accumulate IQ data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'pagenum',1, ...
                   'scaleFactor', 1.0, ...
                   'regionnum', 1,...
                   'threadSync',1), 1, P.numRays + m*P.dopPRIs);
% - Set specific ReconInfo attributes.
%   - ReconInfos for 2D frame.
ReconInfo(1).Pre = 'clearInterBuf';
for j = 1:P.numRays
    ReconInfo(j).txnum = j;
    ReconInfo(j).rcvnum = 2*j-1; % 1st receive of each interleaved sample pair
    ReconInfo(j).regionnum = j;
end
ReconInfo(P.numRays).Post = 'IQ2IntensityImageBuf';

%  - ReconInfos for Doppler ensembles.
ReconInfo(P.numRays+1).Pre = 'clearInterBuf';
for n = 1:m  % for each Doppler beam
    k = P.numRays + P.dopPRIs*(n-1);
    for j = 1:P.dopPRIs
        ReconInfo(k+j).mode = 'replaceIQ_normalize';
        ReconInfo(k+j).txnum = 2*P.numRays + n;
        ReconInfo(k+j).rcvnum = 2*P.numRays + m*(j-1) + n;
        ReconInfo(k+j).regionnum = n;
        ReconInfo(k+j).pagenum = j;
        ReconInfo(k+j).scaleFactor = 1;
    end
end

% Specify Process structure arrays.
DopState = 'freq';

Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',3.0,...            % pgain is image processing gain
                         'reject',2,...      % reject level 
                         'persistMethod','simple',...
                         'persistLevel',20,...
                         'interpMethod','4pt',...  %method of interp. (1=4pt)
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
Process(2).Parameters = {'IntBufSrc',[2,1],...          % number,frame of IQ (Inter)buffer to process.
                         'SrcPages',[3,P.dopPRIs-2],... % start page number, num. pages in IQ buffer                         
                         'ImgBufDest',[2,-1],... 
                         'pdatanum',2,...             % number of PData structure
                         'prf',P.dopPRF,...           % Doppler PRF in Hz
                         'wallFilter','regression',...  % 1 -> quadratic regression
                         'pwrThreshold',P.pwrThres,...
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
                         'persistMethod','simple',...
                         'persistLevel',P.persfreq,...
                         'interpMethod','4pt',...  %method of interp.
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
% TTNI (Time To Next Interleave pulse) is determined by the image depth. 
% In this case, the Trans.frequency is 30MHz, and the
% doppler transmit frequency is 25 MHz. P.endDepth is 200 wavelengths (~10mm), 
% so the round trip travel time is 2*200/25 = 16 us. 
% Therefore, 20 us TTNI should be enough and the Doppler PRF won't be 
% depenent on the depth if the maximum depth is ~15mm.

TTNI = 40; % Time To Next Interleave pulse
TTNW = 100 - TTNI; % Time To Next WibeBeam pulse after interleave
% - set receive profile for 2D
SeqControl(1).command = 'setRcvProfile';        
SeqControl(1).argument = 1;  % 
% - time between wide beam interleave acquisitions
SeqControl(2).command = 'timeToNextAcq'; 
SeqControl(2).argument = TTNI;
% - time to next wide beam acquisition after interleave 
SeqControl(3).command = 'timeToNextAcq'; 
SeqControl(3).argument = TTNW;
% -- Time between 2D acquisition and Doppler ensemble. Set to allow time for profile change.
SeqControl(4).command = 'timeToNextAcq';
SeqControl(4).argument = 2500; % time in usec
% -- Change to TPC Profile 2 (Doppler)
SeqControl(5).command = 'setTPCProfile';
SeqControl(5).condition = 'next';
SeqControl(5).argument = 2;
% - set Receive profile for Doppler
SeqControl(6).command = 'setRcvProfile';        
SeqControl(6).argument = 2;  % 
% - time between interleaved Doppler beam pulses. All Doppler beams are
%   acquired at the PRF rate.  e.g. if dopPRF = 2000 and no. of beams, m,
%   equals 4, the time between beams would be (1/2000)/4.
SeqControl(7).command = 'timeToNextAcq';
SeqControl(7).argument = round(1e6/(P.dopPRF*m)); 
% -- Change to Profile 1 (2D)
SeqControl(8).command = 'setTPCProfile';
SeqControl(8).condition = 'next';
SeqControl(8).argument = 1;
% - time between Doppler and 2D (15000 usec = 15msec)
% time from the first 2D pulse to the last Doppler pulse
timeElapse = 2*(P.numRays)*(TTNW+TTNI)+SeqControl(4).argument + P.dopPRIs*SeqControl(7).argument;
SeqControl(9).command = 'timeToNextAcq';
SeqControl(9).argument = 39000; %timeElapse + 2500;
FrameRate2D = 1e6/SeqControl(9).argument;
% - return to Matlab
SeqControl(10).command = 'returnToMatlab';
% - Jump back to start.
SeqControl(11).command = 'jump';
SeqControl(11).argument = 2;
nsc = 12; % nsc is count of SeqControl objects

% Specify Event structure arrays.
n = 1;

Event(n).info = 'ROI plot';
Event(n).tx = 0;         % no transmit
Event(n).rcv = 0;        % no rcv
Event(n).recon = 0;      % no reconstruction
Event(n).process = 4;    % outline doppler region
Event(n).seqControl = 0; % no seqCntrl
n = n+1;
    
for i = 1:Resource.RcvBuffer(1).numFrames
    rcvInd = (2*P.numRays + m*P.dopPRIs)*(i-1); % rcvInd keeps track of Receive index increment per frame.

    Event(n).info = 'Set RcvProfile for 2D';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 0;      % no reconstruction
    Event(n).process = 0;    % no processing
    Event(n).seqControl = 1;
    n = n+1;    
    
    % Acquire 2D frame
    for j = 1:P.numRays             
        Event(n).info = 'acquire 1st interleave pair';
        Event(n).tx = j;   % use 1st set of TX structures.
        Event(n).rcv = rcvInd+2*j-1;   
        Event(n).recon = 0;      % no reconstruction.
        Event(n).process = 0;    % no processing
        Event(n).seqControl = 2; % TTNI     
        n = n+1;
        Event(n).info = 'acquire 2nd interleave pair';
        Event(n).tx = P.numRays+j; % use 2nd set of TX structures.
        Event(n).rcv = rcvInd+2*j;   
        Event(n).recon = 0;      % no reconstruction.
        Event(n).process = 0;    % no processing
        Event(n).seqControl = 3;      
        n = n+1;        
    end
    Event(n-1).seqControl = [4,5];   % replace last 2D acquisition Event's seqControl
    
%     Event(n).info = 'Set RcvProfile for Doppler';
%     Event(n).tx = 0;         % no transmit
%     Event(n).rcv = 0;        % no rcv
%     Event(n).recon = 0;      % no reconstruction
%     Event(n).process = 0;    % no processing
%     Event(n).seqControl = 6;
%     n = n+1;  
    
    % Acquire ensembles for all Doppler beams.
    for j = 1:P.dopPRIs
        for k = 1:m
            Event(n).info = 'Acquire Doppler ensembles';
            Event(n).tx = 2*P.numRays+k;      % use next TX structure after 2D.
            Event(n).rcv = rcvInd+(2*P.numRays+(j-1)*m+k);
            Event(n).recon = 0;      % no reconstruction.
            Event(n).process = 0;    % no processing
            Event(n).seqControl = 7; 
            n = n+1;
        end
    end
    Event(n-1).seqControl = [8,9,nsc]; % replace last Doppler acquisition Event's seqControl
      SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
      nsc = nsc+1;

    Event(n).info = 'recons and 2D process'; 
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = [1,2];  % reconstruction for both 2D and Doppler
    Event(n).process = 1;    % process 2D
    Event(n).seqControl = 0; % no seqCntrl
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
        Event(n).seqControl = 10;        
    end
    n = n+1;   
    
end

Event(n).info = 'Jump back';
Event(n).tx = 0;        % no TX
Event(n).rcv = 0;       % no Rcv
Event(n).recon = 0;     % no Recon
Event(n).process = 0; 
Event(n).seqControl = 11;


%% User specified UI Control Elements
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
PRFmin = 1000; PRFmax = 5000; stepDiff = PRFmax-PRFmin;
UI(7).Control = VsSliderControl('LocationCode','UserC3','Label','PRF (Hz)',...
                  'SliderMinMaxVal',[PRFmin,PRFmax,P.dopPRF],...
                  'SliderStep',[50/stepDiff,500/stepDiff],'ValueFormat','%3.0f', ...
                  'Callback', @PRFChange);

% - Steer Angle adjustment
UI(8).Control = VsSliderControl('LocationCode','UserC2','Label','Doppler Angle',...
                  'SliderMinMaxVal',[-20,20,round(P.dopAngle*180/pi)],...
                  'SliderStep',[1/40,5/40],'ValueFormat','%3.0f', ...
                  'Callback', @SteerAngle);

% - Range Change
MinMaxmm = [2,16]; % min max in mm
AxesUnit = 'wls';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
        AxesUnit = 'mm';
        MinMaxVal = [MinMaxmm, P.endDepth * (Resource.Parameters.speedOfSound/1000/Trans.frequency)];
    else
        MinMaxVal = [MinMaxmm * (Trans.frequency/Resource.Parameters.speedOfSound*1000), P.endDepth];
    end
end
stepBase = MinMaxVal(2)-MinMaxVal(1);
UI(9).Control = VsSliderControl('LocationCode','UserA1','Label',['Range (',AxesUnit,')'],...
                  'SliderMinMaxVal',MinMaxVal,...
                  'SliderStep',[1/stepBase,2/stepBase],'ValueFormat','%3.1f', ...
                  'Callback', @RangeChangeCallback);

% - Transmit frequency
UI(10).Control = VsSliderControl('LocationCode','UserC1','Label','Dop Transmit Freq',...
                  'SliderMinMaxVal',[15,40,26.7858],...
                  'SliderStep',[0.025,0.1],'ValueFormat','%2.3f', ...
                  'Callback', @TransmitFreq);
              
% - Power compression factor
UI(11).Control = VsSliderControl('LocationCode','UserB5','Label','pwrCompression',...
                  'SliderMinMaxVal',[0.1,1.0,0.5],...
                  'SliderStep',[1/100,0.1],'ValueFormat','%1.2f', ...
                  'Callback', @PwrCompression);
             
% - External function for ROIplot
EF(1).Function = vsv.seq.function.ExFunctionDef('ROIplot', @ROIplot);

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 3;    

% Save all the structures to a .mat file.
filename = fullfile( vsv.file.getVSXDir, 'MatFiles/SetUpL39_21gDWideBeamDoppler_QSApp.mat');

save(filename);
% filename = 'L39-21gDWideBeamDoppler'; %VSX


% **** Callback functions used by UIs. ****
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

function DopplerModeCallback(~, ~, UIState)
    P = evalin('base','P');
    Control = repmat(struct('Command','set&Run','Parameters',[]),1,4);
    Process = evalin('base','Process');
    Resource = evalin('base','Resource');
    hDisplay = Resource.DisplayWindow(1).figureHandle;
    currentMap = get(hDisplay,'Colormap');

    switch UIState
       case 1  % Velocity mode
          newMap = grayscaleCFImap;
          newMap(1:128,:) = currentMap(1:128,:);      
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
          newMap(1:128,:) = currentMap(1:128,:);      
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

function DopPowerThreshold(~, ~, UIValue)
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

function ColorPriorityLevel(~, ~, UIValue)
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

function ColorPersistenceLevel(~, ~, UIValue)
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

function ReplotROI(~, ~, UIState)
    Process = evalin('base','Process');
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
    return
end

function PRFChange(~, ~, UIValue)
    P = evalin('base','P');
    TTNI = evalin('base','TTNI');
    Process = evalin('base','Process');
    SeqControl = evalin('base','SeqControl');

    P.dopPRF = UIValue;
    m = P.dopNumRays;
    SeqControl(6).argument = round(1e6/(P.dopPRF*m)-TTNI); 

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

function SteerAngle(~, ~, UIValue)
    P = evalin('base','P');
    TX = evalin('base','TX');
    PData = evalin('base','PData');

    P.dopAngle = UIValue * pi/180;
    assignin('base','P',P);

    for n = 1:P.dopNumRays
        PData(2).Region(n).Shape.angle = P.dopAngle;
    end
    PData(3).Region.Shape.angle = P.dopAngle;
    PData(2).Region = computeRegions(PData(2));
    PData(3).Region = computeRegions(PData(3));
    assignin('base','PData',PData);

    for n = 1:P.dopNumRays    
        TX(2*P.numRays+n).Steer = [P.dopAngle,0.0];
        TX(2*P.numRays+n).Delay = computeTXDelays(TX(2*P.numRays+n));
        TX(2*P.numRays+n).TXPD = computeTXPD(TX(2*P.numRays+n),PData(2));
    end
    assignin('base','TX',TX);

    Control.Command = 'update&Run';
    Control.Parameters = {'PData','InterBuffer','ImageBuffer','Recon','TX'};
    assignin('base','Control',Control);
end

function RangeChangeCallback(hObject, ~, UIValue)
    simMode = evalin('base','Resource.Parameters.simulateMode');
    % No range change if in simulate mode 2.
    if simMode == 2
        set(hObject,'Value',evalin('base','P.endDepth'));
        return
    end

    % P.endDepth
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
    P.dopEndDepth = P.endDepth;
    P.txFocus = 5*P.dopEndDepth;

    % PData
    PData = evalin('base','PData');
    m = P.dopNumRays;

    PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3)); % rows
    for i = 1:size(PData(1).Region,2)
        PData(1).Region(i).Shape.height = P.endDepth;
    end
    PData(2) = PData(1); % default is the same as PData(1), but the user may want to change the PDelta
    % - Doppler PData structure for TXPD calculation
    PData(2) = PData(1);
    PData(2).PDelta(1) = Trans.spacing/2;
    PData(2).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(2).PDelta(1)); % cols
    PData(2).Region = repmat(struct('Shape',struct(...
        'Name','Parallelogram',...
        'Position',[0,0,P.dopStartDepth],... %wil be cahnged later
        'width', 22*Trans.spacing,...
        'height', P.dopEndDepth-P.dopStartDepth,...
        'angle',P.dopAngle)),1,m);

    TxDopOrgX = evalin('base','TxDopOrgX');
    for n = 1:m    
        PData(2).Region(n).Shape.Position(1) = TxDopOrgX(n);
    end

    % - PData(3) is used for doppler processing and outline the doppler area
    PData(3) = PData(2); PData(3).Region = [];
    PData(3).Region.Shape = struct(...
        'Name','Parallelogram',...
        'Position',[TxDopOrgX(1)-PData(2).Region(1).Shape.width/2+44*Trans.spacing,0,P.dopStartDepth],...
        'width', 88*Trans.spacing,...
        'height', P.dopEndDepth-P.dopStartDepth,...
        'angle',P.dopAngle);
    for n = 1:3
        PData(n).Region = computeRegions(PData(n));
    end
    assignin('base','PData',PData);
    evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);');

    % Update TXPD data of TX structures.
    TX = evalin('base','TX');
    h = waitbar(0,'Program TX parameters, please wait!');
    for n = 1:P.numRays
        TX(n).focus = P.txFocus;
        TX(n+P.numRays).focus = P.txFocus;
        TX(n+P.numRays).Delay = computeTXDelays(TX(n+P.numRays));
        TX(n).Delay = TX(n+P.numRays).Delay + 0.5*Trans.frequency/62.5; % 0.24 wls of Trans.frequency=30 equivalent to 8 nsec.    
        TX(n).TXPD = computeTXPD(TX(n),PData(1));
        TX(n+P.numRays).TXPD = computeTXPD(TX(n+P.numRays),PData(1));    
        waitbar(n/(P.numRays+m))
    end
    k = 2*P.numRays;
    for n = 1:m
        TX(n+k).TXPD = computeTXPD(TX(n+k),PData(2));
        waitbar((P.numRays+n)/(P.numRays+m))
    end
    close(h)
    assignin('base','TX',TX);

    % Receive
    Receive = evalin('base', 'Receive');
    maxAcqLngth2D = sqrt(P.endDepth^2 + (Trans.numelements*Trans.spacing)^2);
    maxAcqLngthDop =  sqrt((P.dopEndDepth/cos(P.dopAngle))^2 + (22*Trans.spacing)^2);
    for i = 1:Resource.RcvBuffer(1).numFrames
        k = (2*P.numRays + m*P.dopPRIs)*(i-1); % k keeps track of Receive index increment per frame.
        for j = 1:2*P.numRays  % acquisitions for 2D        
            Receive(j+k).endDepth = maxAcqLngth2D; 
        end
        % Doppler acquisitions
        n = 2*P.numRays;
        for j = 1:P.dopPRIs*m
            Receive(k+n+j).endDepth = maxAcqLngthDop;        
        end
    end
    assignin('base','Receive',Receive);

    % TGC
    TGC = evalin('base','TGC');
    TGC(1).rangeMax = P.endDepth;
    TGC(1).Waveform = computeTGCWaveform(TGC(1));
    TGC(2).rangeMax = P.dopEndDepth;
    TGC(2).Waveform = computeTGCWaveform(TGC(2));
    assignin('base','TGC',TGC);
    assignin('base','nTGC',get(findobj('tag','TGCnum'),'Value'));
    evalin('base','if VDAS==1, Result = loadTgcWaveform(nTGC); end');

    Control(1).Command = 'set&Run';
    Control(1).Parameters = {'Process',2,'prf',P.dopPRF};
    Control(2).Command = 'update&Run';
    Control(2).Parameters = {'PData','InterBuffer','ImageBuffer','DisplayWindow','SeqControl','TX','Receive','Recon','TGC'};
    assignin('base','Control', Control);
    assignin('base', 'action', 'displayChange');
end

function TransmitFreq(~, ~, UIValue)
    P = evalin('base','P');
    P.numRays = evalin('base','P.numRays');
    P.dopPRIs = evalin('base','P.dopPRIs');
    TW = evalin('base','TW');
    Receive = evalin('base','Receive');
    Resource = evalin('base','Resource');

    m = P.dopNumRays;
    TW(2).Parameters(1) = UIValue;

    for i = 1:Resource.RcvBuffer(1).numFrames
        k = 2*(P.numRays + m*P.dopPRIs)*(i-1); % k keeps track of Receive index increment per frame.
        for j = (2*P.numRays+1):2*(P.numRays+P.dopPRIs*m)
                % Doppler acquisition
                Receive(j+k).demodFrequency = TW(2).Parameters(1);
        end
    end

    assignin('base','TW',TW);
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'TW','Receive'};
end

function PwrCompression(~, ~, UIValue)
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
end


%% **** Callback function used by External function definition (EF) ****
function ROIplot(varargin)
% ROIplot External Processing function

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

    return
end
