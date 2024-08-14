%% Copyright 2021 Verasonics, Inc. Verasonics Registered U.S. Patent and Trademark Office.
% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: SetUpVermon1024_DivergingSpiralDoppler_2system.m

% Description:  - Doppler sequence based off of the paper:
% "3-D Ultrafast Doppler Imaging Applied to the Noninvasive and
% Quantitative Imaging of Blood Vessels in Vivo"
% by J. Provost et al
% EEE Trans Ultrason Ferroelectr Freq Control. 2015 August
% 62(8): 1467–1472. doi:10.1109/TUFFC.2015.007032
%
%  Angles: -2, 0, 2 degrees in both directions (9 total)
%  compoounding all 9 angles ==> 1 volume
%  PRF = 10304 Hz;
%  2000 volumes acquired
%

clear all

Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.
filename = 'Vermon512_DivergingSpiralDoppler_IQ_2system';

%=== Settings ===%
XYdepth = 42;            % depth of C plane imaging plane
% Bmode
P(1).numTx = 1;
P(1).startDepth = 4;   % Acquisition depth in wavelengths
P(1).endDepth = 64;   % This should preferrably be a multipled of 128 samples.
P(1).viewAngle = 20 * pi/180;                    % angle between z-axis and surface of cone.
P(1).curved = 0; % set to 0 for flat 1 for defocused at the reconstruction apex
P(1).numFrames = 10;
% Doppler
P(2).numTx = 1;
P(2).startDepth = 4;   % Acquisition depth in wavelengths
P(2).endDepth = 64;   % This should preferrably be a multipled of 128 samples.
P(2).viewAngle = 20 * pi/180;                    % angle between z-axis and surface of cone.
P(2).curved = 0; % set to 0 for flat 1 for defocused at the reconstruction apex
P(2).numFrames = 10;

% Set Bmode parameters
% - Bmode properties - %
reject = 0;
pgain = 20;              % digital gain
persist = 30;            % persistance level
senscutoff = 0.6;        % sensitivity cutoff
CompressionFactor = 75;  % compression Factor
na = 10;      % Set na = number of angles that would be use to rotate the imaging plane around 4*pi.

% - Set Doppler parameters - %
nDopAngles = 9; % (+-2 in both directions + 0)
nVol = 32;%128;       % Set nVol = number of volume acquisitions in Doppler ensemble. 200x9
dopAngle = 2 * pi/180;
dopPRF = 10304; % Doppler PRF in Hz.
acqPeriod = nDopAngles*nVol/dopPRF; %in seconds
acqFPD = 1/acqPeriod;
pwrThres = 0.02;
cpt = 255;
persf = 80;
dopFilter = 'FIRLow';


% colorflow Settings
DopState='computeCFIFreqEst';
dopData = 'signedColor3D';
dopColorMap = grayscaleCFImap;
% Power Doppler Setting
%DopState= 'computeCFIPowerEst';
%dopData = 'unsignedColor3D';
%dopColorMap = grayscaleCPAmap;

RcvProfile(1).AntiAliasCutoff = 10;
RcvProfile(1).LnaGain = 18; % Profile used for imaging
RcvProfile(2).LnaGain = 24; % Profile used for Doppler
RcvProfile(2).LnaZinSel = 31; % Force high-Z state for best Doppler sensitivity
frameRateFactor = 5;% Specify factor for converting sequenceRate to frameRate.

%%===================== RDMA SETUP ======================================%%
usingMultiSys = 1;
[RDMAconfig, RDMAsystem]=vsv.multi.generateConfig(2);
%%=======================================================================%%

% Define system parameters.
Resource.Parameters.numTransmit = 256;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 256;   % number of receive channels.
Resource.Parameters.speedCorrectionFactor = 1.0;
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 0;
Resource.Parameters.initializeOnly = 0;

%% Specify Trans structure array.
Trans.name = 'Matrix1024-3';
Trans.units = 'wavelengths';  % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);  % computeTrans is used for known transducers.
Trans.id = hex2dec('0D0100'); % need this still
Trans.numelements = 512;
Trans.ElementPos = Trans.ElementPos(1:Trans.numelements,:);
Trans.ConnectorES=zeros(1,Trans.numelements)';
Trans.ConnectorES(RDMAconfig.nodeIndex*256+[1:256])=[128:-1:65, 1:1:64, 192:-1:129, 193:1:256];

%Intermediate Variables
waveLength = (Resource.Parameters.speedOfSound/1000)/Trans.frequency;
if strcmp(Trans.units,'mm')
    Trans.ElementPosMm = Trans.ElementPos;
    Trans.ElementPosWL = Trans.ElementPos./waveLength;
else
    Trans.ElementPosMm = Trans.ElementPos.*waveLength;
    Trans.ElementPosWL = Trans.ElementPos;
end

%% PData
PData(1).PDelta = [1,1,1]*1;
if ~exist('extent','var'), extent = max(max(Trans.ElementPosWL(:,1)),max(Trans.ElementPosWL(:,2))); end
zApex = -extent/tan(P(1).viewAngle); % in WL

PData(1).Size(1) = ceil(2.0*(P(1).endDepth-zApex)*tan(P(1).viewAngle)/PData(1).PDelta(2));  if mod(PData(1).Size(1),2)==0, PData(1).Size(1) =  PData(1).Size(1)+1; end %rows
PData(1).Size(2) = PData(1).Size(1);%columns
PData(1).Size(3) = ceil((P(1).endDepth)/PData(1).PDelta(3)); %sections
PData(1).Origin = [-((PData(1).Size(2)-1)/2)*PData(1).PDelta(1), ((PData(1).Size(1)-1)/2)*PData(1).PDelta(1), 0];

PData(1).Region(1) = struct(...
    'Shape',struct('Name','Pyramid',...
    'Position',[0,0,zApex],...
    'angle',P(1).viewAngle,...
    'z1',P(1).startDepth,...
    'z2',P(1).endDepth));
PData(1).Region(2) = struct(...
    'Shape',struct('Name','Slice',...
    'Orientation','xz',...
    'oPAIntersect', 0, 'andWithPrev',1));
PData(1).Region(3) = struct(...
    'Shape',struct('Name','Slice',...
    'Orientation','yz',...
    'oPAIntersect', 0, 'andWithPrev',1));
PData(1).Region(4) = struct(...
    'Shape',struct('Name','Slice',...
    'Orientation','xy',...
    'oPAIntersect', XYdepth, 'andWithPrev',1));
PData(1).Region = computeRegions(PData(1));

PData.Region(5).PixelsLA = unique([PData(1).Region(2).PixelsLA;PData(1).Region(3).PixelsLA;PData(1).Region(4).PixelsLA]);
PData.Region(5).Shape.Name = 'Custom';
PData.Region(5).numPixels = length(PData(1).Region(5).PixelsLA);
PData(2)=PData(1);

%% Specify Media.  Use point targets in middle of PData.
for i=1:61
    Media.MP(i,:)= [0,-30+i,50,0.1];
    Media.MP(61+i,:)= [-30+i,0,72,0.1];
end
Media.MP(size(Media.MP,1)+1,:)=[18,32,82,1];
Media.numPoints = size(Media.MP,1);

%% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 697600;
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = P.numFrames;  % buffer for primary IQ
Resource.InterBuffer(1).numFrames = 1;
Resource.InterBuffer(2).numFrames = 1;
Resource.InterBuffer(2).pagesPerFrame = nVol;
if strcmp(RDMAconfig.rdmaRole,'primary')
    %--- Interbuffer---%
    % Bmode
    Resource.InterBuffer(3).numFrames = 1;  %
    % Doppler
    Resource.InterBuffer(4).numFrames = 1;  %
    Resource.InterBuffer(4).pagesPerFrame = nVol;
    
    %---Destintation InterBuffer for Just Doppler Processing
    Resource.InterBuffer(5).numFrames = 1;  %
    Resource.InterBuffer(5).pagesPerFrame = nVol;

    %--- ImageBuffer ---%
    Resource.ImageBuffer(1).numFrames = P(1).numFrames; %Bmode
    Resource.ImageBuffer(2).numFrames = P(2).numFrames; %Doppler
    %-- Display Windows --%
    Resource.DisplayWindow(1).Type = 'Verasonics';
    Resource.DisplayWindow(1).Title = '3D Flash Image - XZ plane';
    Resource.DisplayWindow(1).pdelta = 0.25;
    Resource.DisplayWindow(1).Position = [0,580, ...
        ceil(PData(1).Size(2)*PData(1).PDelta(2)/Resource.DisplayWindow(1).pdelta), ... % width
        ceil(PData(1).Size(3)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta)];    % height
    Resource.DisplayWindow(1).Orientation = 'xz';
    Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0.0,0.0];
    Resource.DisplayWindow(1).Colormap = dopColorMap;
    Resource.DisplayWindow(1).splitPalette = 1;
    Resource.DisplayWindow(1).AxesUnits = 'mm';
    %-------------------------------------------------%
    Resource.DisplayWindow(2).Type = 'Verasonics';
    Resource.DisplayWindow(2).Title = '3D Flash Image - YZ plane';
    Resource.DisplayWindow(2).pdelta = Resource.DisplayWindow(1).pdelta;
    Resource.DisplayWindow(2).Position = [430,580, ...
        ceil(PData(1).Size(1)*PData(1).PDelta(1)/Resource.DisplayWindow(2).pdelta), ... % width
        ceil(PData(1).Size(3)*PData(1).PDelta(3)/Resource.DisplayWindow(2).pdelta)];    % height
    Resource.DisplayWindow(2).Orientation = 'yz';
    Resource.DisplayWindow(2).ReferencePt = [0,-PData(1).Origin(2),0];
    Resource.DisplayWindow(2).Colormap = dopColorMap;
    Resource.DisplayWindow(2).splitPalette = 1;
    Resource.DisplayWindow(2).AxesUnits = 'mm';
    %-------------------------------------------------%
    Resource.DisplayWindow(3).Type = 'Verasonics';
    Resource.DisplayWindow(3).Title = '3D Flash Image - XY plane';
    Resource.DisplayWindow(3).pdelta = Resource.DisplayWindow(1).pdelta;
    Resource.DisplayWindow(3).Position = [0,40, ...
        ceil(PData(1).Size(2)*PData(1).PDelta(2)/Resource.DisplayWindow(3).pdelta), ... % width
        ceil(PData(1).Size(1)*PData(1).PDelta(1)/Resource.DisplayWindow(3).pdelta)];    % height
    Resource.DisplayWindow(3).Orientation = 'xy';
    Resource.DisplayWindow(3).ReferencePt = [PData(1).Origin(1),...
        -PData(1).Origin(2), XYdepth];%PData(1).Region(end).Shape.oPAIntersect];
    Resource.DisplayWindow(3).Colormap = dopColorMap;
    Resource.DisplayWindow(3).splitPalette = 1;
    Resource.DisplayWindow(3).AxesUnits = 'mm';
end

% Specify Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,.67,2,1];  % A, B, C, D

TW(2).type = 'parametric';
TW(2).Parameters = [Trans.frequency,.67,8,1];  % A, B, C, D

% Specify TPC structure.
TPC(1).hv = 30;

% Specify TX structure array.
rotang=linspace(0,4*pi,na); % na=9 creates a plane every 45 degrees (0, 45, 90, 135, 180, 225, 270, 315 and 360).

r=min(abs([Trans.ElementPosWL(1,1) Trans.ElementPosWL(1,2)])); % In WL

R=linspace(r,3*r/4,round(na/2));
R=[R(1:end-1) linspace(3*r/4,r/4,na-round(na/2)) 0];

AZ= [atan(R.*cos(rotang)/zApex) ]; % This add an extra value for the for AZ and EL to make a non steer diverging wave
EL= [atan(R.*sin(rotang)/zApex) ];

TX = repmat(struct('waveform', 1, ...
    'Origin', [0.0,0.0,0.0], ...
    'focus', P(2).curved*zApex, ...
    'Steer', [0.0,0.0], ...
    'Apod', ones(1, Trans.numelements), ...
    'Delay', zeros(1, Trans.numelements), ...
    'TXPD', [], ...
    'peakCutOff', 2,...
    'peakBLMax', 20), 1, na+nDopAngles+1); % na TXs for 2D + 1 for Doppler + 1 dummy


% - Set event specific TX attributes.
for n = 1:na   % na transmit events for 2D
    TX(n).Steer = [AZ(n), EL(n)];
    TX(n).Delay = computeTXDelays(TX(n));
end
% -- 4 TX struct needed for Doppler
%DopAngles = [0 0; -2 0; 2 0; 0 -2; 0 2; -2 -2; 2 2; -2 2; 2 -2];
DopAngles = [AZ;EL]';
for n = 1:nDopAngles
    TX(na+n).waveform = 2;
    TX(na+n).focus = 0;
    TX(na+n).Steer = [DopAngles(n,1),DopAngles(n,2)];
    TX(na+n).Delay = computeTXDelays(TX(na+n));
    %TX(na+n).peakCutOff= 2;
    %TX(na+n).peakBLMax = 10;
    %TX(na+n).TXPD = computeTXPD(TX(na+1),PData(2));
end
%- Dummy Transmit -%
TX(na+nDopAngles+1)=TX(na);
TX(na+nDopAngles+1).Apod(:)=0;

%% Specify TGC Waveform structure.

TGC(1).CntrlPts = [0 256 358 512 665 767 870 900];
TGC(1).rangeMax = 128;
TGC(1).Waveform = computeTGCWaveform(TGC(1));

%TGC(2).CntrlPts =[ 0 272 662 662 662 662 662 662];
TGC(2).CntrlPts = [600 700 800 850 900 900 900 900];
TGC(2).rangeMax = 128;
TGC(2).Waveform = computeTGCWaveform(TGC(1));

temp(1) = (P(1).endDepth-zApex)*tan(P(1).viewAngle)+ max([max(Trans.ElementPosWL(:,1)), max(Trans.ElementPosWL(:,2))]);
temp(2) = (P(2).endDepth-zApex)*tan(P(2).viewAngle)+ max([max(Trans.ElementPosWL(:,1)), max(Trans.ElementPosWL(:,2))]);
maxAcqLength2D = sqrt(P(1).endDepth^2 + temp(1)^2);
maxAcqLengthDop = sqrt(P(2).endDepth^2 + temp(2)^2);
wl4sPer128 = 128/(4*2);  % wavelengths in a 128 sample block for 4 smpls per wave round trip.
wl2sPer128 = 128/(2*2);  % wavelengths in a 128 sample block for 2 smpls per wave round trip.

%Receive

Receive = repmat(struct(...
    'Apod', zeros(1, Trans.numelements), ...
    'startDepth', 0, ...
    'endDepth', wl4sPer128*ceil(maxAcqLength2D/(wl4sPer128)), ...
    'TGC', 1, ...
    'mode', 0, ...
    'bufnum', 1, ...
    'framenum', 1, ...
    'acqNum', 1, ...
    'sampleMode', 'NS200BW', ...
    'demodFrequency',TW(1).Parameters(1),...
    'callMediaFunc', 0), 1, (nVol*nDopAngles+na)*P(1).numFrames);

%Set event specific Receive attributes.
for i = 1:P(1).numFrames
    k = (nVol*nDopAngles+na)*(i-1); % k keeps track of Receive index increment per frame.
    % - Set attributes for each frame.
    Receive(k+1).callMediaFunc = 1;
    for j = 1:na
        Receive(k+j).Apod(RDMAconfig.nodeIndex*256+[1:256]) = 1;
        Receive(k+j).framenum = i;
        Receive(k+j).acqNum = j;
    end
    for j = (na+1):(na+nVol*nDopAngles)
        % Doppler acquisition
        Receive(k+j).Apod(RDMAconfig.nodeIndex*256+[1:256])=1;
        Receive(k+j).startDepth = P(2).startDepth;
        Receive(k+j).endDepth = P(2).startDepth + wl2sPer128*ceil(maxAcqLengthDop/wl2sPer128);
        Receive(k+j).sampleMode = 'BS100BW';
        Receive(k+j).demodFrequency = TW(2).Parameters(1);
        Receive(k+j).TGC = 2;
        Receive(k+j).framenum = i;
        Receive(k+j).acqNum = j;        % Doppler acqNums continue after 2D
        Receive(k+j).InputFilter =...
            [+0.00085 +0.00000 -0.00345 +0.00000 +0.00839 +0.00000 -0.01669 ...
            +0.00000 +0.02881 +0.00000 -0.04434 +0.00000 +0.06207 +0.00000 ...
            -0.07993 +0.00000 +0.09531 +0.00000 -0.10574 +0.00000 +0.10944];
    end
end

%% Reconstruction
%Recon 1 is preview Recon 1 flash image from ensemble
Recon(1) = struct('senscutoff', senscutoff, ...
    'pdatanum', 1, ...
    'rcvBufFrame', -1, ...
    'IntBufDest', [1,1], ...
    'RINums', 1:na);

if strcmp(RDMAconfig.rdmaRole,'primary')% primary uses buffer #5
    Recon(2) = struct('senscutoff', senscutoff, ...
        'pdatanum', 2, ...
        'rcvBufFrame', -1, ...
        'IntBufDest', [5,1],...
        'RINums', na+[1:nVol*nDopAngles]);
else %each secondary just uses the second
    Recon(2) = struct('senscutoff', senscutoff, ...
        'pdatanum', 2, ...
        'rcvBufFrame', -1, ...
        'IntBufDest', [2,1],...
        'RINums', na + [1:nVol*nDopAngles]);
end

% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 'accumIQ', ...
    'Pre',[],...
    'Post', [],...
    'txnum', 1, ...
    'rcvnum', 1, ...
    'pagenum',1, ...
    'scaleFactor', 1, ...
    'regionnum', 5), 1, na+nVol*nDopAngles);
ReconInfo(1).Pre = 'clearInterBuf';
for i = 1:na %go through each bmode angle
    ReconInfo(i).txnum = i;
    ReconInfo(i).rcvnum = i;
    ReconInfo(i).regionnum = 5;
end
for i = 1:nVol %go through ensemble
    for j = 1:nDopAngles
        if j == 1
            ReconInfo(na+(i-1)*nDopAngles+j).mode = 'replaceIQ';
        else
            ReconInfo(na+(i-1)*nDopAngles+j).mode = 'accumIQ';
        end
        
        ReconInfo(na+(i-1)*nDopAngles+j).txnum = na + j;
        ReconInfo(na+(i-1)*nDopAngles+j).rcvnum = na + (i-1)*nDopAngles + j;
        ReconInfo(na+(i-1)*nDopAngles+j).pagenum = i;
        ReconInfo(na+(i-1)*nDopAngles+j).regionnum = 5;
    end
end


%% Specify Process structure array.
if strcmp(RDMAconfig.rdmaRole,'primary')
    Process(1).classname = 'Image';
    Process(1).method = 'imageDisplay';
    Process(1).Parameters = {'imgbufnum',1,...
        'framenum',-1,...
        'pdatanum',1,...
        'srcData','intensity3D',...
        'pgain', pgain,...
        'reject', reject,...
        'persistMethod','simple',...
        'persistLevel',persist,...
        'interpMethod','4pt',...
        'grainRemoval','none',...
        'processMethod','none',...
        'averageMethod','none',...
        'compressMethod','power',...
        'compressFactor', CompressionFactor,...
        'mappingMethod','lowerHalf',...
        'display',0,...      % don't display image after processing
        'displayWindow',1};
    
    Process(2).classname = 'Image';
    Process(2).method = 'imageDisplay';
    Process(2).Parameters = {'imgbufnum',1,...   % number of buffer to process.
        'framenum',-1,...    % (-1 => lastFrame)
        'pdatanum',1,...     % number of PData structure to use
        'srcData','intensity3D',...
        'pgain',pgain,...    % pgain is image processing gain
        'reject', reject,... % reject level
        'persistMethod','simple',...
        'persistLevel',persist,...
        'interpMethod','4pt',...
        'grainRemoval','none',...
        'processMethod','none',...
        'averageMethod','none',...
        'compressMethod','power',...
        'compressFactor', CompressionFactor,...
        'mappingMethod','lowerHalf',...
        'display',0,...      % don't display image after processing
        'displayWindow',2};
    Process(3).classname = 'Image';
    Process(3).method = 'imageDisplay';
    Process(3).Parameters = {'imgbufnum',1,...   % number of buffer to process.
        'framenum',-1,...    % (-1 => lastFrame)
        'pdatanum',1,...     % number of PData structure to use
        'srcData','intensity3D',...
        'pgain',pgain,...    % pgain is image processing gain
        'reject', reject,... % reject level
        'persistMethod','simple',...
        'persistLevel',persist,...
        'interpMethod','4pt',...
        'grainRemoval','none',...
        'processMethod','none',...
        'averageMethod','none',...
        'compressMethod','power',...
        'compressFactor', CompressionFactor,...
        'mappingMethod','lowerHalf',...
        'display',0,...      % don't display image after processing
        'displayWindow',3};
    
    %---------------------------------------------------------------------%
    Process(4).classname = 'Process';
    Process(4).method = 'sumIQNormMag';  %combine buffers and put it into image buffer
    Process(4).Parameters = {'IntBufSrc',[1,1,3,1],... %Bmode Primary, Bmode Secondary
        'ImgBufDest',[1,-1]};
    
    Process(5).classname = 'Process';
    Process(5).method = 'sumIQNorm';  %combine buffers and put it into inter buffer
    Process(5).Parameters = {'IntBufSrc',[2,1,4,1],... %Doppler Primary, Doppler Secondary
        'IntBufDest',[5,1]};
    %---------------------------------------------------------------------%
    Process(6).classname = 'Doppler'; % process structure for 1st Doppler ensemble
    Process(6).method = DopState;
    Process(6).Parameters = {'IntBufSrc', [5, 1],...      % buffer and frame num of interbuffer
        'SrcPages', [3, nVol-2],...   % Start and last pagenum
        'ImgBufDest', [2, 1],...   % Buffer and frame num of image desitation
        'pdatanum', 2,...           % Number of PData structure
        'prf', dopPRF,...           % Doppler PRF in Hz
        'wallFilter', dopFilter,...  % 1 -> quadratic regression
        'pwrThreshold', pwrThres,...
        'maxPower', 50,...
        'postFilter', 1};
    %---------------------------------------------------------------------%
    Process(7).classname = 'Image';
    Process(7).method = 'imageDisplay';
    Process(7).Parameters = {'imgbufnum',2,...   % number of buffer to process.
        'framenum',1,...   % (-1 => lastFrame)
        'srcData',dopData,... % type of data to display.
        'pdatanum',2,...    % number of PData structure to use
        'pgain',1.0,...     % pgain is image processing gain
        'reject',2,...      % reject level
        'persistMethod','none',...
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
    
    Process(8).classname = 'Image';
    Process(8).method = 'imageDisplay';
    Process(8).Parameters = {'imgbufnum',2,...   % number of buffer to process.
        'framenum',1,...   % (-1 => lastFrame)
        'srcData',dopData,... % type of data to display.
        'pdatanum',2,...    % number of PData structure to use
        'pgain',1.0,...     % pgain is image processing gain
        'reject',2,...      % reject level
        'persistMethod','none',...
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
        'displayWindow',2};
    
    Process(9).classname = 'Image';
    Process(9).method = 'imageDisplay';
    Process(9).Parameters = {'imgbufnum',2,...   % number of buffer to process.
        'framenum',1,...   % (-1 => lastFrame)
        'srcData',dopData,... % type of data to display.
        'pdatanum',2,...    % number of PData structure to use
        'pgain',1.0,...     % pgain is image processing gain
        'reject',2,...      % reject level
        'persistMethod','none',...
        'persistLevel',persf,...
        'interpMethod','4pt',...
        'grainRemoval','none',...
        'processMethod','none',...
        'averageMethod','none',...
        'compressMethod','power',...
        'compressFactor', 40,...
        'mappingMethod','upperHalf',...
        'threshold', cpt,...
        'display', 1,...      % don't display image after processing
        'displayWindow', 3};
end

%% Specify SeqControl structure arrays.
% Jump to start of sequence
JUMP = 1;
SeqControl(JUMP).command = 'jump'; % jump back to start
SeqControl(JUMP).argument = 1;
% Time to next acquisition
TTNAFLASH = 2;
SeqControl(TTNAFLASH).command = 'timeToNextAcq';  % time between bmode pulses
SeqControl(TTNAFLASH).argument = 250;  % 250 usec
% Time to next acquisition
TTNAPRF = 3;
SeqControl(TTNAPRF).command = 'timeToNextAcq';  % time between doppler pulses
SeqControl(TTNAPRF).argument = round((1/dopPRF)/1e-6);
%
TTNAFRAME = 4;
SeqControl(TTNAFRAME).command = 'timeToNextAcq';  % time between frames
SeqControl(TTNAFRAME).argument = 10;%100000 - SeqControl(TTNAFLASH).argument*na - SeqControl(TTNAPRF).argument*(nVol*nDopAngles-1);  % 20 msec
% Return to Matlab
RETMAT = 5;
SeqControl(RETMAT).command = 'returnToMatlab';
% Multisystem Sync seqcontrol
MSSYNC = 6;
SeqControl(MSSYNC).command = 'multiSysSync';
SeqControl(MSSYNC).condition = 'normal';       % 'normal' uses internal trigger; other options are 'BNC_Rising' and 'BNC_Falling'
SeqControl(MSSYNC).argument = 1000;  % 1 sec TIMEOUT
%
SYNC = 7;
SeqControl(SYNC).command = 'sync'; % synchronize SW and HW sequencers
SeqControl(SYNC).argument = 20e3;     % 20msec TIMEOUT
%
TRIGO = 8;
SeqControl(TRIGO).command = 'triggerOut';
%
RDMASYNC = 9;
SeqControl(RDMASYNC).command = 'rdmaSync';
%
TTNADUMMY = 10;
SeqControl(TTNADUMMY).command = 'timeToNextAcq';  % time between dummy TX and start of Frame
SeqControl(TTNADUMMY).argument = 10; % 10us minimum time

nsc = length(SeqControl)+1; % nsc is count of SeqControl objects

%% Event Sequence
n = 1; % n is count of Events
for i = 1:Resource.RcvBuffer(1).numFrames
    if strcmp(RDMAconfig.rdmaRole,'primary')
        Event(n).info = 'DUMMY TX';
        Event(n).tx = na+nDopAngles+1;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = TTNADUMMY;
        n=n+1;
    end
    
    for j = 1:na  %acquire Bmode
        Event(n).info = 'Full aperature.';
        Event(n).tx = j;
        Event(n).rcv = (na+nVol*nDopAngles)*(i-1)+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = TTNAFLASH;
        if (j == 1)
            Event(n).seqControl = [TRIGO, MSSYNC, TTNAFLASH] ;  %just sync on the 1st Doppler Bmode TX
        end
        n = n + 1;
    end
    
    for j = 1:nVol  %acquire ensemble
        for k = 1:nDopAngles
            Event(n).info = 'Full aperature.';
            Event(n).tx = (na+k);
            Event(n).rcv = (na+nVol*nDopAngles)*(i-1) + na + (j-1)*nDopAngles + k;
            Event(n).recon = 0;
            Event(n).process = 0;
            Event(n).seqControl = TTNAPRF;
            if (j==nVol)&&(k==nDopAngles)
                if strcmp(RDMAconfig.rdmaRole,'primary')
                    Event(n).seqControl = TTNAFRAME;  %the last secondary TX shouldn't have a TTNA, it will sync with the primary on the next frame
                else
                    Event(n).seqControl = 0;
                end
            end
            n = n + 1;
        end
    end
    
    Event(n).info = 'Transfer to Host';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = [nsc, nsc+1];
    SeqControl(nsc).command = 'transferToHost';
    SeqControl(nsc+1).command = 'waitForTransferComplete';
    SeqControl(nsc+1).argument = nsc;
    nsc = nsc + 2;
    n = n+1;
    
    Event(n).info = 'Reconstruct';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = [1, 2]; %reconstruct Bmode & Doppler
    Event(n).process = 0;
    Event(n).seqControl = 0;
    n = n+1;
    
    if strcmp(RDMAconfig.rdmaRole,'secondary')    %--- secondary ---%
        %--- this part will be hidden in future ---%
        RDMAcontrol(1).bufferType = 'inter';
        RDMAcontrol(1).primaryDstBufNum = 3;
        RDMAcontrol(1).primaryDstFrameNum = 1;
        RDMAcontrol(1).secondarySrcBufNum = 1;
        RDMAcontrol(1).secondarySrcFrameNum = 1;
        RDMAcontrol(1).secondaryIndex = RDMAconfig.nodeIndex;
        
        RDMAcontrol(2).bufferType = 'inter';
        RDMAcontrol(2).primaryDstBufNum = 4;
        RDMAcontrol(2).primaryDstFrameNum = 1;
        RDMAcontrol(2).secondarySrcBufNum = 2;
        RDMAcontrol(2).secondarySrcFrameNum  = 1;
        RDMAcontrol(2).secondaryIndex = RDMAconfig.nodeIndex;
        
        Event(n).info = 'RDMA transfer';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = [nsc, nsc+1]; %RDMA Synchronous write 1 Bmode, 1 for Doppler
        SeqControl(nsc).command = 'rdmaWrite';
        SeqControl(nsc).argument = 1; % link it to RDMA structure
        SeqControl(nsc+1).command = 'rdmaWrite';
        SeqControl(nsc+1).argument = 2; % link it to RDMA structure
        nsc = nsc + 2;
        n = n + 1;

        Event(n).info = 'MultiSystem SW Sync';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = RDMASYNC;
        n = n+1;
        
    else
        %--- primary ---%
        RDMAcontrol(1).bufferType = 'inter';
        RDMAcontrol(1).primaryDstBufNum = 3;
        RDMAcontrol(1).primaryDstFrameNum = 1;
        RDMAcontrol(1).secondarySrcBufNum = 1;
        RDMAcontrol(1).secondarySrcFrameNum  = 1;
        RDMAcontrol(1).secondaryIndex = 1;
        
        RDMAcontrol(2).bufferType = 'inter';
        RDMAcontrol(2).primaryDstBufNum = 4;
        RDMAcontrol(2).primaryDstFrameNum =1;
        RDMAcontrol(2).secondarySrcBufNum = 2;
        RDMAcontrol(2).secondarySrcFrameNum  = 1;
        RDMAcontrol(2).secondaryIndex = 1;
        
        Event(n).info = 'MultiSystem SW Sync';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = RDMASYNC;
        n = n+1;

        %--------------------%
        Event(n).info = 'Combine Bmode IQ buffers';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 4;
        Event(n).seqControl = 0;
        n = n+1;
        
        Event(n).info = 'Combine Doppler IQ buffers';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 5;
        Event(n).seqControl = 0;
        n = n+1;
        
        Event(n).info = ['Doppler Process ' num2str(i)];
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 6;
        Event(n).seqControl = 0;
        n = n+1;
        
        Event(n).info = ['Diplay Bmode Process XZ - Frame ' num2str(i)];
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 1;
        Event(n).seqControl = 0;
        n = n+1;
        
        Event(n).info = ['Diplay Bmode Process YZ - Frame ' num2str(i)];
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 2;
        Event(n).seqControl = 0;
        n = n+1;
        
        Event(n).info = ['Diplay Bmode Process XY - Frame ' num2str(i)];
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 3;
        Event(n).seqControl = 0;
        n = n+1;
        
        Event(n).info = ['Diplay Doppler Process XZ - Frame ' num2str(i)];
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 7;
        Event(n).seqControl = 0;
        n = n+1;
        
        Event(n).info = ['Diplay Doppler Process YZ - Frame ' num2str(i)];
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 8;
        Event(n).seqControl = 0;
        n = n+1;
        
        Event(n).info = ['Diplay Doppler Process XY - Frame ' num2str(i)];
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 9;
        Event(n).seqControl = 0;
        n = n+1;
        
    end
    
    Event(n).info = 'Exit to Matlab';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    if floor(i/frameRateFactor) == i/frameRateFactor     % Exit to Matlab every 5th frame
        Event(n).seqControl = RETMAT;
    else
        Event(n).seqControl = RDMASYNC;
    end
    n = n+1;
    
    if strcmp(RDMAconfig.rdmaRole,'primary')
        Event(n).info = 'Sync';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = SYNC;
        n = n + 1;
    end
    
end
Event(n).info = 'Jump back to first event';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = JUMP;
n = n+1;

% User specified UI Control Elements
% - Sensitivity Cutoff
UI(1).Control =  {'UserB7','Style','VsSlider','Label','Sens. Cutoff',...
    'SliderMinMaxVal',[0, 1.0, Recon(1).senscutoff],...
    'SliderStep',[0.025,0.1],'ValueFormat','%1.3f'};
UI(1).Callback = text2cell('%SensCutoff');

if strcmp(RDMAconfig.rdmaRole,'primary')
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
end
% Save all the structures to a .mat file.
save(['/tmp/' filename]);
return

% **** Callback routines to be encoded by text2cell function. ****
%SensCutoff - Sensitivity cutoff change
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
%SensCutoff

%-UI#2Callback - Doppler mode change
Control = repmat(struct('Command','set&Run','Parameters',[]),1,4);
Process = evalin('base','Process');
Resource = evalin('base','Resource');

switch UIState
    case 1  % Velocity mode
        newMap = grayscaleCFImap;
        newMap(1:128,:) = Resource.DisplayWindow(1).Colormap(1:128,:);
        Resource.DisplayWindow(1).Colormap = newMap;
        Resource.DisplayWindow(2).Colormap = newMap;
        Resource.DisplayWindow(3).Colormap = newMap;
        assignin('base','persp',get(findobj('Tag','UserB1Slider'),'Value'));
        persf = evalin('base','persf'); persValue = persf(1);
        Control(1).Parameters = {'Process',6,'method','computeCFIFreqEst'};% process 6 is doppler processing
        Control(2).Parameters = {'Process',7,'srcData','signedColor3D','persistMethod','dynamic','persistLevel',persValue};
        Control(3).Parameters = {'Process',8,'srcData','signedColor3D','persistMethod','dynamic','persistLevel',persValue};
        Control(4).Parameters = {'Process',9,'srcData','signedColor3D','persistMethod','dynamic','persistLevel',persValue};
        Control(5).Parameters = {'DisplayWindow',1,'colormap',newMap};
        Control(6).Parameters = {'DisplayWindow',2,'colormap',newMap};
        Control(7).Parameters = {'DisplayWindow',3,'colormap',newMap};
        Control(8).Parameters = {'ImageBuffer',1,'lastFrame',0};
        set(findobj('tag','UserB1Edit'),'String',num2str(persValue,'%3.0f'));
        set(findobj('tag','UserB1Slider'),'Value',persValue);
        assignin('base','DopState','freq');
        % Set modified Process attributes in base Matlab environment.
        Process(6).method = 'computeCFIFreqEst';
        for k = 1:2:length(Process(7).Parameters)
            if strcmp(Process(7).Parameters{k},'srcData')
                Process(7).Parameters{k+1} = 'signedColor3D';
                Process(8).Parameters{k+1} = 'signedColor3D';
                Process(9).Parameters{k+1} = 'signedColor3D';
            elseif strcmp(Process(7).Parameters{k},'persistMethod')
                Process(7).Parameters{k+1} = 'dynamic';
                Process(8).Parameters{k+1} = 'dynamic';
                Process(9).Parameters{k+1} = 'dynamic';
            elseif strcmp(Process(7).Parameters{k},'persistLevel')
                Process(7).Parameters{k+1} = persf;
                Process(8).Parameters{k+1} = persf;
                Process(9).Parameters{k+1} = persf;
            end
        end
    case 2  % Power mode
        newMap = grayscaleCPAmap;
        newMap(1:128,:) = Resource.DisplayWindow(1).Colormap(1:128,:);
        Resource.DisplayWindow(1).Colormap = newMap;
        Resource.DisplayWindow(2).Colormap = newMap;
        Resource.DisplayWindow(3).Colormap = newMap;
        for k = 1:2:length(Process(7).Parameters)
            if strcmp(Process(7).Parameters{k},'persistLevel'), persf = Process(7).Parameters{k+1}; end
        end
        assignin('base','persf',persf);
        persValue = evalin('base','persp');
        Control(1).Parameters = {'Process',6,'method','computeCFIPowerEst'};% process 6 is doppler processing
        Control(2).Parameters = {'Process',7,'srcData','unsignedColor3D','persistMethod','simple','persistLevel',persValue};
        Control(3).Parameters = {'Process',8,'srcData','unsignedColor3D','persistMethod','simple','persistLevel',persValue};
        Control(4).Parameters = {'Process',9,'srcData','unsignedColor3D','persistMethod','simple','persistLevel',persValue};
        Control(5).Parameters = {'DisplayWindow',1,'colormap',newMap};
        Control(6).Parameters = {'DisplayWindow',2,'colormap',newMap};
        Control(7).Parameters = {'DisplayWindow',3,'colormap',newMap};
        Control(8).Parameters = {'ImageBuffer',1,'lastFrame',0};
        set(findobj('tag','UserB1Edit'),'String',num2str(persValue,'%3.0f'));
        set(findobj('tag','UserB1Slider'),'Value',persValue);
        assignin('base','DopState','power');
        Process(2).method = 'computeCFIPowerEst';
        for k = 1:2:length(Process(7).Parameters)
            if strcmp(Process(7).Parameters{k},'srcData')
                Process(7).Parameters{k+1} = 'unsignedColor3D';
                Process(8).Parameters{k+1} = 'unsignedColor3D';
                Process(9).Parameters{k+1} = 'unsignedColor3D';
            elseif strcmp(Process(7).Parameters{k},'persistMethod')
                Process(7).Parameters{k+1} = 'simple';
                Process(8).Parameters{k+1} = 'simple';
                Process(9).Parameters{k+1} = 'simple';
            elseif strcmp(Process(7).Parameters{k},'persistLevel')
                Process(7).Parameters{k+1} = persValue;
                Process(8).Parameters{k+1} = persValue;
                Process(9).Parameters{k+1} = persValue;
            end
        end
end

assignin('base','newMap',newMap);
evalin('base','Resource.DisplayWindow(1).Colormap = newMap;');
evalin('base','Resource.DisplayWindow(2).Colormap = newMap;');
evalin('base','Resource.DisplayWindow(3).Colormap = newMap;');
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
if ishandle(hCMTool),
    delete(hCMTool);
    set(findobj('tag','toolsMenu'),'Value',1); % set tools selection back to none
end

return
%-UI#2Callback

%-UI#3Callback - Doppler Power change
Process = evalin('base','Process');
for k = 1:2:length(Process(6).Parameters)
    if strcmp(Process(6).Parameters{k},'pwrThreshold'), Process(6).Parameters{k+1} = UIValue; end
end
assignin('base','Process',Process);
% Set Control.Command to set Doppler threshold.
Control = evalin('base','Control');
Control.Command = 'set&Run';
Control.Parameters = {'Process',6,'pwrThreshold',UIValue};
assignin('base','Control', Control);
%-UI#3Callback

%-UI#4Callback - Color Priority change
% Set the value in the Process structure for use in cineloop playback.
Process = evalin('base','Process');
for k = 1:2:length(Process(7).Parameters)
    if strcmp(Process(7).Parameters{k},'threshold'), Process(7).Parameters{k+1} = UIValue; end
    if strcmp(Process(8).Parameters{k},'threshold'), Process(8).Parameters{k+1} = UIValue; end
    if strcmp(Process(9).Parameters{k},'threshold'), Process(9).Parameters{k+1} = UIValue; end
end
assignin('base','Process',Process);
% Set Control.Command to set Image.threshold.
Control = evalin('base','Control');
Control(1).Command = 'set&Run';
Control(1).Parameters = {'Process',7,'threshold',UIValue};
Control(2).Command = 'set&Run';
Control(2).Parameters = {'Process',8,'threshold',UIValue};
Control(3).Command = 'set&Run';
Control(3).Parameters = {'Process',9,'threshold',UIValue};
assignin('base','Control', Control);
%-UI#4Callback

%-UI#5Callback - Color Persistence change
% Set the value in the Process structure for use in cineloop playback.
Process = evalin('base','Process');
for k = 1:2:length(Process(7).Parameters)
    if strcmp(Process(7).Parameters{k},'persistLevel'), Process(7).Parameters{k+1} = UIValue; end
    if strcmp(Process(8).Parameters{k},'persistLevel'), Process(8).Parameters{k+1} = UIValue; end
    if strcmp(Process(9).Parameters{k},'persistLevel'), Process(9).Parameters{k+1} = UIValue; end
end
assignin('base','Process',Process);
% Set Control.Command to set Image.persistLevel.
Control = evalin('base','Control');
Control(1).Command = 'set&Run';
Control(1).Parameters = {'Process',7,'persistLevel',UIValue};
Control(2).Command = 'set&Run';
Control(2).Parameters = {'Process',8,'persistLevel',UIValue};
Control(3).Command = 'set&Run';
Control(3).Parameters = {'Process',9,'persistLevel',UIValue};
assignin('base','Control', Control);

% If PTool window is open, adjust persistLevel1 in Process(7)
hPTool = findobj('tag','ProcessTool');
if ishandle(hPTool),
    hPNum = findobj('tag','processNum');
    if isequal(get(findobj('tag','processNum'),'Value'),7)
        set(findobj('tag','persistSlider1'),'Value',UIValue);
        set(findobj('tag','persistValue1'),'String',num2str(UIValue));
    end
end
return
%-UI#5Callback
