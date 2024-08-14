% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use.
%
% File name: SetUpC5_2vWideBeamHISC_64LE.m - Example of curved array imaging with
%            wide beam transmits and 3 angle spatial compounding
%
% Description:
%   Sequence program for C5-2v using wide beam transmit scanning with
%   pulse inversion harmonic imaging and spatial compounding.
%   For thewide beam transmits, P.numTx elements are used, with floor(P.numTx/2)
%   transmitters on each side of the center element (where possible). Transmit
%   origins are spaced equally across the array, starting with the element closest to
%   P.numTx/10 and ending with the element closest to 128-P.numTx/10. (By not centering
%   transmits on the end elements, the image quality at the edges of the scan is
%   improved.) The transmit focus is set to below the bottom of the scan to generate
%   a wide beam that is as wide at the bottom as a single PData region.

%   Three RcvBuffers are used to store wide beam scan frames acquired with three different steering
%   angles.  The different scan angles (steered left, no steering, steered right) are
%   acquired in series using P.numRays scan lines in each frame.  Each wide beam is
%   acquired with two complementary transmits with the receive data summed in local
%   memory on the acquisition modules.
%
%   The most recent frame in each buffer is processed and displayed sequentially.
%   Frames are averaged after image processing using a running average of three
%   frames that updates every frame. The transmit aperture size can be set by the
%   variable - numTx. The receive aperture always covers the full 128 element
%   aperture.
%
% Last update:
%   08/26/2021 - fixed index off by one in ReconInfo.regionnum

%clear[ 	]+all

P.numTx = 60;   % no. of elements in TX aperture. Should no larger than 63
P.numRays = 64; % no. of rays in frame
P.txFocusMm = 600; % focus in mm
P.startDepthMm = 2;  % startDepth in mm
P.endDepthMm = 120;
P.maxDepthMm = 160;  % maxDepth for RangeChange and RcvBuffer
P.dtheta = 8*(pi/180); % steering angle for beams (in radians)

% Specify system parameters.
Resource.Parameters.numTransmit = 128;
Resource.Parameters.speedOfSound = 1540;
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.
Resource.System.SoftwareVersion = [4 0 0]; % Minimum software release for this script.
Resource.System.UTA = '260-S'; % This script requires the 260-S UTA.

% Specify Trans structure array.
Trans.name = 'C5-2v';
Trans.units = 'wavelengths'; % required in Gen3 to prevent default to mm units
Trans = computeTrans(Trans);
radius = Trans.radius;
theta = -63.5*Trans.spacing/radius; % angle to element 1 from centerline
thetaLeft = -(63.5 - round(P.numTx/8))*Trans.spacing/radius; % angle to lft edge of scan
scanangle = 2*(-thetaLeft);
rayDelta = scanangle/(P.numRays-1); % angle between rays
Angle = thetaLeft:rayDelta:(-thetaLeft);

% Convert mm to wavelength
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);
P.txFocus = P.txFocusMm*scaleToWvl;
P.startDepth = P.startDepthMm*scaleToWvl;  % startDepth in wavelength
P.endDepth = P.endDepthMm*scaleToWvl;
P.maxDepth = P.maxDepthMm*scaleToWvl;
maxBufLength = ceil(sqrt((P.maxDepth+radius)^2 + radius^2 - ...
                     2*(P.maxDepth+radius)*radius*cos(scanangle)));
maxBufSizePerAcq = 128*ceil(maxBufLength*8/128);

% Specify PData structure array.
PData.PDelta = [1.0, 0, 0.5];  % x, y and z pdeltas
sizeRows = 10 + ceil((P.endDepth + radius - (radius * cos(scanangle/2)))/PData.PDelta(3));
sizeCols = 10 + ceil(2*(P.endDepth + radius)*sin(scanangle/2)/PData.PDelta(1));
PData.Size = [sizeRows,sizeCols,1];     % size, origin and pdelta set region of interest.
PData.Origin(1,1) = (P.endDepth+radius)*sin(-scanangle/2) - 5;
PData.Origin(1,2) = 0;
PData.Origin(1,3) = ceil(radius * cos(thetaLeft)) - radius - 5;
% Define the first PData region as the entire scan region to use as a mask.
PData.Region(1) = struct(...
    'Shape',struct('Name','Sector',...
                   'Position',[0,0,-radius],...
                   'r1',radius+P.startDepth,...
                   'r2',radius+P.endDepth,...
                   'angle',scanangle,...
                   'steer',0));
m = 1;
% Define PData Regions for P.numRays scanlines (steered straight ahead)
for n = 1:P.numRays
    PData.Region(n+m) = struct(...
        'Shape',struct('Name','Sector',...
                       'Position',[0,0,-radius],...
                       'r1',radius+P.startDepth,...
                       'r2',radius+P.endDepth,...
                       'angle',rayDelta*6,...
                       'steer',Angle(n),...
                       'andWithPrev',1));
end
% Define steered left regions
m = m + P.numRays;
for n = 1:P.numRays
    if n<10
        steer = -((n-1)/10)*P.dtheta;
    elseif n>(P.numRays-10)
        steer = -((P.numRays-n)/10)*P.dtheta;
    else
        steer = -P.dtheta;
    end
    d = radius*tan(-steer);
    b = sqrt(radius*radius + d*d);
    c = sqrt(2*radius*radius + d*d - 2*radius*b*cos(-steer));
    PData.Region(n+m) = struct(...
        'Shape',struct('Name','Sector',...
                       'Position',[c,0,-radius],...
                       'r1',b+P.startDepth,...
                       'r2',b+P.endDepth,...
                       'angle',rayDelta*6,...
                       'steer',(Angle(n)+steer),...
                       'andWithPrev',1));
end
m = m + P.numRays;
for n = 1:P.numRays
    if n<=10
        steer = ((n-1)/10)*P.dtheta;
    elseif n>(P.numRays-10)
        steer = ((P.numRays-n)/10)*P.dtheta;
    else
        steer = P.dtheta;
    end
    d = radius*tan(steer);
    b = sqrt(radius*radius + d*d);
    c = sqrt(2*radius*radius + d*d - 2*radius*b*cos(steer));
    PData.Region(n+m) = struct(...
        'Shape',struct('Name','Sector',...
                       'Position',[-c,0,-radius],...
                       'r1',b+P.startDepth,...
                       'r2',b+P.endDepth,...
                       'angle',rayDelta*6,...
                       'steer',(Angle(n)+steer),...
                       'andWithPrev',1));
end
PData.Region = computeRegions(PData);

%  Media points for curved array.
% - Uncomment for speckle
% Media.numPoints = (20000);
% Media.MP = rand(Media.numPoints,4);
% Media.MP(:,2) = 0;
% Media.MP(:,4) = 0.01 + 0.04*Media.MP(:,4);  % Random low amplitude
% RandR = P.endDepth*scaleToWvl*Media.MP(:,1)+radius;
% RandTheta = scanangle*(Media.MP(:,3)-0.5);
% Media.MP(:,1) = RandR.*sin(RandTheta);
% Media.MP(:,3) = RandR.*cos(RandTheta)-radius;
% - Define points
%Media.MP(1,:) = [0,0,radius+70,1.0];
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
Media.attenuation = -0.5;

rowsPerFrame = P.numRays*maxBufSizePerAcq;
% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 2*P.numRays*maxBufSizePerAcq; % this size allows for all rays, with range of up to 400 wvlngths
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numTransmit;
Resource.RcvBuffer(1).numFrames = numRcvFrames;
Resource.RcvBuffer(2).datatype = 'int16';
Resource.RcvBuffer(2).rowsPerFrame = 2*P.numRays*maxBufSizePerAcq; % this size allows for all rays, with range of up to 400 wvlngths
Resource.RcvBuffer(2).colsPerFrame = Resource.Parameters.numTransmit;
Resource.RcvBuffer(2).numFrames = numRcvFrames;  % all RcvBuffers should have same no. of frames
Resource.RcvBuffer(3).datatype = 'int16';
Resource.RcvBuffer(3).rowsPerFrame = 2*P.numRays*maxBufSizePerAcq; % this size allows for all rays, with range of up to 400 wvlngths
Resource.RcvBuffer(3).colsPerFrame = Resource.Parameters.numTransmit;
Resource.RcvBuffer(3).numFrames = numRcvFrames;
Resource.InterBuffer.datatype = 'complex double';
Resource.InterBuffer.numFrames = 1;
Resource.ImageBuffer.datatype = 'double';
Resource.ImageBuffer.numFrames = 20;
Resource.DisplayWindow(1).Title = 'C5-2vWideBeamHISC_64LE';
Resource.DisplayWindow(1).pdelta = 0.45;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).numFrames = numImageFrames;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow.Colormap = gray(256);

% Specify Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,.67,2,-1];
TW(2).type = 'parametric';
TW(2).Parameters = [Trans.frequency,.67,2,1];

% Specify TX structure array.
% - We need 128 transmit specifications.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'focus', P.txFocus, ...
                   'Steer', [0.0,0.0], ...
                   'Apod', zeros(1,Resource.Parameters.numTransmit), ...
                   'Delay', zeros(1,Resource.Parameters.numTransmit), ...
                   'TXPD', [], ...
                   'peakCutOff', 0.6, ...
                   'peakBLMax', 10.0), 1, 6*P.numRays);

% - Set event specific TX attributes.
%    floor(P.numTx/2) is the number of elements to include on each side of the
%    transmit origin center element, for the specified focus and sensitivity cutoff.
if floor(P.numTx/2) == P.numTx/2, L = P.numTx+1; else, L = P.numTx; end
W = kaiser(L,2.0);
Ce = zeros(1,P.numRays);
for n = 1:P.numRays   % P.numRays transmit events
    % Set transmit Origins.
    TX(n).Origin = [radius*sin(Angle(n)), 0.0, radius*cos(Angle(n))-radius];
    Ce(n) = round(1+(Trans.numelements-1)*(Angle(n) - theta)/(-2*theta));
    % Set transmit Apodization so that a maximum of numTx transmitters are active.
    w1 = 1;
    w2 = L;
    lft = Ce(n) - floor(P.numTx/2);
    if lft < 1, w1 = -lft+2; lft = 1; end
    rt = Ce(n) + floor(P.numTx/2);
    if rt > Trans.numelements, w2 = L - (rt - Trans.numelements); rt = Trans.numelements; end
    TX(n).Apod(lft:rt) = W(w1:w2);
    TX(n).Delay = computeTXDelays(TX(n));
    % Inverted pulse TXs
    TX(n+P.numRays).waveform = 2;
    TX(n+P.numRays).Origin = TX(n).Origin;
    TX(n+P.numRays).Apod = TX(n).Apod;
    TX(n+P.numRays).Delay = TX(n).Delay;
end
m = 2*P.numRays;
for n = 1:P.numRays   % P.numRays transmit events
    TX(n+m).Origin = TX(n).Origin;
    TX(n+m).Apod = TX(n).Apod;
    if n<=10
        TX(n+m).Steer = [-((n-1)/10)*P.dtheta,0.0];
    elseif n>(P.numRays-10)
        TX(n+m).Steer = [-((P.numRays-n)/10)*P.dtheta,0.0];
    else
        TX(n+m).Steer = [-P.dtheta,0.0];
    end
    TX(n+m).Delay = computeTXDelays(TX(n+m));
    % Inverted pulse TXs
    TX(n+m+P.numRays).waveform = 2;
    TX(n+m+P.numRays).Steer = TX(n+m).Steer;
    TX(n+m+P.numRays).Origin = TX(n+m).Origin;
    TX(n+m+P.numRays).Apod = TX(n+m).Apod;
    TX(n+m+P.numRays).Delay = TX(n+m).Delay;
end
m = m + 2*P.numRays;
for n = 1:P.numRays   % P.numRays transmit events
    TX(n+m).Origin = TX(n).Origin;
    TX(n+m).Apod = TX(n).Apod;
    if n<=10
        TX(n+m).Steer = [((n-1)/10)*P.dtheta,0.0];
    elseif n>(P.numRays-10)
        TX(n+m).Steer = [((P.numRays-n)/10)*P.dtheta,0.0];
    else
        TX(n+m).Steer = [P.dtheta,0.0];
    end
    TX(n+m).Delay = computeTXDelays(TX(n+m));
    % Inverted pulse TXs
    TX(n+m+P.numRays).waveform = 2;
    TX(n+m+P.numRays).Origin = TX(n+m).Origin;
    TX(n+m+P.numRays).Steer = TX(n+m).Steer;
    TX(n+m+P.numRays).Apod = TX(n+m).Apod;
    TX(n+m+P.numRays).Delay = TX(n+m).Delay;
end

% calculate TXPD and TX delay
h = waitbar(0,'Program TX parameters, please wait!');
steps = 3*P.numRays;
for i = 1:3
    for j = 1:P.numRays
        ind = P.numRays*2*(i-1)+j;
        TX(ind).TXPD = computeTXPD(TX(ind),PData);
        TX(ind+P.numRays).TXPD = TX(ind).TXPD;
        waitbar((P.numRays*(i-1)+j)/steps)
    end
end
close(h)

% Specify TGC Waveform structure.
TGC.CntrlPts = [0,250,538,658,770,830,890,950];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);


BPF = ...
[-0.00058 +0.00024 +0.00070 +0.00208 -0.00366 -0.00076 -0.00290 ...
 +0.01392 -0.00729 +0.00186 -0.02802 +0.03461 -0.00436 +0.03033 ...
 -0.08264 +0.02850 +0.00638 +0.13232 -0.13400 -0.24313 +0.51282];

% Specify Receive structure arrays.
% -- Compute the maximum receive path length.
maxAcqLength = ceil(sqrt((P.endDepth+radius)^2 + radius^2 - ...
                     2*(P.endDepth+radius)*radius*cos(scanangle)));
Receive = repmat(struct('Apod', zeros(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, 12*P.numRays*Resource.RcvBuffer(1).numFrames);
% Set event specific Receive attributes for each frame.
% - Receives are defined for each of the three RcvBuffers on a frame by frame basis:
%   [Buf1_Frm1_Receives, Buf1_Frm2_Receives, Buf1_Frm3_Receives, ..., Buf1_FrmN_Receives,
%    Buf2_Frm1_Receives, Buf2_Frm2_Receives, Buf2_Frm3_Receives, ..., Buf2_FrmN_Receives,
%    Buf3_Frm1_Receives, Buf3_Frm2_Receives, Buf3_Frm3_Receives, ..., Buf3_FrmN_Receives]
% - Each frame in a buffer consists of 4*P.numRays acquisitions, where each widebeam
%   Ray has 4 Receives - 2 Receives for synthetic aperture, and 2 for
%   pulse inversion harmonic imaging. Each buffer then has 4*P.numRays*numFrames
%   acquisitions.
nr = P.numRays;
for nb = 1:3  % for each of 3 buffers (steering directions)
    m = 4*P.numRays*Resource.RcvBuffer(1).numFrames*(nb-1); % m sets starting index for buffer
    for i = 1:Resource.RcvBuffer(1).numFrames
        k = 4*nr*(i-1);  % k sets starting Receive for a frame
        Receive(k+1).callMediaFunc = 1;
        for j = 1:4:4*nr  % j counts acquisitions for each frame.
            % Receives for a Ray
            %   Synthetic aperture 1, normal & inverted transmits
            Receive(m+k+j).bufnum = nb;
            Receive(m+k+j).Apod(1:64) = 1.0;
            Receive(m+k+j).framenum = i;
            Receive(m+k+j).acqNum = (j+1)/2;
            Receive(m+k+j+1).bufnum = nb;
            Receive(m+k+j+1).Apod(1:64) = 1.0;
            Receive(m+k+j+1).framenum = i;
            Receive(m+k+j+1).acqNum = (j+1)/2;
            Receive(m+k+j+1).mode = 1;
            %   Synthetic aperture 2, normal and inverted transmits
            Receive(m+k+j+2).bufnum = nb;
            Receive(m+k+j+2).Apod(65:128) = 1.0;
            Receive(m+k+j+2).framenum = i;
            Receive(m+k+j+2).acqNum = (j+1)/2 + 1;
            Receive(m+k+j+3).bufnum = nb;
            Receive(m+k+j+3).Apod(65:128) = 1.0;
            Receive(m+k+j+3).framenum = i;
            Receive(m+k+j+3).acqNum = (j+1)/2 + 1;
            Receive(m+k+j+3).mode = 1;
        end
    end
end


% Specify Recon structure arrays.
Recon = repmat(struct('senscutoff', 0.5, ...
               'pdatanum', 1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...
               'rcvBufFrame',-1, ...
               'RINums',1:2*P.numRays), 1, 3);
% - Set specific Recon attributes.
Recon(2).RINums = (2*P.numRays+1):(4*P.numRays);
Recon(3).RINums = (4*P.numRays+1):(6*P.numRays);

% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 'accumIQ', ...  % default is to accumulate IQ data.
                   'Pre',[],...
                   'Post',[],...
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'scaleFactor', .25, ...
                   'regionnum', 0), 1, 6*P.numRays);
% - Set specific ReconInfo attributes.
%   No steering
ReconInfo(1).Pre = 'clearInterBuf';
for j = 1:2:2*P.numRays
    ReconInfo(j).txnum = (j+1)/2;
    ReconInfo(j).rcvnum = 2*j-1;
    ReconInfo(j).regionnum = 1+(j+1)/2;
    ReconInfo(j+1).txnum = (j+1)/2;
    ReconInfo(j+1).rcvnum = 2*j+1;
    ReconInfo(j+1).regionnum = 1+(j+1)/2;
end
ReconInfo(2*P.numRays).Post = 'IQ2IntensityImageBuf';
%   Steered left
m = 4*P.numRays*Resource.RcvBuffer(1).numFrames;
k = 2*P.numRays;
ReconInfo(k+1).Pre = 'clearInterBuf';
for j = 1:2:2*P.numRays
    ReconInfo(j+k).txnum = (j+1)/2 + 2*P.numRays;
    ReconInfo(j+k).rcvnum = m+2*j-1;
    ReconInfo(j+k).regionnum = 1+(j+1)/2 + P.numRays;
    ReconInfo(j+k+1).txnum = (j+1)/2 + 2*P.numRays;
    ReconInfo(j+k+1).rcvnum = m+2*j+1;
    ReconInfo(j+k+1).regionnum = 1+(j+1)/2 + P.numRays;
end
ReconInfo(k+2*P.numRays).Post = 'IQ2IntensityImageBuf';
%   Steered right
m = m + 4*P.numRays*Resource.RcvBuffer(1).numFrames;
k = k+2*P.numRays;
ReconInfo(k+1).Pre = 'clearInterBuf';
for j = 1:2:2*P.numRays
    ReconInfo(j+k).txnum = (j+1)/2 + 4*P.numRays;
    ReconInfo(j+k).rcvnum = m+2*j-1;
    ReconInfo(j+k).regionnum = 1+(j+1)/2 + 2*P.numRays;
    ReconInfo(j+k+1).txnum = (j+1)/2 + 4*P.numRays;
    ReconInfo(j+k+1).rcvnum = m+2*j+1;
    ReconInfo(j+k+1).regionnum = 1+(j+1)/2 + 2*P.numRays;
end
ReconInfo(k+2*P.numRays).Post = 'IQ2IntensityImageBuf';

% Specify Process structure array.
pers = 20;
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,... % frame number in src buffer (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure (defines output figure).
                         'norm',1,...        % normalization method(1 means fixed)
                         'pgain',2.5,...            % pgain is image processing gain
                         'reject',2,...
                         'persistMethod','simple',...
                         'persistLevel',pers,...
                         'interpMethod','4pt',...
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','runAverage3',...
                         'compressMethod','power',...
                         'compressFactor',45,...
                         'mappingMethod','full',...
                         'display',1,...      % display image after processing
                         'displayWindow',1};

% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; %  jump back to start.
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq'; % time between rays
SeqControl(2).argument = 260; % 200usec
SeqControl(3).command = 'timeToNextAcq';  % optional time between frames
SeqControl(3).argument = 20000;  % 5000 usec = 5msec time between frames
SeqControl(4).command = 'returnToMatlab';
nsc = 5;

% Specify Event structure arrays.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    % Acquire frame with straight ahead wide beams
    k = 4*P.numRays*(i-1);
    for j = 1:4:4*P.numRays       % Acquire frame
        Event(n).info = '1st half of acquire aperture.';
        Event(n).tx = (j+3)/4;
        Event(n).rcv = k+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;

        Event(n).info = '1st half of acquire beam, inverted TX.';
        Event(n).tx = (j+3)/4+ P.numRays;
        Event(n).rcv = k+j+1;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;

        Event(n).info = '2nd half of aperture, normal TX.';
        Event(n).tx = (j+3)/4;
        Event(n).rcv = k+j+2;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;

        Event(n).info = '2nd half of aperture, inverted TX.';
        Event(n).tx = (j+3)/4 + P.numRays;
        Event(n).rcv = k+j+3;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;
    end
    Event(n-1).seqControl = [3,nsc]; % modify last acquisition Event's seqControl
      SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
      nsc = nsc+1;

    Event(n).info = 'recon and process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 1;
    Event(n).seqControl = 0;
    n = n+1;

    % Acquire frame with steered left wide beams
    % Acquire frame with steered left wide beams
    m = 4*P.numRays*Resource.RcvBuffer(1).numFrames;
    for j = 1:4:4*P.numRays        % Acquire frame
        Event(n).info = '1st half of aperture, normal TX.';
        Event(n).tx = (j+3)/4 + 2*P.numRays;
        Event(n).rcv = m+k+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;
        Event(n).info = '1st half of aperture, inverted TX.';
        Event(n).tx = (j+3)/4 + 3*P.numRays;
        Event(n).rcv = m+k+j+1;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;

        Event(n).info = '2nd half of aperture, normal TX.';
        Event(n).tx = (j+3)/4 + 2*P.numRays;
        Event(n).rcv = m+k+j+2;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;
        Event(n).info = '2nd half of aperture, inverted TX.';
        Event(n).tx = (j+3)/4 + 3*P.numRays;
        Event(n).rcv = m+k+j+3;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;
    end
    Event(n-1).seqControl = [3,nsc]; % modify last acquisition Event's seqControl
      SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
      nsc = nsc+1;

    Event(n).info = 'recon and process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 2;
    Event(n).process = 1;
    Event(n).seqControl = 0;
    n = n+1;

    % Acquire frame with steered right wide beams
    m = m + 4*P.numRays*Resource.RcvBuffer(1).numFrames;
    for j = 1:4:4*P.numRays        % Acquire frame
        Event(n).info = '1st half of aperture, normal TX.';
        Event(n).tx = (j+3)/4 + 4*P.numRays;
        Event(n).rcv = m+k+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;
        Event(n).info = '1st half of aperture, inverted TX.';
        Event(n).tx = (j+3)/4 + 5*P.numRays;
        Event(n).rcv = m+k+j+1;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;

        Event(n).info = '2nd half of aperture, normal TX.';
        Event(n).tx = (j+3)/4 + 4*P.numRays;
        Event(n).rcv = m+k+j+2;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;
        Event(n).info = '2nd half of aperture, inverted TX.';
        Event(n).tx = (j+3)/4 + 5*P.numRays;
        Event(n).rcv = m+k+j+3;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;
    end
    Event(n-1).seqControl = [3,nsc]; % modify last acquisition Event's seqControl
      SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
      nsc = nsc+1;

    Event(n).info = 'recon and process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 3;
    Event(n).process = 1;
    Event(n).seqControl = 0;
    if floor(i/3) == i/3     % Exit to Matlab every 3rd frame
        Event(n).seqControl = 4;
    end
    n = n+1;
end

Event(n).info = 'Jump back';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 1;

% User specified UI Control Elements
import vsv.seq.uicontrol.VsSliderControl

% - Sensitivity Cutoff
UI(1).Control = VsSliderControl('LocationCode','UserB7','Label','Sens. Cutoff',...
                   'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                   'SliderStep',[0.025,0.1],'ValueFormat','%1.3f',...
                   'Callback', @SensCutoffCallback);

% - Range Change
MinMaxMm = [20,P.maxDepthMm]; % min max in mm
if isfield(Resource.DisplayWindow(1),'AxesUnits') && strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
    AxesUnit = 'mm';
    MinMaxVal = [MinMaxMm, P.endDepth/scaleToWvl];
else
    AxesUnit = 'wls';
    MinMaxVal = [MinMaxMm*scaleToWvl, P.endDepth];
end
stepBase = MinMaxVal(2)-MinMaxVal(1);
UI(2).Control = VsSliderControl('LocationCode','UserA1','Label',['Range (',AxesUnit,')'],...
                  'SliderMinMaxVal',MinMaxVal,...
                  'SliderStep',[5/stepBase,10/stepBase],'ValueFormat','%3.0f', ...
                  'Callback', @RangeChangeCallback );

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 6;

% Save all the structures to a .mat file.
filename = fullfile( vsv.file.getVSXDir, 'MatFiles/SetUpC5_2vWideBeamHISC_64LE_QSApp.mat');

save(filename);



% **** Callback routines used by UI Controls ****
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

function RangeChangeCallback(hObject,~,UIValue)
    simMode = evalin('base','Resource.Parameters.simulateMode');
    % No range change if in simulate mode 2.
    if simMode == 2
        set(hObject,'Value',evalin('base','P.endDepth'));
        return
    end
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
    assignin('base','P',P);

    % Modify PData for new range
    PData = evalin('base','PData');
    radius = evalin('base','radius');
    scanangle = evalin('base','scanangle');

    sizeRows = 10 + ceil((P.endDepth + radius - (radius * cos(scanangle/2)))/PData(1).PDelta(3));
    sizeCols = 10 + ceil(2*(P.endDepth + radius)*sin(scanangle/2)/PData(1).PDelta(1));
    PData(1).Size = [sizeRows,sizeCols,1];     % size, origin and pdelta set region of interest.
    PData(1).Origin(1,1) = (P.endDepth+radius)*sin(-scanangle/2) - 5;
    % Define PData Regions for P.numRays scanlines
    PData.Region(1).r2 = radius+P.endDepth;
    m = 1;
    for n = 1:P.numRays
        PData.Region(n+m).Shape.r2 = radius+P.endDepth;
    end
    % Define steered left regions
    m = m+P.numRays;
    for n = 1:P.numRays
        if n<10
            steer = -((n-1)/10)*P.dtheta;
        elseif n>(P.numRays-10)
            steer = -((P.numRays-n)/10)*P.dtheta;
        else
            steer = -P.dtheta;
        end
        d = radius*tan(-steer);
        b = sqrt(radius*radius + d*d);
        c = sqrt(2*radius*radius + d*d - 2*radius*b*cos(-steer));
        PData.Region(n+m).Position(1) = c;
        PData.Region(n+m).Shape.r2 = b+P.endDepth;
    end
    m = m + P.numRays;
    for n = 1:P.numRays
        if n<=10
            steer = ((n-1)/10)*P.dtheta;
        elseif n>(P.numRays-10)
            steer = ((P.numRays-n)/10)*P.dtheta;
        else
            steer = P.dtheta;
        end
        d = radius*tan(steer);
        b = sqrt(radius*radius + d*d);
        c = sqrt(2*radius*radius + d*d - 2*radius*b*cos(steer));
        PData.Region(n+m).Position(1) = -c;
        PData.Region(n+m).Shape.r2 = b+P.endDepth;
    end
    PData.Region = computeRegions(PData);

    assignin('base','PData',PData);

    evalin('base','Resource.DisplayWindow(1).Position(3) = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);');
    evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);');
    evalin('base','Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];');

    % Update TXPD data of TX structures.
    TX = evalin('base','TX');
    h = waitbar(0,'Program TX parameters, please wait!');
    steps = 3*P.numRays;
    for i = 1:3
        for j = 1:P.numRays
            ind = P.numRays*2*(i-1)+j;
            TX(ind).TXPD = computeTXPD(TX(ind),PData);
            TX(ind+P.numRays).TXPD = TX(ind).TXPD;
            waitbar((P.numRays*(i-1)+j)/steps)
        end
    end
    close(h)
    assignin('base','TX',TX);
    % Update Receive structures
    Receive = evalin('base', 'Receive');
    maxAcqLength = ceil(sqrt((P.endDepth+radius)^2 + radius^2 - ...
                         2*(P.endDepth+radius)*radius*cos(scanangle)));
    for i = 1:size(Receive,2)
        Receive(i).endDepth = maxAcqLength;
    end
    assignin('base','Receive',Receive);
    evalin('base','TGC.rangeMax = P.endDepth;');
    evalin('base','TGC.Waveform = computeTGCWaveform(TGC);');
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'PData','InterBuffer','ImageBuffer','TX','Receive','TGC','Recon','DisplayWindow'};
    assignin('base','Control', Control);
    assignin('base', 'action', 'displayChange');
end