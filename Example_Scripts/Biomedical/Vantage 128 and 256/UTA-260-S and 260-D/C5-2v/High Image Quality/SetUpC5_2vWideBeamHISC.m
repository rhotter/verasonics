% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use.
%
% File name: SetUpC5_2vWideBeamHISC.m - Example of curved array imaging with
%            wide beam transmits, spatial compounding and harmonic imaging.
%
% Description:
%   Sequence programming file for C5-2v curved array using wide beam transmit scanning with
%   harmonic imaging and spatial compounding. Of the 128 transmit channels, only
%   P.numTx are used, with P.numTx transmitters on each side of the center element
%   (where possible). All 128 receive channels are used, although the
%   element sensitivity cutoff will limit the useful aperture. The receive acquisitions
%   use 200% bandwidth. Reconstruction and Processing are asynchronous with
%   respect to acquisition.
%
%   The PData regions are generated to provide an overlap of 6 regions at every pixel.
%   The transmit focus is set to below the bottom of the PData region to generate a wide
%   beam that is as wide at the bottom of a single region.
%
%   Three wide beam scans are performed with P.numRays scan lines in each
%   frame using pulse inversion averaging for each beam direction. Each scan
%   is made with different steering directions (steered left, no steering,
%   steered right).  The frames are averaged after image processing using a
%   running average that updates each frame. This version acquires all three
%   steering directions into a 'super' frame to allow asynchronous acquisition
%   and processing.
%
% Last update:
% 09/11/2020 - Update to SW 4.3 format for new user UIControls and External function definitions (VTS-1691).
%   More info:(?/Example_Scripts/Vantage_Features/New UI Scheme/SetUpL11_5vFlash_NewUI)

clear all

demodFreq = 4.5;

P.numTx = 60;   % no. of elements in TX aperture.
P.numRays = 64; % no. of rays in frame
P.txFocusMm = 600; % focus in mm
P.startDepthMm = 2;  % startDepth in mm
P.endDepthMm = 120;
P.maxDepthMm = 160;  % maxDepth for RangeChange and RcvBuffer
P.dtheta = 8*(pi/180); % steering angle for beams

% Specify system parameters.
Resource.Parameters.numTransmit = 128;  % number of transmit channels.
Resource.Parameters.speedOfSound = 1540;
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Specify Trans structure array.
Trans.name = 'C5-2v';
Trans.units = 'wavelengths'; % required in Gen3 to prevent default to mm units
Trans = computeTrans(Trans);  % C5-2v transducer is 'known' transducer so we can use computeTrans.
Trans.maxHighVoltage = 60;  % set maximum high voltage limit for pulser supply.
radius = Trans.radius;
scanangle = Trans.numelements*Trans.spacing/radius;
rayDelta = scanangle/P.numRays; % angle between rays
theta = -(scanangle/2) + 0.5*rayDelta; % angle to left edge from centerline
Angle = theta:rayDelta:(-theta);

% Convert mm to wavelength
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);
P.txFocus = P.txFocusMm*scaleToWvl;
P.startDepth = P.startDepthMm*scaleToWvl;  % startDepth in wavelength
P.endDepth = P.endDepthMm*scaleToWvl;
P.maxDepth = P.maxDepthMm*scaleToWvl;
maxBufLength = ceil(sqrt((P.maxDepth+radius)^2 + radius^2 - ...
                     2*(P.maxDepth+radius)*radius*cos(scanangle)));
maxBufSizePerAcq = 128*ceil(maxBufLength*8*(demodFreq/Trans.frequency)/128); % demodulation frequency

% Specify PData structure array.
PData.PDelta = [1.0, 0, 0.5];  % x, y and z pdeltas
sizeRows = 10 + ceil((P.endDepth + radius - (radius * cos(scanangle/2)))/PData.PDelta(3));
sizeCols = 10 + ceil(2*(P.endDepth + radius)*sin(scanangle/2)/PData.PDelta(1));
PData.Size = [sizeRows,sizeCols,1];     % size, origin and pdelta set region of interest.
PData.Origin(1,1) = (P.endDepth+radius)*sin(-scanangle/2) - 5;
PData.Origin(1,2) = 0;
PData.Origin(1,3) = ceil(radius * cos(theta)) - radius - 5;
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

% Specify Resources.
rowsPerFrame = P.numRays*maxBufSizePerAcq;
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = rowsPerFrame; % this size allows for all rays, with range of up to 400 wvlngths
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numTransmit;
Resource.RcvBuffer(1).numFrames = 10;
Resource.RcvBuffer(2).datatype = 'int16';
Resource.RcvBuffer(2).rowsPerFrame = rowsPerFrame; % this size allows for all rays, with range of up to 400 wvlngths
Resource.RcvBuffer(2).colsPerFrame = Resource.Parameters.numTransmit;
Resource.RcvBuffer(2).numFrames = 10;  % all RcvBuffers should have same no. of frames
Resource.RcvBuffer(3).datatype = 'int16';
Resource.RcvBuffer(3).rowsPerFrame = rowsPerFrame; % this size allows for all rays, with range of up to 400 wvlngths
Resource.RcvBuffer(3).colsPerFrame = Resource.Parameters.numTransmit;
Resource.RcvBuffer(3).numFrames = 10;
Resource.InterBuffer.datatype = 'complex double';
Resource.InterBuffer.numFrames = 1;
Resource.ImageBuffer.datatype = 'double';
Resource.ImageBuffer.numFrames = 20;
Resource.DisplayWindow(1).Title = 'C5-2vWideBeamHISC';
Resource.DisplayWindow(1).pdelta = 0.45;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).numFrames = 50;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow.Colormap = gray(256);

% Specify Transmit waveform structure.
TW(1) = struct('type','parametric','Parameters',[demodFreq/2,.67,2,-1]);
TW(2) = struct('type','parametric','Parameters',[demodFreq/2,.67,2,1]);

% Specify TX structure array.
% - We need 128 transmit specifications.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'focus', P.txFocus, ...
                   'Steer', [0.0,0.0], ...
                   'Apod', zeros(1,Trans.numelements), ...
                   'Delay', zeros(1,Trans.numelements), ...
                   'TXPD', [], ...
                   'peakCutOff', 0.6, ...
                   'peakBLMax', 15.0), 1, 6*P.numRays);

% - Set event specific TX attributes.
%    P.numTx/2 is the number of elements to include on each side of the
%    center element, for the specified focus and sensitivity cutoff.
%    Thus the full transmit aperture will be P.numTx + 1 elements.
for n = 1:P.numRays   % P.numRays transmit events
    % Set transmit Origins.
    TX(n).Origin = [radius*sin(Angle(n)), 0.0, radius*cos(Angle(n))-radius];
    ce = round(1+(Trans.numelements-1)*(Angle(n) - theta)/(-2*theta));
    % Set transmit Apodization so that a maximum of numTx + 1 transmitters are active.
    L = P.numTx+1;
    W = blackman(L);
    w1 = 1;
    w2 = L;
    lft = round(ce - P.numTx/2);
    if lft < 1, w1 = -lft+2; lft = 1; end
    rt = round(ce + P.numTx/2);
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
    TX(n+m+P.numRays).Apod = TX(n+m).Apod;
    TX(n+m+P.numRays).Delay = TX(n+m).Delay;
end

% Calculating TXPD
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

BPF = ...
[-0.00058 +0.00024 +0.00070 +0.00208 -0.00366 -0.00076 -0.00290 ...
 +0.01392 -0.00729 +0.00186 -0.02802 +0.03461 -0.00436 +0.03033 ...
 -0.08264 +0.02850 +0.00638 +0.13232 -0.13400 -0.24313 +0.51282];

% [-0.00046 +0.00000 +0.00247 +0.00000 -0.00696 +0.00000 +0.01279 ...
%     +0.00000 -0.01477 +0.00000 +0.00409 +0.00000 +0.02704 +0.00000 ...
%     -0.07898 +0.00000 +0.14020 +0.00000 -0.19037 +0.00000 +0.20990], ...

% Specify Receive structure arrays.
% -- Compute the maximum receive path length, using the law of cosines.
maxAcqLength = ceil(sqrt((P.endDepth+radius)^2 + radius^2 - ...
                     2*(P.endDepth+radius)*radius*cos(scanangle)));
Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'demodFrequency', demodFreq, ...
                        'mode', 0, ...
                        'InputFilter', BPF,...
                        'callMediaFunc', 0), 1, 6*P.numRays*Resource.RcvBuffer(1).numFrames);

m = P.numRays;
% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    ind = 6*m*(i-1);
    Receive(ind+1).callMediaFunc = 1;
    for j = 1:m
        for k = 1:3 % for three buffers
        Receive(ind+j+(2*k-2)*m).bufnum = k;
        Receive(ind+j+(2*k-2)*m).framenum = i;
        Receive(ind+j+(2*k-2)*m).acqNum = j;
        Receive(ind+j+(2*k-1)*m).bufnum = k;
        Receive(ind+j+(2*k-1)*m).framenum = i;
        Receive(ind+j+(2*k-1)*m).acqNum = j;
        Receive(ind+j+(2*k-1)*m).mode = 1;
        end
    end
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [0,250,590,710,770,830,890,950];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Recon structure arrays.
Recon = repmat(struct('senscutoff', 0.5, ...
               'pdatanum', 1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...
               'rcvBufFrame',-1, ...
               'RINums',1:P.numRays), 1, 3);
% - Set specific Recon attributes.
Recon(2).RINums = (P.numRays+1):(2*P.numRays);
Recon(3).RINums = (2*P.numRays+1):(3*P.numRays);

% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 'accumIQ', ...  % default is to accumulate IQ data.
                   'Pre',[],...
                   'Post',[],...
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'scaleFactor', 0.33, ...
                   'regionnum', 0), 1, 3*P.numRays);
% - Set specific ReconInfo attributes.
ReconInfo(1).Pre = 'clearInterBuf';
for j = 1:P.numRays
    ReconInfo(j).txnum = j;
    ReconInfo(j).rcvnum = j;
    ReconInfo(j).regionnum = j+1;
end
ReconInfo(P.numRays).Post = 'IQ2IntensityImageBuf';
k = P.numRays;
ReconInfo(k+1).Pre = 'clearInterBuf';
for j = 1:P.numRays
    ReconInfo(j+k).txnum = 2*P.numRays+j;
    ReconInfo(j+k).rcvnum = 2*P.numRays+j;
    ReconInfo(j+k).regionnum = j+k+1;
end
ReconInfo(k+P.numRays).Post = 'IQ2IntensityImageBuf';
k = k + P.numRays;
ReconInfo(k+1).Pre = 'clearInterBuf';
for j = 1:P.numRays
    ReconInfo(j+k).txnum = 4*P.numRays+j;
    ReconInfo(j+k).rcvnum = 4*P.numRays+j;
    ReconInfo(j+k).regionnum = j+k+1;
end
ReconInfo(k+P.numRays).Post = 'IQ2IntensityImageBuf';

% Specify Process structure array.
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,... % frame number in src buffer (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure (defines output figure).
                         'norm',1,...        % normalization method(1 means fixed)
                         'pgain',5,...            % pgain is image processing gain
                         'reject',2,...
                         'persistMethod','simple',...
                         'persistLevel',20,...
                         'interpMethod','4pt',...
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','runAverage3',...
                         'compressMethod','power',...      % X^0.5 normalized to output word size
                         'compressFactor',45,...
                         'mappingMethod','full',...
                         'display',1,...      % display image after processing
                         'displayWindow',1};

% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; %  jump back to start.
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq'; % time between rays
SeqControl(2).argument = 250; % 250usec
SeqControl(3).command = 'timeToNextAcq';  % optional time between frames (not needed for synchronous operation)
SeqControl(3).argument = 5000;
SeqControl(4).command = 'returnToMatlab';
nsc = 5;

% Specify Event structure arrays.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    % Acquire frame with straight ahead wide beams
    k = 6*m*(i-1);
    for j = 1:m     % Acquire frame
        Event(n).info = 'Acquire beam, normal TX.';
        Event(n).tx = j;
        Event(n).rcv = j+k;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;

        Event(n).info = 'Acquire beam, inverted TX.';
        Event(n).tx = j+m;
        Event(n).rcv = j+k+m;
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
    for j = 1:m        % Acquire frame
        Event(n).info = 'Acquire beam, normal TX.';
        Event(n).tx = j+2*m;
        Event(n).rcv = j+k+2*m;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;

        Event(n).info = 'Acquire beam, inverted TX.';
        Event(n).tx = j+3*m;
        Event(n).rcv = j+k+3*m;
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
    for j = 1:m        % Acquire frame
        Event(n).info = 'Acquire beam, normal TX.';
        Event(n).tx = j+4*m;
        Event(n).rcv = j+k+4*m;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;

        Event(n).info = 'Acquire beam, inverted TX.';
        Event(n).tx = j+5*m;
        Event(n).rcv = j+k+5*m;
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
    Event(n).seqControl = 4;
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
UI(1).Control = VsSliderControl('LocationCode','UserB7',...
    'Label','Sens. Cutoff',...
    'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
    'SliderStep',[0.025,0.1],'ValueFormat','%1.3f',...
    'Callback',@SensCutOffCallback);

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
UI(2).Control = VsSliderControl('LocationCode','UserA1',...
    'Label',['Range (',AxesUnit,')'],'SliderMinMaxVal',MinMaxVal,...
    'SliderStep',[5/stepBase,10/stepBase],'ValueFormat','%3.0f',...
    'Callback',@RangeChangeCallback);

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 3;

% Save all the structures to a .mat file.
save('MatFiles/C5-2vWideBeamHISC');

%% Automatic VSX Execution:
% Uncomment the following line to automatically run VSX every time you run
% this SetUp script (note that if VSX finds the variable 'filename' in the
% Matlab workspace, it will load and run that file instead of prompting the
% user for the file to be used):

% filename = 'C5-2vWideBeamHISC';  VSX;


%% **** Callback routines used by UIControls (UI) ****

function SensCutOffCallback(~, ~, UIValue)
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

function RangeChangeCallback(hObject, ~, UIValue)
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
    for n = 1:P.numRays+1
        PData.Region(n).Shape.r2 = radius+P.endDepth;
    end
    % Define steered left regions
    m = P.numRays+1;  % include first mask region
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
        PData.Region(n+m).Shape.r2 = b+P.endDepth;
    end
    assignin('base','PData',PData);
    evalin('base','PData(1).Region = computeRegions(PData(1));');
    evalin('base','Resource.DisplayWindow(1).Position(3) = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);');
    evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);');
    evalin('base','Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];');

    % Update TXPD data of TX structures.
    TX = evalin('base','TX');
    PData = evalin('base','PData');
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