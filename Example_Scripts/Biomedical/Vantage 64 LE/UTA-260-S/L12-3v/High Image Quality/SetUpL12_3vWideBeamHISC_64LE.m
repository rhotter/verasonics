% Notice: This file is provided by Verasonics to end users as a programming
%    example for the Verasonics Vantage Research Ultrasound System.
%    Verasonics makes no claims as to the functionality or intended
%    application of this program and the user assumes all responsibility for
%    its use.
%
% File Name: SetUpL12_3vWideBeamHISC_64LE.m - Wide beam, Harmonic Imaging, Spatial Compounding example
%
% Description: The L12-3v 38mm Linear array is used to acquire frames using wide beam pulse inversion transmits.
%    The wide beam transmits are steered between successive frames - left, straight ahead, right.
%    The reconstructed intensity data from each frame are averaged using a running average of 3
%    frames after image processing.
%    Wide beam pulse inversion transmits and receive acquisitions are used with a 128 element aperture that is
%    moved across the 192 elements of the array.  The active transmit aperture is a subset of the
%    full 128 element aperture.
%
%    Example apertures: In the following aperture examples, each space represents 2 elements. The
%    number of transmitting elements is given by numTx.
%       Aperture 1 (numTx=48):
%    tttttttttttt---------------------------------------------------\--------------------------------
%    rrrrrrrrrrrrrrrrrrrrrrrrrr----------------------------------------------------------------------
%       Aperture for almost middle wide beam 24 (numTx=48, numRays=48):
%    ---------------/-------------------tttttttttttttttttttttttt-------------------\-----------------
%    -----------------------------rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr------------------------------
%       Aperture for last wide beam (numTx=48, numRays=48):
%    --------------------------------/---------------------------------------------------tttttttttttt
%    -------------------------------------------------------------------rrrrrrrrrrrrrrrrrrrrrrrrrrrrr
%
%   The receive data from each aperture for each wide beam are stored under different acqNums
%   in the Receive buffer. The reconstruction sums the IQ data from the aquisitions for all beams
%   that contribute to a pixel and computes intensity values to produce the full frame.
%   This version does synchronous acquisition and processing.
%
% Last update:
% 08/22/2016 - tested with software 3.2

clear all

% key parameters
P.numTx = 48; % Number of Transmit Elements in the aperture (should be even).
P.numRays = 48;  % Number of wide beams in a frame.
P.txFocusMm = 200;   % transmit focal pt in wavelengths
P.dtheta = 10*(pi/180);  % angle delta for beams steered left or right
P.startDepthMm = 0.5;  % startDepth in mm
P.endDepthMm = 50;  % endDepth in mm
P.maxDepthMm = 50;  % maxDepth for RangeChange and RcvBuffer

% Define system parameters.
Resource.Parameters.numTransmit = 128;  % number of transmit channels.
Resource.Parameters.numRcvChannels = 64;  % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Specify Trans structure array.
Trans.name = 'L12-3v';
Trans.units = 'wavelengths'; % required in Gen3 to prevent default to mm units
Trans = computeTrans(Trans);  % L12-3v transducer is 'known' transducer so we can use computeTrans.
% - nominal center frequency from computeTrans is 7.813 MHz

% Convert mm to wavelength
demodFreq = 9; % demodulation frequency
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);
P.txFocus = P.txFocusMm*scaleToWvl;
P.startDepth = P.startDepthMm*scaleToWvl;  % startDepth in wavelength
P.endDepth = P.endDepthMm*scaleToWvl;
P.maxDepth = P.maxDepthMm*scaleToWvl;
maxBufLength = ceil(sqrt(P.maxDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
maxBufSizePerAcq = 128*ceil(maxBufLength*8*(demodFreq/Trans.frequency)/128);

% Specify PData structure arrays.
% - PData contains rectangular regions for no beam steering.
PData.PDelta = [1.0, 0, 0.5];
PData.Size(1) = ceil((P.endDepth-P.startDepth)/PData.PDelta(3)); % P.startDepth, P.endDepth and pdelta set PData.Size.
PData.Size(2) = ceil((Trans.numelements*Trans.spacing)/PData.PDelta(1));
PData.Size(3) = 1;      % single image page
PData.Origin = [-Trans.spacing*(Trans.numelements/2),0,P.startDepth]; % x,y,z of upper lft crnr.
% Define numRays rectangular regions, positioned around the origin of the wide beams.
TxOrgX = (-95.5*Trans.spacing):(191*Trans.spacing/(P.numRays-1)):(95.5*Trans.spacing); % x coords of beam centers
for n = 1:P.numRays
    PData.Region(n) = struct('Shape',struct('Name','Rectangle','Position',[TxOrgX(n),0,0],'width',32,'height',P.endDepth));
end
m = P.numRays;
% Define numRays steered left parallelogram regions, positioned around the origin of the wide beams.
for n = 1:P.numRays
    if n<=10
        angle = -((n-1)/10)*P.dtheta;
    elseif n>(P.numRays-10)
        angle = -((P.numRays-n)/10)*P.dtheta;
    else
        angle = -P.dtheta;
    end
    PData.Region(n+m) = struct('Shape',struct('Name','Parallelogram','Position',[TxOrgX(n),0,0],'width',32,'height',P.endDepth/cos(P.dtheta),'angle',angle));
end
m = m + P.numRays;
% Define numRays steered right parallelogram regions, positioned around the origin of the wide beams.
for n = 1:P.numRays
    if n<=10
        angle = ((n-1)/10)*P.dtheta;
    elseif n>(P.numRays-10)
        angle = ((P.numRays-n)/10)*P.dtheta;
    else
        angle = P.dtheta;
    end
    PData.Region(n+m) = struct('Shape',struct('Name','Parallelogram','Position',[TxOrgX(n),0,0],'width',32,'height',P.endDepth/cos(P.dtheta),'angle',angle));
end
PData.Region = computeRegions(PData);

% Specify Media object. 'pt1.m' script defines array of point targets.
pt1;
Media.attenuation = -0.5;
Media.function = 'movePoints';

% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = P.numRays*maxBufSizePerAcq; % this size allows for all rays, with maximum range
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 12;
Resource.RcvBuffer(2).datatype = 'int16';
Resource.RcvBuffer(2).rowsPerFrame = P.numRays*maxBufSizePerAcq; % this size allows for all rays, with maximum range
Resource.RcvBuffer(2).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(2).numFrames = 12;  % all RcvBuffers should have same no. of frames
Resource.RcvBuffer(3).datatype = 'int16';
Resource.RcvBuffer(3).rowsPerFrame = P.numRays*maxBufSizePerAcq; % this size allows for all rays, with maximum range
Resource.RcvBuffer(3).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(3).numFrames = 12;
Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).numFrames = 1;  % one intermediate buffer needed.
Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).numFrames = 10;
Resource.DisplayWindow.Title = 'L12-3vWideBeamHISC';
Resource.DisplayWindow(1).pdelta = 0.4;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).numFrames = 50;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = gray(256);

% Specify Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [4.5,0.67,1,1];

TW(2).type = 'parametric';
TW(2).Parameters = [4.5,0.67,1,-1];

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'aperture', 1, ...
                   'Apod', zeros(1,Resource.Parameters.numTransmit), ...
                   'focus', P.txFocus, ...
                   'Steer', [0.0,0.0], ...
                   'TXPD', [], ...
                   'peakCutOff', 1.0, ...
                   'peakBLMax', 15.0, ...
                   'Delay', zeros(1,Resource.Parameters.numTransmit)), 1, 6*P.numRays);

scaleToWvl = 1;
if strcmp(Trans.units, 'mm')
    scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);
end

% - Set event specific TX attributes.
for n = 1:P.numRays
    TX(n).Origin(1) = TxOrgX(n);
    % Compute transmit aperture number (first element of 128 element aperture)
    [Dummy,ce(n)] = min(abs(scaleToWvl*Trans.ElementPos(:,1)-TxOrgX(n)));
    lft = ce(n) - 64;
    if lft < 1, lft = 1; end
    if lft > 65, lft = 65; end
    TX(n).aperture = lft;
    % Compute TX.Apod within mux aperture.
    lft = ce(n) - P.numTx/2;
    if lft < 1, lft = 1; end
    rt = ce(n) + P.numTx/2;
    if rt > 192, rt = 192; end
    TX(n).Apod((lft-(TX(n).aperture-1)):(rt-(TX(n).aperture-1))) = 1;
%     % Apply apodization function.
%     [RIndices,CIndices,V] = find(TX(n).Apod);
%     V = kaiser(size(V,2),1);
%     TX(n).Apod(CIndices) = V;
    % Compute transmit delays
    TX(n).Delay = computeTXDelays(TX(n));
    % Define parameters for inverted transmit
    TX(n+P.numRays).waveform = 2;
    TX(n+P.numRays).Origin(1) = TxOrgX(n);
    TX(n+P.numRays).aperture = TX(n).aperture;
    TX(n+P.numRays).Apod = TX(n).Apod;
    TX(n+P.numRays).Steer = TX(n).Steer;
    TX(n+P.numRays).Delay = TX(n).Delay;
end
m = 2*P.numRays;
for n = 1:P.numRays
    TX(n+m).Origin(1) = TX(n).Origin(1);
    TX(n+m).aperture = TX(n).aperture;
    TX(n+m).Apod = TX(n).Apod;
    if n<=10
        TX(n+m).Steer = [-((n-1)/10)*P.dtheta,0.0];
    elseif n>(P.numRays-10)
        TX(n+m).Steer = [-((P.numRays-n)/10)*P.dtheta,0.0];
    else
        TX(n+m).Steer = [-P.dtheta,0.0];
    end
    % Compute transmit delays
    TX(n+m).Delay = computeTXDelays(TX(n+m));
    % Define parameters for inverted transmit
    TX(n+m+P.numRays).waveform = 2;
    TX(n+m+P.numRays).Origin(1) = TxOrgX(n);
    TX(n+m+P.numRays).aperture = TX(n).aperture;
    TX(n+m+P.numRays).Apod = TX(n).Apod;
    TX(n+m+P.numRays).Steer = TX(n+m).Steer;
    TX(n+m+P.numRays).Delay = TX(n+m).Delay;
end
m = 4*P.numRays;
for n = 1:P.numRays
    TX(n+m).Origin(1) = TX(n).Origin(1);
    TX(n+m).aperture = TX(n).aperture;
    TX(n+m).Apod = TX(n).Apod;
    if n<=10
        TX(n+m).Steer = [((n-1)/10)*P.dtheta,0.0];
    elseif n>(P.numRays-10)
        TX(n+m).Steer = [((P.numRays-n)/10)*P.dtheta,0.0];
    else
        TX(n+m).Steer = [P.dtheta,0.0];
    end
    % Compute transmit delays
    TX(n+m).Delay = computeTXDelays(TX(n+m));
    % Define parameters for inverted transmit
    TX(n+m+P.numRays).waveform = 2;
    TX(n+m+P.numRays).Origin(1) = TxOrgX(n);
    TX(n+m+P.numRays).aperture = TX(n).aperture;
    TX(n+m+P.numRays).Apod = TX(n).Apod;
    TX(n+m+P.numRays).Steer = TX(n+m).Steer;
    TX(n+m+P.numRays).Delay = TX(n+m).Delay;
end

% calculate TXPD
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
TGC.CntrlPts = [0,250,548,710,820,880,940,990];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

BPF = ...
[-0.00058 +0.00024 +0.00070 +0.00208 -0.00366 -0.00076 -0.00290 ...
 +0.01392 -0.00729 +0.00186 -0.02802 +0.03461 -0.00436 +0.03033 ...
 -0.08264 +0.02850 +0.00638 +0.13232 -0.13400 -0.24313 +0.51282];

% Specify Receive structure arrays
%   P.endDepth - add additional acquisition depth to account for some channels
%              having longer path lengths.
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
wlsPer128 = 128/(4*2); % wavelengths in 128 samples for 4 samplesPerWave
Receive = repmat(struct('aperture', 1, ...
                        'Apod', zeros(1,128), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', P.startDepth + wlsPer128*ceil(maxAcqLength/wlsPer128), ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'demodFrequency', demodFreq, ...
                        'InputFilter', BPF, ...
                        'mode', 0, ...
                        'callMediaFunc', 0),1,6*P.numRays*Resource.RcvBuffer(1).numFrames);

% - Set event specific Receive attributes for each frame.
m = P.numRays;
% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    ind = 6*m*(i-1);
    Receive(ind+1).callMediaFunc = 1;
    for j = 1:m
        for k = 1:3 % for three buffers
            Receive(ind+j+(2*k-2)*m).aperture = TX(j).aperture;
            lft = ce(j) - 32;
            if lft < 1, lft = 1; end
            rt = ce(j) + 33;
            if rt > 192, rt = 192; end
            Receive(ind+j+(2*k-2)*m).Apod((lft-(TX(j).aperture-1)):(rt-(TX(j).aperture-1))) = 1;
            Receive(ind+j+(2*k-2)*m).bufnum = k;
            Receive(ind+j+(2*k-2)*m).framenum = i;
            Receive(ind+j+(2*k-2)*m).acqNum = j;

            Receive(ind+j+(2*k-1)*m).aperture = TX(j).aperture;
            Receive(ind+j+(2*k-1)*m).Apod = Receive(ind+j+(2*k-2)*m).Apod;
            Receive(ind+j+(2*k-1)*m).bufnum = k;
            Receive(ind+j+(2*k-1)*m).framenum = i;
            Receive(ind+j+(2*k-1)*m).acqNum = j;
            Receive(ind+j+(2*k-1)*m).mode = 1;
        end
    end
end

% Specify Recon structure arrays.
% - We need three Recon structures, one for each steering direction.
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
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'scaleFactor', 0.5, ...
                   'regionnum', 0), 1, 3*P.numRays);
% - Set specific ReconInfo attributes.
ReconInfo(1).Pre = 'clearInterBuf';
for j = 1:P.numRays
    ReconInfo(j).txnum = j;
    ReconInfo(j).rcvnum = j;
    ReconInfo(j).regionnum = j;
end
ReconInfo(P.numRays).Post = 'IQ2IntensityImageBuf';
k = P.numRays;
ReconInfo(k+1).Pre = 'clearInterBuf';
for j = 1:P.numRays
    ReconInfo(j+k).txnum = 2*P.numRays+j;
    ReconInfo(j+k).rcvnum = 2*P.numRays+j;
    ReconInfo(j+k).regionnum = j+k;
end
ReconInfo(k+P.numRays).Post = 'IQ2IntensityImageBuf';
k = k + P.numRays;
ReconInfo(k+1).Pre = 'clearInterBuf';
for j = 1:P.numRays
    ReconInfo(j+k).txnum = 4*P.numRays+j;
    ReconInfo(j+k).rcvnum = 4*P.numRays+j;
    ReconInfo(j+k).regionnum = j+k;
end
ReconInfo(k+P.numRays).Post = 'IQ2IntensityImageBuf';

% Specify Process structure array.
pers = 20;
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',3.0,...            % pgain is image processing gain
                         'reject',2,...      % reject level
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
SeqControl(2).argument = 300; % 300usec
SeqControl(3).command = 'timeToNextAcq';  % optional time between frames (not needed for synchronous operation)
SeqControl(3).argument = 25000;  % 25000 usec = 25msec time between frames
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
scaleToWvl = 1;
AxesUnit = 'mm';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if ~strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
        AxesUnit = 'wls';
        scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);
    end
end
UI(2).Control = VsSliderControl('LocationCode','UserA1',...
    'Label',['Range (',AxesUnit,')'],'SliderMinMaxVal',[5,P.maxDepthMm,P.endDepthMm]*scaleToWvl,...
    'SliderStep',[0.1,0.2],'ValueFormat','%3.0f',...
    'Callback',@RangeChangeCallback);

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 1;
% Save all the structures to a .mat file.
save('MatFiles/L12-3vWideBeamHISC_64LE');


% **** Callback functions that used by UI controls. ****
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

    PData = evalin('base','PData');
    PData.Size(1) = ceil((P.endDepth-P.startDepth)/PData.PDelta(3)); % P.startDepth, P.endDepth and pdelta set PData.Size(1).
    for i = 1:size(PData.Region,2)
        PData.Region(i).Shape.height = P.endDepth;
    end
    PData.Region = computeRegions(PData);
    assignin('base','PData',PData);

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
    evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData.Size(1)*PData.PDelta(3)/Resource.DisplayWindow(1).pdelta);');
    Receive = evalin('base', 'Receive');
    wlsPer128 = evalin('base','wlsPer128');
    maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
    for i = 1:size(Receive,2)
        Receive(i).endDepth = P.startDepth + wlsPer128*ceil(maxAcqLength/wlsPer128);
    end
    assignin('base','Receive',Receive);

    evalin('base','TGC.rangeMax = P.endDepth;');
    evalin('base','TGC.Waveform = computeTGCWaveform(TGC);');
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'PData','InterBuffer','ImageBuffer','DisplayWindow','TX','Receive','TGC','Recon'};
    assignin('base','Control', Control);
    assignin('base', 'action', 'displayChange');
end
