% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: SetUpL38_22vWideBeamISC.m - Example of wide beam imaging
%                                        with spatial compounding.
% Description:
%    Transmit/Receive is performed with a 128 element mux aperture that is
%    moved across the 256 element aperture. The Tx aperture is moved with each
%    ray line, but the mux aperture is moved only for ray lines in the
%    central portion of the transducer aperture where the full 128 channel
%    aperture can be centered around the Tx aperture.
%
%    In the following aperture examples, each space represents 4 elements.
%    Aperture 1 (P.numTx=12):
%      ttt-------------------------------------------------------------
%      rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr--------------------------------
%    Aperture for middle beam 128 (TX.aperture=65, P.numTx=24, numRay=128):
%      -----------------------------tttttt-----------------------------
%      ----------------rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr----------------
%    Aperture for wide beam (TX.aperture=129, P.numTx=12, numRay=256):
%      -------------------------------------------------------------ttt
%      --------------------------------rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
%
%   The receive data from each aperture are stored under different acqNums in
%   the Receive buffer.
%
%   This version does asynchronous acquisition and processing.
%
% Last update:
%   11/01/2019

%clear[ 	]+all

P.startDepth = 5;  % Define P.startDepth and P.endDepth at top for use in defining other parameters.
P.endDepth = 200;
P.numRays = 85;
P.numTx = 40;    % number of active transmitters in TX aperture.
P.txFocus = P.endDepth;   % transmit focal pt in wavelengths
P.dtheta = 10*(pi/180);  % angle delta for beams steered left or right

% Define system parameters.
Resource.Parameters.numTransmit = 128;  % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;  % number of receive channels.
Resource.Parameters.speedOfSound = 1540;  % speed of sound in m/s
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Specify Trans structure array.
Trans.name = 'L38-22v';
Trans.units = 'wavelengths'; % required in Gen3 to prevent default to mm units
Trans = computeTrans (Trans);

mm2wl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);

% Specify PData structure array.
PData.PDelta = [Trans.spacing, 0, 0.5];
PData.Size(1) = ceil((P.endDepth-P.startDepth)/PData.PDelta(3)); % P.startDepth, P.endDepth and pdelta set PData.Size.
PData.Size(2) = ceil((Trans.numelements*Trans.spacing)/PData.PDelta(1));
PData.Size(3) = 1;      % single image page
PData.Origin = [-((Trans.numelements-1)/2)*Trans.spacing,0,P.startDepth]; % x,y,z of upper lft crnr.

% Define P.numRays rectangular regions, positioned around the origin of the wide beams.
txdx = 255*Trans.spacing/(P.numRays-1);
TxOrgX = (-127.5*Trans.spacing):txdx:(127.5*Trans.spacing); % x coords of beam centers
for n = 1:P.numRays
    PData.Region(n) = struct('Shape',struct('Name','Trapezoid','Position',[TxOrgX(n),0,0],...
                             'top',P.numTx*Trans.spacing,'bottom',txdx,'height',P.endDepth));
end
m = P.numRays;
% Define P.numRays steered left parallelogram regions, positioned around the origin of the wide beams.
for n = 1:P.numRays
    if n<=10
        angle = -((n-1)/10)*P.dtheta;
    elseif n>(P.numRays-10)
        angle = -((P.numRays-n)/10)*P.dtheta;
    else
        angle = -P.dtheta;
    end
    PData.Region(n+m) = struct('Shape',struct('Name','Trapezoid','Position',[TxOrgX(n),0,0],...
                               'top',P.numTx*Trans.spacing,'bottom',txdx,'height',P.endDepth/cos(P.dtheta),'steer',angle));
end
m = m + P.numRays;
% Define P.numRays steered right parallelogram regions, positioned around the origin of the wide beams.
for n = 1:P.numRays
    if n<=10
        angle = ((n-1)/10)*P.dtheta;
    elseif n>(P.numRays-10)
        angle = ((P.numRays-n)/10)*P.dtheta;
    else
        angle = P.dtheta;
    end
    PData.Region(n+m) = struct('Shape',struct('Name','Trapezoid','Position',[TxOrgX(n),0,0],...
                               'top',P.numTx*Trans.spacing,'bottom',txdx,'height',P.endDepth/cos(P.dtheta),'steer',angle));
end
PData.Region = computeRegions(PData);

% Specify Media object. 'pt1.m' script defines array of point targets.
pt1;
Media.attenuation = -0.5;
Media.function = 'movePoints';

% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 2*P.numRays*400*8; % this size allows for all rays, with range of up to 400 wvlngths
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = numRcvFrames;
Resource.RcvBuffer(2).datatype = 'int16';
Resource.RcvBuffer(2).rowsPerFrame = 2*P.numRays*400*8; % this size allows for all rays, with range of up to 400 wvlngths
Resource.RcvBuffer(2).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(2).numFrames = numRcvFrames;  % all RcvBuffers should have same no. of frames
Resource.RcvBuffer(3).datatype = 'int16';
Resource.RcvBuffer(3).rowsPerFrame = 2*P.numRays*400*8; % this size allows for all rays, with range of up to 400 wvlngths
Resource.RcvBuffer(3).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(3).numFrames = numRcvFrames;
Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).numFrames = 1;  % one intermediate buffer needed.
Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).numFrames = 10;
Resource.DisplayWindow(1).Title = 'L38-22vWideBeamISC';
Resource.DisplayWindow(1).pdelta = 0.4;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).numFrames = numImageFrames;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = gray(256);

% Specify Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,0.67,2,1];

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'aperture', 1, ...
                   'Apod', zeros(1,Resource.Parameters.numTransmit), ...
                   'focus', P.txFocus, ...
                   'Steer', [0.0,0.0], ...
                   'Delay', zeros(1,Resource.Parameters.numTransmit), ...
                   'TXPD', [], ...
                   'peakCutOff', 1.0, ...
                   'peakBLMax', 10.0), 1, 6*P.numRays);

scaleToWvl = 1;
if strcmp(Trans.units, 'mm')
    scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);
end

% - Set event specific TX attributes.
for n = 1:P.numRays  % specify P.numRays transmit events
    % Set transmit Origins.
    TX(n).Origin(1) = TxOrgX(n);
    % Compute available transmit mux aperture
    [Dummy,ce] = min(abs(scaleToWvl*Trans.ElementPos(:,1)-TxOrgX(n))); % ce is closest element to center of aperture.
    lft = round(ce - 64);
    if lft < 1, lft = 1; end
    if lft > 129, lft = 129; end
    TX(n).aperture = lft;
    % Compute TX.Apod within mux aperture.
    lft = round(ce - P.numTx/2);
    if lft < 1, lft = 1; end
    rt = round(ce + P.numTx/2);
    if rt > 256, rt = 256; end
    TX(n).Apod((lft-(TX(n).aperture-1)):(rt-(TX(n).aperture-1))) = 1;
    % Apply apodization function.
    [RIndices,CIndices,V] = find(TX(n).Apod);
    V = kaiser(size(V,2),1);
    TX(n).Apod(CIndices) = V;
    % Define 2nd TXs for interleave.
    TX(n+P.numRays) = TX(n);
    % Compute transmit delays for normal and interleave
    TX(n+P.numRays).Delay = computeTXDelays(TX(n+P.numRays));
    % Add the delay offset representing one half of A/D sample period (1/62.5MHz = 16ns) to first acq's delay.
    TX(n).Delay = TX(n+P.numRays).Delay + 0.5*Trans.frequency/62.5; % 0.24 wls of Trans.frequency=30 equivalent to 8 nsec.
    % Compute transmit pixel data only for the first TX in the interleave pair.
    TX(n).TXPD = computeTXPD(TX(n),PData);
    TX(n+P.numRays).TXPD = TX(n).TXPD;
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
    % Define 2nd TXs for interleave.
    TX(n+m+P.numRays) = TX(n+m);
    % Compute transmit delays for normal and interleave
    TX(n+m+P.numRays).Delay = computeTXDelays(TX(n+m+P.numRays));
    % Add the delay offset representing one half of A/D sample period (1/62.5MHz = 16ns) to first acq's delay.
    TX(n+m).Delay = TX(n+m+P.numRays).Delay + 0.5*Trans.frequency/62.5; % 0.24 wls of Trans.frequency=30 equivalent to 8 nsec.
    % Compute transmit pixel data
    TX(n+m).TXPD = computeTXPD(TX(n+m),PData);
    TX(n+m+P.numRays).TXPD = TX(n+m).TXPD;
end
m = m + 2*P.numRays;
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
    % Define 2nd TXs for interleave.
    TX(n+m+P.numRays) = TX(n+m);
    % Compute transmit delays for normal and interleave
    TX(n+m+P.numRays).Delay = computeTXDelays(TX(n+m+P.numRays));
    % Add the delay offset representing one half of A/D sample period (1/62.5MHz = 16ns) to first acq's delay.
    TX(n+m).Delay = TX(n+m+P.numRays).Delay + 0.5*Trans.frequency/62.5; % 0.24 wls of Trans.frequency=30 equivalent to 8 nsec.
    % Compute transmit pixel data
    TX(n+m).TXPD = computeTXPD(TX(n+m),PData);
    TX(n+m+P.numRays).TXPD = TX(n+m).TXPD;
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [150,468,644,788,925,965,1000,1023];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Receive structure arrays
RcvProfile.LnaZinSel = 31;
%   P.endDepth - add additional acquisition depth to account for some channels
%              having longer path lengths.
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
Receive = repmat(struct('aperture', 1, ...
                        'Apod', zeros(1,128), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BWI', ...
                        'mode', 0, ...
                        'callMediaFunc', 0),1,6*P.numRays*Resource.RcvBuffer(1).numFrames);

% - Set event specific Receive attributes for each frame.
for i = 1:Resource.RcvBuffer(1).numFrames
    k = 6*P.numRays*(i-1);
    Receive(k+1).callMediaFunc = 1;
    for j = 1:P.numRays
        Receive(k+2*j-1).aperture = TX(j).aperture; % mux aperture same as transmit
        Receive(k+2*j-1).Apod(1:128) = 1.0;
        Receive(k+2*j-1).framenum = i;
        Receive(k+2*j-1).acqNum = 2*j-1;
        Receive(k+2*j).aperture = TX(P.numRays+j).aperture; % mux aperture same as transmit
        Receive(k+2*j).Apod(1:128) = 1.0;
        Receive(k+2*j).framenum = i;
        Receive(k+2*j).acqNum = 2*j;
        % Receives for 2nd angle
        m = k + 2*P.numRays;
        Receive(m+2*j-1).aperture = TX(2*P.numRays+j).aperture; % mux aperture same as transmit
        Receive(m+2*j-1).Apod(1:128) = 1.0;
        Receive(m+2*j-1).bufnum = 2;
        Receive(m+2*j-1).framenum = i;
        Receive(m+2*j-1).acqNum = 2*j-1;
        Receive(m+2*j).aperture = TX(3*P.numRays+j).aperture; % mux aperture same as transmit
        Receive(m+2*j).Apod(1:128) = 1.0;
        Receive(m+2*j).bufnum = 2;
        Receive(m+2*j).framenum = i;
        Receive(m+2*j).acqNum = 2*j;
        % Receives for 3rd angle
        m = m + 2*P.numRays;
        Receive(m+2*j-1).aperture = TX(4*P.numRays+j).aperture; % mux aperture same as transmit
        Receive(m+2*j-1).Apod(1:128) = 1.0;
        Receive(m+2*j-1).bufnum = 3;
        Receive(m+2*j-1).framenum = i;
        Receive(m+2*j-1).acqNum = 2*j-1;
        Receive(m+2*j).aperture = TX(5*P.numRays+j).aperture; % mux aperture same as transmit
        Receive(m+2*j).Apod(1:128) = 1.0;
        Receive(m+2*j).bufnum = 3;
        Receive(m+2*j).framenum = i;
        Receive(m+2*j).acqNum = 2*j;
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
                   'scaleFactor', 1.0, ...
                   'regionnum', 0), 1, 3*P.numRays);
% - Set specific ReconInfo attributes.
ReconInfo(1).Pre = 'clearInterBuf';
for j = 1:P.numRays
    ReconInfo(j).txnum = j;
    ReconInfo(j).rcvnum = 2*j-1;
    ReconInfo(j).regionnum = j;
end
ReconInfo(P.numRays).Post = 'IQ2IntensityImageBuf';
k = P.numRays;
ReconInfo(k+1).Pre = 'clearInterBuf';
for j = 1:P.numRays
    ReconInfo(j+k).txnum = 2*P.numRays+j;
    ReconInfo(j+k).rcvnum = 2*P.numRays+2*j-1;
    ReconInfo(j+k).regionnum = P.numRays+j;
end
ReconInfo(k+P.numRays).Post = 'IQ2IntensityImageBuf';
k = k + P.numRays;
ReconInfo(k+1).Pre = 'clearInterBuf';
for j = 1:P.numRays
    ReconInfo(j+k).txnum = 4*P.numRays+j;
    ReconInfo(j+k).rcvnum = 4*P.numRays+2*j-1;
    ReconInfo(j+k).regionnum = 2*P.numRays+j;
end
ReconInfo(k+P.numRays).Post = 'IQ2IntensityImageBuf';

% Specify Process structure array.
pers = 20;
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',1.0,...            % pgain is image processing gain
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
%  - Jump back to start.
SeqControl(1).command = 'jump';
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(2).argument = 50;  % 50 usec
SeqControl(3).command = 'timeToNextAcq';  % time to next wide beam acquisition
SeqControl(3).argument = 100;  % 100 usec
SeqControl(4).command = 'timeToNextAcq';
SeqControl(4).argument = 20000; %20000 usec = 20msec time between frames
SeqControl(5).command = 'returnToMatlab';
nsc = 6; % nsc is count of SeqControl objects

% Specify Event structure arrays.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    % Acquire frame with straight ahead wide beams
    k = 6*P.numRays*(i-1);
    for j = 1:P.numRays        % Acquire frame
        Event(n).info = 'acquire aperture.';
        Event(n).tx = j;
        Event(n).rcv = k+2*j-1;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;
        Event(n).info = 'acquire interleave aperture.';
        Event(n).tx = j+P.numRays;
        Event(n).rcv = k+2*j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 3;
        n = n+1;
    end
    Event(n-1).seqControl = [4,nsc]; % modify last acquisition Event's seqControl
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
    for j = 1:P.numRays        % Acquire frame
        Event(n).info = 'acquire aperture.';
        Event(n).tx = j+2*P.numRays;
        Event(n).rcv = k+2*P.numRays+2*j-1;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;
        Event(n).info = 'acquire interleave aperture.';
        Event(n).tx = j+3*P.numRays;
        Event(n).rcv = k+2*P.numRays+2*j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 3;
        n = n+1;
    end
    Event(n-1).seqControl = [4,nsc]; % modify last acquisition Event's seqControl
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
    for j = 1:P.numRays        % Acquire frame
        Event(n).info = 'acquire aperture.';
        Event(n).tx = j+4*P.numRays;
        Event(n).rcv = k+4*P.numRays+2*j-1;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;
        Event(n).info = 'acquire interleave aperture.';
        Event(n).tx = j+5*P.numRays;
        Event(n).rcv = k+4*P.numRays+2*j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 3;
        n = n+1;
    end
    Event(n-1).seqControl = [4,nsc]; % modify last acquisition Event's seqControl
      SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
      nsc = nsc+1;

    Event(n).info = 'recon and process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 3;
    Event(n).process = 1;
    Event(n).seqControl = 0;
    if i ~= Resource.RcvBuffer(1).numFrames  % Exit to Matlab every 3rd frame
        Event(n).seqControl = 5;
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
% - Sensitivity Cutoff
UI(1).Control =  {'UserB7','Style','VsSlider','Label','Sens. Cutoff',...
                  'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                  'SliderStep',[0.025,0.1],'ValueFormat','%1.3f'};
UI(1).Callback = text2cell('%SensCutOffCallback');
% - Range Change
wls2mm = 1;
AxesUnit = 'wls';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
        AxesUnit = 'mm';
        wls2mm = Resource.Parameters.speedOfSound/1000/Trans.frequency;
    end
end
UI(2).Control = {'UserA1','Style','VsSlider','Label',['Range (',AxesUnit,')'],...
                 'SliderMinMaxVal',[64,200,P.endDepth]*wls2mm,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f'};
UI(2).Callback = text2cell('%RangeChangeCallback');

% - TX Aperture change
UI(3).Control = {'UserB3','Style','VsSlider','Label','TX Aper',...
                 'SliderMinMaxVal',[1,128,P.numTx],'SliderStep',[1/128,1/64],'ValueFormat','%3.0f'};
UI(3).Callback = text2cell('%ApertureCallback');

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 3;
% Save all the structures to a .mat file.
filename = fullfile( vsv.file.getVSXDir, 'MatFiles/SetUpL38_22vWideBeamISC_QSApp.mat');

save(filename);

% Automatic VSX Execution:
% Uncomment the following line to automatically run VSX every time you run
% this SetUp script (note that if VSX finds the variable 'filename' in the
% Matlab workspace, it will load and run that file instead of prompting the
% user for the file to be used):

% filename = 'L38-22vWideBeamISC';  VSX;

return

% **** Callback routines to be converted by text2cell function. ****
%SensCutOffCallback - Sensitivity cutoff change
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
%SensCutOffCallback

%RangeChangeCallback - Range change
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
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
        P.endDepth = UIValue*scaleToWvl;
    end
end
P.txFocus = P.endDepth;
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
for n = 1:P.numRays
    TX(n).focus = P.txFocus;
    TX(n+P.numRays).focus = P.txFocus;
    TX(n+2*P.numRays).focus = P.txFocus;
    TX(n+3*P.numRays).focus = P.txFocus;
    TX(n+4*P.numRays).focus = P.txFocus;
    TX(n+5*P.numRays).focus = P.txFocus;
    % Compute transmit delays for normal and interleave
    TX(n+P.numRays).Delay = computeTXDelays(TX(n+P.numRays));
    TX(n+3*P.numRays).Delay = computeTXDelays(TX(n+3*P.numRays));
    TX(n+5*P.numRays).Delay = computeTXDelays(TX(n+5*P.numRays));
    % Add the delay offset representing one half of A/D sample period (1/62.5MHz = 16ns) to first acq's delay.
    TX(n).Delay = TX(n+P.numRays).Delay + 0.5*Trans.frequency/62.5; % 0.24 wls of Trans.frequency=30 equivalent to 8 nsec.
    TX(n+2*P.numRays).Delay = TX(n+3*P.numRays).Delay + 0.5*Trans.frequency/62.5; % 0.24 wls of Trans.frequency=30 equivalent to 8 nsec.
    TX(n+4*P.numRays).Delay = TX(n+5*P.numRays).Delay + 0.5*Trans.frequency/62.5; % 0.24 wls of Trans.frequency=30 equivalent to 8 nsec.
    TX(n).TXPD = computeTXPD(TX(n),PData);
    TX(n+2*P.numRays).TXPD = computeTXPD(TX(n+2*P.numRays),PData);
    TX(n+4*P.numRays).TXPD = computeTXPD(TX(n+4*P.numRays),PData);
    waitbar(n/P.numRays)
end
close(h)
assignin('base','TX',TX);
evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData.Size(1)*PData.PDelta(3)/Resource.DisplayWindow(1).pdelta);');
Receive = evalin('base', 'Receive');
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
for i = 1:size(Receive,2)
    Receive(i).endDepth = maxAcqLength;
end
assignin('base','Receive',Receive);
evalin('base','TGC.rangeMax = P.endDepth;');
evalin('base','TGC.Waveform = computeTGCWaveform(TGC);');
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'PData','InterBuffer','ImageBuffer','DisplayWindow','TX','Receive','TGC','Recon'};
assignin('base','Control', Control);
assignin('base', 'action', 'displayChange');
return
%RangeChangeCallback

%ApertureCallback
simMode = evalin('base','Resource.Parameters.simulateMode');
% No focus change if in simulate mode 2.
if simMode == 2
    set(hObject,'Value',evalin('base','P.numTx'));
    return
end
P = evalin('base','P');
P.numTx = UIValue;
assignin('base','P',P);

Trans = evalin('base', 'Trans');
scaleToWvl = evalin('base','scaleToWvl');
PData = evalin('base','PData');
for n = 1:size(PData.Region,2)
    PData.Region(n).Shape.top = P.numTx*Trans.spacing;
end
PData.Region = computeRegions(PData);
assignin('base','PData',PData);

% - Redefine event specific TX attributes for the new aperture.
TX = evalin('base', 'TX');
TxOrgX = evalin('base','TxOrgX');
h = waitbar(0,'Program TX parameters, please wait!');
steps = 3*P.numRays;
for n = 1:P.numRays  % specify P.numRays transmit events
    TX(n).Apod = zeros(1,128);
    TX(n).Origin(1) = TxOrgX(n);
    % Compute available transmit mux aperture
    [Dummy,ce] = min(abs(scaleToWvl*Trans.ElementPos(:,1)-TxOrgX(n))); % ce is closest ele to cntr of aper.
    lft = round(ce - 64);
    if lft < 1, lft = 1; end
    if lft > 129, lft = 129; end
    TX(n).aperture = lft;
    % Compute TX.Apod within mux aperture.
    lft = round(ce - P.numTx/2);
    if lft < 1, lft = 1; end
    rt = round(ce + P.numTx/2);
    if rt > 256, rt = 256; end
    TX(n).Apod((lft-(TX(n).aperture-1)):(rt-(TX(n).aperture-1))) = 1;
    % Apply apodization function.
    [RIndices,CIndices,V] = find(TX(n).Apod);
    V = kaiser(size(V,2),1);
    TX(n+P.numRays) = TX(n);
    % Compute transmit delays for normal and interleave
    TX(n+P.numRays).Delay = computeTXDelays(TX(n+P.numRays));
    % Add the delay offset representing one half of A/D sample period (1/62.5MHz = 16ns) to first acq's delay.
    TX(n).Delay = TX(n+P.numRays).Delay + 0.5*Trans.frequency/62.5; % 0.24 wls of Trans.frequency=30 equivalent to 8 nsec.
    % Compute transmit pixel data only for the first TX in the interleave pair.
    TX(n).TXPD = computeTXPD(TX(n),PData);
    waitbar(n/steps);
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
    % Define 2nd TXs for interleave.
    TX(n+m+P.numRays) = TX(n+m);
    % Compute transmit delays for normal and interleave
    TX(n+m+P.numRays).Delay = computeTXDelays(TX(n+m+P.numRays));
    % Add the delay offset representing one half of A/D sample period (1/62.5MHz = 16ns) to first acq's delay.
    TX(n+m).Delay = TX(n+m+P.numRays).Delay + 0.5*Trans.frequency/62.5; % 0.24 wls of Trans.frequency=30 equivalent to 8 nsec.
    % Compute transmit pixel data
    TX(n+m).TXPD = computeTXPD(TX(n+m),PData);
    TX(n+m+P.numRays).TXPD = TX(n+m).TXPD;
    waitbar((n+P.numRays)/steps);
end
m = m + 2*P.numRays;
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
    % Define 2nd TXs for interleave.
    TX(n+m+P.numRays) = TX(n+m);
    % Compute transmit delays for normal and interleave
    TX(n+m+P.numRays).Delay = computeTXDelays(TX(n+m+P.numRays));
    % Add the delay offset representing one half of A/D sample period (1/62.5MHz = 16ns) to first acq's delay.
    TX(n+m).Delay = TX(n+m+P.numRays).Delay + 0.5*Trans.frequency/62.5; % 0.24 wls of Trans.frequency=30 equivalent to 8 nsec.
    % Compute transmit pixel data
    TX(n+m).TXPD = computeTXPD(TX(n+m),PData);
    TX(n+m+P.numRays).TXPD = TX(n+m).TXPD;
    waitbar((n+2*P.numRays)/steps);
end
close(h)
assignin('base','TX', TX);
% Set Control command to update TX
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'PData','TX','Recon'};
assignin('base','Control', Control);
return
%ApertureCallback
