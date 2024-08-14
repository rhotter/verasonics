% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: SetUpL12_5_50mm256RyLnsMultiFocal.m - Example of 256 ray line imaging
%                                        with multi focal zone.
% Description:
%    Generate .mat file for L12-5 50mm Linear array for Vantage system.
%    Transmit/Receive is performed with a 128 element mux aperture that is
%    moved across the 256 element aperture. The Tx aperture is moved with each
%    ray line, but the mux aperture is moved only for ray lines in the
%    central portion of the transducer aperture where the full 128 channel
%    aperture can be centered around the Tx aperture.
%
% In the following aperture examples, each space represents 2 elements.
% Aperture 1 (numTx=12):
%   tttttt--------------------------------------------------------\-------------------------------------------------------------
%   rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr--------------------------------------------------------------
% Aperture for almost middle beam 128 (TX.aperture=65, numTx=24, numRay=128):
%   ----------------------------------------------------------tttttttttttt------------------------------------------------------
%   --------------------------------rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr----------------------------
% Aperture for last wide beam (TX.aperture=129, numTx=12, numRay=256):
%   -----------------------------------------------------------/----------------------------------------------------------tttttt
%   ------------------------------------------------------------rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
%
%   The receive data from each aperture are stored under different acqNums in the Receive buffer.
%   This version does asynchronous acquisition and processing.
%
% Testing: Tested with software release 2.10.x on Vantage 128.
%
% Last update:
% 12/06/2015 - modified for SW 3.0

clear all

P.numRays = 256;
P.startDepth = 2;   % Acquisition depth in wavelengths
P.endDepth = 192;   % This should preferrably be a multiple of 128 samples.

% Number of Focal Zones
fz = 3;

% Transmit focus in wavelengths
TXFocus(1) = 120;
TXFocus(2) = 95;
TXFocus(3) = 70;

% Receive zone endDepth in wavelengths
RcvZone(1) = 192;
RcvZone(2) = 100;
RcvZone(3) = 70;

% Define system parameters.
Resource.Parameters.numTransmit = 128;  % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;  % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Specify Trans structure array.
Trans.name = 'L12-5 50mm';
Trans.units = 'wavelengths'; % required in Gen3 to prevent default to mm units
% note nominal center frequency in computeTrans is 7.813 MHz
Trans = computeTrans(Trans);  % L12-5 transducer is 'known' transducer so we can use computeTrans.

% Specify PData structure array.
PData(1).PDelta = [Trans.spacing, 0, 0.5];  % x, y, z pdeltas
PData.Size(1) = ceil((P.endDepth-P.startDepth)/PData.PDelta(3)); % startDepth, endDepth and pdelta set PData.Size.
PData.Size(2) = ceil((Trans.numelements*Trans.spacing)/PData.PDelta(1));
PData.Size(3) = 1;      % single image page
PData.Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of upper lft crnr.
% Define numRays rectangular regions centered on TX beam origins.
PData(1).Region = repmat(struct('Shape',struct('Name','Rectangle',...
                                               'Position',[0,0,P.startDepth],...
                                               'width',Trans.spacing,...
                                               'height',P.endDepth-P.startDepth)),1,P.numRays);
firstRayLocX = -((Trans.numelements-1)/2)*Trans.spacing;
for i = 1:P.numRays
    PData.Region(i).Shape.Position(1) = firstRayLocX + (i-1)*Trans.spacing;
end
PData.Region = computeRegions(PData);

% Specify Media object. 'pt1.m' script defines array of point targets.
pt1;
Media.attenuation = -0.5;
Media.function = 'movePoints';

% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = P.numRays*400*8; % this size allows for 192 rays, with range of 384 wvlngths
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 10;        % 10 frames used for RF cineloop.
Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).numFrames = 1;  % one intermediate buffer needed.
Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).numFrames = 10;
Resource.DisplayWindow(1).Title = 'L12-5_50mm256RyLnsMultiFocal';
Resource.DisplayWindow(1).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).numFrames = 20;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = gray(256);

% Specify Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [8.929,0.67,2,1];

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'aperture', 1, ...
                   'Apod', zeros(1,Resource.Parameters.numTransmit), ...
                   'focus', TXFocus(1), ...
                   'Steer', [0.0,0.0], ...
                   'Delay', zeros(1,Resource.Parameters.numTransmit)), 1, fz*P.numRays);

% Determine TX aperture based on focal point and desired f number.
txFNum = 4.5;  % set to desired f-number value for transmit (range: 1.0 - 20)
txNumEl = zeros(1,fz);
for j = 1:fz
    txNumEl(j)=round((TXFocus(j)/txFNum)/Trans.spacing/2); % no. of elements in 1/2 aperture.
    if txNumEl(j) > (Trans.numelements/2 - 1), txNumEl(j) = floor(Trans.numelements/2 - 1); end
end

scaleToWvl = 1;
if strcmp(Trans.units, 'mm')
    scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);
end

% - Set event specific TX attributes.
for j = 1:P.numRays  % specify P.numRays transmit events
    k = fz*(j-1);
    for n = 1:fz
        % Set transmit Origins to positions of elements.
        TX(n+k).Origin = PData.Region(j).Shape.Position;
        % Compute available transmit mux aperture
        [Dummy,ce] = min(abs(scaleToWvl*Trans.ElementPos(:,1)-TX(n+k).Origin(1))); % ce is closest ele to cntr of aper.
        lft = round(ce - 64);
        if lft < 1, lft = 1; end
        if lft > 129, lft = 129; end
        TX(n+k).aperture = lft;
        % Compute TX.Apod within mux aperture.
        lft = round(ce - txNumEl(n));
        if lft < 1, lft = 1; end
        rt = round(ce + txNumEl(n));
        if rt > 256, rt = 256; end
        TX(n+k).Apod((lft-(TX(n+k).aperture-1)):(rt-(TX(n+k).aperture-1))) = 1;
        % Compute transmit delays
        TX(n+k).focus = TXFocus(n);
        TX(n+k).Delay = computeTXDelays(TX(n+k));
    end
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [395,535,650,710,770,830,890,950];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Receive structure arrays
%   endDepth - add additional acquisition depth to account for some channels
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
                        'sampleMode', 'NS200BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 0),1,fz*P.numRays*Resource.RcvBuffer(1).numFrames);

% - Set event specific Receive attributes for each frame.
for i = 1:Resource.RcvBuffer(1).numFrames
    k = P.numRays*fz*(i-1);
    Receive(k+1).callMediaFunc = 1;
    for j = 1:P.numRays
        w = fz*(j-1);
        for z = 1:fz
            %First half of aperture - Set max acq length for first TX in each set of focal zones
            if z == 1
                Receive(k+w+z).endDepth = maxAcqLength;
            else
                Receive(k+w+z).endDepth = RcvZone(z);
            end
            Receive(k+w+z).aperture = TX(w+z).aperture; % mux aperture same as transmit
            Receive(k+w+z).Apod(1:128) = 1.0;
            Receive(k+w+z).framenum = i;
            Receive(k+w+z).acqNum = j;
        end
    end
end

% Specify Recon structure arrays.
% - We need one Recon structure which will be used for each frame.
Recon = struct('senscutoff', 0.5, ...
               'pdatanum', 1, ...
               'rcvBufFrame',-1, ...
               'ImgBufDest', [1,-1], ...
               'RINums',(1:P.numRays)');

% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 0, ...  % default is to accumulate IQ data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 1), 1, P.numRays);
% - Set specific ReconInfo attributes.
for j = 1:P.numRays
    ReconInfo(j).txnum = 1+fz*(j-1);
    ReconInfo(j).rcvnum = 1+fz*(j-1);
    ReconInfo(j).regionnum = j;
end

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
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',40,...
                         'mappingMethod','full',...
                         'display',1,...      % display image after processing
                         'displayWindow',1};

% Specify SeqControl structure arrays.
%  - Jump back to start.
SeqControl(1).command = 'jump';
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between ray line acquisitions
SeqControl(2).argument = 160;  % 160 usec
SeqControl(3).command = 'timeToNextAcq';
SeqControl(3).argument = 80000;  % 15000 usec = 15msec time between frames
SeqControl(4).command = 'returnToMatlab';
nsc = 5; % nsc is count of SeqControl objects

% Specify Event structure arrays.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    k = P.numRays*fz*(i-1);
    for j = 1:P.numRays                 % Acquire all ray lines for frame
        w = fz*(j-1);
        for z = 1:fz
            Event(n).info = 'Acquire ray line';
            Event(n).tx = z+w;
            Event(n).rcv = k+z+w;
            Event(n).recon = 0;
            Event(n).process = 0;
            Event(n).seqControl = 2;
            n = n+1;
        end
    end
    Event(n-1).seqControl = [3,nsc]; % modify last acquisition Event's seqControl
    SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
    nsc = nsc+1;

    Event(n).info = 'recon and process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
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
% - Sensitivity Cutoff
UI(1).Control =  {'UserB7','Style','VsSlider','Label','Sens. Cutoff',...
                  'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                  'SliderStep',[0.025,0.1],'ValueFormat','%1.3f'};
UI(1).Callback = text2cell('%SensCutOffCallback');

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 1;
% Save all the structures to a .mat file.
save('MatFiles/L12-5_50mm_256RyLnsMultiFocal');

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
