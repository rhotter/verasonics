% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: SetUpL22_14vLF_FlashAngles_32LE.m - Example of plane wave imaging with multiple steered angle transmits
%
% Description:
%   Sequence programming file for L22-14vLF Linear array, using plane wave
%   transmits with multiple steering angles and 4-1 synthetic aperture
%   (three acquisitions per line to acquire data from all 128 receive
%   elements). 64 transmit channels and 32 receive channels
%   are active and positioned as follows (each char represents 2 elements)
%   for each of the 4 synthetic apertures. This script uses 4X sampling
%   with A/D sample rate of 62.5 MHz for a 15.625 MHz processing center
%   frequency.  Transmit is at 17.8 MHz and receive bandpass filter has
%   been shifted to 18 MHz center frequency, 13.9MHz -3 dB bandwidth to
%   support the 12 MHz bandwidth of the L22-14vLF (18MHz center frequency,
%   67% bandwidth). Processing is asynchronous with respect to acquisition.
%
%
%   Element Nos.                                                               1
%                               3   4           6         8    9               2
%               1               3   2           5         6    6               8
%   Aperture 1: |               |   |           |         |    |               |
%               tttttttttttttttttttttttttttttttt--------------------------------
%               rrrrrrrrrrrrrrrr------------------------------------------------
%               |               |    |          |         |    |               |
%   Aperture 2: |               |    |          |         |    |               |
%               --------tttttttttttttttttttttttttttttttt------------------------
%               ----------------rrrrrrrrrrrrrrrr--------------------------------
%               |               |    |          |         |    |               |
%   Aperture 3: |               |    |          |         |    |               |
%               ------------------------ttttttttttttttttttttttttttttttt---------
%               --------------------------------rrrrrrrrrrrrrrr-----------------
%               |               |    |          |         |    |               |
%   Aperture 4: |               |    |          |         |    |               |
%               --------------------------------tttttttttttttttttttttttttttttttt
%               ------------------------------------------------rrrrrrrrrrrrrrrr
%               |               |    |          |         |    |               |
%
%   The receive data from each of these apertures are stored under
%   different acqNums in the Receive buffer. The reconstruction sums the
%   IQ data from the 4 aquisitions and multiple angles and computes intensity values
%   to produce the full frame. Processing is asynchronous with respect to acquisition.
%
% Last update:
% 09/03/2018 - modified for 32LE system

clear all
P.startDepth = 5;   % Acquisition depth in wavelengths
P.endDepth = 192;   % This should preferrably be a multiple of 128 samples.

RcvProfile.LnaZinSel = 31;

na = 7;      % Set na = number of angles.
if (na > 1)
    dtheta = (36*pi/180)/(na-1);
    P.startAngle = -36*pi/180/2;
else
    dtheta = 0;
    P.startAngle=0;
end % set dtheta to range over +/- 18 degrees.

% Define system parameters.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.
Resource.System.SoftwareVersion = [3 5 0]; % Minimum software release for this script.
Resource.System.UTA = '260-MUX'; % This script requires the 260-MUX UTA.
Resource.System.transmitChannels = 64;
Resource.System.receiveChannels = 32;

% Specify Trans structure array.
Trans.name = 'L22-14v LF';
Trans.units = 'wavelengths';
Trans = computeTrans(Trans);
Trans.maxHighVoltage = 25; % mfr data sheet lists 30 Volt limit
Trans = computeUTAMux64(Trans); % Add HVMux field for use with UTA 260-Mux

% Specify PData structure array.
PData(1).PDelta = [0.4, 0, 0.25];
PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3)); % startDepth, endDepth and pdelta set PData(1).Size.
PData(1).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(3) = 1;      % single image page
PData(1).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of upper lft crnr.
% No PData.Region specified, so a default Region for the entire PData array will be created by computeRegions.

% Specify Media object. 'pt1.m' script defines array of point targets.
pt1;
Media.attenuation = -0.5;
Media.function = 'movePoints';

% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 4*na*4096; % this size allows for maximum range
Resource.RcvBuffer(1).colsPerFrame = Resource.System.receiveChannels;
Resource.RcvBuffer(1).numFrames = 30;    % 30 frames stored in RcvBuffer.
Resource.InterBuffer(1).numFrames = 1;   % one intermediate buffer needed.
Resource.ImageBuffer(1).numFrames = 10;
Resource.DisplayWindow(1).Title = 'L22-14vLF_FlashAngles_32LE';
Resource.DisplayWindow(1).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = gray(256);

% Specify Transmit waveform structure.
TW.type = 'parametric';
TW.Parameters = [18, 0.67, 3, 1];   % 18 MHz center frequency, 67% pulsewidth 1.5 cycle burst

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'aperture',[],...
                   'Apod', ones(1,Resource.System.transmitChannels), ...
                   'focus', 0.0, ...
                   'Steer', [0.0,0.0], ...
                   'Delay', zeros(1,Resource.System.transmitChannels)), 1, 4*na);
% - Set event specific TX attributes.
if fix(na/2) == na/2       % if na even
    P.startAngle = (-(fix(na/2) - 1) - 0.5)*dtheta;
else
    P.startAngle = -fix(na/2)*dtheta;
end
angle = P.startAngle;
for j = 1:4:4*na
    TX(j).aperture = 1; % use computeMuxAperture function to get aperture number
    TX(j).Steer = [angle,0.0];
    TX(j).Delay = computeTXDelays(TX(j));
    TX(j+1).aperture = 17;
    TX(j+1).Steer = [angle,0.0];
    TX(j+1).Delay = computeTXDelays(TX(j+1));
    TX(j+2).aperture = 49;
    TX(j+2).Steer = [angle,0.0];
    TX(j+2).Delay = computeTXDelays(TX(j+2));
    TX(j+3).aperture = 65; % Use the tx aperture that starts at element 65.
    TX(j+3).Steer = [angle,0.0];
    TX(j+3).Delay = computeTXDelays(TX(j+3));
    angle = angle + dtheta;
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [330 560 780 1010 1023 1023 1023 1023]; % [0,511,716,920,1023,1023,1023,1023];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Receive structure arrays.
% - We need 4*na Receives for every frame.

% sampling center frequency is 15.625, but we want the bandpass filter
% centered on the actual transducer center frequency of 18 MHz with 67%
% bandwidth, or 12 to 24 MHz.  Coefficients below were set using
% "filterTool" with normalized cf=1.15 (18 MHz), bw=0.85,
% xsn wdth=0.41 resulting in -3 dB 0.71 to 1.6 (11.1 to 25 MHz), and
% -20 dB 0.57 to 1.74 (8.9 to 27.2 MHz)
%
BPF1 = [ -0.00009 -0.00128 +0.00104 +0.00085 +0.00159 +0.00244 -0.00955 ...
         +0.00079 -0.00476 +0.01108 +0.02103 -0.01892 +0.00281 -0.05206 ...
         +0.01358 +0.06165 +0.00735 +0.09698 -0.27612 -0.10144 +0.48608 ];

maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
Receive = repmat(struct('Apod', zeros(1,Resource.System.transmitChannels), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength,...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'InputFilter', BPF1, ...
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, 4*na*Resource.RcvBuffer(1).numFrames);

% - Set event specific Receive attributes for each frame.
for i = 1:Resource.RcvBuffer(1).numFrames
    k = 4*na*(i-1);
    Receive(k+1).callMediaFunc = 1;
    for j = 1:4:4*na
        % -- 1st synthetic aperture acquisition for acq 1.
        Receive(k+j).Apod(1:32) = 1.0;
        Receive(k+j).aperture = 1;
        Receive(k+j).framenum = i;
        Receive(k+j).acqNum = j;
        % -- 2nd synthetic aperture acquisition for acq 2.
        Receive(k+j+1).Apod(33:64) = 1.0;
        Receive(k+j+1).aperture = 1;
        Receive(k+j+1).framenum = i;
        Receive(k+j+1).acqNum = j+1;
        % -- 3rd synthetic aperture acquisition for acq 3.
        Receive(k+j+2).Apod(1:32) = 1.0;
        Receive(k+j+2).aperture = 65;
        Receive(k+j+2).framenum = i;
        Receive(k+j+2).acqNum = j+2;
        % -- 4rd synthetic aperture acquisition for acq 4.
        Receive(k+j+3).Apod(33:64) = 1.0;
        Receive(k+j+3).aperture = 65;
        Receive(k+j+3).framenum = i;
        Receive(k+j+3).acqNum = j+3;
    end
end

% Specify Recon structure arrays.
% - We need one Recon structure which will be used for each frame.
Recon = struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame',-1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...
               'RINums', 1:4*na);

% Define ReconInfo structures.
% We need na ReconInfo structures for na steering angles.
ReconInfo = repmat(struct('mode', 'accumIQ', ...  % default is to accumulate IQ data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 1), 1, 4*na);
% - Set specific ReconInfo attributes.
ReconInfo(1).mode = 'replaceIQ'; % replace IQ data.
    for j = 1:4:4*na  % For each row in the column
        ReconInfo(j).txnum = j;
        ReconInfo(j).rcvnum = j;
        ReconInfo(j+1).txnum = j+1;
        ReconInfo(j+1).rcvnum = j+1;
        ReconInfo(j+2).txnum = j+2;
        ReconInfo(j+2).rcvnum = j+2;
        ReconInfo(j+3).txnum = j+3;
        ReconInfo(j+3).rcvnum = j+3;
    end
ReconInfo(4*na).mode = 'accumIQ_replaceIntensity'; % accum and detect

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
SeqControl(1).command = 'jump'; % jump back to start
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(2).argument = 160;  % 160 usec
SeqControl(3).command = 'timeToNextAcq';  % time between frames
SeqControl(3).argument = 20000 - (na-1)*160;  % 20 msec
SeqControl(4).command = 'returnToMatlab';
nsc = 5; % nsc is count of SeqControl objects

% Specify Event structure arrays.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    k = 4*na*(i-1);
    for j = 1:4:4*na                      % Acquire frame
        Event(n).info = '1st aperture.';
        Event(n).tx = j;
        Event(n).rcv = k+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;

        Event(n).info = '2nd aperture.';
        Event(n).tx = j+1;
        Event(n).rcv = k+j+1;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;

        Event(n).info = '3rd aperture.';
        Event(n).tx = j+2;
        Event(n).rcv = k+j+2;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n=n+1;

        Event(n).info = '4th aperture.';
        Event(n).tx = j+2;
        Event(n).rcv = k+j+3;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n=n+1;
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
    if floor(i/5) == i/5     % Exit to Matlab every 5th frame
        Event(n).seqControl = 4;
    else
        Event(n).seqControl = 0;
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
UI(1).Callback = text2cell('%SensCutoffCallback');

% - Range Change
MinMaxVal = [64,300,P.endDepth]; % default unit is wavelength
AxesUnit = 'wls';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
        AxesUnit = 'mm';
        MinMaxVal = MinMaxVal * (Resource.Parameters.speedOfSound/1000/Trans.frequency);
    end
end
UI(2).Control = {'UserA1','Style','VsSlider','Label',['Range (',AxesUnit,')'],...
                 'SliderMinMaxVal',MinMaxVal,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f'};
UI(2).Callback = text2cell('%RangeChangeCallback');

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 5;

% Save all the structures to a .mat file.
save('MatFiles/L22-14vLF_FlashAngles_32LE');
return

% **** Callback routines to be converted by text2cell function. ****
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
return
%SensCutoffCallback

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
assignin('base','P',P);

evalin('base','PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));');
evalin('base','PData(1).Region = computeRegions(PData(1));');
evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);');
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
Control.Parameters = {'PData','InterBuffer','ImageBuffer','DisplayWindow','Receive','TGC','Recon'};
assignin('base','Control', Control);
assignin('base', 'action', 'displayChange');
return
%RangeChangeCallback
