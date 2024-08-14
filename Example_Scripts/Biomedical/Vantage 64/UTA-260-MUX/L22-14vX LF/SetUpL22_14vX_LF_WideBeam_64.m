% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: SetUpL22_14vX_LF_WideBeam_64.m - Example of imaging with partially focused "wide beam" transmits
%
% Description:
%   Sequence programming file for L22-14vX LF Linear array, using numrays widebeam
%   transmit and receive acquisitions. The transmit aperture size is set
%   by the variable - P.numTx, and the receive aperture always covers the full 64 element aperture.
%   This script uses 4X sampling with A/D sample rate of
%   62.5 MHz for a 15.625 MHz processing center frequency.  Transmit is at
%   17.8 MHz and receive bandpass filter has been shifted to 18 MHz center
%   frequency, 13.9MHz -3 dB bandwidth to support the 12 MHz bandwidth of
%   the L22-14vX LF (18MHz center frequency, 67% bandwidth). Processing is
%   asynchronous with respect to acquisition.
%
%   The 260-MUX UTA connects two elements with each system channel through 
%   a 2-1 multiplexer, allowing sweeping sub-apertures of up to 64 elements 
%   across the 128-element transducer.
%
% Last update:
%   04/03/2018 - modified for SW 3.0
%   08/22/2018 - scirpt improvement
%   04/12/2021 - updated for 4.5 release


clear all

% key parameters
P.numTx = 48;     % Number of Transmit Elements in the aperture.
P.numRays = 48;   % Number of wide beams in a frame.
P.startDepth = 5;
P.endDepth = 192;

% Define system parameters.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.
Resource.System.Product = 'Vantage64';
Resource.System.SoftwareVersion = [3 5 0]; % Minimum software release for this script.
Resource.System.UTA = '260-MUX'; % This script requires the 260-MUX UTA.
Resource.System.transmitChannels = 64;
Resource.System.receiveChannels = 64;

% Specify Trans structure array.
Trans.name = 'L22-14vX LF';
Trans.frequency = 62.5/4;
Trans.units = 'wavelengths'; % required in Gen3 to prevent default to mm units
Trans = computeTrans(Trans);  % L22-14vX LF transducer is 'known' transducer so we can use computeTrans.
Trans = computeUTAMux64(Trans); % Add HVMux field for use with UTA 260-Mux
Trans.elBias =  0; % To disable elBias due to the UTA 260 MUX used.
Trans.maxHighVoltage = 25;

% - nominal center frequency from computeTrans is 15.625 MHz

% Specify PData structure arrays.
% - PData contains rectangular regions for no beam steering.
PData.PDelta = [0.4*Trans.spacing, 0, 0.5];
PData.Size(1) = ceil((P.endDepth-P.startDepth)/PData.PDelta(3)); % P.startDepth, P.endDepth and pdelta set PData.Size.
PData.Size(2) = ceil((Trans.numelements*Trans.spacing)/PData.PDelta(1));
PData.Size(3) = 1;      % single image page
PData.Origin = [-Trans.spacing*(Trans.numelements/2),0,P.startDepth]; % x,y,z of upper lft crnr.
% Define P.numRays rectangular regions, positioned around the origin of the wide beams.
TxOrgX = (-63.5*Trans.spacing):(127*Trans.spacing/(P.numRays-1)):(63.5*Trans.spacing); % x coords of beam centers
for n = 1:P.numRays
    PData.Region(n).Shape = struct('Name','Rectangle','Position',[TxOrgX(n),0,0],'width',32,'height',P.endDepth);
end

% Specify Media object.
pt1;
Media.attenuation = -0.5;
Media.function = 'movePoints';

% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = P.numRays*4096; % this size allows for maximum range
Resource.RcvBuffer(1).colsPerFrame = Resource.System.receiveChannels;
Resource.RcvBuffer(1).numFrames = 10;    % 10 frames stored in RcvBuffer.
Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).numFrames = 1;   % one intermediate buffer needed.
Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).numFrames = 10;
Resource.DisplayWindow(1).Title = 'L22-14vX LF Wide Beam, 4X sampling at 62.5 MHz';
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
TW.type = 'parametric';
TW.Parameters = [18, 0.67, 4, 1];

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'Apod', zeros(1,Resource.System.transmitChannels), ...
                   'focus', 3*P.endDepth, ...
                   'Steer', [0.0,0.0], ...
                   'Delay', zeros(1,Resource.System.transmitChannels), ...
                   'TXPD', [], ...
                   'peakCutOff', 1.0, ...
                   'peakBLMax', 4.0), 1, P.numRays);

scaleToWvl = 1;
if strcmp(Trans.units, 'mm')
    scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);
end

% - Set event specific TX attributes.
TxOrgX = (-63.5*Trans.spacing):(127*Trans.spacing/(P.numTx-1)):(63.5*Trans.spacing);
Ce = fix(1:(127/(P.numRays-1)):128); % calculate center element numbers.
for n = 1:P.numRays
    TX(n).Origin(1) = TxOrgX(n);

    % Set transmit Apodization.
    lft = Ce(n) - floor(P.numTx/2);
    if lft < 1, lft = 1; end
    rt = Ce(n) + floor(P.numTx/2);
    if rt > Trans.numelements, rt = Trans.numelements; end
    TX(n).aperture = max(1, rt-63);
    TX(n).Apod((lft:rt)-TX(n).aperture + 1) = 1;
    [RIndices,CIndices,V] = find(TX(n).Apod);
    % Add window on aperture apodization
    V = kaiser(size(V,2),1);
    TX(n).Apod(CIndices) = V;

    % Add small steering angle on sides
    if n < 4
        TX(n).Steer(1) = - (4 - n)*(pi/180);
    end
    if n > (P.numRays-3)
        TX(n).Steer(1) = (4 - (P.numRays - n))*(pi/180);
    end
    % Compute transmit delays
    TX(n).Delay = computeTXDelays(TX(n));
    % Compute transmit pixel data
    TX(n).TXPD = computeTXPD(TX(n),PData);
end

% Specify Receive structure arrays.
% - We need numRay Receives for every frame.

% sampling center frequency is 15.625, but we want the bandpass filter
% centered on the actual transducer center frequency of 18 MHz with 67%
% bandwidth, or 12 to 24 MHz.  Coefficients below were set using
% "G3_BPFdevelopment" with normalized cf=1.15 (18 MHz), bw=0.85,
% xsn wdth=0.41 resulting in -3 dB 0.71 to 1.6 (11.1 to 25 MHz), and
% -20 dB 0.57 to 1.74 (8.9 to 27.2 MHz)
%
BPF1 = [ -0.00009 -0.00128 +0.00104 +0.00085 +0.00159 +0.00244 -0.00955 ...
         +0.00079 -0.00476 +0.01108 +0.02103 -0.01892 +0.00281 -0.05206 ...
         +0.01358 +0.06165 +0.00735 +0.09698 -0.27612 -0.10144 +0.48608 ];


maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
Receive = repmat(struct('Apod', ones(1,Resource.System.receiveChannels), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'SampleMode', 'NS200BW', ...
                        'mode', 0, ...
                        'InputFilter', BPF1, ...
                        'callMediaFunc', 0), 1, P.numRays*Resource.RcvBuffer(1).numFrames);
% - Set event specific Receive attributes for each frame.
for i = 1:Resource.RcvBuffer(1).numFrames
    Receive(P.numRays*(i-1)+1).callMediaFunc = 1;
    for j = 1:P.numRays
        Receive(P.numRays*(i-1)+j).aperture = TX(j).aperture;
        Receive(P.numRays*(i-1)+j).framenum = i;
        Receive(P.numRays*(i-1)+j).acqNum = j;
    end
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [0,511,716,920,1023,1023,1023,1023]; %[400,550,650,710,770,830,890,950];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Recon structure arrays.
% - We need one Recon structure which will be used for each frame.
Recon = struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame',-1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...
               'RINums',1:P.numRays);

% Define ReconInfo structures.
% We need numRays ReconInfo structures.
ReconInfo = repmat(struct('mode', 'accumIQ', ...  % default is to accumulate IQ data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'scaleFactor', 0.2, ...
                   'regionnum', 1), 1, P.numRays);

ReconInfo(1).Pre = 'clearInterBuf';
for j = 1:P.numRays
    ReconInfo(j).txnum = j;
    ReconInfo(j).rcvnum = j;
    ReconInfo(j).regionnum = j;
end
ReconInfo(P.numRays).Post = 'IQ2IntensityImageBuf';

% Specify Process structure array.
pers = 20;
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',2.0,...            % pgain is image processing gain
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
SeqControl(2).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(2).argument = 220;  % 220 usec
SeqControl(3).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(3).argument = 35000 - (P.numRays-1)*220;  % 35000 usec = 35msec time between frames
SeqControl(4).command = 'returnToMatlab';
nsc = 5; % nsc is count of SeqControl objects

% Specify Event structure arrays.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    for j = 1:P.numRays        % Acquire frame
        Event(n).info = 'acquire beam';
        Event(n).tx = j;
        Event(n).rcv = P.numRays*(i-1)+j;
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
import vsv.seq.uicontrol.VsSliderControl;

% - Sensitivity Cutoff
UI(1).Control = VsSliderControl('LocationCode','UserB7','Label','Sens. Cutoff',...
                  'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                  'SliderStep',[0.025,0.1],'ValueFormat','%1.3f', ...
                  'Callback', @SensCutoffCallback ); 

% - Range Change
wls2mm = 1;
AxesUnit = 'wls';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
        AxesUnit = 'mm';
        wls2mm = Resource.Parameters.speedOfSound/1000/Trans.frequency;
    end
end
UI(2).Control = VsSliderControl('LocationCode','UserA1','Label',['Range (',AxesUnit,')'],...
                 'SliderMinMaxVal',[64,300,P.endDepth]*wls2mm,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f', ...
                 'Callback', @RangeChangeCallback);  


% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 1;

% Save all the structures to a .mat file.
save('MatFiles/L22-14vX_LF_WideBeam_64');


%% **** Callback routines used by UIControls (UI) ****

%SensCutoffCallback - Sensitivity cutoff change
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

    PData = evalin('base','PData');
    PData.Size(1) = ceil((P.endDepth-P.startDepth)/PData.PDelta(3));
    TxOrgX = (-63.5*Trans.spacing):(127*Trans.spacing/(P.numRays-1)):(63.5*Trans.spacing);

    for n = 1:P.numRays
        PData.Region(n).Shape = struct('Name','Rectangle','Position',[TxOrgX(n),0,0],'width',32,'height',P.endDepth);
    end

    PData.Region = computeRegions(PData);
    assignin('base','PData',PData);
    % Update TXPD data of TX structures.
    TX = evalin('base','TX');
    for i = 1:size(TX,2)
        TX(i).TXPD = computeTXPD(TX(i),PData);
    end
    assignin('base','TX',TX);
    % Update Receive structures.
    Receive = evalin('base', 'Receive');
    maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
    for i = 1:size(Receive,2)
        Receive(i).endDepth = maxAcqLength;
    end
    assignin('base','Receive',Receive);
    evalin('base','TGC.rangeMax = P.endDepth;');
    evalin('base','TGC.Waveform = computeTGCWaveform(TGC);');
    evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData.Size(1)*PData.PDelta(3)/Resource.DisplayWindow(1).pdelta);');
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'PData','InterBuffer','ImageBuffer','TX','Receive','TGC','Recon','DisplayWindow'};
    assignin('base','Control', Control);
    assignin('base', 'action', 'displayChange');
end