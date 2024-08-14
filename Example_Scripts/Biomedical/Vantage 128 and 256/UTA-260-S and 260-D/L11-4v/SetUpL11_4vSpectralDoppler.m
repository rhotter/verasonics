% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: SetUpL11_4vSpectralDoppler.m - Spectral doppler with plane wave transmit.
%
% Description:
%   Sequence programming file for L11-4v Linear array, using plane wave
%   transmits for 2D and doppler imaging.
%   - for the 2D image, a single plane wave transmit is used. All transmit
%     and receive channels are active for each acquisition.
%   - for the doppler ensemble, multiple pulses are transmitted at a
%     steering angle of dopAngle radians. 128 transmit and 64 receive
%     channels are active for each acquisition. The Doppler measurement
%     region is the region under the central 64 transducer elements.
%   - processing is asynchronous with respect to acquisition
%   - An acquisition 'frame' consists of the number of transmit/receive events
%     that equals 10msec of time. Doppler and 2D acquisitions are interleaved.
%   - Doppler PRF will be quantized such that the acquisition PRF (twice
%     the Doppler PRF) has a PRI of integer microseconds.
%   - For the Doppler measurement, the sample volume in the measurement
%     region is selected by clicking the left mouse button in the DisplayWindow.
%     The sample volume is selected as a single I,Q pixel whose location is
%     indicated by a red circle. This action also switches the processing mode
%     to Doppler only. Clicking the right mouse button cancels the Doppler
%     processing and returns to 2D only processing. While the processing
%     operates in duplex mode - either 2D or Doppler, acquisition is always
%     the same, allowing Doppler operation in RcvData Loop playback mode.
%
% Last update:
% 3/27/2020 - support change in Mouse Click tool @Author: Dr. Daniel
% Rohrbach, needed for new UI scheeme to support callback functions.
% 05/12/2020 - Update to SW 4.3 format for new user UIControls and External function definitions (VTS 1691)
%   (~/Example_Scripts/Vantage_Features/New UI Scheme/SetUpL11_5vFlash_NewUI)

clear all

% - 2D P(1) structure
P(1).startDepth = 2;          % Acquisition depth in wavelengths
P(1).endDepth = 128;          % This should preferrably be a multiple of 128 samples.
% - Doppler P(2) structure (This is the pixel region for possible Doppler reconstructions.)
P(2).startDepth = 0;          % Acquisition depth in wavelengths
P(2).endDepth = 128;          % This should preferrably equate to a multiple of 128 samples.

nfrms = 10; % no. of frames to acquire in RcvBuffer (sets size of RF cineloop).
dopPRFnominal = 2000;
framePeriodnominal_uSec = 40000; % nominal frame period for collections of PRIs.
[psched ]   = spectralDoppler('restrictAcqPRF',framePeriodnominal_uSec,dopPRFnominal);
[pschedAll ]= spectralDoppler('restrictAcqPRF',framePeriodnominal_uSec,[]);
disp(' '),
disp([mfilename,': Allowed PRF values: '])
disp(sprintf('%4.4f Hz\n',pschedAll.PRFDop))
framePeriod = psched.FrameTime_uSec*1e-6;
dopPRF = psched.PRFDop;
ttna_microsec = psched.PRIAcq_uSec;
acqPRF = psched.PRFAcq;
maxAcqPRF = psched.maxPRFAcq;
maxDopPRF = psched.maxPRFDop;
minAcqPRF = psched.minPRFAcq;
minDopPRF = psched.minPRFDop;
NWindow = round(dopPRF*0.014); %samples in FFT window (not including zero padding)
SDDespeckle = 0.4;
dopWFCutoffNorm = 0.025; %allfilter normalized cutoff frequency
dopSweepTimeIndex = 2;%1-4 avail for 20msec frame time
SDDynRangeDB = 13;% not including FFT processing gain
sgp= spectralDoppler('spectrogramParameters',framePeriod,dopSweepTimeIndex);
sweepTime = sgp.THist;
baseLineShiftNorm = -0.375;

nPRIs = psched.pulsesPerFrameAcq;  % needs to start as even .default no. of acquisition PRIs in a 'frame'
nDopPRIs =  psched.pulsesPerFrameDop;
dopAngle = 14 * pi/180; % angle for Doppler flash transmits.

maxPRIs = psched.pulsesPerFrameAcqMax ; % max no. of acquisition PRIs in a frame .

% - Specify system parameters
Resource.Parameters.numTransmit = 128;          % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;        % number of receive channels.
Resource.Parameters.speedOfSound = 1540;
Resource.Parameters.speedCorrectionFactor = 1.0;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 %forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 %stops sequence and processes RcvData continuously.

% Specify Trans structure array.
Trans.name = 'L11-4v';
Trans.units = 'wavelengths'; % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);  % computeTrans is used for known transducers.
Trans.maxHighVoltage = 50;    % set a reasonable high voltage limit.

% Specify PData structure arrays.
% - 2D PData structure
PData(1).PDelta = [1.0,0,0.5];
PData(1).Size(1,1) = ceil(P(1).endDepth/PData(1).PDelta(3)); % rows
PData(1).Size(1,2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1)); % cols
PData(1).Size(1,3) = 1;      % single image page
PData(1).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P(1).startDepth]; % x,y,z of upr lft crnr w respect to xducer origin.
% - Doppler PData structure
PData(2).PDelta = [1.0,0,1.0];
PData(2).Size(1,1) = 4;             % 4 IQ samples for Doppler
PData(2).Size(1,2) = 1;             % single column for Doppler
PData(2).Size(1,3) = 1;             % single image page
PData(2).Origin = [0,0,ceil(P(2).endDepth/2)];
% - Region Outline for Region Constrain of Doppler Processing
DopRegionLim.x = [-32,32]; % in wavelength
DopRegionLim.y = [4,128];
PData(3) = PData(1);
PData(3).Region = struct('Shape',struct('Name','Rectangle',...
    'Position',[0,0,DopRegionLim.y(1)],...
    'width',diff(DopRegionLim.x),...
    'height',diff(DopRegionLim.y)));

% Specify Media object.
pt1;
Media.function = 'movePoints';

% - RcvBuffer(1) is for 2D, and Doppler steering angle.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 2048*maxPRIs;
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = nfrms;
% - InterBuffer(1) is for 2D data.
Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).numFrames = 1;          % one intermediate frame needed for 2D.
% - InterBuffer(2) is for Doppler data.
Resource.InterBuffer(2).datatype = 'complex';
Resource.InterBuffer(2).numFrames = 1;          % two intermediate frames needed for color
Resource.InterBuffer(2).pagesPerFrame = maxPRIs/2; % no. of pages set by maxPRIs
% - ImageBuffer(1) is for 2D image.
Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).numFrames = 10;
% - DisplayWindow for 2D image
Resource.DisplayWindow(1).Title = 'L11-4vSpectralDoppler';
Resource.DisplayWindow(1).pdelta = 0.4;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];  % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).numFrames = 20;
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = gray(256);

% Specify Transmit waveform structure.
% - 2D transmit waveform
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,0.67,2,1];
% - Doppler transmit waveform
TW(2).type = 'parametric';
TW(2).Parameters = [Trans.frequency,0.67,6,1];

% Specify TX structure arrays.
% We need 1 TX specification for 2D and 1 TX specification for Doppler (1 steering angle).
TX = repmat(struct('waveform', 1, ...
    'Origin', [0.0,0.0,0.0], ...
    'focus', 0.0, ...
    'Steer', [0.0,0.0], ...
    'Apod', ones(1,Trans.numelements), ...
    'Delay', zeros(1,Trans.numelements)), 1, 2); %
% - Set event specific TX attributes.
TX(1).Delay = computeTXDelays(TX(1));
% -- Second TX struct needed for Doppler
TX(2).waveform = 2;
TX(2).Steer = [dopAngle,0.0];
TX(2).Delay = computeTXDelays(TX(2));

% Specify TGC Waveform structures.
% 2D TGC waveform
TGC(1).CntrlPts = [51,331,430,515,603,690,777,887];
TGC(1).rangeMax = P(1).endDepth;
TGC(1).Waveform = computeTGCWaveform(TGC(1));
% Doppler TGC waveform
TGC(2).CntrlPts = [400,500,640,710,770,830,890,950];
TGC(2).rangeMax = P(1).endDepth;
TGC(2).Waveform = computeTGCWaveform(TGC(2));

% Specify Receive structure arrays -
%   Define enough Receives to handle the maxPRI case.  For lower PRIs, some of the Receives will go
%      unused in the Event list.
%   endDepth - add additional acquisition depth to account for some channels
%              having longer path lengths.
maxAcqLength = ceil(sqrt(P(1).endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
Receive = repmat(struct('Apod', zeros(1,Trans.numelements), ...
                'startDepth', P(1).startDepth, ...
                'endDepth', maxAcqLength, ...
                'TGC', 1, ...
                'bufnum', 1, ...
                'framenum', 1, ...
                'acqNum', 1, ...
                'sampleMode', 'BS100BW',...
                'mode', 0, ...
                'callMediaFunc', 0),1,maxPRIs*nfrms);

% - Set event specific Receive attributes.
for i = 1:nfrms
    k = maxPRIs*(i-1);
    for j = 1:2:maxPRIs
        % -- 2D acquisition.
        Receive(k+j).Apod(1:128) = 1.0;
        Receive(k+j).framenum = i;
        Receive(k+j).acqNum = j;
        if j==1, Receive(k+j).callMediaFunc = 1; end
        % -- Doppler.
        Receive(k+j+1).Apod(33:96) = 1.0; % start Dop aperture at element 33
        Receive(k+j+1).TGC = 2;
        Receive(k+j+1).framenum = i;
        Receive(k+j+1).acqNum = j+1;
    end
end

% Specify Recon structure arrays.
% - We need a Recon structure for 2D and one for Doppler.  The 2D frame will only use
%   the first two 2D acquisitions of the acquisition 'frame'.
% - We only need one Recon structure for the Doppler line.  The Doppler reconstruction will use
%   nPRIs/2 acquisitions from each 'frame' for reconstruction, with each Doppler PRI reconstructed
%   to a different page of the IQData buffer.
Recon = repmat(struct('senscutoff', 0.6, ...
                      'pdatanum', 1, ...
                      'rcvBufFrame', -1, ...  % use the most recent frame
                      'IntBufDest', [1,1], ...
                      'ImgBufDest', [1,-1], ...
                      'RINums', zeros(1,1)), 1, 2);
% - Set Recon values for 2D frame.
Recon(1).RINums(1,1:2) = (1:2);  % 2 ReconInfos needed for 2-1 synthetic aperture
k = 3;
% - Set Recon values for first Doppler ensemble.
Recon(2).pdatanum = 2;
Recon(2).IntBufDest = [2,1];
Recon(2).ImgBufDest = [];
Recon(2).RINums(1,1:(nPRIs/2)) = (k:(k+(nPRIs/2)-1));   % nPRIs/2 ReconInfos needed for Doppler.
Recon(2).senscutoff = 0.8; %reduce aperture to reduce angle spread
% Define ReconInfo structures.
% - For 2D, we need 2 ReconInfo structures.
% - For Doppler, we only need nPRIs/2 ReconInfo structures for the default case, but we will
%   define maxPRIs/2 structures to allow handling the maxPRI case.
ReconInfo = repmat(struct('mode', 'replaceIQ', ...    % set replace IQ data as default.
    'txnum', 1, ...
    'rcvnum', 1, ...
    'pagenum',1, ...
    'regionnum', 1), 1, 2 + nPRIs/2);
% - Set specific ReconInfo attributes.
%   - ReconInfos for 2D frame.
ReconInfo(1).mode = 'replaceIQ'; % replace IQ data
ReconInfo(1).txnum = 1;
ReconInfo(1).rcvnum = 1;
ReconInfo(2).mode = 'accumIQ_replaceIntensity'; % accum and detect
ReconInfo(2).txnum = 1;
ReconInfo(2).rcvnum = 3;
k = 2;
%  - ReconInfos for Doppler ensemble.
for j = 1:nPRIs/2
    ReconInfo(k+j).mode = 'replaceIQ';
    ReconInfo(k+j).txnum = 2;
    ReconInfo(k+j).rcvnum = 2*(j-1) + 2;
    ReconInfo(k+j).pagenum = j;
    ReconInfo(k+j).regionnum = 1;  % only 1 col..
end


% Specify Process structure arrays.
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

Process(2).classname = 'External';
Process(2).method = 'spectralDoppler';
Process(2).Parameters = {'srcbuffer','inter',... % name of buffer to process.
            'srcbufnum',2,...
            'srcframenum',1,...
            'dstbuffer','none'};

% EF2 is external function for ROI plot
Process(3).classname = 'External';
Process(3).method = 'ROIplot';
Process(3).Parameters = {'srcbuffer','none'};

% Specify SeqControl structure arrays.  Missing fields are set to NULL.
%  - Jump back to start of acquisition.
SeqControl(1).command = 'jump';
SeqControl(1).argument = [];  % will be modified later
SeqControl(2).command = 'timeToNextAcq';
SeqControl(2).argument = ttna_microsec ; % time in usec
% loopCnt/loopTst are used to implement a conditional jump
SeqControl(3).command = 'loopCnt';
SeqControl(3).argument = 0;
SeqControl(4).command = 'loopCnt';
SeqControl(4).argument = 1;

nsc = 5; % nsc is count for next SeqControl object

n = 1; % n is count for next Event

Event(n).info = 'ROI plot';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 3;
Event(n).seqControl = 0;
n = n+1;

Event(n).info = 'Start here for 2D only processing.';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 3;
n = n+1;

Event(n).info = 'Jump to start of acquisitions';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 1;
n = n+1;

Event(n).info = 'Start here for Doppler only processing.';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 4;
n = n+1;

nStartAcq = n;
SeqControl(1).argument = nStartAcq;

for i = 1:nfrms  % Acquire frames
    k = maxPRIs*(i-1);
    for j = 1:2:nPRIs % acquisition repeats every 2 acquisitions
        Event(n).info = '2D acquisition.';
        Event(n).tx = 1;
        Event(n).rcv = k+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;

        Event(n).info = 'Doppler acquisition.';
        Event(n).tx = 2;
        Event(n).rcv = k+j+1;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;
    end

    Event(n).info = 'transfer acquired frame data to host';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = nsc;
    SeqControl(nsc).command = 'transferToHost';
    nsc = nsc+1;
    n = n+1;

    Event(n).info = 'test loopCnt';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = nsc;
    SeqControl(nsc).command = 'loopTst';
    SeqControl(nsc).argument = n+3; % jump to Doppler if nz
    nsc = nsc+1;
    n = n+1;

    Event(n).info = 'recon and process 2D';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 1;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = 'jump over Doppler';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = nsc;
    SeqControl(nsc).command = 'jump';
    SeqControl(nsc).argument = n+3;
    nsc = nsc+1;
    n = n+1;

    Event(n).info = 'recon and process Doppler';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 2;
    Event(n).process = 2;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = 'set loopCnt back to 1';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = 4;
    n = n+1;

    Event(n).info = 'exit to Matlab';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = 0;
    % Exit to Matlab every 3rd frame
    if floor(i/3) == i/3
        Event(n).seqControl = nsc;
        SeqControl(nsc).command = 'returnToMatlab';
        nsc = nsc+1;
    end
    n = n+1;
end

Event(n).info = 'Jump back to the start of acquisition';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 1;


% User specified UI Control Elements
import vsv.seq.uicontrol.VsSliderControl; 

% - Sensitivity Cutoff
UI(1).Control = VsSliderControl( 'LocationCode',    'UserB7', ...
                             'Label',           'Sens. Cutoff',   ...
                             'SliderMinMaxVal', [ 0,1.0,Recon(1).senscutoff], ...
                             'SliderStep',      [0.025,0.1], ...
                             'ValueFormat',     '%1.3f', ...
                             'Callback',        @sensitivityCutoffCallback );


% - Range Change
MinMaxVal = [64,300,P(1).endDepth ]; % default unit is wavelength
AxesUnit = 'wls';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
        AxesUnit = 'mm';
        MinMaxVal = MinMaxVal * (Resource.Parameters.speedOfSound/1000/Trans.frequency);
    end
end
UI(2).Control = VsSliderControl( 'LocationCode', 'UserA1','Label',['Range (',AxesUnit,')'],...
                 'SliderMinMaxVal',MinMaxVal,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f', ...
                 'Callback', @rangeChangeCallback);

% - PRF control
UI(3).Control =  VsSliderControl( 'LocationCode', 'UserB4','Label','Doppler PRF',...
    'SliderMinMaxVal',[minDopPRF,maxDopPRF,dopPRF],...
    'SliderStep',[100/ (maxDopPRF-minDopPRF),200/ (maxDopPRF-minDopPRF)],'ValueFormat','%4.0f Hz', ...
    'Callback', @dopplerPRFCallback);

% - Wall Filter control
UI(4).Control =  VsSliderControl( 'LocationCode','UserB3','Label','Doppler WF',...
                 'SliderMinMaxVal',[8,0.1*maxDopPRF,dopWFCutoffNorm*dopPRF],...
                 'SliderStep',[0.02,0.2],'ValueFormat','%4.0f Hz', ...
                 'Callback', @dopplerWallFilterCallback);

% - Sweep Time control
UI(5).Control = VsSliderControl('LocationCode', 'UserB2','Label','Sweep Time',...
                 'SliderMinMaxVal',[min(sgp.THistAllow),max(sgp.THistAllow),sweepTime],...
                 'SliderStep',[0.25,0.25],'ValueFormat','%1.0f s', ...
                 'Callback', @sweepTimeCallback);

% - Baseline control
UI(6).Control =  VsSliderControl('LocationCode', 'UserB1','Label','Baseline',...
                 'SliderMinMaxVal',[-0.5,0.5,baseLineShiftNorm],...
                 'SliderStep',[0.05,0.1],'ValueFormat','%1.2f', ...
                 'Callback', @baselineCallback );

% - Window ButtonDown callback for setting sample volume position.
UI(7).Callback = @wbdCallback; 
UI(7).Statement = @()setCallback(Resource.DisplayWindow(1).Type);

% external functions
import vsv.seq.function.ExFunctionDef;
EF(1).Function = ExFunctionDef('setCallback', @setCallback);
EF(2).Function = ExFunctionDef('ROIplot', @ROIplot ); 


% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 3;

% Save all the structures to a .mat file.
save('MatFiles/L11-4vSpectralDoppler');


%% **** Callback routines used by UIControls (UI) ****

function sensitivityCutoffCallback(~, ~, UIValue)

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

function rangeChangeCallback(hObject, ~, UIValue )
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
    P(1).endDepth = UIValue;
    if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
        if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
            P(1).endDepth = UIValue*scaleToWvl;
        end
    end
    assignin('base','P',P);

    evalin('base','PData(1).Size(1) = ceil((P(1).endDepth-P(1).startDepth)/PData(1).PDelta(3));');
    evalin('base','PData(1).Region = computeRegions(PData(1));');
    evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);');
    Receive = evalin('base', 'Receive');

    maxAcqLength = ceil(sqrt(P(1).endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
    for i = 1:size(Receive,2)
        Receive(i).endDepth = maxAcqLength;
    end
    assignin('base','Receive',Receive);
    evalin('base','TGC(1).rangeMax = P(1).endDepth;');
    evalin('base','TGC(1).Waveform = computeTGCWaveform(TGC(1));');

    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'PData','InterBuffer','ImageBuffer','Receive','TGC','Recon','DisplayWindow'};
    assignin('base','Control', Control);
    assignin('base', 'action', 'displayChange');
end

function dopplerPRFCallback( hObject, ~, UIValue)

    dopPRFnom = UIValue;
    framePeriod = evalin('base','framePeriod');
    psched = spectralDoppler('restrictAcqPRF',round(framePeriod*1e6),dopPRFnom);
    acqPRF = psched.PRFAcq ;
    dopPRF = psched.PRFDop ;
    ttna_microsec =   psched.PRIAcq_uSec ;
%     ttnaNom_microsec = psched.PRIAcq_uSec ;
    framePeriod = psched.FrameTime_uSec*1e-6 ;
    nPRIs= psched.pulsesPerFrameAcq ;
    nDopPRIs = psched.pulsesPerFrameDop ;
%     maxPRIs = psched.pulsesPerFrameAcqMax  ;

    set(hObject,'Value',dopPRF);
    assignin('base', 'dopPRF', dopPRF);
    assignin('base', 'acqPRF', acqPRF);
    assignin('base', 'framePeriod', framePeriod);

    nfrms = evalin('base','nfrms');
    maxPRIs = evalin('base','maxPRIs');
    assignin('base','nDopPRIs',nDopPRIs);
    assignin('base','nPRIs',nPRIs);
    % Update the Doppler Recon to specify nDopPRIs ReconInfo structures.
    Recon = evalin('base','Recon');
    Recon(2).RINums(1,1:nDopPRIs) = (3:(3+(nDopPRIs)-1));
    assignin('base', 'Recon', Recon);

    ReconInfo = evalin('base','ReconInfo');
    %  - ReconInfos for Doppler ensemble.
    k = 2;
    for j = 1:nPRIs/2
        ReconInfo(k+j).mode = 'replaceIQ';
        ReconInfo(k+j).txnum = 2;
        ReconInfo(k+j).rcvnum = 2*(j-1) + 2;
        ReconInfo(k+j).pagenum = j;
        ReconInfo(k+j).regionnum = 1;  % only 1 col..
    end
    assignin('base','ReconInfo',ReconInfo);

    % Re-compute Event structures.
    Event = evalin('base','Event');
    nStartAcq = evalin('base','nStartAcq');
    n = nStartAcq;
    Event = Event(1:n-1); % throw out all but preamble
    SeqControl = evalin('base','SeqControl');
    nsc = 5;
    SeqControl = SeqControl(1:nsc-1); % throw out all but those prior to nsc
    SeqControl(2).argument = ttna_microsec;
    % Re-build the Event list.
    for i = 1:nfrms  % Acquire frames
        k = maxPRIs*(i-1);
        for j = 1:2:nPRIs % acquisition repeats every 2 acquisitions
            Event(n)=struct('info','2D acquisition.',...
               'tx',1,'rcv',k+j,'recon',0,'process',0,'seqControl',2);
               n = n+1;
            Event(n)=struct('info','Doppler acquisition.',...
               'tx',2,'rcv',k+j+1,'recon',0,'process',0,'seqControl',2);
               n = n+1;
        end
        Event(n)=struct('info','transfer acquired frame data to host',...
           'tx',0,'rcv',0,'recon',0,'process',0,'seqControl',nsc);
           SeqControl(nsc).command = 'transferToHost';
           nsc = nsc+1;
           n = n+1;
        Event(n)=struct('info','test loopCnt',...
           'tx',0,'rcv',0,'recon',0,'process',0,'seqControl',nsc);
           SeqControl(nsc).command = 'loopTst';
           SeqControl(nsc).argument = n+3; % jump to Doppler if nz
           nsc = nsc+1;
           n = n+1;
        Event(n)=struct('info','recon and process 2D',...
           'tx',0,'rcv',0,'recon',1,'process',1,'seqControl',0);
           n = n+1;
        Event(n)=struct('info','jump over Doppler',...
           'tx',0,'rcv',0,'recon',0,'process',0,'seqControl',nsc);
           SeqControl(nsc).command = 'jump';
           SeqControl(nsc).argument = n+3;
           nsc = nsc+1;
           n = n+1;
        Event(n)=struct('info','recon and process Doppler',...
           'tx',0,'rcv',0,'recon',2,'process',2,'seqControl',0);
           n = n+1;
        Event(n)=struct('info','set loopCnt back to 1',...
           'tx',0,'rcv',0,'recon',0,'process',0,'seqControl',4);
           n = n+1;
        Event(n)=struct('info','exit to Matlab',...
           'tx',0,'rcv',0,'recon',0,'process',0,'seqControl',0);
           % Exit to Matlab every 3rd frame
           if floor(i/3) == i/3
               Event(n).seqControl = nsc;
               SeqControl(nsc).command = 'returnToMatlab';
               nsc = nsc+1;
           end
           n = n+1;
    end
    Event(n)=struct('info','Jump back to the start of acquisition',...
       'tx',0,'rcv',0,'recon',0,'process',0,'seqControl',1);

    assignin('base', 'Event', Event);
    assignin('base', 'SeqControl', SeqControl);
    % Update VDAS parameters for Receive objects
    evalin('base','VsUpdate(''Receive'')');
    % Set Control.Command to update changed structures in runAcq.
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'Recon','Event','SeqControl'};
    assignin('base','Control', Control);
    % Update the WF control:
    hwf = findobj('tag','UserB3Slider');
    dopWFCutoffNorm = evalin('base' ,'dopWFCutoffNorm');
    set(hwf,'Value',dopPRF*dopWFCutoffNorm);
    hwf = findobj('tag','UserB3Edit');
    set(hwf,'String',num2str(dopPRF*dopWFCutoffNorm,'%4.0f Hz'));
end 

function dopplerWallFilterCallback(~, ~, UIValue)
    dopPRF = evalin('base','dopPRF');
    dopWFCutoffNorm = UIValue/dopPRF;
    assignin('base', 'dopWFCutoffNorm', dopWFCutoffNorm);
end

function sweepTimeCallback( hObject, ~, UIValue)
%-UI#5Callback - Sweep Time
    sweepTime = round(UIValue);
    TFrame=evalin('base','framePeriod');
    sgp=spectralDoppler('spectrogramParameters',TFrame,1);
    [~,swInd] = min(abs(sgp.THistAllow-sweepTime));
    %sweepTimeAllowed = sgp.THistAllow(swInd);

    dopSweepTimeIndex = swInd;

    assignin('base', 'dopSweepTimeIndex', dopSweepTimeIndex);
    set(hObject,'Value',sweepTime);
    assignin('base', 'sweepTime', sweepTime);
end

function baselineCallback( hObject, ~, UIValue)
%-UI#6Callback - Baseline
    baseLineShiftNorm = UIValue;
    baseLineShiftNorm = round(baseLineShiftNorm*20)/20;
    set(hObject,'Value',baseLineShiftNorm);
    assignin('base', 'baseLineShiftNorm', baseLineShiftNorm);
end

function wbdCallback(varargin)
% Window ButtonDown callback for sample volume positioning.
    persistent init wbFig wbAxes elypHndl
    Trans = evalin('base','Trans');
    Resource = evalin('base','Resource');
    DopRegionLim = evalin('base','DopRegionLim');

    if isempty(init)
        init = 1;
        elypHndl = 0;
    end

    % Get pos ofmouse click and ensure the unit is wavelength for spectral doppler
    VsType = Resource.DisplayWindow(1).Type;
    if strcmp(VsType,'Matlab')
        hObject = varargin{1};
        wbFig = evalin('base','Resource.DisplayWindow(1).figureHandle');
        wbAxes = get(wbFig,'CurrentAxes');
        cp = get(wbAxes,'CurrentPoint');
        if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
            if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
                scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);
%                 wl2mm = (Resource.Parameters.speedOfSound/1000)/Trans.frequency;
                cp = cp*scaleToWvl;
                unitMm = 1;
            end
        else
            unitMm = 0;
        end
        % Check for sample volume outside of measurement region.
        x = min(max(cp(1,1),DopRegionLim.x(1)),DopRegionLim.x(2));
        z = min(max(cp(1,2),DopRegionLim.y(1)),DopRegionLim.y(2));
    else
        mct = evalin('base','mct');
        mouseButton = mct.buttonName;
        switch string(mouseButton)
            case 'Start'
                return
            case 'Left'
                hObject.SelectionType = 'normal';
                imageViewer = evalin('base','imageViewer');
                imageCanvas = imageViewer(1).getImageCanvas();
                srcImage = imageCanvas.getSrcImage;
                x = mct.xImgPos/srcImage.pixelUnitScaling;
                z = mct.yImgPos/srcImage.pixelUnitScaling;
            case 'Right'
                hObject.SelectionType = 'alt';
        end
    end

    % Change Event Sequence if left mouse button is pressed
    if strcmp(hObject.SelectionType,'normal')
        if elypHndl == 0
            % set startEvent for Doppler only.
            Control = evalin('base','Control');
            Control(2).Command = 'set&Run';
            Control(2).Parameters = {'Parameters',1,'startEvent',4};
            evalin('base','Resource.Parameters.startEvent = 4;');
        else
            delete(elypHndl); % delete previous annotation circle
        end

        if strcmp(VsType,'Matlab')
            % set annotation in the Matlab Viewer (Verasonics Viewer is done in MouseClickTool)
            abox = [x-1,z-1,2,2];
            axun = get(wbAxes,'Units'); % save axis units
            set(wbAxes,'Units','normalized'); % make axes units normalized
            axpos = get(wbAxes,'Position'); % get axes position
            axlim = axis(wbAxes); % get the axis limits [xlim ylim (zlim)]
            if unitMm, axlim = axlim*scaleToWvl; end
            axwidth = diff(axlim(1:2));
            axheight = diff(axlim(3:4));
            % Transform from data space coordinates to normalized figure coordinates
            bbox(1) = (abox(1)-axlim(1))/axwidth * axpos(3) + axpos(1);
            bbox(2) = (axlim(4)-abox(2))/axheight * axpos(4) + axpos(2);
            bbox(3) = abox(3) * axpos(3)/axwidth;
            bbox(4) = abox(4) * axpos(4)/axheight;
            set(wbAxes,'Units',axun); % restore axis units
            elypHndl = annotation('ellipse',bbox,'EdgeColor',[1 .2 .2],'LineWidth',1.0);
        end

        assignin('base','SVPos',[x,z-2])
        %Modify structures for new SV position.
        evalin('base','PData(2).Origin = [SVPos(1),0,SVPos(2)];');
        evalin('base','PData(2).Region=computeRegions(PData(2));');

        Control(1).Command = 'update&Run';
        Control(1).Parameters = {'PData'};
        assignin('base','Control', Control);
        % if right mouse button ...
    elseif strcmp(hObject.SelectionType,'alt')
        % delete annotation circle
        if elypHndl ~= 0
            delete(elypHndl);
            elypHndl = 0;
        end
        % set startEvent for 2D only
        Control.Command = 'set&Run';
        Control.Parameters = {'Parameters',1,'startEvent',1};
        evalin('base','Resource.Parameters.startEvent = 1;');
        assignin('base','Control',Control);
    end
end

%% **** Callback routines used by External function definition (EF) ****

function setCallback(VsType)
    switch VsType
        case 'Verasonics'
            
            vantageWindow = evalin('base','vantageWindow');
            imageViewer = evalin('base','imageViewer');
            DopRegionLim = evalin('base','DopRegionLim');
            mct  = javaObjectEDT( 'com.verasonics.viewer.tools.mouseclicktool.MouseClickTool', vantageWindow, imageViewer(1));
            mcth = javaObjectEDT( handle(mct, 'callbackproperties' ) );
            mct.setToolEnabled(true);
            mct.lowerLimitX = DopRegionLim.x(1);
            mct.upperLimitX = DopRegionLim.x(2);
            mct.lowerLimitY = DopRegionLim.y(1);
            mct.upperLimitY = DopRegionLim.y(2);
            mcth.ClickedEventCallback = @wbdCallback;
            assignin('base','mct',mct);
        case 'Matlab'
            evalin('base','set(Resource.DisplayWindow(1).figureHandle,''WindowButtonDownFcn'',@wbdCallback);');
    end
end

function ROIplot(varargin)
%ROIplot
    persistent drawROI

    % drawRegionOutline(WinNum,PDataNum,RegionNum) creates an outline from
    % PData(PDataNum).Region(RegionNum) on displayWindow(WinNum) with default color - white.
    if isempty(drawROI)
        drawROI = 1;
        evalin('base','hROI = drawRegionOutline(1,3,1);')
    else
        evalin('base','drawRegionOutline(hROI,1,3,1);')
    end
end