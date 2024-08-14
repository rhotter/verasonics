% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: SetUpL11_5vFlashDoppler_GPU_IQ.m - Example of external CDI
% processing on the GPU
%
% Description:
%  - This example is derived from the original Doppler example script:
%    "SetUpL11_5vFlashDoppler.m"
%
%  - The purpose is to illustrate to the user how to pass IQ data to the GPU
%    for parallel processing.
%
%  - There are 4 examples of processing Doppler data externally included:
%       -- 'externalDoppler_CPU_M';
%           External Doppler processing on the CPU using a matlab script
%       -- 'externalDoppler_GPU_PCT'
%           External Doppler processing on the GPU using Matlab's parallel
%           computing toolbox
%       -- 'externalDoppler_CPU_MEX'
%           External Doppler processing on the CPU using c code compiled
%           into a .mex file
%       -- 'externalDoppler_CPU_MULTI_MEX'
%           External Doppler multi-threaded processing on the CPU using c 
%           code compiled into a .mex file
%       -- 'externalDoppler_GPU_MEXCUDA';
%           External Doppler processing on the GPU using cuda code
%           compiled into a .mex file
%
%   - The script "compileExternalDopplerIQ.m" must be run prior to this
%     script to generate the mex files.
%
%   - A Doppler "string phantom" is available in the Media structure to be
%     used for testing the code in simulate mode
%
%  NOTE:  GPU specific code is commented with @gpu




clear all

% P(1) is used for Bmode, P(2) is used for Doppler
P(1).startDepth = 5;   % Acquisition depth in wavelengths
P(1).endDepth = 192;   % This should preferrably be a multiple of 128 samples.
P(2).startDepth = 0;   % Acquisition depth in wavelengths
P(2).endDepth = 128;   % This should preferrably be a multiple of 128 samples.

% Set 2D parameters
na = 7;      % Set na = number of flash angles for 2D.
if (na > 1)
    dtheta2D = (30*pi/180)/(na-1);
    startAngle = -30*pi/180/2;
else
    dtheta2D = 0;
    startAngle=0;
end % set dtheta2D to range over +/- 15 degrees.

% Set Doppler parameters
ne = 14;     % Set ne = number of acquisitions in Doppler ensemble.
drops = 2;  %number of doppler pulses to drop
Np = ne - drops;
dopAngle = 12 * pi/180;
dopPRF = 4.0e+03; % Doppler PRF in Hz.
pwrThres = 0.65;
DopState = 'computeCFIFreqEst';
cpt = 28;  % define here so we can use in UIControl below
persf = 80;
persp = 90;
threshParams = struct('threshold', pwrThres, 'PowerEstMode', strcmp(DopState,'computeCFIPowerEst'), 'Npri', Np);

% Specify system parameters.
Resource.Parameters.numTransmit = 128;          % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;        % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Specify Trans structure array.
Trans.name = 'L11-5v';
Trans.units = 'wavelengths'; % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);  % computeTrans is used for known transducers.
Trans.maxHighVoltage = 50;    % set a reasonable high voltage limit.

% Specify PData structure arrays.
% - 2D PData structure
PData(1).PDelta = [1, 0, 0.5];
PData(1).Size(1) = ceil((P(1).endDepth-P(1).startDepth)/PData(1).PDelta(3)); % rows
PData(1).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1)); % cols
PData(1).Size(3) = 1;             % single image page
PData(1).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P(1).startDepth]; % x,y,z of uppr lft crnr.
% - Doppler PData structure
PData(2).PDelta = [1, 0, 0.5];
PData(2).Size(1) = ceil((P(2).endDepth-P(2).startDepth)/PData(2).PDelta(3)); % flow window rows
PData(2).Size(2) = ceil((128*Trans.spacing)/PData(2).PDelta(1));  % flow window columns
PData(2).Size(3) = 1;             % single image page
PData(2).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P(2).startDepth]; % x,y,z of upper lft crnr.

% Specify Media object.
Media = vsv.gpu.media.createDopplerMedia();

% Specify Resources.
% - RcvBuffer(1) is for both 2D and Doppler acquisitions.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 2048*(na + ne);
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 20;           % 20 frames allocated for RF acqusitions.
% InterBuffer(1) is for 2D reconstructions.
Resource.InterBuffer(1).numFrames = 1;          % one intermediate frame needed for 2D.
% InterBuffer(2) is for Doppler reconstructions.
Resource.InterBuffer(2).numFrames = 1;          % one intermediate frame needed for Doppler.
Resource.InterBuffer(2).pagesPerFrame = ne;     % ne pages per ensemble
Resource.InterBuffer(2).gpuMemoryType = 'CPU';
% ImageBuffer(1) is for 2D image.
Resource.ImageBuffer(1).numFrames = 20;
% ImageBuffer(2) is for Doppler image.
Resource.ImageBuffer(2).numFrames = 20;
Resource.ImageBuffer(2).gpuMemoryType = 'CPU';
% DisplayWindow is for 2D combined with Doppler
Resource.DisplayWindow(1).Title = 'L11-5vFlashDoppler';
Resource.DisplayWindow(1).pdelta = 0.4;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
    DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)]; % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).numFrames = 20;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = grayscaleCFImap;
Resource.DisplayWindow(1).splitPalette = 1;

% ------Specify structures used in Events------
% Specify Transmit waveform structures.
% - 2D transmit waveform
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,0.67,2,1];
% - Doppler transmit waveform, transmit frequency should be equivalent to
% supported demodFrequency
TW(2).type = 'parametric';
TW(2).Parameters = [6.25,0.67,6,1];

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
    'Origin', [0.0,0.0,0.0], ...
    'focus', 0.0, ...
    'Steer', [0.0,0.0], ...
    'Apod', ones(1,Trans.numelements), ...
    'Delay', zeros(1,Trans.numelements)), 1, na+1); % na TXs for 2D + 1 for Doppler
% - Set event specific TX attributes.
for n = 1:na   % na transmit events for 2D
    TX(n).Steer = [(startAngle+(n-1)*dtheta2D),0.0];
    TX(n).Delay = computeTXDelays(TX(n));
end
% -- only one TX struct needed for Doppler
TX(na+1).waveform = 2;
TX(na+1).Steer = [dopAngle,0.0];
TX(na+1).Delay = computeTXDelays(TX(na+1));

% Specify TPC structures.
TPC(1).name = '2D';
TPC(1).maxHighVoltage = 50;
TPC(2).name = 'Doppler';
TPC(2).maxHighVoltage = 35;

RcvProfile(1).LnaGain = 18; % Profile used for imaging
RcvProfile(2).LnaGain = 24; % Profile used for Doppler
RcvProfile(2).LnaZinSel = 31; % Force high-Z state for best Doppler sensitivity

% Specify TGC Waveform structures.
% - 2D TGC
TGC.CntrlPts = [0,141,275,404,510,603,702,782];
TGC(1).rangeMax = P(1).endDepth;
TGC(1).Waveform = computeTGCWaveform(TGC(1));
% - Doppler TGC
TGC(2).CntrlPts =[ 0 272 662 662 662 662 662 662];
TGC(2).rangeMax = P(2).endDepth;
TGC(2).Waveform = computeTGCWaveform(TGC(2));

% Specify Receive structure arrays.
%   We need to acquire all the 2D and Doppler data within a single RcvBuffer frame.  This allows
%   the transfer-to-host DMA after each frame to transfer a large amount of data, improving throughput.
% - We need na Receives for a 2D frame and ne Receives for a Doppler frame.
maxAcqLngth2D = sqrt(P(1).endDepth^2 + (Trans.numelements*Trans.spacing)^2) - P(1).startDepth;
maxAcqLngthDop =  sqrt(P(2).endDepth^2 + (96*Trans.spacing)^2) - P(2).startDepth;
wl4sPer128 = 128/(4*2);  % wavelengths in a 128 sample block for 4 smpls per wave round trip.
wl2sPer128 = 128/(2*2);  % wavelengths in a 128 sample block for 2 smpls per wave round trip.
Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
    'startDepth', P(1).startDepth, ...
    'endDepth', P(1).startDepth + wl4sPer128*ceil(maxAcqLngth2D/wl4sPer128), ...
    'TGC', 1, ...
    'bufnum', 1, ...
    'framenum', 1, ...
    'acqNum', 1, ...
    'sampleMode', 'NS200BW', ... % 200% Bandwidth for 2D
    'mode', 0, ...
    'callMediaFunc', 0), 1, (na+ne)*Resource.RcvBuffer(1).numFrames);
%   'InputFilter', [-0.0141,0.0420,-0.0825,0.1267,-0.1612,0.1783], ...  % bandpass filter for Doppler (~ 20% BW)
% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    k = (na + ne)*(i-1); % k keeps track of Receive index increment per frame.
    % - Set attributes for each frame.
    %Receive(k+1).callMediaFunc = 1;
    for j = 1:na  % acquisitions for 2D
        Receive(j+k).callMediaFunc = 1;
        Receive(j+k).framenum = i;
        Receive(j+k).acqNum = j;
    end
    for j = (na+1):(na+ne)
        % Doppler acquisition
        Receive(j+k).callMediaFunc = 1;
        Receive(j+k).startDepth = P(2).startDepth;
        Receive(j+k).endDepth = P(2).startDepth + wl2sPer128*ceil(maxAcqLngthDop/wl2sPer128);
        Receive(j+k).sampleMode = 'BS100BW';
        Receive(j+k).demodFrequency = TW(2).Parameters(1);
        Receive(j+k).TGC = 2;
        Receive(j+k).framenum = i;
        Receive(j+k).acqNum = j;        % Doppler acqNums continue after 2D
    end
end

% Specify Recon structure arrays.
% - We need two Recon structures, one for 2D, one for Doppler. These will be referenced in the same
%   event, so that they will use the same (most recent) acquisition frame.
Recon = repmat(struct('senscutoff', 0.85, ...
    'pdatanum', 1, ...
    'rcvBufFrame', -1, ...
    'IntBufDest', [1,1], ...
    'ImgBufDest', [1,-1], ...
    'RINums', zeros(1,1)), 1, 2);
% - Set Recon values for 2D frame.
Recon(1).RINums(1,1:na) = (1:na);  % na ReconInfos needed for na angles
k = na + 1;
% - Set Recon values for Doppler ensemble.
Recon(2).pdatanum = 2;
Recon(2).IntBufDest = [2,1];
Recon(2).ImgBufDest = [0,0];
Recon(2).RINums(1,1:(ne-drops)) = (k:(k+ne-1-drops));   % ne ReconInfos needed for Doppler ensemble.

% Define ReconInfo structures.
% - For 2D, we need na ReconInfo structures for na steering angles.
% - For Doppler, we need ne ReconInfo structures.
ReconInfo = repmat(struct('mode', 'accumIQ', ...    % accumulate IQ data.
    'txnum', 1, ...
    'rcvnum', 1, ...
    'pagenum',1, ...
    'regionnum', 1), 1, na + (ne-drops));
% - Set specific ReconInfo attributes.
%   - ReconInfos for 2D frame.
if na>1
    ReconInfo(1).mode = 'replaceIQ'; % replace IQ data
    for j = 1:na
        ReconInfo(j).txnum = j;
        ReconInfo(j).rcvnum = j;
    end
    ReconInfo(na).mode = 'accumIQ_replaceIntensity'; % accum and detect
else
    ReconInfo(1).mode = 'replaceIntensity';
end

%  - ReconInfos for Doppler ensemble.
k = na;

for j = 1:(ne-drops)
    ReconInfo(k+j).mode = 'replaceIQ';
    ReconInfo(k+j).txnum = na + 1;
    ReconInfo(k+j).rcvnum = drops + na + j;
    ReconInfo(k+j).pagenum = j;
end

% Specify Process structure arrays.
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
    'framenum',-1,...   % (-1 => lastFrame)
    'pdatanum',1,...    % number of PData structure to use
    'pgain',1.0,...            % pgain is image processing gain
    'reject',2,...      % reject level
    'persistMethod','simple',...
    'persistLevel',20,...
    'interpMethod','4pt',...
    'grainRemoval','none',...
    'processMethod','none',...
    'averageMethod','none',...
    'compressMethod','power',...
    'compressFactor',40,...
    'mappingMethod','lowerHalf',...
    'display',0,...      % don't display image after processing
    'displayWindow',1};

%-- @gpu: This process is used for accessing all of the external GPU processing functions --%
Process(2).classname = 'External';
Process(2).method = 'externalDoppler_CPU_M';  %start with external .m function
Process(2).Parameters = {'srcbuffer','inter', ...
    'srcbufnum', 2, ... % no. of buffer to process.
    'srcframenum', 1, ... % starting frame no.
    'dstbuffer', 'image', ...
    'dstbufnum', 2, ...
    'dstframenum', -2}; %<--- NOTE the '-2' is used.  See Verasonics sequence progrmaming manual for more details

%-- @gpu: This process is used for accessing all of the external GPU processing functions --%
Process(3).classname = 'External';
Process(3).method = 'externalDoppler_GPU_PCT';  %start with external .m function
Process(3).Parameters = {'srcbuffer','inter', ...
    'srcbufnum', 2, ... % no. of buffer to process.
    'srcframenum', 1, ... % starting frame no.
    'dstbuffer', 'image', ...
    'dstbufnum', 2, ...
    'dstframenum', -2}; %<--- NOTE the '-2' is used.  See Verasonics sequence progrmaming manual for more details


Process(4).classname = 'External';
Process(4).method = 'externalDoppler_CPU_MEX';
Process(4).Parameters = {'srcbuffer','inter', ...
    'srcbufnum', 2, ... % no. of buffer to process.
    'srcframenum', 1, ... % starting frame no.
    'dstbuffer', 'image', ...
    'dstbufnum', 2, ...
    'dstframenum', -2}; %<--- NOTE the '-2' is used.  See Verasonics sequence progrmaming manual for more details

Process(5).classname = 'External';
Process(5).method = 'externalDoppler_CPU_MULTI_MEX';
Process(5).Parameters = {'srcbuffer','inter', ...
    'srcbufnum', 2, ... % no. of buffer to process.
    'srcframenum', 1, ... % starting frame no.
    'dstbuffer', 'image', ...
    'dstbufnum', 2, ...
    'dstframenum', -2}; %<--- NOTE the '-2' is used.  See Verasonics sequence progrmaming manual for more details


Process(6).classname = 'External';
Process(6).method = 'externalDoppler_GPU_MEXCUDA';  %start with external .m function
Process(6).Parameters = {'srcbuffer','inter', ...
    'srcbufnum', 2, ... % no. of buffer to process.
    'srcframenum', 1, ... % starting frame no.
    'dstbuffer', 'image', ...
    'dstbufnum', 2, ...
    'dstframenum', -2}; %<--- NOTE the '-2' is used.  See Verasonics sequence progrmaming manual for more details

IMGPROC = 7;
Process(IMGPROC).classname = 'Image';
Process(IMGPROC).method = 'imageDisplay';
Process(IMGPROC).Parameters = {'imgbufnum',2,...   % number of buffer to process.
    'framenum',-1,...   % (-1 => lastFrame)
    'srcData','signedColor',... % type of data to display.
    'pdatanum',2,...    % number of PData structure to use
    'pgain',1.0,...            % pgain is image processing gain
    'reject',2,...      % reject level
    'persistMethod','dynamic',...
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

% Specify SeqControl structure arrays.
% -- Time between 2D flash angle acquisitions
SeqControl(1).command = 'timeToNextAcq';
SeqControl(1).argument = 200;
% -- Change to Profile 2 (Doppler)
SeqControl(2).command = 'setTPCProfile';
SeqControl(2).condition = 'next';
SeqControl(2).argument = 2;
% -- Time between 2D acquisition and Doppler ensemble. Set to allow time for profile change.
SeqControl(3).command = 'timeToNextAcq';
SeqControl(3).argument = 7000; % time in usec
% -- PRF for Doppler ensemble
SeqControl(4).command = 'timeToNextAcq';
SeqControl(4).argument = round(1/(dopPRF*1e-06)); % (for 3KHz dopPRF & 14 ensemble = 4.7 msecs)
% -- Change to Profile 1 (2D)
SeqControl(5).command = 'setTPCProfile';
SeqControl(5).condition = 'next';
SeqControl(5).argument = 1;
% -- Time between Doppler and next 2D acquisition. Set to allow time for profile change.
SeqControl(6).command = 'timeToNextAcq';
SeqControl(6).argument = 7000; % time in usec
% -- Jump back to start.
SeqControl(7).command = 'jump';
SeqControl(7).argument = 1;
% set receive profile
SeqControl(8).command = 'setRcvProfile';
SeqControl(8).argument = 1;
SeqControl(9).command = 'setRcvProfile';
SeqControl(9).argument = 2;
% - The remainder of the SeqControl structures are defined dynamically in the sequence events.
%   The variable nsc keeps track of the next SeqControl structure index.
nsc = 10;  % next SeqControl number

% Specify Event structure arrays.
externalProcessEvents=[];
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    % Acquire 2D frame
    for j = 1:na
        Event(n).info = 'Acquire 2D flash angle';
        Event(n).tx = j;
        Event(n).rcv = (na+ne)*(i-1)+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 1;
        if j == 1
            Event(n).seqControl = [1,8];
        end
        n = n+1;
    end
    Event(n-1).seqControl = [2,3];   % replace last 2D acquisition Event's seqControl
    % Acquire Doppler ensemble.
    for j = (na+1):(na+ne)
        Event(n).info = 'Acquire Doppler ensemble';
        Event(n).tx = na+1;
        Event(n).rcv = (na+ne)*(i-1)+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 4;
        if j == na+1
            Event(n).seqControl = [4,9];
        end
        n = n+1;
    end
    Event(n-1).seqControl = [5,6,nsc]; % replace last Doppler acquisition Event's seqControl
    SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
    nsc = nsc+1;
    
    Event(n).info = 'recons and 2D process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = [1,2];
    Event(n).process = 1;
    Event(n).seqControl = 0;
    n = n+1;
    
    externalProcessEvents(end+1) = n;  %@gpu: keep track of which events contain the external process for switching from GUI
    Event(n).info = 'External Doppler processing';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 2; % @gpu: call the external Doppler processing Process here
    Event(n).seqControl = 0;
    n = n+1;
    
    Event(n).info = 'Doppler image display';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = IMGPROC;
    Event(n).seqControl = 0;
    if floor(i/3) == i/3     % Exit to Matlab every 3rd frame
        Event(n).seqControl = nsc;
        SeqControl(nsc).command = 'returnToMatlab';
        nsc = nsc+1;
    end
    n = n+1;
end

Event(n).info = 'Jump back';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 7;


% User specified UI Control Elements
import vsv.seq.uicontrol.VsButtonGroupControl;
import vsv.seq.uicontrol.VsSliderControl

% - External Process Selection - @GPU:  a UI to allow the user to quickly
% switch between different external processing methods (using update&run)
% without having to reload the sequence.
UI(1).Control = VsButtonGroupControl('LocationCode','UserB7',...
    'Title','External Process', ...
    'NumButtons', 5, ...
    'PossibleCases', {'CPU .m','GPU PCT','CPU .mex', 'CPU Multi .mex', 'GPU .mexcuda'}, ...
    'Callback', @ExternalProcessMethodCallback);

% Doppler Method selection
UI(2).Control = vsv.seq.uicontrol.VsButtonGroupControl(...
    'LocationCode','UserB5',...
    'Title','Doppler Mode',...
    'PossibleCases',{'Velocity','Power'},...
    'Callback',@velocityPowerCallback);

% Doppler Power Threshold Slider
UI(3).Control = VsSliderControl('LocationCode','UserB4',...
    'Label','DopPwerThres',...
    'SliderMinMaxVal',[0,1.0,pwrThres],...
    'SliderStep',[0.025 0.1],'ValueFormat','%1.3f',...
    'Callback',@powerThresCallback);

% Color Priority Threshold Slider
UI(4).Control = VsSliderControl('LocationCode','UserB3',...
    'Label','Color Priority',...
    'SliderMinMaxVal',[0,255,cpt],...
    'SliderStep',[0.025 0.1],'ValueFormat','%3.0f',...
    'Callback',@ColorPriorityCallback);

% - Color Persistence Slider
UI(5).Control = VsSliderControl('LocationCode','UserB2',...
    'Label','Color Persistence','SliderMinMaxVal',[0,100,persf],...
    'SliderStep',[1/100,0.1],'ValueFormat','%3.0f',...
    'Callback', @ColorPersistenceCallback);


% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 3;

% Save all the structures to a .mat file.
save('MatFiles/L11-5vFlashDoppler_GPU');

return



function ExternalProcessMethodCallback(~, ~, UIState)
%@GPU: This is the callback for the UI to quickly switch between different
%external functions to compare their framerates without having to restart
%the sequence each time.
externalIndex = [];
switch UIState
    case 1
        disp('Switch to externalDoppler_CPU_M')
        externalIndex = 2;
    case 2
        v = ver;
        if (license('test','Distrib_Computing_Toolbox')  && any(strcmp('Parallel Computing Toolbox', {v.Name})) && (gpuDeviceCount()>0))
            disp('Switch to externalDoppler_GPU_PCT')
            externalIndex = 3;
        else
            warning('A compatible GPU Device & Matlab Parallel Computing Toolbox need to be installed to be able to test this function')
        end
    case 3
        if ~isempty(which('externalDoppler_CPU_MEX'))
            disp('Switch to externalDoppler_CPU_MEX')
            externalIndex = 4;
        else
            warning('externalDoppler_CPU_MEX.mex does not exist or is not in the path.  Please type ''compileExampleGPUCode'' from the command line to generate the mex file')
        end
        
    case 4
        if ~isempty(which('externalDoppler_CPU_MULTI_MEX'))
            disp('Switch to externalDoppler_CPU_MULTI_MEX')
            externalIndex = 5;
        else
            warning('externalDoppler_CPU_MEX.mex does not exist or is not in the path.  Please type ''compileExampleGPUCode'' from the command line to generate the mex file')
        end        
        
    case 5
        if ~isempty(which('externalDoppler_GPU_MEXCUDA'))
            disp('Switch to externalDoppler_GPU_MEXCUDA')
            externalIndex = 6;
        else
            warning('externalDoppler_CPU_MEX.mex does not exist or is not in the path.  Please type ''compileExampleGPUCode'' from the command line to generate the mex file')
        end
end

if ~isempty(externalIndex)
    Event=evalin('base','Event');
    externalProcessEvents = evalin('base','externalProcessEvents');
    for i = 1:length(externalProcessEvents)
        Event(externalProcessEvents(i)).process = externalIndex;  %replace external events with new process index
    end
    assignin('base','Event', Event);
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'Event'};
    assignin('base','Control', Control);
end

end

function velocityPowerCallback(~, ~, UIValue)
    IMGPROC = evalin('base','IMGPROC');
    Control = repmat(struct('Command','set&Run','Parameters',[]),1,3);
    Process = evalin('base','Process');
    Resource = evalin('base','Resource');
    
    switch UIValue
        case 1 % Velocity mode
            newMap = grayscaleCFImap;
            newMap(1:128,:) = Resource.DisplayWindow(1).Colormap(1:128,:);
            Resource.DisplayWindow(1).Colormap = newMap;
            assignin('base','persp',get(findobj('Tag','UserB2Slider'),'Value'));
            persf = evalin('base','persf'); persValue = persf(1);
            Control(1).Parameters = {'Process',IMGPROC,'srcData',...
                'signedColor','persistMethod',...
                'dynamic','persistLevel',persValue};
            Control(2).Parameters = {'DisplayWindow',1,'colormap',newMap};
            Control(3).Parameters = {'ImageBuffer',1,'lastFrame',0};
            set(findobj('tag','UserB2Edit'),'String',num2str(persValue,'%3.0f'));
            set(findobj('tag','UserB2Slider'),'Value',persValue);
            assignin('base','DopState','computeCFIFreqEst');
            % Set modified Process attributes in base Matlab environment.
            for k = 1:2:length(Process(IMGPROC).Parameters)
                if strcmp(Process(IMGPROC).Parameters{k},'srcData')
                    Process(IMGPROC).Parameters{k+1} = 'signedColor';
                elseif strcmp(Process(IMGPROC).Parameters{k},'persistMethod')
                    Process(IMGPROC).Parameters{k+1} = 'dynamic';
                elseif strcmp(Process(IMGPROC).Parameters{k},'persistLevel')
                    Process(IMGPROC).Parameters{k+1} = persf;
                end
            end
        case 2 % Power mode
            newMap = grayscaleCPAmap;
            newMap(1:128,:) = Resource.DisplayWindow(1).Colormap(1:128,:);
            Resource.DisplayWindow(1).Colormap = newMap;
            for k = 1:2:length(Process(IMGPROC).Parameters)
                if strcmp(Process(IMGPROC).Parameters{k},'persistLevel')
                    persf = Process(IMGPROC).Parameters{k+1};
                end
            end
            assignin('base','persf',persf);
            persValue = evalin('base','persp');
            Control(1).Parameters = {'Process',IMGPROC,'srcData',...
                'unsignedColor','persistMethod',...
                'simple','persistLevel',persValue};
            Control(2).Parameters = {'DisplayWindow',1,'colormap',newMap};
            Control(3).Parameters = {'ImageBuffer',1,'lastFrame',0};
            set(findobj('tag','UserB2Edit'),'String',num2str(persValue,'%3.0f'));
            set(findobj('tag','UserB2Slider'),'Value',persValue);
            assignin('base','DopState','computeCFIPowerEst');
            for k = 1:2:length(Process(IMGPROC).Parameters)
                if strcmp(Process(IMGPROC).Parameters{k},'srcData')
                    Process(IMGPROC).Parameters{k+1} = 'unsignedColor';
                elseif strcmp(Process(IMGPROC).Parameters{k},'persistMethod')
                    Process(IMGPROC).Parameters{k+1} = 'simple';
                elseif strcmp(Process(IMGPROC).Parameters{k},'persistLevel')
                    Process(IMGPROC).Parameters{k+1} = persValue;
                end
            end
    end
    assignin('base','newMap',newMap);
    evalin('base','Resource.DisplayWindow(1).Colormap = newMap;');
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
    if ishandle(hCMTool)
        delete(hCMTool);
        set(findobj('tag','toolsMenu'),'Value',1); % set tools selection back to none
    end

    evalin('base', 'threshParams = struct(''threshold'', pwrThres, ''PowerEstMode'', strcmp(DopState,''computeCFIPowerEst''), ''Npri'', Np);');
end
function powerThresCallback(~,~,UIValue)
    assignin('base','pwrThres', UIValue);
    evalin('base', 'threshParams = struct(''threshold'', pwrThres, ''PowerEstMode'', strcmp(DopState,''computeCFIPowerEst''), ''Npri'', Np);');
end

function ColorPriorityCallback(~, ~, UIValue)
    % Set the value in the Process structure for use in cineloop playback.
    IMGPROC = evalin('base','IMGPROC');
    Process = evalin('base','Process');
    for k = 1:2:length(Process(IMGPROC).Parameters)
        if strcmp(Process(IMGPROC).Parameters{k},'threshold')
            Process(IMGPROC).Parameters{k+1} = UIValue;
        end
    end
    assignin('base','Process',Process);
    % Set Control.Command to set Image.threshold.
    Control = evalin('base','Control');
    Control.Command = 'set&Run';
    Control.Parameters = {'Process', IMGPROC,'threshold',UIValue};
    assignin('base','Control', Control);
end

function ColorPersistenceCallback(~, ~, UIValue)
    % Set the value in the Process structure for use in cineloop playback.
    IMGPROC = evalin('base','IMGPROC');
    Process = evalin('base','Process');
    for k = 1:2:length(Process(IMGPROC).Parameters)
        if strcmp(Process(IMGPROC).Parameters{k},'persistLevel')
            Process(IMGPROC).Parameters{k+1} = UIValue;
        end
    end
    assignin('base','Process',Process);
    % Set Control.Command to set Image.persistLevel.
    Control = evalin('base','Control');
    Control.Command = 'set&Run';
    Control.Parameters = {'Process', IMGPROC, 'persistLevel', UIValue};
    assignin('base','Control', Control);

    % If PTool window is open, adjust persistLevel1 in Process(IMGPROC)
    hPTool = findobj('tag','ProcessTool');
    if ishandle(hPTool)
        if isequal(get(findobj('tag','processNum'),'Value'), IMGPROC)
            set(findobj('tag','persistSlider1'),'Value',UIValue);
            set(findobj('tag','persistValue1'),'String',num2str(UIValue));
        end
    end
end