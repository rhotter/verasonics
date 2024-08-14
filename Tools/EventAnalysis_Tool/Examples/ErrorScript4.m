% ErrorScript4 has a bug about Receive.endDepth
% The endDepth of all Receives used for one Recon should be identical,
% otherwise the image will not be display correctly because of incorrect
% data source for image reconstruction

clear all
% P(1) is used for Bmode, P(2) is used for Doppler
P(1).startDepth = 5;   % Acquisition depth in wavelengths
P(1).endDepth = 9*16;  % This should preferrably be a multiple of 128 samples (4 samples/wv).
P(2) = P(1);           % Default Doppler depth is the same as Bmode

% Set 2D parameters
na = 7;      % Set na = number of flash angles for 2D.
if (na > 1), dtheta = (30*pi/180)/(na-1); startAngle = -30*pi/180/2; else dtheta = 0; startAngle=0; end % set dtheta to range over +/- 15 degrees.

% Set Doppler parameters
ne = 14;     % Set ne = number of acquisitions in Doppler ensemble.
dopAngle = 12 * pi/180;
dopPRF = 3.0e+03; % Doppler PRF in Hz.
pwrThres = 0.5;

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
Trans.name = 'L12-3v';
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
PData(2) = PData(1); % default is the same as PData(1), but the user may want to change the PDelta
regionH = P(2).endDepth-P(2).startDepth;
regionW = ceil((128*Trans.spacing)/PData(2).PDelta(1)); % 128 elements for Doppler
PData(2).Region.Shape = struct(...
    'Name','Parallelogram',...
    'Position',[PData(2).Origin(1)+PData(2).Size(2)/2,0,P(2).startDepth],...
    'width', regionW,...
    'height',regionH,...
    'angle',dopAngle);
PData(2).Region = computeRegions(PData(2));

% Specify Media object.
pt1;
Media.attenuation = -0.5;
Media.function = 'movePoints';

% Specify Resources.
% - RcvBuffer(1) is for both 2D and Doppler acquisitions.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 4096*(2*na + ne);
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 20;           % 20 frames allocated for RF acqusitions.
% InterBuffer(1) is for 2D reconstructions.
Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).numFrames = 1;          % one intermediate frame needed for 2D.
% InterBuffer(2) is for Doppler reconstructions.
Resource.InterBuffer(2).datatype = 'complex';
Resource.InterBuffer(2).numFrames = 1;          % one intermediate frame needed for Doppler.
Resource.InterBuffer(2).pagesPerFrame = ne;     % ne pages per ensemble
% ImageBuffer(1) is for 2D image.
Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).numFrames = 20;
% ImageBuffer(2) is for Doppler image.
Resource.ImageBuffer(2).datatype = 'double';    % image buffer for Doppler
Resource.ImageBuffer(2).numFrames = 20;
% DisplayWindow is for 2D combined with Doppler
Resource.DisplayWindow(1).Title = 'Error4';
Resource.DisplayWindow(1).pdelta = 0.4;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)]; % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).numFrames = 20;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = grayscaleCFImap;
Resource.DisplayWindow(1).splitPalette = 1;

% ------Specify structures used in Events------
% Specify Transmit waveform structures.
% - 2D transmit waveform
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,0.67,2,1];   % A, B, C, D
% - Doppler transmit waveform, transmit frequency should be equivalent to
% supported demodFrequency
TW(2).type = 'parametric';
TW(2).Parameters = [7.813,0.67,6,1];   % A, B, C, D

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'aperture', 1, ...
                   'Apod', ones(1,Resource.Parameters.numTransmit), ...
                   'focus', 0.0, ...
                   'Steer', [0.0,0.0], ...
                   'Delay', zeros(1,Resource.Parameters.numTransmit)), 1, 2*na+1);
% - Set event specific TX attributes.
for n = 1:2:2*na   % 2*na transmit events
    angle = startAngle+((n+1)/2-1)*dtheta;
    TX(n).Steer = [angle,0.0];
    TX(n).aperture = 1; % Use the tx aperture that starts at element 1.
    TX(n).Delay = computeTXDelays(TX(n),'TOAE'); % use 'TransmitOnAllElements' flag
    TX(n+1).Steer = [angle,0.0];
    TX(n+1).aperture = 65; % Use the tx aperture that starts at element 65.
    TX(n+1).Delay = computeTXDelays(TX(n+1),'TOAE');
end
% -- only one TX struct needed for Doppler, middle 96 elements are used for
% Doppler beam transmission
TX(2*na+1).waveform = 2;
TX(2*na+1).aperture = 33;
TX(2*na+1).Steer = [dopAngle,0.0];
TX(2*na+1).Delay = computeTXDelays(TX(2*na+1));

% Specify TPC structures.
TPC(1).name = '2D';
TPC(1).maxHighVoltage = 50;
TPC(2).name = 'Doppler';
TPC(2).maxHighVoltage = 35;

RcvProfile(1).LnaGain = 18; % Profile used for imaging
RcvProfile(2).LnaGain = 24; % Profile used for Doppler
RcvProfile(2).LnaZinSel = 31; % Force high-Z state for best Doppler sensitivity

% Specify TGC Waveform structure.
% - 2D TGC
TGC(1).CntrlPts = [139,535,650,710,770,932,992,1012];
TGC(1).rangeMax = P(1).endDepth;
TGC(1).Waveform = computeTGCWaveform(TGC(1));
% - Doppler TGC
TGC(2).CntrlPts =[ 0 272 662 662 662 662 662 662];
TGC(2).rangeMax = P(2).endDepth;
TGC(2).Waveform = computeTGCWaveform(TGC(2));

% Specify Receive structure arrays -
%   endDepth - add additional acquisition depth to account for some channels
%              having longer path lengths.
%   InputFilter - The same coefficients are used for all channels. The
%              coefficients below give a broad bandwidth bandpass filter.
maxAcqLngth2D = sqrt(P(1).endDepth^2 + (Trans.numelements*Trans.spacing)^2) - P(1).startDepth;
maxAcqLngthDop =  sqrt(P(2).endDepth^2 + (128*Trans.spacing)^2) - P(2).startDepth;
samplesPerWaveDop = 2*TW(2).Parameters(1)/Trans.frequency;

wl4sPer128 = 128/(4*2);  % wavelengths in a 128 sample block for 4 smpls per wave round trip.
wl2sPer128 = 128/(2*samplesPerWaveDop);  % wavelengths in a 128 sample block for 2 smpls per wave round trip.

Receive = repmat(struct('Apod', zeros(1,Resource.Parameters.numTransmit), ...
                        'startDepth', P(1).startDepth, ...
                        'endDepth', P(1).startDepth + wl4sPer128*ceil(maxAcqLngth2D/wl4sPer128), ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ... % 200% Bandwidth for 2D
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, (2*na+ne)*Resource.RcvBuffer(1).numFrames);
% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames  % 2*na acquisitions per frame
    k = (2*na+ne)*(i-1);
    Receive(k+1).callMediaFunc = 1;
    % -- 1st synthetic aperture acquisition for aperture 1.
    for j = 1:2:2*na; % for each aperture acquire all angles
        Receive(k+j).Apod(1:96) = 1.0;
        Receive(k+j).aperture = 1;
        Receive(k+j).framenum = i;
        Receive(k+j).acqNum = j;
    % -- 2nd synthetic aperture acquisition for aperture 65.
        Receive(k+j+1).Apod(33:128) = 1.0;
        Receive(k+j+1).aperture = 65;
        Receive(k+j+1).framenum = i;
        Receive(k+j+1).acqNum = j+1;
    end
    for j = (2*na+1):(2*na+ne)
        Receive(j+k).Apod(1:128) = 1;
        Receive(j+k).aperture = 33;
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
Recon = repmat(struct('senscutoff', 0.7, ...
               'pdatanum', 1, ...
               'rcvBufFrame', -1, ...     % use most recently transferred frame
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...  % auto-increment ImageBuffer each recon
               'RINums', zeros(1,1)), 1, 2);
% - Set Recon values for 2D frame.
Recon(1).RINums(1,1:2*na) = 1:2*na;   % 3*na ReconInfos needed for na angles
k = 2*na + 1;
% - Set Recon values for Doppler ensemble.
Recon(2).pdatanum = 2;
Recon(2).IntBufDest = [2,1];
Recon(2).ImgBufDest = [0,0];
Recon(2).RINums(1,1:ne) = k:(k+ne-1);   % ne ReconInfos needed for Doppler ensemble.

% Define ReconInfo structures.
% - For 2D, we need na ReconInfo structures for na steering angles.
% - For Doppler, we need ne ReconInfo structures.
ReconInfo = repmat(struct('mode', 'accumIQ', ...    % accumulate IQ data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'pagenum',1, ...
                   'regionnum', 1), 1, 2*na + ne);
% - Set specific ReconInfo attributes.
%   - ReconInfos for 2D frame.
if na>1
    ReconInfo(1).mode = 'replaceIQ'; % replace IQ data
    for j = 1:2*na
        ReconInfo(j).txnum = j;
        ReconInfo(j).rcvnum = j;
    end
    ReconInfo(2*na).mode = 'accumIQ_replaceIntensity'; % accum and detect
else
    ReconInfo(1).mode = 'replaceIntensity';
end

%  - ReconInfos for Doppler ensemble.
k = 2*na;
for j = 1:ne
    ReconInfo(k+j).mode = 'replaceIQ';
    ReconInfo(k+j).txnum = 2*na + 1;
    ReconInfo(k+j).rcvnum = 2*na + j;
    ReconInfo(k+j).pagenum = j;
end

% Specify Process structure arrays.
cpt = 28;  % define here so we can use in UIControl below
persf = 80;
persp = 90;
DopState = 'freq';

Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',1.0,...            % pgain is image processing gain
                         'reject',2,...      % reject level
                         'persistMethod','simple',...
                         'persistLevel',20,...
                         'interpMethod','4pt',...  %method of interp. (1=4pt)
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',40,...
                         'mappingMethod','lowerHalf',...
                         'display',0,...      % don't display image after processing
                         'displayWindow',1};

Process(2).classname = 'Doppler';                   % process structure for 1st Doppler ensemble
Process(2).method = 'computeCFIFreqEst';
Process(2).Parameters = {'IntBufSrc',[2,1],...          % number of (Inter?)buffer to process.
                         'SrcPages',[3,ne-2],...        % start frame number in source buffer
                         'ImgBufDest',[2,-1],...
                         'pdatanum',2,...           % number of PData structure
                         'prf',dopPRF,...           % Doppler PRF in Hz
                         'wallFilter','regression',...  % 1 -> quadratic regression
                         'pwrThreshold',pwrThres,...
                         'maxPower',50,...
                         'postFilter',1};

Process(3).classname = 'Image';
Process(3).method = 'imageDisplay';
Process(3).Parameters = {'imgbufnum',2,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'srcData','signedColor',... % type of data to display.
                         'pdatanum',2,...    % number of PData structure to use
                         'pgain',1.0,...            % pgain is image processing gain
                         'reject',2,...      % reject level
                         'persistMethod','dynamic',...
                         'persistLevel',persf,...
                         'interpMethod','4pt',...  %method of interp. (1=4pt)
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',40,...
                         'mappingMethod','upperHalf',...
                         'threshold',cpt,...
                         'display',1,...      % display image after processing
                         'displayWindow',1};

% EF1 is external function for ROI plot
Process(4).classname = 'External';
Process(4).method = 'ROIplot';
Process(4).Parameters = {'srcbuffer','none'};

dopDuration = round(1/(dopPRF*1e-06))*ne;
switchTime = (20000 - dopDuration - 200*na*2)/2; % Used to keep 50 frames/second

% Specify SeqControl structure arrays.
% -- Time between 2D flash angle acquisitions
SeqControl(1).command = 'timeToNextAcq';
SeqControl(1).argument = 200;  % time in usec (2.8 msecs for na = 7)
% -- Change to Profile 2 (Doppler)
SeqControl(2).command = 'setTPCProfile';
SeqControl(2).condition = 'next';
SeqControl(2).argument = 2;
% -- Time between 2D acquisition and Doppler ensemble. Set to allow time for profile change.
SeqControl(3).command = 'timeToNextAcq';
SeqControl(3).argument = switchTime; % time in usec
% -- PRF for Doppler ensemble
SeqControl(4).command = 'timeToNextAcq';
SeqControl(4).argument = round(1/(dopPRF*1e-06)); % (for 3KHz dopPRF & 14 ensemble = 4.7 msecs)
% -- Change to Profile 1 (2D)
SeqControl(5).command = 'setTPCProfile';
SeqControl(5).condition = 'next';
SeqControl(5).argument = 1;
% -- Time between Doppler and next 2D acquisition. Set to allow time for profile change.
SeqControl(6).command = 'timeToNextAcq';
SeqControl(6).argument = switchTime; % time in usec
% -- Jump back to start.
SeqControl(7).command = 'jump';
SeqControl(7).argument = 2;
% set receive profile
SeqControl(8).command = 'setRcvProfile';
SeqControl(8).argument = 1;  %
SeqControl(9).command = 'setRcvProfile';
SeqControl(9).argument = 2;  %
% - The remainder of the SeqControl structures are defined dynamically in the sequence events.
%   The variable nsc keeps track of the next SeqControl structure index.
nsc = 10;  % next SeqControl number

% Specify Event structure arrays.
n = 1;

Event(n).info = 'ROI plot';
Event(n).tx = 0;         % no transmit
Event(n).rcv = 0;        % no rcv
Event(n).recon = 0;      % no reconstruction
Event(n).process = 4;    % outline doppler region
Event(n).seqControl = 0; % no seqCntrl
n = n+1;

for i = 1:Resource.RcvBuffer(1).numFrames
    k = (2*na+ne)*(i-1);
    for j = 1:2:2*na
        Event(n).info = '1st half of aperture.';
        Event(n).tx = j;         % use TX structure of 1st aperture.
        Event(n).rcv = k+j;      % use Rcv structure of 1st aperture.
        Event(n).recon = 0;      % no reconstruction.
        Event(n).process = 0;    % no processing
        Event(n).seqControl = 1; % time between flash angles
        if j == 1
            Event(n).seqControl = [1,8]; % Set correct RcvProfile at the 1st Doppler transmit
        end
        n = n+1;

        Event(n).info = '2nd half of aperture.';
        Event(n).tx = j+1;       % use TX structure of 2nd aperture
        Event(n).rcv = k+j+1;    % use Rcv structure of 2nd aperture
        Event(n).recon = 0;      % no reconstruction.
        Event(n).process = 0;    % no processing
        Event(n).seqControl = 1; % time between syn. aper. acqs.
        n = n+1;
    end
    Event(n-1).seqControl = [2,3];   % replace last 2D acquisition Event's seqControl
    % Acquire Doppler ensemble.
    for j = (2*na+1):(2*na+ne)
        Event(n).info = 'Acquire Doppler ensemble';
        Event(n).tx = 2*na+1;      % use next TX structure after 2D.
        Event(n).rcv = k+j;
        Event(n).recon = 0;      % no reconstruction.
        Event(n).process = 0;    % no processing
        Event(n).seqControl = 4; % TTNA for Doppler
        if j == 2*na+1
            Event(n).seqControl = [4,9]; % Set correct RcvProfile at the 1st Doppler transmit
        end
        n = n+1;
    end
    Event(n-1).seqControl = [5,6,nsc]; % replace last Doppler acquisition Event's seqControl
      SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
      nsc = nsc+1;

    Event(n).info = 'recons and 2D process';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = [1,2];  % reconstruction for both 2D and Doppler
    Event(n).process = 1;    % process 2D
    Event(n).seqControl = 0; % no seqCntrl
    n = n+1;

    Event(n).info = 'Doppler processing';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 0;      % no reconstruction
    Event(n).process = 2;    % IQ to estimate processing
    Event(n).seqControl = 0; % no seqCntrl
    n = n+1;

    Event(n).info = 'Doppler image display';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 0;      % no reconstruction
    Event(n).process = 3;    % display processing
    Event(n).seqControl = 0; % no seqCntrl
    if floor(i/3) == i/3     % Exit to Matlab every 3rd frame
        Event(n).seqControl = nsc;
        SeqControl(nsc).command = 'returnToMatlab';
        nsc = nsc+1;
    end
    n = n+1;

end

Event(n).info = 'Jump back';
Event(n).tx = 0;        % no TX
Event(n).rcv = 0;       % no Rcv
Event(n).recon = 0;     % no Recon
Event(n).process = 0;
Event(n).seqControl = 7;


% User specified UI Control Elements
% - Sensitivity Cutoff
UI(1).Control =  {'UserB7','Style','VsSlider','Label','Sens. Cutoff',...
                  'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                  'SliderStep',[0.025,0.1],'ValueFormat','%1.3f'};
UI(1).Callback = text2cell('%SensCutoffCallback');

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

% - Replot ROI for doppler
UI(6).Control =  {'UserC4','Style','VsButtonGroup','Title','Doppler display',...
    'NumButtons',2,'Labels',{'on','off'}};
UI(6).Callback = text2cell('%-UI#6Callback');

% - PRF adjustment
PRFmin = 250; PRFmax = 4250; stepDiff = PRFmax-PRFmin;
UI(7).Control = {'UserC3','Style','VsSlider','Label','PRF (Hz)','SliderMinMaxVal',[PRFmin,PRFmax,dopPRF],...
                 'SliderStep',[50/stepDiff,500/stepDiff],'ValueFormat','%3.0f'};
UI(7).Callback = text2cell('%-UI#7Callback');

% - Steer Angle adjustment
UI(8).Control = {'UserC2','Style','VsSlider','Label','Doppler Angle','SliderMinMaxVal',[-20,20,round(dopAngle*180/pi)],...
                 'SliderStep',[1/40,5/40],'ValueFormat','%3.0f'};
UI(8).Callback = text2cell('%-UI#8Callback');

% - Range Change
MinMaxmm = [5,85]; % min max in mm
AxesUnit = 'wls';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
        AxesUnit = 'mm';
        MinMaxVal = [MinMaxmm, P(2).endDepth * (Resource.Parameters.speedOfSound/1000/Trans.frequency)];
    else
        MinMaxVal = [MinMaxmm * (Trans.frequency/Resource.Parameters.speedOfSound*1000), P(2).endDepth];
    end
end
stepBase = MinMaxVal(2)-MinMaxVal(1);
UI(9).Control = {'UserA1','Style','VsSlider','Label',['Range (',AxesUnit,')'],'SliderMinMaxVal',MinMaxVal,...
                 'SliderStep',[5/stepBase,10/stepBase],'ValueFormat','%3.0f'};
UI(9).Callback = text2cell('%-UI#9Callback');

% - External function for ROIplot
EF(1).Function = vsv.seq.function.ExFunctionDef('ROIplot', @ROIplot);
% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 3;

% Save all the structures to a .mat file.
filename = 'MatFiles/Error4';
save(filename);
% VSX

return

% **** Callback routines to be encoded by text2cell function. ****
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

%-UI#2Callback - Doppler mode change
Control = repmat(struct('Command','set&Run','Parameters',[]),1,4);
Process = evalin('base','Process');
Resource = evalin('base','Resource');
hDisplay = Resource.DisplayWindow(1).figureHandle;
currentMap = get(hDisplay,'Colormap');

switch UIState
   case 1  % Velocity mode
      newMap = grayscaleCFImap;
      newMap(1:128,:) = currentMap(1:128,:);
      assignin('base','persp',get(findobj('Tag','UserB1Slider'),'Value'));
      persf = evalin('base','persf'); persValue = persf(1);
      Control(1).Parameters = {'Process',2,'method','computeCFIFreqEst'};
      Control(2).Parameters = {'Process',3,'srcData','signedColor','persistMethod','dynamic','persistLevel',persValue};
      Control(3).Parameters = {'DisplayWindow',1,'colormap',newMap};
      Control(4).Parameters = {'ImageBuffer',1,'lastFrame',0};
      set(findobj('tag','UserB1Edit'),'String',num2str(persValue,'%3.0f'));
      set(findobj('tag','UserB1Slider'),'Value',persValue);
      assignin('base','DopState','freq');
      % Set modified Process attributes in base Matlab environment.
      Process(2).method = 'computeCFIFreqEst';
      for k = 1:2:length(Process(3).Parameters)
          if strcmp(Process(3).Parameters{k},'srcData'), Process(3).Parameters{k+1} = 'signedColor';
          elseif strcmp(Process(3).Parameters{k},'persistMethod'), Process(3).Parameters{k+1} = 'dynamic';
          elseif strcmp(Process(3).Parameters{k},'persistLevel'), Process(3).Parameters{k+1} = persf;
          end
      end
      assignin('base','Process',Process);
   case 2  % Power mode
      newMap = grayscaleCPAmap;
      newMap(1:128,:) = currentMap(1:128,:);
      for k = 1:2:length(Process(3).Parameters)
          if strcmp(Process(3).Parameters{k},'persistLevel'), persf = Process(3).Parameters{k+1}; end
      end
      assignin('base','persf',persf);
      persValue = evalin('base','persp');
      Control(1).Parameters = {'Process',2,'method','computeCFIPowerEst'};
      Control(2).Parameters = {'Process',3,'srcData','unsignedColor','persistMethod','simple','persistLevel',persValue};
      Control(3).Parameters = {'DisplayWindow',1,'colormap',newMap};
      Control(4).Parameters = {'ImageBuffer',1,'lastFrame',0};
      set(findobj('tag','UserB1Edit'),'String',num2str(persValue,'%3.0f'));
      set(findobj('tag','UserB1Slider'),'Value',persValue);
      assignin('base','DopState','power');
      Process(2).method = 'computeCFIPowerEst';
      for k = 1:2:length(Process(3).Parameters)
          if strcmp(Process(3).Parameters{k},'srcData'), Process(3).Parameters{k+1} = 'unsignedColor';
          elseif strcmp(Process(3).Parameters{k},'persistMethod'), Process(3).Parameters{k+1} = 'simple';
          elseif strcmp(Process(3).Parameters{k},'persistLevel'), Process(3).Parameters{k+1} = persValue;
          end
      end
      assignin('base','Process',Process);
end
assignin('base','Control', Control);

% If PTool window is open, adjust all uicontrols
hPTool = findobj('tag','ProcessTool');
if ishandle(hPTool),
    posPTool = get(hPTool,'position');
    PTool;
    set(findobj('tag','ProcessTool'),'position',posPTool);
end
return
%-UI#2Callback

%-UI#3Callback - Doppler Power change
Process = evalin('base','Process');
for k = 1:2:length(Process(2).Parameters)
    if strcmp(Process(2).Parameters{k},'pwrThreshold'), Process(2).Parameters{k+1} = UIValue; end
end
assignin('base','Process',Process);
% Set Control.Command to set Doppler threshold.
Control = evalin('base','Control');
Control.Command = 'set&Run';
Control.Parameters = {'Process',2,'pwrThreshold',UIValue};
assignin('base','Control', Control);
%-UI#3Callback

%-UI#4Callback - Color Priority change
% Set the value in the Process structure for use in cineloop playback.
Process = evalin('base','Process');
for k = 1:2:length(Process(3).Parameters)
    if strcmp(Process(3).Parameters{k},'threshold'), Process(3).Parameters{k+1} = UIValue; end
end
assignin('base','Process',Process);
% Set Control.Command to set Image.threshold.
Control = evalin('base','Control');
Control.Command = 'set&Run';
Control.Parameters = {'Process',3,'threshold',UIValue};
assignin('base','Control', Control);
%-UI#4Callback

%-UI#5Callback - Color Persistence change
% Set the value in the Process structure for use in cineloop playback.
Process = evalin('base','Process');
for k = 1:2:length(Process(3).Parameters)
    if strcmp(Process(3).Parameters{k},'persistLevel'), Process(3).Parameters{k+1} = UIValue; end
end
assignin('base','Process',Process);
% Set Control.Command to set Image.persistLevel.
Control = evalin('base','Control');
Control.Command = 'set&Run';
Control.Parameters = {'Process',3,'persistLevel',UIValue};
assignin('base','Control', Control);

% If PTool window is open, adjust persistLevel1 in Process(3)
hPTool = findobj('tag','ProcessTool');
if ishandle(hPTool),
    hPNum = findobj('tag','processNum');
    if isequal(get(findobj('tag','processNum'),'Value'),3)
        set(findobj('tag','persistSlider1'),'Value',UIValue);
        set(findobj('tag','persistValue1'),'String',num2str(UIValue));
    end
end
return
%-UI#5Callback

%-UI#6Callback
Process = evalin('base','Process');
if UIState == 1
    evalin('base','set(ROIHandle,''Visible'',''on'');');
    for k = 1:2:length(Process(1).Parameters)
        if strcmp(Process(1).Parameters{k},'display'), Process(1).Parameters{k+1} = 0; end
    end
    for k = 1:2:length(Process(3).Parameters)
        if strcmp(Process(3).Parameters{k},'display'), Process(3).Parameters{k+1} = 1; end
    end
    % Set Control.Command to set Image.persistLevel.
    Control = evalin('base','Control');
    Control(1).Command = 'set&Run';
    Control(2).Command = 'set&Run';
    Control(1).Parameters = {'Process',1,'display',0};
    Control(2).Parameters = {'Process',3,'display',1};
else
    evalin('base','set(ROIHandle,''Visible'',''off'');');
    for k = 1:2:length(Process(1).Parameters)
        if strcmp(Process(1).Parameters{k},'display'), Process(1).Parameters{k+1} = 1; end
    end
    for k = 1:2:length(Process(3).Parameters)
        if strcmp(Process(3).Parameters{k},'display'), Process(3).Parameters{k+1} = 0; end
    end
    % Set Control.Command to set Image.persistLevel.
    Control = evalin('base','Control');
    Control(1).Command = 'set&Run';
    Control(2).Command = 'set&Run';
    Control(1).Parameters = {'Process',1,'display',1};
    Control(2).Parameters = {'Process',3,'display',0};
end
assignin('base','Control', Control);
return
%-UI#6Callback

%-UI#7Callback - PRF
dopPRF = UIValue;
na = evalin('base','na');
ne = evalin('base','ne');
dopDuration = round(1/(dopPRF*1e-06))*ne;
switchTime = (20000 - dopDuration - 200*na*2)/2; % Used to keep 50 frames/second

SeqControl = evalin('base','SeqControl');
SeqControl(3).argument = switchTime;
SeqControl(4).argument = round(1/(dopPRF*1e-06));
SeqControl(6).argument = switchTime;

assignin('base','dopPRF',dopPRF);
assignin('base','SeqControl',SeqControl);

Control.Command = 'update&Run';
Control.Parameters = {'SeqControl'};
assignin('base','Control',Control);

return
%-UI#7Callback

%-UI#8Callback - Dopler Angle will change TX beam and PData
dopAngle = UIValue * pi/180;
assignin('base','dopAngle',dopAngle);
PData = evalin('base','PData');
PData(2).Region.Shape.angle = dopAngle;
PData(2).Region = computeRegions(PData(2));
assignin('base','PData',PData);

TX = evalin('base','TX');
na = evalin('base','na');
TX(2*na+1).Steer = [dopAngle,0.0];
TX(2*na+1).Delay = computeTXDelays(TX(2*na+1));
assignin('base','TX',TX);

Control.Command = 'update&Run';
Control.Parameters = {'PData','ImageBuffer','Recon','TX'};
assignin('base','Control',Control);

return
%-UI#8Callback

%-UI#9Callback
simMode = evalin('base','Resource.Parameters.simulateMode');
% No range change if in simulate mode 2.
if simMode == 2
    set(hObject,'Value',evalin('base','P.endDepth'));
    return
end

% P.endDepth
Trans = evalin('base','Trans');
Resource = evalin('base','Resource');
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);
P = evalin('base','P');
P(1).endDepth = UIValue;
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
        P(1).endDepth = UIValue*scaleToWvl;
    end
end
P(2) = P(1);
assignin('base','P',P);

% PData
PData = evalin('base','PData');
dopAngle = evalin('base','dopAngle');
PData(1).Size(1) = ceil((P(1).endDepth-P(1).startDepth)/PData(1).PDelta(3)); % rows
PData(2) = PData(1); % default is the same as PData(1), but the user may want to change the PDelta
regionH = P(2).endDepth-P(2).startDepth;
regionW = ceil((128*Trans.spacing)/PData(2).PDelta(1)); % 128 elements for Doppler
PData(2).Region.Shape = struct(...
    'Name','Parallelogram',...
    'Position',[PData(2).Origin(1)+PData(2).Size(2)/2,0,P(2).startDepth],...
    'width', regionW,...
    'height',regionH,...
    'angle',dopAngle);
PData(1).Region = computeRegions(PData(1));
PData(2).Region = computeRegions(PData(2));
assignin('base','PData',PData);
evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);');

% Receive
Receive = evalin('base', 'Receive');
maxAcqLngth2D = sqrt(P(1).endDepth^2 + (Trans.numelements*Trans.spacing)^2) - P(1).startDepth;
maxAcqLngthDop =  sqrt(P(2).endDepth^2 + (96*Trans.spacing)^2) - P(2).startDepth;
wl4sPer128 = 128/(4*2);  % wavelengths in a 128 sample block for 4 smpls per wave round trip.
wl2sPer128 = 128/(2*2);  % wavelengths in a 128 sample block for 2 smpls per wave round trip.
na = evalin('base','na');
ne = evalin('base','ne');
for i = 1:Resource.RcvBuffer(1).numFrames
    k = (na + ne)*(i-1); % k keeps track of Receive index increment per frame.
    for j = 1:2:2*na  % acquisitions for 2D
        Receive(j+k).endDepth = P(1).startDepth + wl4sPer128*ceil(maxAcqLngth2D/wl4sPer128);
        Receive(j+k+1).endDepth = P(1).startDepth + wl4sPer128*ceil(maxAcqLngth2D/wl4sPer128);
    end
    for j = (2*na+1):(2*na+ne); % Doppler acquisition
        Receive(j+k).endDepth = P(2).startDepth + wl2sPer128*ceil(maxAcqLngthDop/wl2sPer128);
    end
end
assignin('base','Receive',Receive);

% TGC
TGC = evalin('base','TGC');
TGC(1).rangeMax = P(1).endDepth;
TGC(1).Waveform = computeTGCWaveform(TGC(1));
TGC(2).rangeMax = P(2).endDepth;
TGC(2).Waveform = computeTGCWaveform(TGC(2));
assignin('base','TGC',TGC);

Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'PData','InterBuffer','ImageBuffer','DisplayWindow','Receive','TGC','Recon'};
assignin('base','Control', Control);
assignin('base', 'action', 'displayChange');

return
%-UI#9Callback

% External Function
function ROIplot(varargin)

    persistent ROIHandle

    P = evalin('base','P');
    UI = evalin('base','UI');
    UIPos = evalin('base','UIPos');
    Trans = evalin('base','Trans');
    PData = evalin('base','PData');
    dopAngle = evalin('base','dopAngle');
    Resource = evalin('base','Resource');

    % change the txt for voltage slider
    hv1txt = findobj('tag','hv1txt');
    hv2txt = findobj('tag','hv2txt');
    set(hv1txt,'String','Bmode Voltage');
    set(hv2txt,'String','Doppler Voltage');

    % Outline the PData region for doppler processing
    X(1) = PData(2).Region.Shape.Position(1)-PData(2).Region.Shape.width/2;
    Y(1) = PData(2).Region.Shape.Position(3);
    if Y(1) < P(1).startDepth, Y(1) = P(1).startDepth; end
    X(2) = X(1) + PData(2).Region.Shape.height*tan(dopAngle);
    X(3) = X(2) + PData(2).Region.Shape.width;
    Y(2) = Y(1) + PData(2).Region.Shape.height;
    Y(3) = Y(2);
    X(4) = X(1) + PData(2).Region.Shape.width;
    Y(4) = Y(1);
    if dopAngle > 0
        if X(3) > PData(1).Origin(1)+PData(1).Size(2)
            X(5) = X(4); Y(5) = Y(4);
            X(3) = PData(1).Origin(1)+PData(1).Size(2);
            X(4) = X(3);
            Y(4) = 32*Trans.spacing/tan(dopAngle);
        end
    else
        if X(2) < PData(1).Origin(1)
            X(5) = X(4); Y(5) = Y(4);
            X(4) = X(3); Y(4) = Y(3);
            X(3) = PData(1).Origin(1);
            X(2) = X(3);
            Y(2) = 32*Trans.spacing/tan(-dopAngle);
        end
    end

    ROIX = [X X(1)];
    ROIY = [Y Y(1)];

    % scale to mm if the axes unit of the displaywindow is mm
    scaleToMm = (Resource.Parameters.speedOfSound/1000)/Trans.frequency;
    if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
        if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
            ROIX = ROIX * scaleToMm; ROIY = ROIY * scaleToMm;
        end
    end
    bmodeFigHandle = evalin('base','Resource.DisplayWindow(1).figureHandle');
    if ishandle(bmodeFigHandle)
        if  isempty(ROIHandle) || ~ishandle(ROIHandle)
            figure(bmodeFigHandle), hold on,
            ROIHandle = plot(ROIX,ROIY,'Color',[0.7 0.7 0.7],'LineWidth',1.5,'LineStyle','--','Visible','off');
            if get(UI(6).handle(2),'Value') == 1 % 'on' is selected
                set(ROIHandle,'Visible','on');
            end
            assignin('base','ROIHandle',ROIHandle);
        else
            set(ROIHandle,'XData',ROIX);
            set(ROIHandle,'YData',ROIY);
        end
    end
end

