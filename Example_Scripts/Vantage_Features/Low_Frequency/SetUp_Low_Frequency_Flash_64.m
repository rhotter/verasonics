% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use.
%
% File name: SetUp_Low_Freq_Flash - Example of imaging using low
% frequency system configuration with a fictitious 64 element array.
% Frequencies can be selected from low (50kHz and 500kHz) and standard
% (3 MHz), demonstrating the configuration differences between low
% and standard frequency systems.
%
% Description:
%   Sequence programming file for a fictitious 64 element array, single
%   flash angle transmission. For each acquisition, all channels are
%   active for transmit and receive. Reconstruction and processing can
%   be set to synchronous with respect to acquisition with
%   Resource.Parameters.waitForProcessing = 1.
%
% Last update:
% 11/10/2020 - Merged low frequency (50 kHz and 500 kHz) and standard
%   frequency (3 MHz) scripts into a single script.
% 05/06/2020 - Update to SW 4.3 format for new user UIControls and External function definitions (VTS 1691).
%   More info:(?/Example_Scripts/Vantage_Features/New UI Scheme/SetUpL11_5vFlash_NewUI)

clear all

% Frequency selection: use one of the following options
% 'low50k'  for  50 kHz low frequency
% 'low500k' for 500 kHz low frequency
% 'std3m'   for   3 MHz standard frequency
frequency_selection = 'low500k';

P.startDepth = 1;   % Acquisition depth in wavelengths
P.endDepth = 80;    % This should preferrably be a multiple of 128 samples.
P.HV = 25;          % preset voltage
P.Elnum = 64;       % num of el in transmit aperture, must be no more than Trans.numelements

input_filter_on = 1; % digital filter. 1 - default band-pass, 0 - all-pass.
TPCProfile = 5;      % Longer pulse length may necesitate use of TPCProfile 5 at low frequency

% Specify analogue filters for low-frequency acquisition
switch lower(frequency_selection)
    case 'low50k'
        freq_c = 0.05;
        Resource.System.Frequency = 'Low';
        RcvProfile.LnaHPF = 0; % Disable LNA pass high pass filter
        RcvProfile.PgaHPF = 0; % Disable PGA high pass filter (other permissible value is 80 for 80kHz breakpoint)
    case 'low500k'
        freq_c = 0.5;
        Resource.System.Frequency = 'Low';
        RcvProfile.LnaHPF = 100; % LNA high pass filter, 100kHz break frequency
    case 'std3m'
        freq_c = 3;
        Resource.System.Frequency = 'Standard';
        RcvProfile.LnaHPF=200;   % LNA high pass filter, 200kHz break frequency
    otherwise
        fprintf("unrecognized frequency selection: %s\n", frequency_selection)
        return
end

if freq_c < 1
    freq_c_disp = sprintf('%dkHz', freq_c*1e3);
else
    freq_c_disp = sprintf('%dMHz', freq_c);
end

% Specify system parameters.
Resource.Parameters.numTransmit = P.Elnum;     % number of transmit channels.
Resource.Parameters.numRcvChannels = P.Elnum;  % number of receive channels.
Resource.Parameters.speedOfSound = 1540;
Resource.Parameters.waitForProcessing = 1;     % DMA transfers wait for processing.
Resource.Parameters.Connector = 1;
Resource.Parameters.fakeScanhead = 1;          % allows system to run without querying user if no probe is connected
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Specify Trans structure array.
Trans.name = sprintf('Fictitious-%s-64el', freq_c_disp);
Trans.units = 'mm';                                 % required in Gen3 to prevent default to mm units
Trans.numelements = 64;
Trans.frequency = freq_c;                           % nominal frequency in MHz
freq_bw = 1/3;
Trans.Bandwidth = [Trans.frequency*(1-freq_bw), Trans.frequency*(1+freq_bw)];
Trans.type = 0;                                     % Array geometry is linear (x values only).
Trans.id = -1;                                      % -1 means no ID; skip probe ID checks
Trans.connType = -1;                                % -1 means adapt to the connector and UTA that is present (if possible)
Trans.elevationApertureMm = 4;                      % active elevation aperture in mm
Trans.elevationFocusMm = 1000000;                   % nominal elevation focus depth from lens on face of transducer
Trans.spacingMm = Resource.Parameters.speedOfSound/(Trans.frequency*1e3)/2;   % Spacing between elements in mm (half wavelength).
Trans.elementWidth = Trans.spacingMm-50e-3;                          % element width in mm; assumes 50 micron kerf

Trans.ElementPos = zeros(Trans.numelements,4);
Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
% Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
Theta = (-pi/2:pi/100:pi/2);
Theta(51) = 0.0000001;                          % set to almost zero to avoid divide by zero.
eleWidthWl = Trans.elementWidth * Trans.frequency/Resource.Parameters.speedOfSound;
Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
Trans.lensCorrection = 0;
Trans.impedance = 50;
Trans.maxHighVoltage = 50;

% Now convert all units as required, based on Trans.units
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000); % conversion factor from mm to wavelengths
% regardless of units, always provide spacing in wavelengths
Trans.spacing = Trans.spacingMm * scaleToWvl;   % Spacing between elements in wavelengths

P.theta = -pi/4;
P.rayDelta = 2*(-P.theta);
P.aperture = P.Elnum*Trans.spacing; % P.aperture in wavelengths
P.radius = (P.aperture/2)/tan(-P.theta); %f dist. to virt. apex

% Set up PData structure.
PData(1).PDelta = [0.5, 0, 0.5];
PData(1).Size(1) = 10 + ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));
PData(1).Size(2) = 10 + ceil(2*(P.endDepth + P.radius)*sin(-P.theta)/PData(1).PDelta(1));
PData(1).Size(3) = 1;
PData(1).Origin = [-(PData(1).Size(2)/2)*PData(1).PDelta(1), 0, P.startDepth];
PData(1).Region = struct(...
            'Shape',struct('Name','SectorFT', ...
            'Position',[0,0,-P.radius], ...
            'z',P.startDepth, ...
            'r',P.radius+P.endDepth, ...
            'angle',P.rayDelta, ...
            'steer',0));
PData(1).Region = computeRegions(PData(1));

% Specify Media.  Use point targets in middle of PData.
% Set up Media points
Media.MP(1,:) = [-45,0,30,1.0];
Media.MP(2,:) = [-15,0,30,1.0];
Media.MP(3,:) = [15,0,30,1.0];
Media.MP(4,:) = [45,0,30,1.0];
Media.MP(5,:) = [-15,0,60,1.0];
Media.MP(6,:) = [-15,0,90,1.0];
Media.MP(7,:) = [-15,0,120,1.0];
Media.MP(8,:) = [-15,0,150,1.0];
Media.MP(9,:) = [-45,0,120,1.0];
Media.MP(10,:) = [15,0,120,1.0];
Media.numPoints = 10;
Media.function = 'movePoints';

% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 4*2560;  % Sized to accomodate 25 samples per wavelength at 50 kHz fc
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 10;     % 10 frames for RF data.
Resource.InterBuffer(1).numFrames = 1;  % one intermediate buffer needed.
Resource.ImageBuffer(1).numFrames = 10;
Resource.DisplayWindow(1).Title = sprintf('%s frequency configuration: %s single flash angle', Resource.System.Frequency, freq_c_disp);
Resource.DisplayWindow(1).pdelta = PData(1).PDelta(1)/2;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).numFrames = 20;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = gray(256);

% Specify Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,.67,2,1];

ElOffset = floor((Trans.numelements - P.Elnum)/2); if ElOffset<0, error('reset P.Elnum to a value smaller than Trans.numelements.'); end

% Specify TX structure array. Transmit on element n in the array for event n.
TX = repmat(struct('waveform', 1, ...
                   'Origin', zeros(1,3), ...
                   'focus', -P.radius, ...
                   'Steer', [0.0,0.0], ...
                   'Apod', [zeros(1,ElOffset), ones(1,P.Elnum), zeros(1,Trans.numelements-P.Elnum-ElOffset)], ...
                   'Delay', zeros(1,min([Trans.numelements, Resource.Parameters.numTransmit]))), 1, 1);

% - Set event specific TX attributes.
TX(1).Delay = computeTXDelays(TX(1));

% Specify TPC structure.
TPC(1).hv = P.HV;

if TPCProfile==5
    TPC(5).hv = P.HV;
    TPC(5).inUse=1;
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [478, 602, 703, 812, 906, 973, 985, 1022];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Receive structure arrays.
% - Receive structure for each frame.
maxAcqLength = ceil(sqrt(P.aperture^2 + P.endDepth^2 - 2*P.aperture*P.endDepth*cos(P.theta-pi/2)) - P.startDepth);
wlsPer128 = 128/(2*4); % wavelengths in 128 samples for 4 samplesPerWave

Receive = struct('Apod', [zeros(1,ElOffset), ones(1,P.Elnum), zeros(1,Trans.numelements-P.Elnum-ElOffset)], ...
                        'startDepth', P.startDepth, ...
                        'endDepth', P.startDepth + wlsPer128*ceil(maxAcqLength/wlsPer128), ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode','NS200BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 0);

if  input_filter_on == 0
     Receive.InputFilter = [zeros(1,20), 1]; %all-pass
end

Receive = repmat(Receive,1,Resource.RcvBuffer(1).numFrames);


% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    Receive(i).callMediaFunc = 1;
    Receive(i).framenum = i;
    Receive(i).acqNum = 1;
end

% Specify Recon structure arrays.
% - We need one Recon structure.  Each frame will use
%   Single ReconInfo structures, since we are using single flash
%   acquisition.

Recon = struct('senscutoff', 0, ...
               'pdatanum', 1, ...
               'rcvBufFrame', -1, ...     % use most recently transferred frame
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...  % auto-increment ImageBuffer each recon
               'RINums', 1);

% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 'accumIQ', ...  % accumulate IQ data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 1), 1, 1);
% - Set specific ReconInfo attributes.
    ReconInfo(1).mode = 'replaceIntensity';

% Specify Process structure array.
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',1.0,...            % pgain is image processing gain
                         'reject',5,...      % reject level
                         'persistMethod','simple',...
                         'persistLevel',1,...
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
SeqControl(1).command = 'jump'; % jump back to start.
SeqControl(1).argument = 3;
SeqControl(2).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(2).argument = floor(1/30*1e6);  % 30 fps
SeqControl(3).command = 'returnToMatlab';

nsc = 4; % nsc is count of SeqControl objects

% Specify Event structure array.
n = 1; % n is count of Events

Event(n).info = 'select TPC profile';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = nsc;
n = n+1;
  SeqControl(nsc).command = 'setTPCProfile';
  SeqControl(nsc).argument = TPCProfile;
  SeqControl(nsc).condition = 'immediate';
  nsc = nsc + 1;

Event(n).info = 'noop delay';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = nsc;
n = n+1;
  SeqControl(nsc).command = 'noop';
  SeqControl(nsc).argument = 500e3; % 500 msec delay
  nsc = nsc + 1;

for i = 1:Resource.RcvBuffer(1).numFrames
    Event(n).info = 'Acquisition.';
    Event(n).tx = 1;
    Event(n).rcv = i;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = [2 nsc];
    n = n+1;
      SeqControl(nsc).command = 'transferToHost';
      nsc = nsc+1;

    Event(n).info = 'recon and process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 1;
    Event(n).seqControl = 3;
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
    'Callback',@SensCutoffCallback);

% - Range Change
MinMaxVal = [10,150,P.endDepth]; % default unit is wavelength
AxesUnit = 'wls';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
        AxesUnit = 'mm';
        MinMaxVal = MinMaxVal * (Resource.Parameters.speedOfSound/1000/Trans.frequency);
    end
end
UI(2).Control = VsSliderControl('LocationCode','UserA1',...
    'Label',['Range (',AxesUnit,')'],...
    'SliderMinMaxVal',MinMaxVal,'SliderStep',[0.1,0.2],...
    'ValueFormat','%3.0f',...
    'Callback',@RangeChangeCallback);

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 1;

filename = 'Low_Freq_Flash_64';

% Save all the structures to a .mat file.
save(['MatFiles/' filename]);

%% Automatic VSX Execution:
% Uncomment the following line to automatically run VSX every time you run
% this SetUp script (note that if VSX finds the variable 'filename' in the
% Matlab workspace, it will load and run that file instead of prompting the
% user for the file to be used):

% VSX;

clear filename


%% **** Callback routines used by UIControls (UI) ****

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

function RangeChangeCallback(hObject, ~, UIValue)
    simMode = evalin('base','Resource.Parameters.simulateMode');
    % No range change if in simulate mode 2.
    if simMode == 2
        set(hObject,'Value',evalin('base','P.endDepth'));
        return
    end
    Trans = evalin('base','Trans');
    Resource = evalin('base','Resource');
    wlsPer128 = evalin('base','wlsPer128');

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
    PData(1).Size(1) = 10 + ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));
    PData(1).Region = struct(...
                'Shape',struct('Name','SectorFT', ...
                'Position',[0,0,-P.radius], ...
                'z',P.startDepth, ...
                'r',P.radius+P.endDepth, ...
                'angle',P.rayDelta, ...
                'steer',0));
    PData(1).Region = computeRegions(PData(1));
    assignin('base','PData',PData);
    evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta(1));');


    Receive = evalin('base', 'Receive');
    maxAcqLength = ceil(sqrt(P.aperture^2 + P.endDepth^2 - 2*P.aperture*P.endDepth*cos(P.theta-pi/2)) - P.startDepth);
    for i = 1:size(Receive,2)
        Receive(i).endDepth = P.startDepth + wlsPer128*ceil(maxAcqLength/wlsPer128);
    end
    assignin('base','Receive',Receive);
    evalin('base','TGC.rangeMax = P.endDepth;');
    evalin('base','TGC.Waveform = computeTGCWaveform(TGC);');
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'PData','InterBuffer','ImageBuffer','DisplayWindow','Receive','TGC','Recon'};
    assignin('base','Control', Control);
    assignin('base', 'action', 'displayChange');
end
