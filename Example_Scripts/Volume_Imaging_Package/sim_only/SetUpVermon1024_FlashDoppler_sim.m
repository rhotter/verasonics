% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: SetUpVermon1024_FlashDoppler.m - Example of Doppler
% Description:
%
% Last update: 10/31/2018 BWC

clear all


filename = 'MatFiles/Vermon1024_FlashDoppler_sim.mat';
%=== Settings ===%

P(1).numTx = 1;
P(1).startDepth = 0;   % Acquisition depth in wavelengths
P(1).endDepth = 54;   % This should preferrably be a multipled of 128 samples.
P(1).viewAngle = 14 * pi/180;                    % angle between z-axis and surface of cone.
P(1).curved = 1; % set to 0 for flat 1 for defocused at the reconstruction apex
P(1).numFrames = 4;

P(2).numTx = 1;
P(2).startDepth = 0;   % Acquisition depth in wavelengths
P(2).endDepth = 54;   % This should preferrably be a multipled of 128 samples.
P(2).viewAngle = 14 * pi/180;                    % angle between z-axis and surface of cone.
P(2).curved = 1; % set to 0 for flat 1 for defocused at the reconstruction apex
P(2).numFrames = 4;

% Set Bmode parameters
senscutoff = 0.8;
pgain=1;
na = 1;      % Set na = number of flash angles for 2D.
% set dtheta2D to range over +/- 15 degrees.
if (na > 1)
    dtheta2D = (30*pi/180)/(na-1);
    startAngle = -30*pi/180/2;
else dtheta2D = 0;
    startAngle=0;
end

% Set Doppler parameters
ne = 14;     % Set ne = number of acquisitions in Doppler ensemble.
dopAngle = 7 * pi/180;
dopPRF = 3500; % Doppler PRF in Hz.
pwrThres = 0.1;
cpt = 255;
persf = 80;
DopState='computeCFIPowerEst';

RcvProfile(1).AntiAliasCutoff = 10;
RcvProfile(1).LnaGain = 18; % Profile used for imaging
RcvProfile(2).LnaGain = 24; % Profile used for Doppler
RcvProfile(2).LnaZinSel = 31; % Force high-Z state for best Doppler sensitivity

% Define system parameters.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 0;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.numTransmit = 1024;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 1024;    % number of receive channels.
Resource.Parameters.simulateMode = 1; % processes RcvData continuously.
Resource.System.Product= 'SimulateOnly';  % Needed for simulation w/o hardware


%% Specify Trans structure array.
Trans.id = hex2dec('0D0100');
Trans.units = 'mm';
Trans.name = '1024Direct';
%UTA TYPE[ 1  8  1  0 ]
Trans.connType = 8;
Trans.frequency = 3.3;
Trans.Bandwidth = Trans.frequency*[0.7, 1.3];
Trans.type = 2;             % Array geometry is 2D.
Trans.numelements = 1024;
Trans.impedance = 50;
Trans.elementWidth = 0.275; % mm
Trans.spacing = 0.3;        % Spacing between elements in mm.
Trans.maxHighVoltage = 60;
Trans.Connector = [1:1024]';
eleY = Trans.spacing;
eleX = Trans.spacing;
% Define a square array EM (Element Map) that shows the location of
% each element (numbered 1 to 128) on the physical grid of transducer
% elements.  Locations that represent either an element that does not exist
% or an element that is not connected are identified by setting the element
% number to zero.
GSRow = 35;
GSCol = 32;
EM = zeros(GSRow,GSCol);  % define the array and fill it with zeros
% Physical element array order defined in a for loop with three row gaps
for i = 1:GSRow
    if i==1; j=0; end;
    if i==9 || i==18 || i==27
        continue;
    end
    j=j+1;
    EM(i,:) = (j-1)*GSCol+(1:32);
end
% Next use EM array to assign xy position for each element in Trans.ElementPos array
Trans.ElementPos = zeros(1024,5,'double');  % define array and fill with zeros
Trans.ElementPos(:,1:2) = -1;  % fill x,y positions with -1, to allow test for

for i=1:GSRow
    for j=1:GSCol
        if(EM(i,j)>0)   % note positions are relative to center of array
            Trans.ElementPos(EM(i,j),1) = (j-1)*eleX - ((GSCol-1)*eleX)/2;
            Trans.ElementPos(EM(i,j),2) = ((GSRow-1)*eleY)/2 - (i-1)*eleY;
        end
    end
end
clear EM GSCol GSRow


%Intermediate Variables
waveLength = (Resource.Parameters.speedOfSound/1000)/Trans.frequency;
if strcmp(Trans.units,'mm')
    Trans.ElementPosMm = Trans.ElementPos;
    Trans.ElementPosWL = Trans.ElementPos./waveLength;
else
    Trans.ElementPosMm = Trans.ElementPos.*waveLength;
    Trans.ElementPosWL = Trans.ElementPos;
end

% Set element sensitivity function (100 weighting values from -pi/2 to pi/2).
Theta = (-pi/2:pi/100:pi/2);%   Set element sensitivity function (100 weighting values from -pi/2 to pi/2).
Theta(51) = 0.0000001;
Trans.ElementSens = sin(Trans.elementWidth*pi*sin(Theta))./(Trans.elementWidth*pi*sin(Theta));
clear Theta

%% PData
PData(1).PDelta = [1,1,1]*1;
if ~exist('extent','var'), extent = max(max(Trans.ElementPosWL(:,1)),max(Trans.ElementPosWL(:,2))); end
zApex = -extent/tan(P(1).viewAngle);

PData(1).Size(1) = ceil(2.0*(P(1).endDepth-zApex)*tan(P(1).viewAngle)/PData(1).PDelta(2));  if mod(PData(1).Size(1),2)==0, PData(1).Size(1) =  PData(1).Size(1)+1; end %rows
PData(1).Size(2) = PData(1).Size(1);%columns
PData(1).Size(3) = ceil((P(1).endDepth)/PData(1).PDelta(3)); %sections
PData(1).Origin = [-((PData(1).Size(2)-1)/2)*PData(1).PDelta(1), ((PData(1).Size(1)-1)/2)*PData(1).PDelta(1), 0];

PData(1).Region(1) = struct(...
    'Shape',struct('Name','Slice',...
    'Orientation','xz',...
    'oPAIntersect',0));
PData(1).Region(2) = struct(...
    'Shape',struct('Name','Slice',...
    'Orientation','yz',...
    'oPAIntersect',0));
PData(1).Region(3) = struct(...
    'Shape',struct('Name','Slice',...
    'Orientation','xy',...
    'oPAIntersect',40));
PData(1).Region = computeRegions(PData(1));

PData(2)=PData(1);

%% Specify Media.  Use point targets in middle of PData.
% Media.MP(1,:) = [0,0,30,1.0];      % single point.
% Media.program = 'PointTargets3D';   % prog. to execute for creating media pts.
Media.function = 'movePoints';
Media.MP(1,:) = [0,0,30,1.0];
% Media.MP(2,:) = [7,7,30,1.0];
% Media.MP(3,:) = [-7,7,30,1.0];
% Media.MP(4,:) = [7,-7,30,1.0];
% Media.MP(5,:) = [-7,-7,30,1.0];
% Media.MP(6,:) = [0,0,60,1.0];
% Media.MP(7,:) = [14,14,60,1.0];
% Media.MP(8,:) = [-14,14,60,1.0];
% Media.MP(9,:) = [14,-14,60,1.0];
% Media.MP(10,:) = [-14,-14,60,1.0];
% Media.MP(11,:) = [0,0,90,1.0];
% Media.MP(12,:) = [21,21,90,1.0];
% Media.MP(13,:) = [-21,21,90,1.0];
% Media.MP(14,:) = [21,-21,90,1.0];
% Media.MP(15,:) = [-21,-21,90,1.0];
% Media.MP(16,:) = [0,0,120,1.0];
% Media.MP(17,:) = [28,28,120,1.0];
% Media.MP(18,:) = [-28,28,120,1.0];
% Media.MP(19,:) = [28,-28,120,1.0];
% Media.MP(20,:) = [-28,-28,120,1.0];
% Media.MP(21,:) = [0,25,90,1.0];
% Media.MP(22,:) = [0,-25,90,1.0];
% Media.MP(23,:) = [18,0,60,1.0];
% Media.MP(24,:) = [-18,0,60,1.0];
Media.numPoints = size(Media.MP,1);


%% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 1024*(na+ne);
Resource.RcvBuffer(1).colsPerFrame = 1024;
Resource.RcvBuffer(1).numFrames = P.numFrames;
%--- Interbuffer---%
Resource.InterBuffer(1).numFrames = 1;  % one buffer for the primary
Resource.InterBuffer(2).numFrames = 1;  % one buffer for the primary
Resource.InterBuffer(2).pagesPerFrame = ne;     % ne pages per ensemble
%--- ImageBuffer ---%
Resource.ImageBuffer(1).numFrames = 10;
Resource.ImageBuffer(2).numFrames = 10;

Resource.DisplayWindow(1).Type = 'Matlab';
Resource.DisplayWindow(1).Title = '3D Flash Image - XZ plane';
Resource.DisplayWindow(1).pdelta = 0.25;
Resource.DisplayWindow(1).Position = [0,580, ...
    ceil(PData(1).Size(2)*PData(1).PDelta(2)/Resource.DisplayWindow(1).pdelta), ... % width
    ceil(PData(1).Size(3)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta)];    % height
Resource.DisplayWindow(1).Orientation = 'xz';
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,0.0];
Resource.DisplayWindow(1).Colormap = grayscaleCFImap;
Resource.DisplayWindow(1).AxesUnits = 'wavelengths';

Resource.DisplayWindow(2).Type = 'Matlab';
Resource.DisplayWindow(2).Title = '3D Flash Image - YZ plane';
Resource.DisplayWindow(2).pdelta = Resource.DisplayWindow(1).pdelta;
Resource.DisplayWindow(2).Position = [430,580, ...
    ceil(PData(1).Size(1)*PData(1).PDelta(1)/Resource.DisplayWindow(2).pdelta), ... % width
    ceil(PData(1).Size(3)*PData(1).PDelta(3)/Resource.DisplayWindow(2).pdelta)];    % height
Resource.DisplayWindow(2).Orientation = 'yz';
Resource.DisplayWindow(2).ReferencePt = [0,-PData(1).Origin(2),0];
Resource.DisplayWindow(2).Colormap = grayscaleCFImap;
Resource.DisplayWindow(2).AxesUnits = 'wavelengths';

Resource.DisplayWindow(3).Type = 'Matlab';
Resource.DisplayWindow(3).Title = '3D Flash Image - XY plane';
Resource.DisplayWindow(3).pdelta = Resource.DisplayWindow(1).pdelta;
Resource.DisplayWindow(3).Position = [0,40, ...
    ceil(PData(1).Size(2)*PData(1).PDelta(2)/Resource.DisplayWindow(3).pdelta), ... % width
    ceil(PData(1).Size(1)*PData(1).PDelta(1)/Resource.DisplayWindow(3).pdelta)];    % height
Resource.DisplayWindow(3).Orientation = 'xy';
Resource.DisplayWindow(3).ReferencePt = [PData(1).Origin(1),-PData(1).Origin(2),PData(1).Region(end).Shape.oPAIntersect];
Resource.DisplayWindow(3).Colormap = grayscaleCFImap;
Resource.DisplayWindow(3).AxesUnits = 'wavelengths';


% Specify Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,.67,2,1];  % A, B, C, D

TW(2).type = 'parametric';
TW(2).Parameters = [Trans.frequency,.67,6,1];  % A, B, C, D

% Specify TPC structure.
TPC(1).hv = 30.0;
TPC(2).hv = 30.0;


%RcvProfile(1).LnaGain = 18; % Profile used for imaging
RcvProfile(1).LnaGain = 24; % Profile used for Doppler
RcvProfile(1).LnaZinSel = 31; % Force high-Z state for best Doppler sensitivity

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
    'Origin', [0.0,0.0,0.0], ...
    'focus', 0.0, ...
    'Steer', [0.0,0.0], ...
    'Apod', ones(1,Trans.numelements), ...
    'Delay', zeros(1,Trans.numelements),...
    'TXPD', [], ...
    'peakCutOff', 11,...
    'peakBLMax', 7), 1, na+1); % na TXs for 2D + 1 for Doppler

% - Set event specific TX attributes.
for n = 1:na   % na transmit events for 2D
    TX(n).Steer = [(startAngle+(n-1)*dtheta2D),0.0];
    TX(n).Delay = computeTXDelays(TX(n));
end
% -- only one TX struct needed for Doppler
TX(na+1).waveform = 2;
TX(na+1).Steer = [dopAngle,0.0];
TX(na+1).Delay = computeTXDelays(TX(na+1));


%% Specify TGC Waveform structure.

TGC(1).CntrlPts = [600 700 800 850 900 900 900 900];
TGC(1).rangeMax = 128;
TGC(1).Waveform = computeTGCWaveform(TGC(1));

TGC(2).CntrlPts = [600 700 800 850 900 900 900 900];
TGC(2).rangeMax = 128;
TGC(2).Waveform = computeTGCWaveform(TGC(2));

temp = (P(1).endDepth-zApex)*tan(P(1).viewAngle)+ max([max(Trans.ElementPosWL(:,1)), max(Trans.ElementPosWL(:,2))]);
maxAcqLength2D = sqrt(P(1).endDepth^2 + temp^2);
temp = (P(2).endDepth-zApex)*tan(P(2).viewAngle)+ max([max(Trans.ElementPosWL(:,1)), max(Trans.ElementPosWL(:,2))]);
maxAcqLengthDop = sqrt(P(2).endDepth^2 + temp^2);

wl4sPer128 = 128/(4*2);  % wavelengths in a 128 sample block for 4 smpls per wave round trip.
wl2sPer128 = 128/(2*2);  % wavelengths in a 128 sample block for 2 smpls per wave round trip.

%Receive

Receive = repmat(struct(...
    'Apod', ones(1,1024), ...
    'startDepth', 0, ...
    'endDepth', wl4sPer128*ceil(maxAcqLength2D/(wl4sPer128)), ...
    'TGC', 1, ...
    'mode', 0, ...
    'bufnum', 1, ...
    'framenum', 1, ...
    'acqNum', 1, ...
    'sampleMode', 'NS200BW', ...
    'demodFrequency',TW(1).Parameters(1),...
    'callMediaFunc', 0), 1, (na+ne)*P(1).numFrames);

%Set event specific Receive attributes.
for i = 1:P(1).numFrames
    k = (na + ne)*(i-1); % k keeps track of Receive index increment per frame.
    % - Set attributes for each frame.
    Receive(k+1).callMediaFunc = 1;
    for j = 1:na %acquisitions for 2D
        Receive(j+k).framenum = i;
        Receive(j+k).acqNum = j;
    end
    for j = (na+1):(na+ne)
        % Doppler acquisition
        Receive(j+k).startDepth = P(2).startDepth;
        Receive(j+k).endDepth = P(2).startDepth + wl2sPer128*ceil(maxAcqLengthDop/wl2sPer128);
        Receive(j+k).sampleMode = 'BS100BW';
        Receive(j+k).demodFrequency = TW(2).Parameters(1);
        Receive(j+k).TGC = 2;
        Receive(j+k).framenum = i;
        Receive(j+k).acqNum = j;        % Doppler acqNums continue after 2D
    end
end
%%


% Specify Recon structure arrays.
Recon = repmat(struct('senscutoff', senscutoff, ...
    'pdatanum', 1, ...
    'rcvBufFrame', -1, ...
    'IntBufDest', [1,1],...
    'ImgBufDest', [1,-1], ...
    'RINums', zeros(1,1)), 1, 2);

% Define ReconInfo structures.
% - Set Recon values for 2D frame.
Recon(1).RINums(1,1:(3*na)) = (1:(3*na));  % na ReconInfos needed for na angles
k = na*3 +1;
% - Set Recon values for Doppler ensemble.
Recon(2).pdatanum = 2;
Recon(2).IntBufDest = [2,1];
Recon(2).ImgBufDest = [0,0];
Recon(2).RINums(1,1:(3*ne)) = (k:(k+(3*ne)-1));   % ne ReconInfos needed for Doppler ensemble.


% Define ReconInfo structures.
% - For 2D, we need na ReconInfo structures for na steering angles.
% - For Doppler, we need ne ReconInfo structures.
ReconInfo = repmat(struct('mode', 'accumIQ', ...    % accumulate IQ data.
    'txnum', 1, ...
    'rcvnum', 1, ...
    'pagenum',1, ...
    'regionnum', 1), 1, 3*(na + ne));
% - Set specific ReconInfo attributes.
%   - ReconInfos for 2D frame.
if na>1
    ReconInfo(1).mode = 'replaceIQ'; % replace IQ data
    for j = 1:(na*3)
        ReconInfo(j).txnum = j;
        ReconInfo(j).rcvnum = j;
    end
    ReconInfo(na).mode = 'accumIQ_replaceIntensity'; % accum and detect
else
    ReconInfo(1).mode = 'replaceIntensity';
end

%  - ReconInfos for Doppler ensemble.
k = na*3;
for j = 1:ne
    for i = 1:3
        ReconInfo(k+(j-1)*3+i).mode = 'replaceIQ';
        ReconInfo(k+(j-1)*3+i).txnum = na + 1;
        ReconInfo(k+(j-1)*3+i).rcvnum = na + j;
        ReconInfo(k+(j-1)*3+i).pagenum = j;
        ReconInfo(k+(j-1)*3+i).regionnum = i;
    end
end

%% Specify Process structure array.

CompressionFactor = 20;
% Specify Process structure arrays.
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'srcData','intensity3D',...
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

Process(2).classname = 'Image';
Process(2).method = 'imageDisplay';
Process(2).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'srcData','intensity3D',...
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
                         'displayWindow',2};
Process(3).classname = 'Image';
Process(3).method = 'imageDisplay';
Process(3).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'srcData','intensity3D',...
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
                         'displayWindow',3};

%%%---> NOTE: 3D DOPPLER PROCESSING DOES NOT WORK YET.  NEEDS TO BE FIXED
%%%IN RUNACQ
Process(4).classname = 'Doppler';                   % process structure for 1st Doppler ensemble
Process(4).method = DopState;
Process(4).Parameters = {'IntBufSrc',[2,1],...      % buffer and frame num of interbuffer
                         'SrcPages',[3,ne-2],...    % start and last pagenum
                         'ImgBufDest',[2,-1],...    % buffer and frame num of image desitation
                         'pdatanum',2,...           % number of PData structure
                         'prf',dopPRF,...           % Doppler PRF in Hz
                         'wallFilter','regression',...  % 1 -> quadratic regression
                         'pwrThreshold',pwrThres,...
                         'maxPower',50,...
                         'postFilter',1};

Process(5).classname = 'Image';
Process(5).method = 'imageDisplay';
Process(5).Parameters = {'imgbufnum',2,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'srcData','unsignedColor',... % type of data to display.
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

Process(6).classname = 'Image';
Process(6).method = 'imageDisplay';
Process(6).Parameters = {'imgbufnum',2,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'srcData','unsignedColor',... % type of data to display.
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
                         'displayWindow',2};

Process(7).classname = 'Image';
Process(7).method = 'imageDisplay';
Process(7).Parameters = {'imgbufnum',2,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'srcData','unsignedColor',... % type of data to display.
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
                         'displayWindow',3};

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
    Event(n).process = 0;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = 'Doppler processing';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 4;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = 'Bmode Display 1';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 1;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = 'Bmode Display 2';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 2;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = 'Bmode Display 3';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 3;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = 'Doppler image display';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 5;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = 'Doppler image display';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 6;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = 'Doppler image display';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 7;
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



if 0%~HEADLESS
    % User specified UI Control Elements
    % - Sensitivity Cutoff
    UI(1).Control = vsv.seq.uicontrol.VsSliderControl('LocationCode','UserB7','Label','Sens. Cutoff',...
                    'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                    'SliderStep',[0.025,0.1],'ValueFormat','%1.3f',...
                    'Callback', @SensCutoffCallback);

    % - Range Change
    MinMaxVal = [64,300,P.endDepth]; % default unit is wavelength
    AxesUnit = 'wls';

    if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
        if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
            AxesUnit = 'mm';
            MinMaxVal = MinMaxVal * (Resource.Parameters.speedOfSound/1000/Trans.frequency);
        end
    end
    UI(2).Control = vsv.seq.uicontrol.VsSliderControl('LocationCode','UserA1','Label',['Range (',AxesUnit,')'],...
                    'SliderMinMaxVal',MinMaxVal,...
                    'SliderStep',[0.025,0.1],[0.1,0.2],'ValueFormat','%3.0f',...
                    'Callback', @RangeChangeCallback);
    
end

% External function definition.
% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 5;

% Save all the structures to a .mat file.
save(filename);
return

% **** Callback routines to be converted by text2cell function. ****
%SensCutoffCallback - Sensitivity cutoff change
function SensCutoffCallback(hObject,event,UIValue)
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
%SensCutoffCallback

%RangeChangeCallback - Range change
function RangeChangeCallback(hObject,event,UIValue)
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
end
%RangeChangeCallback




