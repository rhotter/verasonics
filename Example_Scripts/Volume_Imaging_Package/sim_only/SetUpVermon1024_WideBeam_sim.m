% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: SetUpVermon1024_RyLns_sim.m - Example of imaging with single
%                                          plane wave transmit
% Description:
%  Single Plane Wave (flash) image with RF transfer (from just 3/4 secondarys,
%  4th secondary is on order).  This version runs faster because it is only
%  reconstructing the planes that are being displayed.
%
% Last update: 10/19/2018 BWC


clear all

Resource.VDAS.dmaTimeout = 30*1000; % (ms) time software sequencer will wait for 'transferToHost'
% set this time long enough to permit launching the Secondary sequence
% and then launching the Primary sequence

filename = 'MatFiles/Vermon1024_RyLns_RF.mat';

%=== Settings ===%
% Bmode properties
pgain = 8; %digital gain
pers = 0; %set persistance level
senscutoff = 0.6;
tn = 7; %MUST BE ODD.  NUMBER OF ELEMENTS USED FOR APODIZATION IN EACH DIRECTION IE 5X5, 3X3,
TXfocus = 54;

P.numTx = tn^2;   % no. of elements in TX aperture.(5x5)
P.numRays = 7*8; % no. of rays in frame (7x8)
P.startDepth = 0;   % start of Acquisition depth range in wavelengths
P.endDepth =80;     %in wavelengths
P.viewAngle = 14 * pi/180;                    % angle between z-axis and surface of cone.
P.curved = 1; % set to 0 for flat 1 for defocused at the reconstruction apex
P.numFrames = 6;


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
Trans.elementWidth = .275; % mm
Trans.spacing = .3;        % Spacing between elements in mm.
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
zApex = -extent/tan(P.viewAngle);

PData(1).Size(1) = ceil(2.0*(P.endDepth-zApex)*tan(P.viewAngle)/PData(1).PDelta(2));  if mod(PData(1).Size(1),2)==0, PData(1).Size(1) =  PData(1).Size(1)+1; end %rows
PData(1).Size(2) = PData(1).Size(1);%columns
PData(1).Size(3) = ceil((P.endDepth)/PData(1).PDelta(3)); %sections
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


%% Specify Media.  Use point targets in middle of PData.
% Media.MP(1,:) = [0,0,30,1.0];      % single point.
% Media.program = 'PointTargets3D';   % prog. to execute for creating media pts.
Media.function = 'movePoints';
Media.MP(1,:) = [0,0,30,1.0];
Media.MP(2,:) = [7,7,30,1.0];
Media.MP(3,:) = [-7,7,30,1.0];
Media.MP(4,:) = [7,-7,30,1.0];
Media.MP(5,:) = [-7,-7,30,1.0];
Media.MP(6,:) = [0,0,60,1.0];
Media.MP(7,:) = [14,14,60,1.0];
Media.MP(8,:) = [-14,14,60,1.0];
Media.MP(9,:) = [14,-14,60,1.0];
Media.MP(10,:) = [-14,-14,60,1.0];
Media.MP(11,:) = [0,0,90,1.0];
Media.MP(12,:) = [21,21,90,1.0];
Media.MP(13,:) = [-21,21,90,1.0];
Media.MP(14,:) = [21,-21,90,1.0];
Media.MP(15,:) = [-21,-21,90,1.0];
Media.MP(16,:) = [0,0,120,1.0];
Media.MP(17,:) = [28,28,120,1.0];
Media.MP(18,:) = [-28,28,120,1.0];
Media.MP(19,:) = [28,-28,120,1.0];
Media.MP(20,:) = [-28,-28,120,1.0];
Media.MP(21,:) = [0,25,90,1.0];
Media.MP(22,:) = [0,-25,90,1.0];
Media.MP(23,:) = [18,0,60,1.0];
Media.MP(24,:) = [-18,0,60,1.0];
Media.numPoints = size(Media.MP,1);


%% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 4096*P.numRays;
Resource.RcvBuffer(1).colsPerFrame = 1024;
Resource.RcvBuffer(1).numFrames = P.numFrames;
Resource.InterBuffer(1).numFrames = 1;  % one buffer for the primary
Resource.ImageBuffer(1).numFrames = 10;

Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).Title = '3D Flash Image - XZ plane';
Resource.DisplayWindow(1).pdelta = 0.3;
Resource.DisplayWindow(1).Position = [0,580, ...
    ceil(PData(1).Size(2)*PData(1).PDelta(2)/Resource.DisplayWindow(1).pdelta), ... % width
    ceil(PData(1).Size(3)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta)];    % height
Resource.DisplayWindow(1).Orientation = 'xz';
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,0.0];
Resource.DisplayWindow(1).Colormap = gray(256);
Resource.DisplayWindow(1).AxesUnits = 'mm';


Resource.DisplayWindow(2).Type = 'Verasonics';
Resource.DisplayWindow(2).Title = '3D Flash Image - YZ plane';
Resource.DisplayWindow(2).pdelta = Resource.DisplayWindow(1).pdelta;
Resource.DisplayWindow(2).Position = [680,580, ...
    ceil(PData(1).Size(1)*PData(1).PDelta(1)/Resource.DisplayWindow(2).pdelta), ... % width
    ceil(PData(1).Size(3)*PData(1).PDelta(3)/Resource.DisplayWindow(2).pdelta)];    % height
Resource.DisplayWindow(2).Orientation = 'yz';
Resource.DisplayWindow(2).ReferencePt = [0,-PData(1).Origin(2),0];
Resource.DisplayWindow(2).Colormap = gray(256);
Resource.DisplayWindow(2).AxesUnits = 'mm';


Resource.DisplayWindow(3).Type = 'Verasonics';
Resource.DisplayWindow(3).Title = '3D Flash Image - XY plane';
Resource.DisplayWindow(3).pdelta = Resource.DisplayWindow(1).pdelta;
Resource.DisplayWindow(3).Position = [0,40, ...
    ceil(PData(1).Size(2)*PData(1).PDelta(2)/Resource.DisplayWindow(3).pdelta), ... % width
    ceil(PData(1).Size(1)*PData(1).PDelta(1)/Resource.DisplayWindow(3).pdelta)];    % height
Resource.DisplayWindow(3).Orientation = 'xy';
Resource.DisplayWindow(3).ReferencePt = [PData(1).Origin(1),-PData(1).Origin(2),PData(1).Region(end).Shape.oPAIntersect];
Resource.DisplayWindow(3).Colormap = gray(256);
Resource.DisplayWindow(3).AxesUnits = 'mm';

% Specify Transmit waveform structure.
TW.type = 'parametric';
TW.Parameters = [Trans.frequency,.67,2,1];  % A, B, C, D

% Specify TPC structure.
TPC(1).hv = 30.0;

%% Specify TX structure
%patterns are based on 5 elements
patternA=[3 8 13 16 20 25 30];  %center elements (center overlap) 7 centers
patternB=[3 6 11 16 21 26 30];  %center elements (outer overlap) 7 centers
%00000000 01111111 11122222 22222333
%12345678_90123456_78901234_56789012
%  ^  ^     ^  ^     ^  ^     ^  ^
patternC=([3 6 11 14 19 22 27 30]-1)*32; %8 centers, 2 per connector set of channels
origElements=repmat( patternA,8,1)+repmat(patternC,7,1)';

TX = repmat(struct('waveform', 1, ...
    'Origin', [0.0,0.0,0.0], ...
    'focus', TXfocus, ...
    'Steer', [0.0,0.0], ...
    'Apod', zeros(1,1024), ...
    'Delay', zeros(1,1024),...
    'peakCutOff', 8.0,...
    'peakBLMax', 8.0),1,P.numRays);

for n = 1:P.numRays
    TX(n).Origin(1) = Trans.ElementPosWL(origElements(n),1);
    TX(n).Origin(2) = Trans.ElementPosWL(origElements(n),2);
    %TX(n).FocalPt = [ TX(n).Origin(1) TX(n).Origin(2) TXfocus];  %focalPt causes errors in computeTXPD
    apods = repmat([-floor(tn/2):floor(tn/2)],tn,1)+repmat(((-floor(tn/2):floor(tn/2))*32)',1,tn)+origElements(n);
    TX(n).Apod(apods(apods>0&apods<1025))=1;
    TX(n).Apod(:)=1; %for debugging
    %--check apodization
    %showA=zeros(32,32);
    %showA(find(TX(n).Apod==1))=1;
    %pcolor(showA)
    %drawnow
    %----%
    TX(n).Delay = computeTXDelays(TX(n));

    %disp(['calculating TXPD ' num2str(n) ' out of ' num2str(length(TX)) ' ...'])
    %TX(n).TXPD = computeTXPD(TX(n),PData(1));
    %disp(['remainging time: ' num2str(round((toc/n)*length(TX)-toc)) ' seconds'])

end

%% Specify TGC Waveform structure.

TGC(1).CntrlPts = [600 700 800 850 900 900 900 900];
TGC(1).rangeMax = 128;
TGC(1).Waveform = computeTGCWaveform(TGC);



%% - Receive
temp = (P.endDepth-zApex)*tan(P.viewAngle)+ max([max(Trans.ElementPosWL(:,1)), max(Trans.ElementPosWL(:,2))]);
maxAcqLength = sqrt(P.endDepth^2 + temp^2);
Receive = repmat(struct(...
    'Apod', ones(1,Resource.Parameters.numRcvChannels), ...
    'startDepth', 0, ...
    'endDepth', 128/(4*2)*ceil(maxAcqLength/(128/(4*2))), ...
    'TGC', 1, ...
    'mode', 0, ...
    'bufnum', 1, ...
    'framenum', 1, ...
    'acqNum', 1, ...
    'sampleMode', 'NS200BW', ...
    'callMediaFunc', 0), 1, P.numRays*P.numFrames);

%Set event specific Receive attributes.
for j = 1:P.numFrames
    Receive(i+(j-1)*P.numRays).callMediaFunc = 1;
    for i = 1:P.numRays
        Receive(i+(j-1)*P.numRays).acqNum = i;
        Receive(i+(j-1)*P.numRays).framenum = j;
    end
end

%%
%% Specify Recon structure arrays.
Recon(1) = struct('senscutoff', senscutoff, ...
    'pdatanum', 1, ...
    'rcvBufFrame', -1, ...
    'IntBufDest', [1,1],...
    'ImgBufDest', [1,-1], ...
    'RINums', 1:P.numRays*3);
%Recon(1).RINums = Recon(1).RINums(rem((1:56*3),3)>0);  %ignore the XY plane for reconstruction
% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 'accumIQ', ...  % default is to accumulate IQ data.
    'Pre',[],...
    'Post',[],...
    'scaleFactor', 1, ...
    'threadSync', 1), 1, P.numRays*3);  %need 3 x number of rays recon, 1 for each imaging plane

ReconInfo(1).Pre = 'clearInterBuf';
for j = 1:P.numRays
    for k = 1:3
        ReconInfo((j-1)*3+k).txnum = j;
        ReconInfo((j-1)*3+k).rcvnum = j;
        ReconInfo((j-1)*3+k).regionnum = k;
    end
end
ReconInfo(end).Post = 'IQ2IntensityImageBuf';


%% Specify Process structure array.

pgain=2;
persist = 10;
CompressionFactor = 20;
% Specify Process structure arrays.
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...
    'framenum',-1,...
    'pdatanum',1,...
    'srcData','intensity3D',...
    'pgain', pgain,...
    'persistMethod','simple',...
    'persistLevel',persist,...
    'interpMethod','4pt',...
    'compressMethod','power',...
    'compressFactor',CompressionFactor,...
    'mappingMethod','full',...
    'display',1,...
    'displayWindow',1};
Process(2).classname = 'Image';
Process(2).method = 'imageDisplay';
Process(2).Parameters = {'imgbufnum',1,...
    'framenum',-1,...
    'pdatanum',1,...
    'srcData','intensity3D',...
    'pgain', pgain,...
    'persistMethod','simple',...
    'persistLevel',persist,...
    'interpMethod','4pt',...
    'compressMethod','power',...
    'compressFactor',CompressionFactor,...
    'mappingMethod','full',...
    'display',1,...
    'displayWindow',2};
Process(3).classname = 'Image';
Process(3).method = 'imageDisplay';
Process(3).Parameters = {'imgbufnum',1,...
    'framenum',-1,...
    'pdatanum',1,...
    'srcData','intensity3D',...
    'pgain', pgain,...
    'persistMethod','simple',...
    'persistLevel',persist,...
    'interpMethod','4pt',...
    'compressMethod','power',...
    'compressFactor',CompressionFactor,...
    'mappingMethod','full',...
    'display',1,...
    'displayWindow',3};


%% Specify SeqControl structure arrays.
% Specify SeqControl structure arrays.
JUMP = 1;
SeqControl(JUMP).command = 'jump'; % jump back to start
SeqControl(JUMP).argument = 1;
%
TTNA = 2;
SeqControl(TTNA).command = 'timeToNextAcq';  % time between acquisitions
%SeqControl(TTNA).argument = round(2*P.endDepth/Trans.frequency); %in usec
SeqControl(TTNA).argument = 10000; %in usec
%
RETMAT = 3;
%-- return to matlab
SeqControl(RETMAT).command = 'returnToMatlab';

nsc = length(SeqControl)+1; % nsc is count of SeqControl objects


%% Event Sequence
n = 1; % n is count of Events
% Acquire all frames defined in RcvBuffer
for j = 1:Resource.RcvBuffer(1).numFrames
    for k = 1:P.numRays
        Event(n).info = 'Transmit/Receive.';
        Event(n).tx = k;   % use 1st TX structure.
        Event(n).rcv = (j-1)*P.numRays+k;  % use 1st Rcv structure of frame.
        Event(n).recon = 0;      % no reconstruction.
        Event(n).process = 0;    % no processing
        Event(n).seqControl = TTNA;
        n = n+1;
    end

    Event(n).info = 'Transfer To Host';
    Event(n).tx = 0;         % use 1st TX structure.
    Event(n).rcv = 0;    % use 1st Rcv structure of frame.
    Event(n).recon = 0;      % no reconstruction.
    Event(n).process = 0;    % no processing
    Event(n).seqControl = [nsc,nsc+1]; % set wait time and transfer data
    SeqControl(nsc).command = 'transferToHost';
    nsc = nsc + 1;
    SeqControl(nsc).command = 'waitForTransferComplete';
    SeqControl(nsc).argument = nsc-1;
    nsc = nsc + 1;
    n = n+1;

    Event(n).info = 'Recon';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 0;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = ['Process XZ - Frame ' num2str(j)];
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 1;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = ['Process YZ - Frame ' num2str(j)];
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 2;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = ['Process XY - Frame ' num2str(j)];
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;%3;
    Event(n).seqControl = 0;
    n = n+1;

end


Event(n).info = 'exit to matlab';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 0;
if floor(j/5) == j/5     % Exit to Matlab every 5th frame
    Event(n).seqControl = RETMAT;
end
n = n+1;

Event(n).info = 'Jump back to first event';
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
MinMaxVal = [25,40,35]; % default unit is wavelength
AxesUnit = 'mm';

if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
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
save(filename);
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
