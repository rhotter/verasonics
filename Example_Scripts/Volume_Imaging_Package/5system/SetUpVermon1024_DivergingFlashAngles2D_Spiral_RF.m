%% Copyright 2021 Verasonics, Inc. Verasonics Registered U.S. Patent and Trademark Office.
% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name:
%
% Description:

%
% Last update: 8/6/2019 BWC


clear all
filename = 'Vermon1024_3DivergingFlashAngles2D_RF';

%=== Settings ===%
% Bmode properties
XYdepth = 10; %depth of C plane imaging plane
TxAngles =[-20 -10 0 10 20];
na1 = length(TxAngles);      % Set na = number of angles.
na2= 9; %this parameter determines the number of planes around 2*pi

pgain = 20;              % digital gain
persist = 30;            % persistance level
senscutoff = 0.6;        % sensitivity cutoff
CompressionFactor = 75;  % compression Factor
P.startDepth = 0;   % Acquisition depth in wavelengths
P.endDepth = 64;   % This should preferrably be a multipled of 128 samples.
P.viewAngle = 20 * pi/180;                    % angle between z-axis and surface of cone.
P.curved = 1; % set to 0 for flat 1 for defocused at the reconstruction apex
P.numFrames = 10;

RcvProfile.AntiAliasCutoff = 15; % 10, 15, 20, 30
% RcvProfile.LnaGain = 18;
% RcvProfile.LnaZinSel = 31;

frameRateFactor = 5;% Specify factor for converting sequenceRate to frameRate.

%%==================== RDMA SETUP ======================================%%
usingMultiSys = 1;
[RDMAconfig, RDMAsystem]=vsv.multi.generateConfig(5);
%%=======================================================================%%

if strcmp(RDMAconfig.rdmaRole,'secondary')
    Resource.Parameters.simulateMode = 0;
else
    Resource.Parameters.simulateMode = 2;
end

% Define system parameters.
if strcmp(RDMAconfig.rdmaRole,'secondary')
    Resource.Parameters.numTransmit = 256;      % number of transmit channels.
    Resource.Parameters.numRcvChannels = 256;    % number of receive channels.
else
    Resource.Parameters.numTransmit = 1024;      % number of transmit channels.
    Resource.Parameters.numRcvChannels = 1024;    % number of receive channels.
end

Resource.Parameters.speedOfSound = 1498;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 0;
Resource.Parameters.initializeOnly = 0;

%% Specify Trans structure array.
Trans.name = 'Matrix1024-3';
Trans.units = 'wavelengths';  % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);  % computeTrans is used for known transducers.
Trans.id = hex2dec('0D0100'); % need this still
if strcmp(RDMAconfig.rdmaRole,'secondary')
    Trans.ConnectorES=zeros(1,1024)';
    Trans.ConnectorES((RDMAconfig.nodeIndex-1)*256+[1:256])=[128:-1:65, 1:1:64, 192:-1:129, 193:1:256];
else
    Trans.connType=0;
    Trans.ConnectorES=[0+[128:-1:65, 1:1:64, 192:-1:129, 193:1:256],...
        256+[128:-1:65, 1:1:64, 192:-1:129, 193:1:256],...
        512+[128:-1:65, 1:1:64, 192:-1:129, 193:1:256],...
        768+[128:-1:65, 1:1:64, 192:-1:129, 193:1:256]]';
end

%Intermediate Variables
waveLength = (Resource.Parameters.speedOfSound/1000)/Trans.frequency;
if strcmp(Trans.units,'mm')
    Trans.ElementPosMm = Trans.ElementPos;
    Trans.ElementPosWL = Trans.ElementPos./waveLength;
else
    Trans.ElementPosMm = Trans.ElementPos.*waveLength;
    Trans.ElementPosWL = Trans.ElementPos;
end

%% PData
PData(1).PDelta = [1,1,1]*1.0;
if ~exist('extent','var'), extent = max(max(Trans.ElementPosWL(:,1)),max(Trans.ElementPosWL(:,2))); end
zApex = -extent/tan(P.viewAngle);

PData(1).Size(1) = ceil(2.0*(P.endDepth-zApex)*tan(P.viewAngle)/PData(1).PDelta(2));  if mod(PData(1).Size(1),2)==0, PData(1).Size(1) =  PData(1).Size(1)+1; end %rows
PData(1).Size(2) = PData(1).Size(1);%columns
PData(1).Size(3) = ceil((P.endDepth)/PData(1).PDelta(3)); %sections
PData(1).Origin = [-((PData(1).Size(2)-1)/2)*PData(1).PDelta(1), ((PData(1).Size(1)-1)/2)*PData(1).PDelta(1), 0]; %sections

PData(1).Region(1) = struct(...
    'Shape',struct('Name','Pyramid',...
    'Position',[0,0,zApex],...
    'angle',P.viewAngle,...
    'z1',P.startDepth,...
    'z2',P.endDepth));
PData(1).Region(2) = struct(...
    'Shape',struct('Name','Slice',...
    'Orientation','xz',...
    'oPAIntersect', 0, 'andWithPrev',1));
PData(1).Region(3) = struct(...
    'Shape',struct('Name','Slice',...
    'Orientation','yz',...
    'oPAIntersect', 0, 'andWithPrev',1));
PData(1).Region(4) = struct(...
    'Shape',struct('Name','Slice',...
    'Orientation','xy',...
    'oPAIntersect', XYdepth, 'andWithPrev',1));
PData(1).Region = computeRegions(PData(1));

PData.Region(5).PixelsLA = unique([PData.Region(2).PixelsLA;PData.Region(3).PixelsLA;PData.Region(4).PixelsLA]);
PData.Region(5).Shape.Name = 'Custom';
PData.Region(5).numPixels = length(PData.Region(5).PixelsLA);


%% Specify Media.  Use point targets in middle of PData.
for i=1:61
    Media.MP(i,:)= [0,-30+i,50,0.1];
    Media.MP(61+i,:)= [-30+i,0,72,0.1];
end
Media.MP(size(Media.MP,1)+1,:)=[18,32,82,1];
Media.numPoints = size(Media.MP,1);

%% Specify Resources.
switch RDMAconfig.rdmaRole
    case 'primary'
        Resource.RcvBuffer(1).datatype = 'int16';
        Resource.RcvBuffer(1).rowsPerFrame = 41344;
        Resource.RcvBuffer(1).colsPerFrame = 1024;
        Resource.RcvBuffer(1).numFrames = P.numFrames;
        Resource.InterBuffer(1).numFrames = 1;  % one buffer for the primary
        Resource.ImageBuffer(1).numFrames = P.numFrames;
        %- XZ Plane
        Resource.DisplayWindow(1).Type = 'Verasonics';
        Resource.DisplayWindow(1).Title = '3D Flash Image - XZ plane';
        Resource.DisplayWindow(1).pdelta = 0.25;
        Resource.DisplayWindow(1).Position = [0,580, ...
            ceil(PData(1).Size(2)*PData(1).PDelta(2)/Resource.DisplayWindow(1).pdelta), ... % width
            ceil(PData(1).Size(3)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta)];    % height
        Resource.DisplayWindow(1).Orientation = 'xz';
        Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0.0,0.0];
        Resource.DisplayWindow(1).Colormap = gray(256);
        Resource.DisplayWindow(1).AxesUnits = 'mm';
        % Resource.DisplayWindow(1).mode = '2d';
        %-- YZ Plane
        Resource.DisplayWindow(2).Type = 'Verasonics';
        Resource.DisplayWindow(2).Title = '3D Flash Image - YZ plane';
        Resource.DisplayWindow(2).pdelta = Resource.DisplayWindow(1).pdelta;
        Resource.DisplayWindow(2).Position = [430,580, ...
            ceil(PData(1).Size(1)*PData(1).PDelta(1)/Resource.DisplayWindow(2).pdelta), ... % width
            ceil(PData(1).Size(3)*PData(1).PDelta(3)/Resource.DisplayWindow(2).pdelta)];    % height
        Resource.DisplayWindow(2).Orientation = 'yz';
        Resource.DisplayWindow(2).ReferencePt = [0,-PData(1).Origin(2),0];
        Resource.DisplayWindow(2).Colormap = gray(256);
        Resource.DisplayWindow(2).AxesUnits = 'mm';
        % Resource.DisplayWindow(2).mode = '2d';
        %-- XY Plane
        Resource.DisplayWindow(3).Type = 'Verasonics';
        Resource.DisplayWindow(3).Title = '3D Flash Image - XY plane';
        Resource.DisplayWindow(3).pdelta = Resource.DisplayWindow(1).pdelta;
        Resource.DisplayWindow(3).Position = [0,40, ...
            ceil(PData(1).Size(2)*PData(1).PDelta(2)/Resource.DisplayWindow(3).pdelta), ... % width
            ceil(PData(1).Size(1)*PData(1).PDelta(1)/Resource.DisplayWindow(3).pdelta)];    % height
        Resource.DisplayWindow(3).Orientation = 'xy';
        Resource.DisplayWindow(3).ReferencePt = [PData(1).Origin(1),...
            -PData(1).Origin(2),85];%PData(1).Region(end).Shape.oPAIntersect];
        Resource.DisplayWindow(3).Colormap = gray(256);
        Resource.DisplayWindow(3).AxesUnits = 'mm';
        % Resource.DisplayWindow(3).mode = '2d';
    case 'secondary'
        Resource.RcvBuffer(1).datatype = 'int16';
        Resource.RcvBuffer(1).rowsPerFrame = 41344;
        Resource.RcvBuffer(1).colsPerFrame = 256;
        Resource.RcvBuffer(1).numFrames = P.numFrames;
end

% Specify Transmit waveform structure.
TW.type = 'parametric';
TW.Parameters = [Trans.frequency,.67,2,1];  % A, B, C, D

% Specify TPC structure.
TPC(1).hv = 7.5;

% Specify TX structure array.

rotang=linspace(0,2*pi,na2); % na=9 creates a plane every 45 degrees (0, 45, 90, 135, 180, 225, 270, 315 and 360).
rotang=rotang(1:end-1); % this eliminates the last angle since 0 = 360

AZ=[];
EL=[];
txAngles=TxAngles.*pi./180;

for i=1:fix(na1/2)
    AZ=[AZ atan(tan(txAngles(i)).*cos(rotang))];
    EL=[EL atan(tan(txAngles(i)).*sin(rotang))];
end

if fix(na1/2) ~= na1/2       % if na is odd, it assumes the script calls for a wave with no steering (plane wave or diverging wave)
    AZ = [AZ 0];
    EL = [EL 0];
end

na=length(AZ);

TX = repmat(struct('waveform', 1, ...
    'Origin', [0.0,0.0,0.0], ...
    'focus', zApex, ...
    'Steer', [0.0,0.0], ...
    'Apod', ones(1,1024), ...
    'Delay', zeros(1,1024),...
    'peakCutOff', 2,...
    'peakBLMax', 20,...
    'TXPD',[]), 1, na);

for n=1:na
    TX(n).Steer = [AZ(n),EL(n) ];
    TX(n).Delay = computeTXDelays(TX(n));
end
%

%- Dummy Transmit -%
TX(na+1)=TX(1);
TX(na+1).Apod(:)=0;

%% Specify TGC Waveform structure.

TGC(1).CntrlPts = [0 256 358 512 665 767 870 900];
TGC(1).rangeMax = 128;
TGC(1).Waveform = computeTGCWaveform(TGC);

%% Receive
temp = (P.endDepth-zApex)*tan(P.viewAngle)+ max([max(Trans.ElementPosWL(:,1)), max(Trans.ElementPosWL(:,2))]);
maxAcqLength = sqrt(P.endDepth^2 + temp^2);
Receive = repmat(struct(...
    'Apod', zeros(1,1024), ...
    'startDepth', 0, ...
    'endDepth', 128/(4*2)*ceil(maxAcqLength/(128/(4*2))), ...
    'TGC', 1, ...
    'mode', 0, ...
    'bufnum', 1, ...
    'framenum', 1, ...
    'acqNum', 1, ...
    'sampleMode', 'NS200BW', ...
    'callMediaFunc', 0), 1, na*P.numFrames);

%Set event specific Receive attributes.

for j = 1:P.numFrames
    Receive(na*(j-1)+1).callMediaFunc = 1;
    for i = 1:na
        if strcmp(RDMAconfig.rdmaRole,'secondary')
            Receive(na*(j-1)+i).Apod((RDMAconfig.nodeIndex-1)*256+[1:256])=1;
        else
            Receive(na*(j-1)+i).Apod(:)=1;
        end
        Receive(na*(j-1)+i).framenum = j;
        Receive(na*(j-1)+i).acqNum = i;
    end
end

%% Reconstruction - Specify Recon structure arrays.
if strcmp(RDMAconfig.rdmaRole,'primary')  %only do this for primary system
    Recon(1) = struct('senscutoff', senscutoff, ...
        'pdatanum', 1, ...
        'rcvBufFrame', -1, ...
        'IntBufDest', [1,1],...
        'ImgBufDest', [1,-1],...
        'RINums', 1:na);

    % Define ReconInfo structures.
    ReconInfo = repmat(struct('mode', 'accumIQ', ...
        'Pre', [],...
        'Post', [],...
        'scaleFactor', 1, ...
        'Aperture', [[128:-1:65, 1:1:64, 192:-1:129, 193:1:256]...
             256+[128:-1:65, 1:1:64, 192:-1:129, 193:1:256]...
             512+[128:-1:65, 1:1:64, 192:-1:129, 193:1:256]...
             768+[128:-1:65, 1:1:64, 192:-1:129, 193:1:256]]), 1, na);

    % - Set specific ReconInfo attributes.
    ReconInfo(1).Pre = 'clearInterBuf';
    for j = 1:na  % For each row in the column
        ReconInfo(j).txnum = j;
        ReconInfo(j).rcvnum = j;
        ReconInfo(j).regionnum = 5; %1 for full 3d Volume 5 for just orthogonal planes
    end
    ReconInfo(end).Post = 'IQ2IntensityImageBuf';

    %% Specify Process structure array.
    Process(1).classname = 'Image';
    Process(1).method = 'imageDisplay';
    Process(1).Parameters = {'imgbufnum',1,...
        'framenum',-1,...
        'pdatanum',1,...
        'srcData','intensity3D',...
        'pgain', pgain,...
        'persistMethod','none',...
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
        'persistMethod','none',...
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
        'persistMethod','none',...
        'persistLevel',persist,...
        'interpMethod','4pt',...
        'compressMethod','power',...
        'compressFactor',CompressionFactor,...
        'mappingMethod','full',...
        'display',1,...
        'displayWindow',3};
end

%% Specify SeqControl structure arrays.
% Jump to start of sequence
JUMP = 1;
SeqControl(JUMP).command = 'jump'; % jump back to start.
SeqControl(JUMP).argument = 1;
% Time to next acquisition
TTNA=2;
SeqControl(TTNA).command = 'timeToNextAcq';  % time between TX's
SeqControl(TTNA).argument = 250;  % 250 usec
TTNAFRAME=3;
SeqControl(TTNAFRAME).command = 'timeToNextAcq';  % time between frames
SeqControl(TTNAFRAME).argument = 200e3;  % 200 msec
% Return to Matlab
RETMAT = 4;
SeqControl(RETMAT).command = 'returnToMatlab';
%-- Multisystem Sync seqcontrol
MSSYNC = 5;
SeqControl(MSSYNC).command = 'multiSysSync';
SeqControl(MSSYNC).condition = 'normal';       % 'normal' uses internal trigger; other options are 'BNC_Rising' and 'BNC_Falling'
SeqControl(MSSYNC).argument = 1000;  % 1s timeout (in ms) 
%
SYNC = 6;
SeqControl(SYNC).command = 'sync'; % synchronize SW and HW sequencers
SeqControl(SYNC).argument = 10e6;     % 5s
%
TRIGO = 7;
SeqControl(TRIGO).command = 'triggerOut';


RDMASYNC = 8;
SeqControl(RDMASYNC).command = 'rdmaSync';
%
TTNADUMMY = 9;
SeqControl(TTNADUMMY).command = 'timeToNextAcq';  % time between dummy TX and start of Frame
SeqControl(TTNADUMMY).argument = 10; % 10us minimum time (in us)

nsc = length(SeqControl)+1; % nsc is count of SeqControl objects


%% Event Sequence
n = 1; % n is count of Events

% Acquire all frames defined in RcvBuffer
for i = 1:Resource.RcvBuffer(1).numFrames
    if (RDMAconfig.nodeIndex==1)
        Event(n).info = 'DUMMY TX';
        Event(n).tx = na+1;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = TTNADUMMY;
        n=n+1;
    end
    for j = 1:na
        Event(n).info = 'Full aperature.';
        Event(n).tx = j;
        Event(n).rcv = na*(i-1)+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = TTNA;
        % -- For Sync between primary and secondary:
        % - sync on 1st TX of each frame
        % - ttna for all systems between each TX
        % - ttnaframe for just primary(when the sync occurs)
        if (j == 1)
            Event(n).seqControl = [TRIGO, MSSYNC, TTNA];
        elseif (j == na)
            if (RDMAconfig.nodeIndex==1)
                Event(n).seqControl = TTNAFRAME;
            else
                Event(n).seqControl = 0;
            end
        else
            Event(n).seqControl = TTNA;
        end
        n = n+1;
    end

    Event(n).info = 'Transfer to Host';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = [nsc, nsc+1];
    SeqControl(nsc).command = 'transferToHost';
    SeqControl(nsc+1).command = 'waitForTransferComplete';
    SeqControl(nsc+1).argument = nsc;
    nsc = nsc + 2;
    n = n+1;

    if strcmp(RDMAconfig.rdmaRole,'secondary') %--- secondary ---%
        %--- this part will be hidden in future ---%
        RDMAcontrol(i).bufferType = 'receive';
        RDMAcontrol(i).primaryDstBufNum = 1;
        RDMAcontrol(i).primaryDstFrameNum = i;
        RDMAcontrol(i).primaryDstStartChNum = 256*(RDMAconfig.nodeIndex-1)+1;
        RDMAcontrol(i).secondarySrcBufNum = 1;
        RDMAcontrol(i).secondarySrcFrameNum  = i;
        RDMAcontrol(i).secondaryIndex  = RDMAconfig.nodeIndex;

        Event(n).info = 'RDMA Transfer';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).processn = 0;
        Event(n).seqControl = nsc;
        SeqControl(nsc).command = 'rdmaWrite';
        SeqControl(nsc).argument = i; % link it to RDMA structure
        SeqControl(nsc).condition = 'asynch'; 
        nsc = nsc + 1;
        n = n+1;

    else %--- primary ---%
        for a = 1:RDMAconfig.numSecondaries
            RDMAcontrol(i,a).bufferType = 'receive';
            RDMAcontrol(i,a).primaryDstBufNum = 1;
            RDMAcontrol(i,a).primaryDstFrameNum = i;
            RDMAcontrol(i,a).primaryDstStartChNum = 256*(a-1)+1;  %only use for receive
            RDMAcontrol(i,a).secondarySrcBufNum = 1;
            RDMAcontrol(i,a).secondarySrcFrameNum  = i;
            RDMAcontrol(i,a).secondaryIndex  = a;
        end

        Event(n).info = 'RDMA Transfer';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = nsc +  [0:3];
        n = n + 1;
        for a = 1:RDMAconfig.numSecondaries
            SeqControl(nsc).command = 'rdmaWrite';
            SeqControl(nsc).argument = (i-1)*RDMAconfig.numSecondaries+a; % link it to RDMA structure
            SeqControl(nsc).condition = 'asynch'; % must match secondary setup
            nsc = nsc + 1;
        end

        Event(n).info = 'Recon';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 1;
        Event(n).process = 0;
        Event(n).seqControl = 0;
        n = n+1;

        Event(n).info = ['Process XZ - Frame ' num2str(i)];
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 1;
        Event(n).seqControl = 0;
        n = n+1;

        Event(n).info = ['Process YZ - Frame ' num2str(i)];
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 2;
        Event(n).seqControl = 0;
        n = n+1;

        Event(n).info = ['Process XY - Frame ' num2str(i)];
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;%3; %debug turn off XY plane for faster processing
        Event(n).seqControl = 0;
        n = n+1;
    end
    
    Event(n).info = 'exit to matlab';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    if floor(i/frameRateFactor) == i/frameRateFactor     % Exit to Matlab every 5th frame
        Event(n).seqControl = RETMAT;
    else
        Event(n).seqControl = RDMASYNC;
    end
    n = n+1;
    
    if RDMAconfig.nodeIndex==1  %HW/SW SYNC on Node 1
        Event(n).info = 'Sync';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = SYNC;
        n = n + 1;
    end
end

Event(n).info = 'Jump back to first event';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = JUMP;

RDMAcontrol=RDMAcontrol';
RDMAcontrol=RDMAcontrol(:);

% User specified UI Control Elements
if strcmp(RDMAconfig.rdmaRole,'primary')
    % - Sensitivity Cutoff
    UI(1).Control =  {'UserB7','Style','VsSlider','Label','Sens. Cutoff',...
        'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
        'SliderStep',[0.025,0.1],'ValueFormat','%1.3f'};
    UI(1).Callback = text2cell('%SensCutoffCallback');

    % - Range Change
    MinMaxVal = [64,300,P.endDepth]; % default unit is wavelength
    AxesUnit = 'wls';

    if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
        if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
            AxesUnit = 'mm';
            MinMaxVal = MinMaxVal * (Resource.Parameters.speedOfSound/1000/Trans.frequency);
        end
    end
    UI(2).Control = {'UserA1','Style','VsSlider','Label',['Range (',AxesUnit,')'],...
        'SliderMinMaxVal',MinMaxVal,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f'};
    UI(2).Callback = text2cell('%RangeChangeCallback');
end


% Save all the structures to a .mat file.
save([tempdir() filename]);
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
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
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

