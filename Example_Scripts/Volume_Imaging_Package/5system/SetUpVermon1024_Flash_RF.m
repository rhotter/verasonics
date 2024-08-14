%% Copyright 2021 Verasonics, Inc. Verasonics Registered U.S. Patent and Trademark Office.
% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: SetUpVermon1024Flash_RF.m - Example of imaging with single plane wave
%                                 transmit
% Description:
%  Single Plane Wave (flash) image with RF transfer
% This version runs faster because it is only reconstructing the planes
% that are being displayed
%
% Last update: 8/6/2019 BWC


clear all

filename = 'Vermon1024_Flash_RF';
%=== Settings ===%
% Bmode properties
XYdepth = 10;            % depth of C plane imaging plane
pgain = 20;              % digital gain
persist = 30;            % persistance level
senscutoff = 0.6;        % sensitivity cutoff
CompressionFactor = 75;  % compression Factor

P.numTx = 1;
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
        Resource.RcvBuffer(1).rowsPerFrame = 2048;
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
        Resource.RcvBuffer(1).rowsPerFrame = 2048;
        Resource.RcvBuffer(1).colsPerFrame = 256;
        Resource.RcvBuffer(1).numFrames = P.numFrames;
end

% Specify Transmit waveform structure.
TW.type = 'parametric';
TW.Parameters = [Trans.frequency,.67,2,1];  % A, B, C, D

% Specify TPC structure.
TPC(1).hv = 7.5;

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
    'Origin', [0.0,0.0,0.0], ...
    'focus', P.curved*zApex, ...
    'Steer', [0.0,0.0], ...
    'Apod', ones(1,1024), ...
    'Delay', zeros(1,1024),...
    'TXPD', [], ...
    'peakCutOff', 2,...
    'peakBLMax', 20), 1,P.numTx);

% Compute transmit delays
TX(1).Delay = computeTXDelays(TX(1));

% Compute transmit pixel data
%TX(1).TXPD = computeTXPD(TX(1),PData(1));
%- Dummy Transmit -%
TX(2)=TX(1);
TX(2).Apod(:)=0;


%-------------------------------------------------------------------------%

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
    'callMediaFunc', 0), 1, P.numFrames);

%Set event specific Receive attributes.
for j = 1:P.numFrames
    Receive(j).callMediaFunc = 1;
    if strcmp(RDMAconfig.rdmaRole,'secondary')
        Receive(j).Apod((RDMAconfig.nodeIndex-1)*256+[1:256])=1;
    else
        Receive(j).Apod(:)=1;
    end
    Receive(j).framenum = j;
    Receive(j).acqNum = 1;
end

%% Reconstruction - Specify Recon structure arrays.
if strcmp(RDMAconfig.rdmaRole,'primary')  %only do this for primary system
    Recon(1) = struct('senscutoff', senscutoff, ...
        'pdatanum', 1, ...
        'rcvBufFrame', -1, ...
        'ImgBufDest', [1,-1],...
        'RINums', 1);

    ReconInfo = struct('mode', 'replaceIntensity', ...
        'txnum', 1, ...
        'rcvnum', 1, ...
        'regionnum', 5, ...
        'Aperture', [[128:-1:65, 1:1:64, 192:-1:129, 193:1:256]...
             256+[128:-1:65, 1:1:64, 192:-1:129, 193:1:256]...
             512+[128:-1:65, 1:1:64, 192:-1:129, 193:1:256]...
             768+[128:-1:65, 1:1:64, 192:-1:129, 193:1:256]]);

    %
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
TTNAFRAME = 2;
SeqControl(TTNAFRAME).command = 'timeToNextAcq';  % time between frames
SeqControl(TTNAFRAME).argument = 200e3; % 200 msec 
% Return to Matlab
RETMAT = 3;
SeqControl(RETMAT).command = 'returnToMatlab';
%-- Multisystem Sync seqcontrol
MSSYNC = 4;
SeqControl(MSSYNC).command = 'multiSysSync';
SeqControl(MSSYNC).condition = 'normal';       % 'normal' uses internal trigger; other options are 'BNC_Rising' and 'BNC_Falling'
SeqControl(MSSYNC).argument = 1000;  % 1 sec
%
SYNC = 5;
SeqControl(SYNC).command = 'sync'; % synchronize SW and HW sequencers
SeqControl(SYNC).argument = 10e6;     % 10 sec
%
TRIGO = 6;
SeqControl(TRIGO).command = 'triggerOut';
%
RDMASYNC = 7;
SeqControl(RDMASYNC).command = 'rdmaSync';
%
TTNADUMMY = 8;
SeqControl(TTNADUMMY).command = 'timeToNextAcq';  % time between dummy TX and start of Frame
SeqControl(TTNADUMMY).argument = 10; % 10us minimum time

nsc = length(SeqControl)+1; % nsc is count of SeqControl objects


%% Event Sequence
n = 1; % n is count of Events
% Acquire all frames defined in RcvBuffer
for j = 1:Resource.RcvBuffer(1).numFrames

    if (RDMAconfig.nodeIndex==1)
        Event(n).info = 'DUMMY TX';
        Event(n).tx = 2;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = TTNADUMMY;
        n=n+1;
    end
    
    Event(n).info = 'HW Sync';
    Event(n).tx = 1;
    Event(n).rcv = j;
    Event(n).recon = 0;
    Event(n).process = 0;
     if RDMAconfig.nodeIndex==1
         Event(n).seqControl = [TRIGO, MSSYNC, TTNAFRAME];  %use only the primary for timing of the frames
     else
         Event(n).seqControl = [TRIGO, MSSYNC];  %secondary systems just wait for primary system to give sync signal
    end
    n = n+1;
    
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
        RDMAcontrol(j).bufferType = 'receive';
        RDMAcontrol(j).primaryDstBufNum = 1;
        RDMAcontrol(j).primaryDstFrameNum = j;
        RDMAcontrol(j).primaryDstStartChNum = 256*(RDMAconfig.nodeIndex-1)+1;
        RDMAcontrol(j).secondarySrcBufNum = 1;
        RDMAcontrol(j).secondarySrcFrameNum  = j;
        RDMAcontrol(j).secondaryIndex  = RDMAconfig.nodeIndex;

        Event(n).info = 'RDMA Transfer';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).processn = 0;
        Event(n).seqControl = nsc;
        SeqControl(nsc).command = 'rdmaWrite';
        SeqControl(nsc).argument = j; % link it to RDMA structure
        SeqControl(nsc).condition = 'asynch';
        nsc = nsc + 1;
        n = n+1;

    else %--- primary ---%
        for a = 1:RDMAconfig.numSecondaries
            RDMAcontrol(j,a).bufferType = 'receive';
            RDMAcontrol(j,a).primaryDstBufNum = 1;
            RDMAcontrol(j,a).primaryDstFrameNum = j;
            RDMAcontrol(j,a).primaryDstStartChNum = 256*(a-1)+1;  %only use for receive
            RDMAcontrol(j,a).secondarySrcBufNum = 1;
            RDMAcontrol(j,a).secondarySrcFrameNum  = j;
            RDMAcontrol(j,a).secondaryIndex  = a;
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
            SeqControl(nsc).argument = (j-1)*RDMAconfig.numSecondaries+a; % link it to RDMA structure
            SeqControl(nsc).condition = 'asynch';
            nsc = nsc + 1;
        end

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
        Event(n).process = 0;%3; %debug turn off XY plane for faster processing
        Event(n).seqControl = 0;
        n = n+1;
    end

    Event(n).info = 'exit to matlab';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    if floor(j/frameRateFactor) == j/frameRateFactor     % Exit to Matlab every 5th frame
        Event(n).seqControl = RETMAT;
    else
        Event(n).seqControl = RDMASYNC;
    end
    n = n+1;
    
    if RDMAconfig.nodeIndex==1
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
