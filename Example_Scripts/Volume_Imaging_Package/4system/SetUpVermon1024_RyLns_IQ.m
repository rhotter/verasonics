%% Copyright 2021 Verasonics, Inc. Verasonics Registered U.S. Patent and Trademark Office.
% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: SetUpVermon1024_RyLns_IQ.m
%
% Description: A volume imaging sequence that transmits 7 x 8 (56) rayLines
%              across the 2D aperture
%

clear all

Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.
filename = 'Vermon1024_RyLns_IQ';

%=== Settings ===%
% Bmode properties
XYdepth = 10;            % depth of C plane imaging plane
pgain = 20;              % digital gain
persist = 30;            % persistance level
senscutoff = 0.6;        % sensitivity cutoff
CompressionFactor = 75;  % compression Factor
tn = 7; %MUST BE ODD.  NUMBER OF ELEMENTS USED FOR APODIZATION IN EACH DIRECTION IE 5X5, 3X3,
P.numTx = tn^2;   % no. of elements in TX aperture.(5x5)
P.numRays = 7*8; % no. of rays in frame (7x8)
P.startDepth = 0;   % start of Acquisition depth range in wavelengths
P.endDepth = 128; %   (128lambda = 8fps)
P.viewAngle = 14 * pi/180;                    % angle between z-axis and surface of cone.
P.curved = 1; % set to 0 for flat 1 for defocused at the reconstruction apex
P.numFrames = 10;
TXfocus = 110;

RcvProfile.AntiAliasCutoff = 15; % 10, 15, 20, 30
% RcvProfile.LnaGain = 18;
% RcvProfile.LnaZinSel = 31;

frameRateFactor = 5;% Specify factor for converting sequenceRate to frameRate.

%%===================== RDMA SETUP ======================================%%
usingMultiSys = 1;
[RDMAconfig, RDMAsystem]=vsv.multi.generateConfig(4);
%%=======================================================================%%


% Define system parameters.
Resource.Parameters.numTransmit = 256;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 256;    % number of receive channels.
Resource.Parameters.speedCorrectionFactor = 1.0;
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 0;
Resource.Parameters.initializeOnly = 0;

%% Specify Trans structure array.
Trans.name = 'Matrix1024-3';
Trans.units = 'wavelengths';  % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);  % computeTrans is used for known transducers.
Trans.id = hex2dec('0D0100'); % need this still
Trans.ConnectorES=zeros(1,1024)';
Trans.ConnectorES(RDMAconfig.nodeIndex*256+[1:256])=[128:-1:65, 1:1:64, 192:-1:129, 193:1:256];

%Intermediate Variables
waveLength = (Resource.Parameters.speedOfSound/1000)/Trans.frequency;
if strcmp(Trans.units,'mm')
    Trans.ElementPosMm = Trans.ElementPos;
    Trans.ElementPosWL = Trans.ElementPos./waveLength;
else
    Trans.ElementPosMm = Trans.ElementPos.*waveLength;
    Trans.ElementPosWL = Trans.ElementPos;
end


%patterns are based on 5 elements
patternA=[3 8 13 16 20 25 30];  %center elements (center overlap) 7 centers
patternB=[3 6 11 16 21 26 30];  %center elements (outer overlap) 7 centers
%00000000 01111111 11122222 22222333
%12345678_90123456_78901234_56789012
%  ^  ^     ^  ^     ^  ^     ^  ^
patternC=([3 6 11 14 19 22 27 30]-1)*32; %8 centers, 2 per connector set of channels
origElements=repmat( patternA,8,1)+repmat(patternC,7,1)';


%% PData
PData(1).PDelta = [1,1,1]*1;
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
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 4096*P.numRays;
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = P.numFrames;
Resource.InterBuffer(1).numFrames = 1;  % one buffer
%--- MultiSys ---%
if strcmp(RDMAconfig.rdmaRole,'primary')
    Resource.InterBuffer(2).numFrames = 1;  % buffer for secondary 1 IQ
    Resource.InterBuffer(3).numFrames = 1;  % buffer for secondary 2 IQ
    Resource.InterBuffer(4).numFrames = 1;  % buffer for secondary 3 IQ
end
Resource.ImageBuffer(1).numFrames = P.numFrames;

if strcmp(RDMAconfig.rdmaRole,'primary')
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
    %-------------------------------------------------%
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
    %-------------------------------------------------%
    Resource.DisplayWindow(3).Type = 'Verasonics';
    Resource.DisplayWindow(3).Title = '3D Flash Image - XY plane';
    Resource.DisplayWindow(3).pdelta = Resource.DisplayWindow(1).pdelta;
    Resource.DisplayWindow(3).Position = [0,40, ...
        ceil(PData(1).Size(2)*PData(1).PDelta(2)/Resource.DisplayWindow(3).pdelta), ... % width
        ceil(PData(1).Size(1)*PData(1).PDelta(1)/Resource.DisplayWindow(3).pdelta)];    % height
    Resource.DisplayWindow(3).Orientation = 'xy';
    Resource.DisplayWindow(3).ReferencePt = [PData(1).Origin(1),...
        -PData(1).Origin(2),XYdepth];%PData(1).Region(end).Shape.oPAIntersect];
    Resource.DisplayWindow(3).Colormap = gray(256);
    Resource.DisplayWindow(3).AxesUnits = 'mm';
    % Resource.DisplayWindow(3).mode = '2d';
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
    'Apod', zeros(1,1024), ...
    'Delay', zeros(1,1024), ...
    'TXPD', [], ...
    'peakCutOff', 2, ...
    'peakBLMax', 20), 1,P.numRays);

tic
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
    %     disp(['calculating TXPD ' num2str(n) ' out of ' num2str(length(TX)) ' ...'])
    %     TX(n).TXPD = computeTXPD(TX(n),PData(1));
    %     disp(['remainging time: ' num2str(round((toc/n)*length(TX)-toc)) ' seconds'])
end

% Dummy TX
TX(P.numRays+1)=TX(1);
TX(P.numRays+1).Apod(:)=0;

%% Specify TGC Waveform structure.

TGC(1).CntrlPts = [0 256 358 512 665 767 870 900];
TGC(1).rangeMax = 128;
TGC(1).Waveform = computeTGCWaveform(TGC);

temp = (P.endDepth-zApex)*tan(P.viewAngle)+ max([max(Trans.ElementPosWL(:,1)), max(Trans.ElementPosWL(:,2))]);
maxAcqLength = sqrt(P.endDepth^2 + temp^2);


%% Receive
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
    'callMediaFunc', 0), 1, P.numRays*P.numFrames);

% Set event specific Receive attributes.

for j = 1:P.numFrames
    Receive(P.numRays*(j-1)+1).callMediaFunc = 1;
    for i = 1:P.numRays
        Receive(P.numRays*(j-1)+i).Apod(RDMAconfig.nodeIndex*256+[1:256]) = 1;
        Receive(P.numRays*(j-1)+i).framenum = j;
        Receive(P.numRays*(j-1)+i).acqNum = i;
    end
end

%% Reconstruction - Specify Recon structure arrays.
Recon(1) = struct('senscutoff', senscutoff, ...
    'pdatanum', 1, ...
    'rcvBufFrame', -1, ...
    'IntBufDest', [1,1],...
    'RINums', 1:P.numRays);

% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 'accumIQ', ...
    'Pre', [],...
    'Post', [],...
    'txnum', 1, ...
    'rcvnum', 1, ...
    'scaleFactor', 1, ...
    'regionnum', 1), 1, P.numRays);

% - Set specific ReconInfo attributes.
ReconInfo(1).Pre = 'clearInterBuf';
for j = 1:P.numRays
    ReconInfo(j).txnum = j;
    ReconInfo(j).rcvnum = j;
    ReconInfo(j).regionnum = 5;
end
%% Specify Process structure array.
if strcmp(RDMAconfig.rdmaRole,'primary')
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
    Process(4).classname = 'Process';
    Process(4).method = 'sumIQNormMag';
    Process(4).Parameters = {'IntBufSrc',[1,1,2,1,3,1,4,1],...
        'ImgBufDest',[1,-1]};
end

%% Specify SeqControl structure arrays.
% Jump to start of sequence
JUMP = 1;
SeqControl(JUMP).command = 'jump'; % jump back to start
SeqControl(JUMP).argument = 1;
%
TTNA = 2;
SeqControl(TTNA).command = 'timeToNextAcq';  % time between acquisitions
SeqControl(TTNA).argument = 250;  % 250 usec
TTNAFRAME = 3;
SeqControl(TTNAFRAME).command = 'timeToNextAcq';  % time between frames
SeqControl(TTNAFRAME).argument = 200000;  % 200 msec
% Return to Matlab
RETMAT = 4;
SeqControl(RETMAT).command = 'returnToMatlab';
% Multisystem Sync seqcontrol
MSSYNC = 5;
SeqControl(MSSYNC).command = 'multiSysSync';
SeqControl(MSSYNC).condition = 'normal';       % 'normal' uses internal trigger; other options are 'BNC_Rising' and 'BNC_Falling'
SeqControl(MSSYNC).argument = 1000;  % 1 sec
%
SYNC = 6;
SeqControl(SYNC).command = 'sync'; % synchronize SW and HW sequencers
SeqControl(SYNC).argument = 20e3;     % 20msec
%
TRIGO = 7;
SeqControl(TRIGO).command = 'triggerOut';
%
RDMASYNC = 8;
SeqControl(RDMASYNC).command = 'rdmaSync';
%
TTNADUMMY = 9;
SeqControl(TTNADUMMY).command = 'timeToNextAcq';  % time between dummy TX and start of Frame
SeqControl(TTNADUMMY).argument = 10; % 10us minimum time

nsc = length(SeqControl)+1; % nsc is count of SeqControl objects


%% Event Sequence
n = 1; % n is count of Events

% Acquire all frames defined in RcvBuffer
for i = 1:Resource.RcvBuffer(1).numFrames
    if strcmp(RDMAconfig.rdmaRole,'primary')
        Event(n).info = 'DUMMY TX';
        Event(n).tx = P.numRays+1;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = TTNADUMMY;
        n=n+1;
    end
    
    for j = 1:P.numRays
        Event(n).info = 'Transmit/Receive.';
        Event(n).tx = j;   % use 1st TX structure.
        Event(n).rcv = (i-1)*P.numRays+j;  % use 1st Rcv structure of frame.
        Event(n).recon = 0;
        Event(n).process = 0;
        % -- For Sync between primary and secondary:
        % - sync on 1st TX of each frame
        % - ttna for all systems between each TX
        % - ttnaframe for just primary(when the sync occurs)
        if (j == 1)
            Event(n).seqControl = [TRIGO, MSSYNC, TTNA];
        elseif (j == P.numRays)
            if strcmp(RDMAconfig.rdmaRole,'primary')
                Event(n).seqControl = TTNAFRAME;
            else
                Event(n).seqControl = 0;
            end
        else
            Event(n).seqControl = TTNA;
        end
        n = n+1;
    end
    
    Event(n).info = 'Transfer To Host';
    Event(n).tx = 0;         % use 1st TX structure.
    Event(n).rcv = 0;    % use 1st Rcv structure of frame.
    Event(n).recon = 0;      % no reconstruction.
    Event(n).process = 0;    % no processing
    Event(n).seqControl = [nsc,nsc+1]; % set wait time and transfer data
    SeqControl(nsc).command = 'transferToHost';
    SeqControl(nsc+1).command = 'waitForTransferComplete';
    SeqControl(nsc+1).argument = nsc;
    nsc = nsc + 2;
    n = n+1;
    
    Event(n).info = 'Reconstruct';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 0;
    Event(n).seqControl = 0;
    n = n+1;
    
    if strcmp(RDMAconfig.rdmaRole,'secondary')    %--- secondary ---%
        %--- this part will be hidden in future ---%
        RDMAcontrol(1).bufferType = 'inter';
        RDMAcontrol(1).primaryDstBufNum = RDMAconfig.nodeIndex+1;
        RDMAcontrol(1).primaryDstFrameNum = 1;
        RDMAcontrol(1).secondarySrcBufNum = 1;
        RDMAcontrol(1).secondarySrcFrameNum  = 1;
        RDMAcontrol(1).secondaryIndex = RDMAconfig.nodeIndex;
        
        Event(n).info = 'RDMA Transfer';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = nsc;
        SeqControl(nsc).command = 'rdmaWrite';
        SeqControl(nsc).argument = 1; % link it to RDMA structure
        nsc = nsc + 1;
        n = n + 1;
        
        Event(n).info = 'Sync';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = RDMASYNC;
        n = n+1;
        
    else %--- primary ---%
        %Generate 3 RDMA structures (3 secondaries)
        for a = 1:RDMAconfig.numSecondaryNodes
            RDMAcontrol(a).bufferType = 'inter';
            RDMAcontrol(a).primaryDstBufNum = a + 1;
            RDMAcontrol(a).primaryDstFrameNum = 1;
            RDMAcontrol(a).secondarySrcBufNum = 1;
            RDMAcontrol(a).secondarySrcFrameNum  = 1;
            RDMAcontrol(a).secondaryIndex = a;
        end
        
        Event(n).info = 'Sync';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = RDMASYNC;
        n = n+1;
        
        Event(n).info = 'Combine IQ buffers';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 4;
        Event(n).seqControl = 0;
        n = n+1;
        
        %--- Only perform the process and display on the primary system ---%
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
        Event(n).process = 3;
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
    
    if strcmp(RDMAconfig.rdmaRole,'primary')
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

% - Sensitivity Cutoff
UI(1).Control =  {'UserB7','Style','VsSlider','Label','Sens. Cutoff',...
    'SliderMinMaxVal',[0, 1.0, Recon(1).senscutoff],...
    'SliderStep',[0.025,0.1],'ValueFormat','%1.3f'};
UI(1).Callback = text2cell('%SensCutoff');


% - Section number change
UI(2).Control = {'UserB1','Style','VsSlider','Label','Z-Section',...
    'SliderMinMaxVal',[1, P.endDepth, XYdepth],...
    'SliderStep',[1/(P.endDepth-1) 10/(P.endDepth-1)],'ValueFormat','%3.0f'};
UI(2).Callback = text2cell('%ZSectionChange');

% Save all the structures to a .mat file.
%NOTE:  you must use '-v7.3' format flag or TX won't save becase it is a
%large structure because of TXPD
save([tempdir() filename],'-v7.3','-nocompression');
return

% **** Callback routines to be converted by text2cell function. ****
%SensCutoff - Sensitivity cutoff change
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
%SensCutoff


%ZSectionChange
if strcmp(evalin('base','RDMAconfig.rdmaRole'),'primary')
    RS = evalin('base','Resource');
    RefPt = RS.DisplayWindow(3).ReferencePt;
    RefPt(3) = UIValue;
    RS.DisplayWindow(3).ReferencePt = RefPt;
    assignin('base','Resource',RS);
    
    PDataT = evalin('base','PData');
    
    PDataT.Region(4).Shape.oPAIntersect = RefPt(3);
    
    PDataT.Region = computeRegions(PDataT);
    PDataT.Region(5).PixelsLA = unique([PDataT.Region(2).PixelsLA;PDataT.Region(3).PixelsLA ;PDataT.Region(4).PixelsLA]);
    PDataT.Region(5).Shape.Name = 'Custom';
    PDataT.Region(5).numPixels = length(PDataT.Region(5).PixelsLA);
    
    assignin('base','PData',PDataT);
    
    Control = repmat(struct('Command',[],'Parameters',[]),1,2);
    Control(1).Command = 'set&Run';
    Control(1).Parameters = {'DisplayWindow',3,'ReferencePt',RefPt};
    Control(2).Command = 'update&Run';
    Control(2).Parameters = {'PData'};
    assignin('base','Control', Control);
else
    PDataT = evalin('base','PData');
    RefPt(3) = UIValue;
    
    PDataT.Region(4).Shape.oPAIntersect = RefPt(3);
    
    PDataT.Region = computeRegions(PDataT(1));
    PDataT.Region(5).PixelsLA = unique([PDataT.Region(2).PixelsLA;PDataT.Region(3).PixelsLA ;PDataT.Region(4).PixelsLA]);
    PDataT.Region(5).Shape.Name = 'Custom';
    PDataT.Region(5).numPixels = length(PDataT.Region(5).PixelsLA);
    
    assignin('base','PData',PDataT);
    
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'PData'};
    assignin('base','Control', Control);
end
%ZSectionChange

