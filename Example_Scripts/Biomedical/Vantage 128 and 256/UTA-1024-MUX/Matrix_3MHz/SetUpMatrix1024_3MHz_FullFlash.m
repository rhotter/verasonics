%  SetUpMatrix1024_3MHz_FullFlash this is a copy of "SetUpMatrix1024_3MHz_Flash",
%  modified to demonstrate use of the dynamic mux programming ability
%  to set HVMux switches to select elements in parallel.  This script does
%  a single flash transmit over all 1024 elements, with each of four
%  256-element synthetic aperture receive events.
%
% Copyright 2001-2022 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.

% Last update:
% 05/12/2020 - Update to SW 4.3 format for new user UIControls and External
% function definitions (VTS 1691)
%   (~/Example_Scripts/Vantage_Features/New UI Scheme/SetUpL11_5vFlash_NewUI)

clear all;clear functions;

P.startDepth = 0;   % start of Acquisition depth range in wavelengths
P.endDepth = 82; %   use longer lines for noise floor testing
P.viewAngle = 20 * pi/180;                    % angle between z-axis and surface of cone.

P.curved = 1; % set to 0 for flat 1 for defocused at the reconstruction apex

P.numFrames = 2;
P.numTx = 1;
P.numSyntheticTx = 1;
P.numSyntheticRcvs = 4;

RcvProfile.AntiAliasCutoff = 15; % 10, 15, 20, 30
% RcvProfile.LnaGain = 18;
% RcvProfile.LnaZinSel = 31;

%% Specify system parameters.
Resource.Parameters.numTransmit = 256;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 256;   % number of receive channels.
Resource.Parameters.verbose = 2; % 1 enables and 0 disables
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.speedCorrectionFactor = 1.0;
Resource.Parameters.fakeScanhead = 0; % allow system HW operation with nothing connected
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.
Resource.Parameters.Connector = 1;

%% Specify Trans structure array.
Trans.name = 'Matrix1024-3';
Trans.units = 'mm';
Trans = computeTrans(Trans);
Trans.HVMux = computeUTAMux1024; % add the HV Mux programming for the UTA 1024-MUX

% Create Papod arrays to be pasted into the TX and Receive structures, with
% the fifth entry enabling all 1024 elements for a single flash transmit
Papod = [ ones(1,256) zeros(1,256) zeros(1,256) zeros(1,256); ...
         zeros(1,256)  ones(1,256) zeros(1,256) zeros(1,256); ...
         zeros(1,256) zeros(1,256)  ones(1,256) zeros(1,256); ...
         zeros(1,256) zeros(1,256) zeros(1,256)  ones(1,256);
          ones(1, 1024)];
% create aperture index for each Papod
for i = 1:size(Papod, 1)
    Paper(i) = computeMuxAperture(Papod(i, :), Trans);
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
PData(1).PDelta = [1,1,1]*1;
if ~exist('extent','var'), extent = max(max(Trans.ElementPosWL(:,1)),max(Trans.ElementPosWL(:,2))); end
zApex = -extent/tan(P.viewAngle);

PData(1).Size(1) = ceil(2.0*(P.endDepth-zApex)*tan(P.viewAngle)/PData(1).PDelta(2));  if mod(PData(1).Size(1),2)==0, PData(1).Size(1) =  PData(1).Size(1)+1; end
PData(1).Size(2) = PData(1).Size(1);
PData(1).Size(3) = ceil((P.endDepth)/PData(1).PDelta(3));
PData(1).Origin = [-((PData(1).Size(2)-1)/2)*PData(1).PDelta(1), ((PData(1).Size(1)-1)/2)*PData(1).PDelta(1), 0];

PData(1).Region(1) = struct('Shape',struct('Name','Pyramid','Position',[0,0,zApex],'angle',P.viewAngle,'z1',P.startDepth,'z2',P.endDepth));

PData(1).Region = computeRegions(PData(1));

%% Specify Media.  Use point targets in middle of PData.
% Media.MP(1,:) = [0,0,30,1.0];      % single point.
% Media.program = 'PointTargets3D';   % prog. to execute for creating media pts.
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


%% Resources
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 4096*P.numTx*P.numSyntheticTx*P.numSyntheticRcvs;
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = P.numFrames;
Resource.ImageBuffer(1).numFrames = 1;
Resource.ImageBuffer(2).numFrames = 1;
Resource.InterBuffer(1).numFrames = 1;
Resource.InterBuffer(2).numFrames = 1;
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).Title = '3D Flash Image - XZ plane';
Resource.DisplayWindow(1).pdelta = 0.25;
Resource.DisplayWindow(1).Position = [0,580, ...
    ceil(PData(1).Size(2)*PData(1).PDelta(2)/Resource.DisplayWindow(1).pdelta), ... % width
    ceil(PData(1).Size(3)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta)];    % height
Resource.DisplayWindow(1).Orientation = 'xz';
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0.0,0.0];
Resource.DisplayWindow(1).Colormap = gray(256);
Resource.DisplayWindow(1).AxesUnits = 'wavelengths';
% Resource.DisplayWindow(1).mode = '2d';
Resource.DisplayWindow(2).Type = 'Verasonics';
Resource.DisplayWindow(2).Title = '3D Flash Image - YZ plane';
Resource.DisplayWindow(2).pdelta = Resource.DisplayWindow(1).pdelta;
Resource.DisplayWindow(2).Position = [430,580, ...
    ceil(PData(1).Size(1)*PData(1).PDelta(1)/Resource.DisplayWindow(2).pdelta), ... % width
    ceil(PData(1).Size(3)*PData(1).PDelta(3)/Resource.DisplayWindow(2).pdelta)];    % height
Resource.DisplayWindow(2).Orientation = 'yz';
Resource.DisplayWindow(2).ReferencePt = [0,-PData(1).Origin(2),0];
Resource.DisplayWindow(2).Colormap = gray(256);
Resource.DisplayWindow(2).AxesUnits = 'wavelengths';
% Resource.DisplayWindow(2).mode = '2d';
Resource.DisplayWindow(3).Type = 'Verasonics';
Resource.DisplayWindow(3).Title = '3D Flash Image - XY plane';
Resource.DisplayWindow(3).pdelta = Resource.DisplayWindow(1).pdelta;
Resource.DisplayWindow(3).Position = [0,40, ...
    ceil(PData(1).Size(2)*PData(1).PDelta(2)/Resource.DisplayWindow(3).pdelta), ... % width
    ceil(PData(1).Size(1)*PData(1).PDelta(1)/Resource.DisplayWindow(3).pdelta)];    % height
Resource.DisplayWindow(3).Orientation = 'xy';
Resource.DisplayWindow(3).ReferencePt = [PData(1).Origin(1),-PData(1).Origin(2),60];%PData.Region(end).Shape.oPAIntersect];
Resource.DisplayWindow(3).Colormap = gray(256);
Resource.DisplayWindow(3).AxesUnits = 'wavelengths';
% Resource.DisplayWindow(3).mode = '2d';

% Specify Transmit waveform structure.
TW.type = 'parametric';
TW.Parameters = [Trans.frequency,.67,2,1];  % A, B, C, D

% Specify TPC structure.
TPC(1).hv = 10.0;

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'focus', P.curved*zApex, ...
                   'Steer', [0.0,0.0], ...
                   'Apod', zeros(1,Trans.numelements), ...
                   'Delay', zeros(1,Trans.numelements),...
                   'peakCutOff', 2,...
                   'peakBLMax', 20,...
                   'aperture',0), 1,P.numTx*P.numSyntheticTx);

for m = 1:P.numTx
    for i=1:P.numSyntheticTx
        TX(i+(m-1)*P.numSyntheticTx).Apod = Papod(5, :);
        TX(i+(m-1)*P.numSyntheticTx).aperture = Paper(5);
        % flash transmit over all 1024 elements with all delays set to zero
        % so don't need to call computeTXDelays
    end
end

%Specify TGC Waveform structure.
TGC(1).CntrlPts = [0 256 358 512 665 767 870 900];
TGC(1).rangeMax = 128;
TGC(1).Waveform = computeTGCWaveform(TGC);

temp = (P.endDepth-zApex)*tan(P.viewAngle)+ max([max(Trans.ElementPosWL(:,1)), max(Trans.ElementPosWL(:,2))]);
maxAcqLength = sqrt(P.endDepth^2 + temp^2);

%Receive
Receive = repmat(struct(...
                'Apod', zeros(1,Trans.numelements), ...
                'startDepth', 0, ...
                'endDepth', 128/(4*2)*ceil(maxAcqLength/(128/(4*2))), ...
                'TGC', 1, ...
                'mode', 0, ...
                'bufnum', 1, ...
                'framenum', 1, ...
                'acqNum', 1, ...
                'aperture',0, ...
                'sampleMode', 'NS200BW', ...
                'callMediaFunc', 0), 1, P.numTx*P.numFrames*P.numSyntheticTx*P.numSyntheticRcvs);

%Set event specific Receive attributes.
for j = 1:P.numFrames
    for i = 1:P.numSyntheticTx*P.numSyntheticRcvs*P.numTx
        Receive(i+(j-1)*P.numSyntheticTx*P.numSyntheticRcvs*P.numTx).acqNum = i;
        Receive(i+(j-1)*P.numSyntheticTx*P.numSyntheticRcvs*P.numTx).framenum = j;
        Receive(i+(j-1)*P.numSyntheticTx*P.numSyntheticRcvs*P.numTx).Apod = Papod(i, :);
        Receive(i+(j-1)*P.numSyntheticTx*P.numSyntheticRcvs*P.numTx).aperture = Paper(i);
    end
end

%%
senscutoff = 0.6;
% Specify Recon structure arrays.
Recon(1) = struct('senscutoff', senscutoff, ...
                   'pdatanum', 1, ...
                   'rcvBufFrame', -1, ...
                   'IntBufDest', [1,1],...
                   'ImgBufDest', [1,-1], ...
                   'RINums', 1:P.numTx*P.numSyntheticTx*P.numSyntheticRcvs);

% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 4, ...
                   'Pre',[],...
                   'Post', [],...
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'scaleFactor', 1, ...
                   'regionnum', 1), 1, P.numTx*P.numSyntheticTx*P.numSyntheticRcvs);

% - Set specific ReconInfo attributes.
% ReconInfo(1).Pre = 'clearImageBuf';
ReconInfo(1).mode = 3;
for m = 1:P.numTx
    for j = 1:P.numSyntheticTx
        for i=1:P.numSyntheticRcvs
            ReconInfo(i+(j-1)*P.numSyntheticRcvs+(m-1)*P.numSyntheticTx*P.numSyntheticRcvs).txnum = j+(m-1)*P.numSyntheticTx;
            ReconInfo(i+(j-1)*P.numSyntheticRcvs+(m-1)*P.numSyntheticTx*P.numSyntheticRcvs).rcvnum = i+(j-1)*P.numSyntheticRcvs+(m-1)*P.numSyntheticTx*P.numSyntheticRcvs;
        end
    end
end
ReconInfo(end).mode = 5;
% ReconInfo(end).Post = 'IQ2IntensityImageBuf';


% %%
pgain=20;
% Specify Process structure arrays.
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...
                         'framenum',-1,...
                         'pdatanum',1,...
                         'srcData','intensity3D',...
                         'pgain', pgain,...
                         'persistMethod','none',...
                         'persistLevel',30,...
                         'interpMethod','4pt',...
                         'compressMethod','power',...
                         'compressFactor',75,...
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
                         'persistLevel',30,...
                         'interpMethod','4pt',...
                         'compressMethod','power',...
                         'compressFactor',75,...
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
                         'persistLevel',30,...
                         'interpMethod','4pt',...
                         'compressMethod','power',...
                         'compressFactor',75,...
                         'mappingMethod','full',...
                         'display',1,...
                         'displayWindow',3};
Process(4).classname = 'External';
Process(4).method = 'volumetricPlot';
Process(4).Parameters = {'srcbuffer','image',...
                         'srcbufnum',1,...
                         'srcframenum',-1,...
                         'dstbuffer','none'};

% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; % jump back to start
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(2).argument = 250;
SeqControl(3).command = 'timeToNextAcq';  % time between frames
SeqControl(3).argument = 20000;  % 20 msec
SeqControl(3).condition = 'ignore'; %Recon processing time
SeqControl(4).command = 'returnToMatlab';
nsc = 5; % nsc is count of SeqControl objects

n = 1;   % start index for Events
for j = 1:Resource.RcvBuffer(1).numFrames
    for m = 1:P.numTx
        for k = 1:P.numSyntheticTx %Loop # of transmits
            for i = 1:P.numSyntheticRcvs %Loop # of apertures
                Event(n).info = 'Transmit/Receive.';
                Event(n).tx = k+(m-1)*P.numSyntheticTx;   % use 1st TX structure.
                Event(n).rcv = i+(k-1)*P.numSyntheticRcvs+(j-1)*P.numSyntheticRcvs*P.numSyntheticTx*P.numTx+(m-1)*P.numSyntheticRcvs*P.numSyntheticTx;  % use 1st Rcv structure of frame.
                Event(n).recon = 0;      % no reconstruction.
                Event(n).process = 0;    % no processing
                Event(n).seqControl = 2; %
                n = n+1;
            end
        end
    end
    Event(n-1).seqControl = 3;

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

    Event(n).info = ['3D Reconstruct Frame ' num2str(i)];
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
    Event(n).process = 3;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = ['3D Process  - Frame' num2str(i)];
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 4;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = 'Return to Matlab';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = 4;
    n = n+1;

    Event(n).info = 'HW/SW Sync';
    Event(n).tx = 0;         % use 1st TX structure.
    Event(n).rcv = 0;      % use ith Rcv structure.
    Event(n).recon = 0;      % no reconstruction.
    Event(n).process = 0;    % no processing
    Event(n).seqControl = nsc;
        SeqControl(nsc).command = 'sync';
        nsc = nsc + 1;
    n = n+1;
end
Event(n).info = 'Jump';
Event(n).tx = 0;         % no TX structure.
Event(n).rcv = 0;        % no Rcv structure.
Event(n).recon = 0;      % no reconstruction.
Event(n).process = 0;    % call processing function
Event(n).seqControl = 1; %


% User specified UI Control Elements
import vsv.seq.uicontrol.VsSliderControl;

% - Sensitivity Cutoff
UI(1).Control = VsSliderControl('LocationCode','UserB7','Label','Sens. Cutoff',...
                  'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                  'SliderStep',[0.025,0.1],'ValueFormat','%1.3f', ...
                  'Callback', @SensCutoffCallback );             
                  
% - Section number change
UI(2).Control = VsSliderControl('LocationCode','UserB1','Label','Z-Section',...
                  'SliderMinMaxVal',[1,P.endDepth,Resource.DisplayWindow(3).ReferencePt(3)],...
                  'SliderStep',[1/(P.endDepth-1) 10/(P.endDepth-1)],'ValueFormat','%3.0f', ...
                  'Callback', @ZSectionChange );

% - Section number change
UI(3).Control = VsSliderControl('LocationCode','UserB2','Label','Y-Section',...
                  'SliderMinMaxVal',[PData(1).Origin(2)-(PData(1).Size(1)-1)*PData(1).PDelta(1), PData(1).Origin(2), Resource.DisplayWindow(1).ReferencePt(2)],...
                  'SliderStep',[1/((PData(1).Size(2)-1)*PData(1).PDelta(2)-1) 10/((PData(1).Size(2)-1)*PData(1).PDelta(2)-1)],'ValueFormat','%3.0f', ...
                  'Callback', @YSectionChange); 

% - Section number change
UI(4).Control = VsSliderControl('LocationCode','UserB3','Label','X-Section',...
                  'SliderMinMaxVal',[PData(1).Origin(1), PData(1).Origin(1)+(PData(1).Size(2)-1)*PData(1).PDelta(2), Resource.DisplayWindow(2).ReferencePt(1)],...
                  'SliderStep',[1/((PData(1).Size(1)-1)*PData(1).PDelta(1)-1) 10/((PData(1).Size(1)-1)*PData(1).PDelta(1)-1)],'ValueFormat','%3.0f', ...
                  'Callback', @XSectionChange); 

% - compression
CompressionFactor=1.3;
UI(5).Control = VsSliderControl('LocationCode','UserA2','Label','CompressionP',...
                  'SliderMinMaxVal',[1,11,CompressionFactor],...
                  'SliderStep',[1/40,4/40],'ValueFormat','%1.1f', ...
                  'Callback', @CompressionP ); 

%external function
EF(1).Function = vsv.seq.function.ExFunctionDef('volumetricPlot', @volumetricPlot);


% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 3;

% Save all the structures to a .mat file.
save('MatFiles/Matrix1024-3MHz-FullFlash.mat'); %filename = 'MatFiles/Matrix1024-8FullFlashDM.mat'; %VSX;


%% **** Callback routines used by UIControls (UI)  ****

function SensCutoffCallback(~,~,UIValue)
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
end


function YSectionChange(~,~,UIValue)
%YSectionChange
    RS = evalin('base','Resource');
    RefPt = RS.DisplayWindow(1).ReferencePt;
    RefPt(2) = UIValue;
    RS.DisplayWindow(1).ReferencePt = RefPt;
    assignin('base','Resource',RS);

    PDataT = evalin('base','PData');
    P = evalin('base','P');

    PDataT(1).Region(1).Shape.oPAIntersect = RefPt(2);

    PDataT(1).Region = computeRegions(PDataT(1));
    assignin('base','PData',PDataT);

    Control = repmat(struct('Command',[],'Parameters',[]),1,2);
    Control(1).Command = 'set&Run';
    Control(1).Parameters = {'DisplayWindow',1,'ReferencePt',RefPt};
    Control(2).Command = 'update&Run';
    Control(2).Parameters = {'PData'};
    assignin('base','Control', Control);
end


function XSectionChange(~,~,UIValue)
%XSectionChange
    RS = evalin('base','Resource');
    RefPt = RS.DisplayWindow(2).ReferencePt;
    RefPt(1) = UIValue;
    RS.DisplayWindow(2).ReferencePt = RefPt;
    assignin('base','Resource',RS);

    PDataT = evalin('base','PData');
    P = evalin('base','P');

    PDataT(1).Region(1).Shape.oPAIntersect = RefPt(1);

    PDataT(1).Region = computeRegions(PDataT(1));
    assignin('base','PData',PDataT);

    Control = repmat(struct('Command',[],'Parameters',[]),1,2);
    Control(1).Command = 'set&Run';
    Control(1).Parameters = {'DisplayWindow',2,'ReferencePt',RefPt};
    Control(2).Command = 'update&Run';
    Control(2).Parameters = {'PData'};
    assignin('base','Control', Control);
end

function ZSectionChange(~,~,UIValue)
%ZSectionChange
    RS = evalin('base','Resource');
    RefPt = RS.DisplayWindow(3).ReferencePt;
    RefPt(3) = UIValue;
    RS.DisplayWindow(3).ReferencePt = RefPt;
    assignin('base','Resource',RS);

    PDataT = evalin('base','PData');
    P = evalin('base','P');

    PDataT(1).Region(1).Shape.oPAIntersect = RefPt(3);

    PDataT(1).Region = computeRegions(PDataT(1));
    assignin('base','PData',PDataT);

    Control = repmat(struct('Command',[],'Parameters',[]),1,2);
    Control(1).Command = 'set&Run';
    Control(1).Parameters = {'DisplayWindow',3,'ReferencePt',RefPt};
    Control(2).Command = 'update&Run';
    Control(2).Parameters = {'PData'};
    assignin('base','Control', Control);
end

function CompressionP(~,~,UIValue)
%CompressionP
    CompressionFactor = UIValue;
    assignin('base','CompressionFactor',CompressionFactor);
end

%% **** Callback routines used by External function definition (EF) ****

function volumetricPlot(ImgData)

    persistent handle3D

    PData = evalin('base','PData');
    Trans = evalin('base','Trans');
    P = evalin('base','P');
    v_depth = 1;
    CFactor = evalin('base','CompressionFactor');
    ImgData=ImgData.^(1/CFactor);
    ImgData(:,:,1:v_depth)=0;
    ImgData = flipdim(ImgData,3);

    if isempty(handle3D) || ~ishandle(handle3D)
        figure('Position', [430, 50, 450, 450]);
        handle3D = axes;
    end

    set(handle3D,'NextPlot','replacechildren');
    vol3d('cdata', ImgData,'texture', '2D','Parent', handle3D);
    grid(handle3D, 'on'); colormap(handle3D,'gray');
    set(handle3D,'NextPlot','add');

    xl=PData(1).Size(2);
    yl=PData(1).Size(1);
    zl=PData(1).Size(3);
    dx=PData(1).PDelta(2);
    dy=PData(1).PDelta(1);
    dz=PData(1).PDelta(3);
    plot3(handle3D,Trans.ElementPosWL(:,1)./dx+(xl-1)/2,Trans.ElementPosWL(:,2)./dy+(yl-1)/2,(Trans.ElementPosWL(:,3)+P.endDepth)./dz,'k.');
    view(handle3D,-45,30);
end

