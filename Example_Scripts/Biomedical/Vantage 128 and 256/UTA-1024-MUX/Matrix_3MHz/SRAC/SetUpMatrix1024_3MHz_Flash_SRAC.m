% SetUpMatrix1024_3MHz_Flash_SRAC.m - Example of Flash imaging for volume
% acquisition using a matrix probe array (3 MHz) and Sparse Random Aperture 
% Compounding (SRAC) technique.
%
% For more information about this technique please refer to the article
% published in Physics in Medicine and Biology in May 2020:
% Phys Med Biol. 2020 May 15. doi: 10.1088/1361-6560/ab9372
%
% Description:
%  This scripts transmits and receives a single diverging wave for volume recontruciton. 
%  The user can select the number of SRAC (P.numSyntheticRcvs), form 1 to 4, that want 
%  to use for compounding. If this number is different than "1", then the script will 
%  transmit the same diverging wave multiple times using the same Random aperture, 
%  however the recieve of the backscatter data will be done with different random 
%  apertures. The random apertures used in this script are complementary
%  to each other (meaning that within a set of apertures, 1-4, no element is repeated). 
%  The script loads a file named "saved100CompRndApod" which has
%  100 sets with 4 complementary random apertures each.
%
% Copyright 2001-2022 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.

%  The script allows the user to choose between reconstructing 3 orthogonal 
%  planes (XZ, YZ and XY) or the whole volume. The three orthogonal planes 
%  were chosen by default in the middle of the probe for the XZ and YZ planes, 
%  and at a depth of 60 wavelengths for the XY plane. 
%  To set the orthogonal plane reconstruction the user must set the variable
%  "ReconRegion" equal to 5. And for the whole volume reconstruction ReconRegion=1. 

% Updates:
% 06/23/2020 - Update to SW 4.3 format for new UIControls and function definitions (VTS 1691). 
%   More info:(.../Example_Scripts/Vantage_Features/New UI Scheme/SetUpL11_5vFlash_NewUI)


clear 

P.startDepth = 0;   % start of Acquisition depth range in wavelengths
P.endDepth = 170; %   use longer lines for noise floor testing
P.viewAngle = 20 * pi/180;                    % angle between z-axis and surface of cone.


P.numFrames = 1;
P.numSyntheticRcvs = 4;

RcvProfile.AntiAliasCutoff = 15; % 10, 15, 20, 30
% RcvProfile.LnaGain = 18;
% RcvProfile.LnaZinSel = 31;

ReconRegion=5; %Set ReconRegion=5 to reconstruct only 3 orthogonal planes, or ReconRegion =1 to reconstruct the whole volume

%% Specify system parameters.
Resource.Parameters.numTransmit = 256;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 256;   % number of receive channels.
Resource.Parameters.verbose = 2; % 1 enables and 0 disables
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.speedCorrectionFactor = 1.0;
Resource.Parameters.fakeScanhead = 0; % allow system HW operation with nothing connected
Resource.Parameters.simulateMode = 0;
Resource.Parameters.Connector = 1;
Resource.VDAS.dmaTimeout=1000;
Resource.VDAS.useNewHal=1;
Resource.VDAS.dmaTimeout=2000;


%% Specify Trans structure array.
Trans.name = 'Matrix1024-3';
Trans.units = 'mm';
Trans = computeTrans(Trans);
Trans.HVMux = computeUTAMux1024; % add the HV Mux programming for the UTA 1024-MUX

% this files contains 100 sets of 4 Rnd Apertures that are complementary. 
% CompRndApod is matrix of 100x4x1024. The sum of every 4 matrices (of the 
% second dimention) account for the full aperture. For this script (Flash)
% only the first set is needed. 

load('saved100CompRndApod'); 
CompRndApod=CompRndApod(1,:,:);

% create aperture index for each Papod
for j=1:4
    Paper(1,j) = computeMuxAperture(squeeze(CompRndApod(1,j,:))', Trans);
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

PData(1).Region(2) = struct('Shape',struct('Name','Slice','Orientation','yz','oPAIntersect',PData.Origin(1)+PData.PDelta(2)*(PData.Size(2)-1)/2,'andWithPrev',1));
PData(1).Region(3) = struct('Shape',struct('Name','Slice','Orientation','xz','oPAIntersect',PData.Origin(2)-PData.PDelta(1)*(PData.Size(1)-1)/2,'andWithPrev',1));
PData(1).Region(4) = struct('Shape',struct('Name','Slice','Orientation','xy','oPAIntersect',60,'andWithPrev',1));

[PData(1).Region] = computeRegions(PData(1));

PData.Region(5).PixelsLA = unique([PData.Region(2).PixelsLA; PData.Region(3).PixelsLA; PData.Region(4).PixelsLA]);
PData.Region(5).Shape.Name = 'Custom';
PData.Region(5).numPixels = length(PData.Region(5).PixelsLA);



%% Specify Media.  Use point targets in middle of PData.

% Media.MP(1,:) = [30,30,90,1.0];      % single point.
% Media.MP(2,:) = [0,10,60,1.0];
% Media.MP(3,:) = [-15,10,30,1.0];
% Media.MP(4,:) = [30,30,60,1.0];
% Media.MP(5,:) = [-30,-30,50,1.0];
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
Resource.RcvBuffer(1).rowsPerFrame = 4096*P.numSyntheticRcvs*2;
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = P.numFrames;
Resource.ImageBuffer(1).numFrames = P.numFrames;
Resource.InterBuffer(1).numFrames = 1;
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
Resource.DisplayWindow(1).numFrames = 20;
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
Resource.DisplayWindow(2).numFrames = 20;
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
Resource.DisplayWindow(3).numFrames = 20;
% Resource.DisplayWindow(3).mode = '2d';

% Specify Transmit waveform structure.
TW.type = 'parametric';
TW.Parameters = [Trans.frequency,.67,2,1];  % A, B, C, D

% Specify TPC structure.
TPC(1).hv = 30.0;

% Specify TX structure array.
TX = struct('waveform', 1, ...
    'Origin', [0.0,0.0,0.0], ...
    'focus', zApex, ...
    'Steer', [0.0,0.0], ...
    'Apod', zeros(1,1024), ...
    'Delay', zeros(1,1024),...
    'peakCutOff', 2,...
    'peakBLMax', 20,...
    'aperture',5,...
    'TXPD',[]);


TX.Apod = squeeze(CompRndApod(1,1,:))';
TX.aperture = Paper(1,1);
TX.Delay = computeTXDelays(TX);

%% allows to visualize the TX in 3D (similar to the EventAnalysisTool
% figure,for i=1:size(TX,2),plot3(Trans.ElementPos(find(TX(i).Apod==1),1),Trans.ElementPos(find(TX(i).Apod==1),2),TX(i).Delay(find(TX(i).Apod==1)),'k*'),zlim([0 10]),view(gca,-12,49),grid on,pause(.5),end

%%
%Specify TGC Waveform structure.
TGC(1).CntrlPts = [0 298 416 595 773 892 1012 1023];
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
    'callMediaFunc', 0), 1, P.numFrames*P.numSyntheticRcvs);


for j = 1:P.numFrames
    k=P.numSyntheticRcvs*(j-1);
    for w=1:P.numSyntheticRcvs
        Receive(w+k).acqNum = w;
        Receive(w+k).framenum = j;
        Receive(w+k).aperture = Paper(1,w);
        Receive(w+k).Apod = squeeze(CompRndApod(1,w,:))';
    end
end



%%
senscutoff = 0.6;
% Specify Recon structure arrays.

Recon = repmat(struct('senscutoff', senscutoff, ...
    'pdatanum', 1, ...
    'rcvBufFrame', -1, ...
    'IntBufDest', [1,1],...
    'ImgBufDest', [1,-1], ...
    'RINums', 1:P.numSyntheticRcvs), 1, 1);


% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 'accumIQ', ...
    'txnum', 1, ...
    'rcvnum', 1, ...
    'scaleFactor', 1, ...
    'regionnum', 1), 1, P.numSyntheticRcvs);

% - Set specific ReconInfo attributes.
ReconInfo(1).Pre = 'clearInterBuf';

for i = 1:P.numSyntheticRcvs
    ReconInfo(i).txnum = 1;
    ReconInfo(i).rcvnum = i;
    ReconInfo(i).regionnum = ReconRegion; %1 for the whole volume, 4 for the slices
end


ReconInfo(i).Post = 'IQ2IntensityImageBuf';

%%
pgain=25;
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
    'compressFactor',55,...
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
    'compressFactor',55,...
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
    'compressFactor',55,...
    'mappingMethod','full',...
    'display',1,...
    'displayWindow',3};

Process(4).classname = 'External';
Process(4).method = 'ModifVSXgui';
Process(4).Parameters = {'srcbuffer','none'};


% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; % jump back to start
SeqControl(1).argument = 2;
SeqControl(2).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(2).argument = 250;
SeqControl(3).command = 'timeToNextAcq';  % time between frames
SeqControl(3).argument = 50000-(250*17*2);%650000;  % 20 msec
% SeqControl(3).condition = 'ignore'; %Recon processing time
SeqControl(4).command = 'returnToMatlab';
nsc = 5; % nsc is count of SeqControl objects

n=1;
Event(n).info = 'Modif VSX gui';
Event(n).tx = 0;   % No TX structure.
Event(n).rcv = 0;  % NoRcv structure of frame.
Event(n).recon = 0;      % no reconstruction.
Event(n).process = 4;    % Processing
Event(n).seqControl = 0; %
n = n+1;


for j = 1:Resource.RcvBuffer(1).numFrames
    v=P.numSyntheticRcvs*(j-1);
    for i = 1:P.numSyntheticRcvs
        Event(n).info = 'Transmit/Receive.';
        Event(n).tx = 1;   % % Loop over the TX structure.
        Event(n).rcv = i+v;  % Rcv structure of frame.
        Event(n).recon = 0;      % no reconstruction.
        Event(n).process = 0;    % no processing
        Event(n).seqControl = 2; %
        n = n+1;
    end
    Event(n-1).seqControl = 3;
    
    Event(n).info = 'Transfer To Host';
    Event(n).tx = 0;         % No TX structure.
    Event(n).rcv = 0;    % No Rcv structure of frame.
    Event(n).recon = 0;      % no reconstruction.
    Event(n).process = 0;    % no processing
    Event(n).seqControl = nsc;%[nsc,nsc+1]; % set wait time and transfer data
    SeqControl(nsc).command = 'transferToHost';
    nsc = nsc + 1;
    n = n+1;
    
    
    Event(n).info = ['3D Reconstruct Frame ' num2str(j)];
    Event(n).tx = 0;    % No TX structure.
    Event(n).rcv = 0;   % No Rcv structure of frame.
    Event(n).recon = 1; % Reconstruction of the slices in real time or the whole volume in "3DRecon" = "on"
    Event(n).process = 0;   % no processing
    Event(n).seqControl = 0;% [nsc];
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
    if floor(j/3) == j/3     % Exit to Matlab every 3rd frame reconstructed
        Event(n).seqControl = 4;
    end
    
    n = n+1;
    
    Event(n).info = ['Process XY - Frame ' num2str(j)];
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 3;
    Event(n).seqControl = 0;
    n = n+1;
    
%     
%             Event(n).info = 'HW/SW Sync';
%     Event(n).tx = 0;         % No TX structure.
%     Event(n).rcv = 0;      % No Rcv structure.
%     Event(n).recon = 0;      % no reconstruction.
%     Event(n).process = 0;    % no processing
%     Event(n).seqControl = nsc;
%     SeqControl(nsc).command = 'sync';
%     nsc = nsc + 1;
%     n = n+1;


end

Event(n).info = 'Return to Matlab';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 4;

n=n+1;

Event(n).info = 'Jump';
Event(n).tx = 0;         % no TX structure.
Event(n).rcv = 0;        % no Rcv structure.
Event(n).recon = 0;      % no reconstruction.
Event(n).process = 0;    % call processing function
Event(n).seqControl = 1; %



%% User specified UI Control Elements
import vsv.seq.uicontrol.VsSliderControl; 


% - Sensitivity Cutoff
UI(1).Control = VsSliderControl('LocationCode','UserB7','Label','Sens. Cutoff',...
                  'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                  'SliderStep',[0.025,0.1],'ValueFormat','%1.3f', ...
                  'Callback', @SensCutoffCallback );  


UI(2).Control = VsSliderControl('LocationCode','UserB1','Label','Z-Section',...
                  'SliderMinMaxVal',[1,P.endDepth,Resource.DisplayWindow(3).ReferencePt(3)],...
                  'SliderStep',[1/(P.endDepth-1) 10/(P.endDepth-1)],'ValueFormat','%3.0f',...
                  'Callback', @ZSectionChangeCallback );  


UI(3).Control = VsSliderControl('LocationCode','UserB2','Label','Y-Section',...
                  'SliderMinMaxVal',[PData(1).Origin(2)-(PData(1).Size(1)-1)*PData(1).PDelta(1), PData(1).Origin(2), Resource.DisplayWindow(1).ReferencePt(2)],...
                  'SliderStep',[1/((PData(1).Size(2)-1)*PData(1).PDelta(2)-1) 10/((PData(1).Size(2)-1)*PData(1).PDelta(2)-1)],'ValueFormat','%3.0f',...
                  'Callback', @YSectionChangeCallback );  


UI(4).Control = VsSliderControl('LocationCode','UserB3','Label','X-Section',...
                  'SliderMinMaxVal',[PData(1).Origin(1), PData(1).Origin(1)+(PData(1).Size(2)-1)*PData(1).PDelta(2), Resource.DisplayWindow(2).ReferencePt(1)],...
                  'SliderStep',[1/((PData(1).Size(1)-1)*PData(1).PDelta(1)-1) 10/((PData(1).Size(1)-1)*PData(1).PDelta(1)-1)],'ValueFormat','%3.0f',...
                  'Callback', @XSectionChangeCallback );  

CompressionFactor=1.3;

UI(5).Control = VsSliderControl('LocationCode','UserA2','Label','CompressionP',...
                  'SliderMinMaxVal',[1,11,CompressionFactor],...
                  'SliderStep',[1/40,4/40],'ValueFormat','%1.1f',...
                  'Callback', @CompressionPCallback );  



UI(6).Control = vsv.seq.uicontrol.VsToggleButtonControl('LocationCode','UserB4','Label','3DRecon off',...
                  'Callback', @SwitchRecon3DCallback );  


EF(1).Function = vsv.seq.function.ExFunctionDef('ModifVSXgui', @ModifVSXgui);

EF(2).Function = vsv.seq.function.ExFunctionDef('CLSliderCallback', @CLSliderCallback);

EF(3).Function = vsv.seq.function.ExFunctionDef('CLValueCallback', @CLValueCallback);


% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 3;


% Save all the structures to a .mat file.
save('MatFiles/Matrix1024-3MHz-Flash_SRAC.mat');
return

%% **** Callback routines to be converted by text2cell function. ****

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

function YSectionChangeCallback(~,~,UIValue)
    %YSectionChange
    RS = evalin('base','Resource');
    RefPt = RS.DisplayWindow(1).ReferencePt;
    RefPt(2) = UIValue;
    RS.DisplayWindow(1).ReferencePt = RefPt;
    assignin('base','Resource',RS);

    PData = evalin('base','PData');
    P = evalin('base','P');


    PData(1).Region(3).Shape.oPAIntersect = RefPt(2);

    PData(1).Region = computeRegions(PData(1));

    PData.Region(5).PixelsLA = unique([PData.Region(2).PixelsLA;PData.Region(3).PixelsLA; PData.Region(4).PixelsLA]);
    PData.Region(5).Shape.Name = 'Custom';
    PData.Region(5).numPixels = length(PData.Region(5).PixelsLA);

    assignin('base','PData',PData);

    Control = repmat(struct('Command',[],'Parameters',[]),1,2);
    Control(1).Command = 'set&Run';
    Control(1).Parameters = {'DisplayWindow',1,'ReferencePt',RefPt};
    Control(2).Command = 'update&Run';
    Control(2).Parameters = {'PData'};
    assignin('base','Control', Control);
end

function XSectionChangeCallback(~,~,UIValue)
    %XSectionChange
    RS = evalin('base','Resource');
    RefPt = RS.DisplayWindow(2).ReferencePt;
    RefPt(1) = UIValue;
    RS.DisplayWindow(2).ReferencePt = RefPt;
    assignin('base','Resource',RS);

    PData = evalin('base','PData');
    P = evalin('base','P');

    PData(1).Region(2).Shape.oPAIntersect = RefPt(1);


    PData(1).Region = computeRegions(PData(1));
    PData.Region(5).PixelsLA = unique([PData.Region(2).PixelsLA;PData.Region(3).PixelsLA ;PData.Region(4).PixelsLA]);
    PData.Region(5).Shape.Name = 'Custom';
    PData.Region(5).numPixels = length(PData.Region(5).PixelsLA);

    assignin('base','PData',PData);

    Control = repmat(struct('Command',[],'Parameters',[]),1,2);
    Control(1).Command = 'set&Run';
    Control(1).Parameters = {'DisplayWindow',2,'ReferencePt',RefPt};
    Control(2).Command = 'update&Run';
    Control(2).Parameters = {'PData'};
    assignin('base','Control', Control);
end

function ZSectionChangeCallback(~,~,UIValue)
    %ZSectionChange
    RS = evalin('base','Resource');
    RefPt = RS.DisplayWindow(3).ReferencePt;
    RefPt(3) = UIValue;
    RS.DisplayWindow(3).ReferencePt = RefPt;
    assignin('base','Resource',RS);

    PData = evalin('base','PData');
    P = evalin('base','P');

    PData(1).Region(4).Shape.oPAIntersect = RefPt(3);

    PData(1).Region = computeRegions(PData(1));
    PData.Region(5).PixelsLA = unique([PData.Region(2).PixelsLA;PData.Region(3).PixelsLA ;PData.Region(4).PixelsLA]);
    PData.Region(5).Shape.Name = 'Custom';
    PData.Region(5).numPixels = length(PData.Region(5).PixelsLA);

    assignin('base','PData',PData);

    Control = repmat(struct('Command',[],'Parameters',[]),1,2);
    Control(1).Command = 'set&Run';
    Control(1).Parameters = {'DisplayWindow',3,'ReferencePt',RefPt};
    Control(2).Command = 'update&Run';
    Control(2).Parameters = {'PData'};
    assignin('base','Control', Control);
end

function CompressionPCallback(~,~,UIValue)
    %CompressionP
    CompressionFactor = UIValue;
    assignin('base','CompressionFactor',CompressionFactor);
end

function SwitchRecon3DCallback(~,~,UIState)
%SwitchRecon3D
    UI = evalin('base','UI');
    ReconInfo=evalin('base','ReconInfo');
    Recon=evalin('base','Recon');
    Resource=evalin('base','Resource');

    if UIState % 1 means pressed
        set(UI(6).handle,'String','3DRecon on');

        for i=1:size(ReconInfo,2)
            ReconInfo(i).regionnum = 1;
        end
        Resource.Parameters.simulateMode = 2;
    else
        set(UI(6).handle,'String','3DRecon off');
        handle3D=evalin('base','handle3D');
        for i=1:size(ReconInfo,2)
            ReconInfo(i).regionnum = 4;
        end
        close(handle3D)
        evalin('base','clear handle3D')
        Resource.Parameters.simulateMode = 0;

    end


    assignin('base','ReconInfo',ReconInfo);
    assignin('base','Resource',Resource);


    Control(1).Command = 'set&Run';
    Control(1).Parameters = {'Parameters',1,'simulateMode',Resource.Parameters.simulateMode,'startEvent',Resource.Parameters.startEvent};

    Control(2).Command = 'update&Run';
    Control(2).Parameters = {'Recon'};
    assignin('base','Control', Control);
end


%% **** Callback routines used by External function definition (EF) ****

function ModifVSXgui()
    h4=findobj('Tag','CLSlider');
    h4.Callback= @CLSliderCallback;

    h4=findobj('Tag','CLValue');
    h4.Callback= @CLValueCallback;


    assignin('base','h4',h4);

end

%-EF#4
function CLSliderCallback(varargin)

    Resource=evalin('base','Resource');
    h4=findobj('Tag','CLSlider');
    cv=h4.Value;

    h5=findobj('Tag','CLValue');
    h5.String=num2str(cv);


    hfreeze=findobj('Tag','FreezeButton');

    if hfreeze.Value==0  % no action if not in freeze
        return
    end

    for i=1:size(Resource.DisplayWindow,2)

        Control.Command = 'cineDisplay';
        Control.Parameters = cell(1,4);
        Control.Parameters{1} = 'displayWindow';
        Control.Parameters{2} = i;
        Control.Parameters{3} = 'frameNumber';
        Control.Parameters{4} = round(cv);
        runAcq(Control);
    end
end

function CLValueCallback(varargin)

    Resource=evalin('base','Resource');
    h4=findobj('Tag','CLSlider');
    h5=findobj('Tag','CLValue');
    cv=str2num(h5.String);

    h4.Value=cv;

    hfreeze=findobj('Tag','FreezeButton');

    if hfreeze.Value==0  % no action if not in freeze
        return
    end

    for i=1:size(Resource.DisplayWindow,2)

        Control.Command = 'cineDisplay';
        Control.Parameters = cell(1,4);
        Control.Parameters{1} = 'displayWindow';
        Control.Parameters{2} = i;
        Control.Parameters{3} = 'frameNumber';
        Control.Parameters{4} = round(cv);
        runAcq(Control);
    end

end
