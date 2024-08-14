% SetUpMatrix1024_3MHz_WideBeam_SRAC.m - Example of using wide beams for
% imaging in volume acquisitions using a matrix probe array (3 MHz) and
% Sparse Random Aperture Compounding (SRAC) technique.
%
% For more information about this technique please refer to the article
% published in Physics in Medicine and Biology in May 2020:
% Phys Med Biol. 2020 May 15. doi: 10.1088/1361-6560/ab9372
%
% Description:
%  This scripts transmits and receives 100 wide beams focused beyong the
%  imaging field which are then compounding in the volume recontruciton. 
%  The scanning patern in this script follows and grid of points at the
%  bottom of the imaging field. The beams are generated from a virtual
%  source at -100 wavelengths (behind the transducer). The parameters
%  P.numRaysX and P.numRaysY control the number of point in the grid in the
%  X and Y axis. The script evenly distributes the centroid of the beams at
%  the P.endDepth and calculates the corresponding steering angles for the 
%  transmits and the Regions that are used in the reconstruction.  
%  The user should keep in mind that a maximum of 255 different MUX
%  apertures can be used. 
%  The user can select the number of SRAC (P.numSyntheticRcvs), form 1 to 4,
%  that want to use for compounding. If this number is different than "1", 
%  then the script will transmit the same beams with differnt angles using 
%  the same Random aperture, however the receive of the backscatter data will
%  be done with different random apertures. The script loads two files, 
%  one with the Random apertures for the Transmit "Saved256RndApodOrigin", with 256 
%  distinct randome apertures. The second file is "saved100CompRndApod" which has 100
%  sets of 4 complementary random apertures to be used in the Receive. However, for 
%  this script only 1 set of the complementary random apertures is used.

%  The script allows the user to choose between reconstructing 3 orthogonal 
%  planes (XZ, YZ and XY) or the whole volume. The three orthogonal planes 
%  were chosen by default in the middle of the probe for the XZ and YZ planes, 
%  and at a depth of 60 wavelengths for the XY plane. 
%  To set the orthogonal plane reconstruction the user must set the variable
%  "ReconRegion" equal to a number different from 1. And to reconstruct the whole 
%  volume, the user must set ReconRegion=1. 

% Updates:
% 06/23/2020 - Update to SW 4.3 format for new UIControls and function definitions (VTS 1691). 
%   More info:(.../Example_Scripts/Vantage_Features/New UI Scheme/SetUpL11_5vFlash_NewUI)


clear 

P.numRaysX = 10;
P.numRaysY = 10;

NumRndApe=P.numRaysX*P.numRaysY; % Determines the number of Random Apertures needed to the transmit

P.startDepth = 0;   % start of Acquisition depth range in wavelengths
P.endDepth = 170; %   
P.viewAngle = 20 * pi/180;                    % angle between z-axis and surface of cone.


P.numFrames = 1;
P.numSyntheticRcvs = 2;

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
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.
Resource.Parameters.Connector = 1;
Resource.VDAS.useNewHal=1;

%% Specify Trans structure array.
Trans.name = 'Matrix1024-3';
Trans.units = 'mm';
Trans = computeTrans(Trans);

Trans.HVMux = computeUTAMux1024;
% add the HV Mux programming for the UTA 1024-MUX

load('Saved256RndApodOrigin'); %loading of the Random Apertures 

% Create Papod arrays to be pasted into the TX and Receive structures
Papod = RndApod(1:NumRndApe,:);

load('saved100CompRndApod'); 
CompRndApod=CompRndApod(1,:,:);

Papod=[Papod; squeeze(CompRndApod)];

% create aperture index for each Papod
for i = 1:size(Papod, 1)
    Paper(i) = computeMuxAperture(Papod(i, :), Trans);
end

waveLength = (Resource.Parameters.speedOfSound/1000)/Trans.frequency; %in Mm

%Intermediate Variables
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
zApex = -extent/tan(P.viewAngle); % in WL

PData(1).Size(1) = ceil(2.0*(P.endDepth-zApex)*tan(P.viewAngle)/PData(1).PDelta(2));
if mod(PData(1).Size(1),2)==0
    PData(1).Size(1) =  PData(1).Size(1)+1;
end

PData(1).Size(2) = PData(1).Size(1);
PData(1).Size(3) = ceil((P.endDepth)/PData(1).PDelta(3));
PData(1).Origin = [-((PData(1).Size(2)-1)/2)*PData(1).PDelta(1), ((PData(1).Size(1)-1)/2)*PData(1).PDelta(1), 0];

% Determining the location of the regions centroids at the P.endDepth
% according to the number of Rays defined in each direction
DistEndDepth=tan(P.viewAngle)*(P.endDepth-zApex)*waveLength;     % in mm
tempTxBmCX2 = -DistEndDepth:(2*DistEndDepth/(P.numRaysX+1)):DistEndDepth;        % x coords of beam centers in mm
TxBmCX2 = tempTxBmCX2(2:end-1);                          %This is to eleminate the edges
tempTxBmCY2 = -DistEndDepth:(2*DistEndDepth/(P.numRaysY+1)):DistEndDepth;        % y coords of beam centers in mm
TxBmCY2 = tempTxBmCY2(2:end-1);                          %This is to eleminate the edges
[BmCX2,BmCY2]=ndgrid(TxBmCX2,TxBmCY2);
bmcx=reshape(BmCX2,1,size(BmCX2,1)*size(BmCX2,2));       %in Mm
bmcy=reshape(BmCY2,1,size(BmCY2,1)*size(BmCY2,2));       %in Mm

AZ=atan(bmcx/((P.endDepth)*waveLength));% in rad from a Mm calculation
EL=atan(bmcy/((P.endDepth)*waveLength));

Rx=zApex*waveLength.*tan(AZ); % in Mm
Ry=zApex*waveLength.*tan(EL); % in Mm

xorig=tan(AZ)*zApex*waveLength; % intersect in the z=0 plane from the zApex to the bmcx and bmcy targets at z=P.endDepth
yorig=tan(EL)*zApex*waveLength;

P.viewAngle=0.26;
PData(1).Region = repmat(struct('Shape',struct('Name','Pyramid',...
    'Position',[0,0,0],...
    'angle',P.viewAngle/P.numRaysX*4,...% The number that multiplies P.viewAngle/P.numRaysX is the number of pixels that would overlap for region (1, would mean no overlap)
    'top',10,...
    'z1',P.startDepth,...
    'z2',P.endDepth,...
    'Steer',[0,0])), 1,NumRndApe);

for n = 1:NumRndApe
    PData(1).Region(n).Shape.Steer = [AZ(n) EL(n)];
end

PData(1).Region(NumRndApe+1) = struct('Shape',struct('Name','Pyramid','Position',[0,0,zApex],'angle',P.viewAngle,'z1',P.startDepth,'z2',P.endDepth,'Steer',[0,0]));

PData(1).Region(NumRndApe+2) = struct('Shape',struct('Name','Slice','Orientation','yz','oPAIntersect',PData.Origin(1)+PData.PDelta(2)*(PData.Size(2)-1)/2,'andWithPrev',NumRndApe+1));
PData(1).Region(NumRndApe+3) = struct('Shape',struct('Name','Slice','Orientation','xz','oPAIntersect',PData.Origin(2)-PData.PDelta(1)*(PData.Size(1)-1)/2,'andWithPrev',NumRndApe+1));
PData(1).Region(NumRndApe+4) = struct('Shape',struct('Name','Slice','Orientation','xy','oPAIntersect',60,'andWithPrev',NumRndApe+1));

[PData(1).Region] = computeRegions(PData(1));

PData.Region(NumRndApe+5).PixelsLA = unique([PData.Region(NumRndApe+2).PixelsLA; PData.Region(NumRndApe+3).PixelsLA; PData.Region(NumRndApe+4).PixelsLA]);
PData.Region(NumRndApe+5).Shape.Name = 'Custom';
PData.Region(NumRndApe+5).numPixels = length(PData.Region(NumRndApe+5).PixelsLA);


%% Specify Media.  Use point targets in middle of PData.


% Media.MP(1,:) = [30,0,90,1.0];      % single point.
% Media.MP(2,:) = [20,0,50,1.0];
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
Resource.RcvBuffer(1).rowsPerFrame = 4096*NumRndApe*P.numSyntheticRcvs;
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
% Resource.DisplayWindow(1).mode = '2d';
Resource.DisplayWindow(1).numFrames = 20;

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
Resource.DisplayWindow(2).numFrames = 20;

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
Resource.DisplayWindow(3).numFrames = 20;

% Specify Transmit waveform structure.
TW.type = 'parametric';
TW.Parameters = [Trans.frequency,.67,2,1];  % A, B, C, D

% Specify TPC structure.
TPC(1).hv = 30.0;

AZtx=atan(bmcx/((P.endDepth)*waveLength));%rad from a Mm calculation
ELtx=atan(bmcy/((P.endDepth)*waveLength));


% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
    'Origin', [0.0,0.0,0.0], ...
    'focus', -100,...
    'Steer', [0.0,0.0], ...
    'Apod', zeros(1,Trans.numelements), ...
    'Delay', zeros(1,Trans.numelements),...
    'peakCutOff', 1,...
    'peakBLMax', 10,...
    'aperture',1), 1,NumRndApe);


for i = 1:NumRndApe
    TX(i).Steer = [AZ(i) EL(i) 0];
    TX(i).Apod = Papod(i, :);
    TX(i).aperture = Paper(i); %Since the synthetic receive is 4, then the rnd apertures start at Paper(5)
    TX(i).Delay = computeTXDelays(TX(i));
    %         TX(i).TXPD = computeTXPD(TX(i),PData(1));
end


%% allows to visualize the TX in 3D (similar to the EventAnalysisTool
% figure,for i=1:size(TX,2),plot3(Trans.ElementPos(find(TX(i).Apod==1),1),Trans.ElementPos(find(TX(i).Apod==1),2),TX(i).Delay(find(TX(i).Apod==1)),'k*'),zlim([0 10]),view(gca,-12,49),grid on,pause(.5),end

%% To plot the Regions and the Tx steering
bmcx2=tan(AZ)*((P.endDepth-zApex)*waveLength);
bmcy2=tan(EL)*((P.endDepth-zApex)*waveLength);

bmcx2tx=tan(AZtx)*((P.endDepth)*waveLength);
bmcy2tx=tan(ELtx)*((P.endDepth)*waveLength);

%% This figure allows the visualization of the centroids and the center line for the beams
% figure,plot3(bmcx,bmcy,ones(1,length(bmcx))*P.endDepth*waveLength,'r*'),grid on
% hold on,plot(Trans.ElementPosMm(:,1),Trans.ElementPosMm(:,2),'*k')
% hold on,plot3(0,0,zApex*waveLength,'b*'),set(gca,'Zdir','reverse')
% for i=1:NumRndApe,hold on,line([Rx(i), bmcx(i)],[Ry(i), bmcy(i)],[zApex*waveLength,P.endDepth*waveLength]),hold off,end
% % % for i=1:NumRndApe,hold on,line([0, bmcx2(i)],[0, bmcy2(i)],[zApex*waveLength,P.endDepth*waveLength],'color','r'),hold off,end

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
    'aperture',5, ...
    'sampleMode', 'NS200BW', ...
    'callMediaFunc', 0), 1, NumRndApe*P.numFrames*P.numSyntheticRcvs);%P.numTx*P.numFrames*P.numSyntheticTx*P.numSyntheticRcvs);


for j = 1:P.numFrames
    k=P.numSyntheticRcvs*NumRndApe*(j-1);
    for i = 1:NumRndApe
        kk=P.numSyntheticRcvs*(i-1);
        for w=1:P.numSyntheticRcvs
            Receive(w+kk+k).acqNum = w+kk;
            Receive(w+kk+k).framenum = j;
            Receive(w+kk+k).aperture = Paper(NumRndApe+w);
            Receive(w+kk+k).Apod = Papod(NumRndApe+w,:);
        end
    end
end


%%

Recon = repmat(struct('senscutoff', 0.5, ...
    'pdatanum', 1, ...
    'IntBufDest', [1,1], ...
    'ImgBufDest', [1,-1], ...
    'rcvBufFrame',-1, ...
    'RINums',1:NumRndApe*P.numSyntheticRcvs), 1, 1);

w=NumRndApe*P.numSyntheticRcvs;


ReconInfo = repmat(struct('mode', 'accumIQ', ...  % default is to accumulate IQ data.
    'txnum', 1, ...
    'rcvnum', 1, ...
    'scaleFactor', 1.0, ...
    'regionnum', 1), 1, NumRndApe*P.numSyntheticRcvs);
% - Set specific ReconInfo attributes.


ReconInfo(1).Pre = 'clearInterBuf';

for w=1:NumRndApe
    k= P.numSyntheticRcvs*(w-1);
    for i = 1:P.numSyntheticRcvs
        ReconInfo(i+k).txnum = w;
        ReconInfo(i+k).rcvnum = i+k;
        if ReconRegion ==1 %ReconRegion = 1 reconstruct the full volume (every region, w), ReconRegion = 5, reconstruct the 3 orthogonal planes
            ReconInfo(i+k).regionnum =  w;
        else
           ReconInfo(i+k).regionnum =  NumRndApe+5;
        end
        
    end
    
end

ReconInfo(i+k).Post = 'IQ2IntensityImageBuf';


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


%
% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; % jump back to start
SeqControl(1).argument = 2;
SeqControl(2).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(2).argument = 250;
SeqControl(3).command = 'timeToNextAcq';  % time between frames
SeqControl(3).argument = 50000-(250*100*1);  % 20 msec
% SeqControl(3).condition = 'ignore'; %Recon processing time
SeqControl(4).command = 'returnToMatlab';
nsc = 5; % nsc is count of SeqControl objects

n=1;
Event(n).info = 'Modif VSX gui';
Event(n).tx = 0;   % No TX structure.
Event(n).rcv = 0;  % No Rcv structure of frame.
Event(n).recon = 0;      % no reconstruction.
Event(n).process = 4;    % no processing
Event(n).seqControl = 0; %
n = n+1;

for j = 1:Resource.RcvBuffer(1).numFrames
    v=P.numSyntheticRcvs*NumRndApe*(j-1);
    for w=1:NumRndApe
        k= P.numSyntheticRcvs*(w-1);
        for i = 1:P.numSyntheticRcvs
            Event(n).info = 'Transmit/Receive.';
            Event(n).tx = w;   % Loop over all the Beams TX structure.
            Event(n).rcv = i+k+v;  % Rcv structure of frame.
            Event(n).recon = 0;      % no reconstruction.
            Event(n).process = 0;    % no processing
            Event(n).seqControl = 2; %
            n = n+1;
        end
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
    if floor(j/3) == j/3     % Exit to Matlab every 3rd frame reconstructed
        Event(n).seqControl = 4;
    end
    
    n = n+1;
    
    Event(n).info = ['Process XY - Frame ' num2str(i)];
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 3;
    Event(n).seqControl = 0;
    n = n+1;
    
if ReconRegion==1
    Event(n).info = 'HW/SW Sync';
    Event(n).tx = 0;         % use 1st TX structure.
    Event(n).rcv = 0;      % use ith Rcv structure.
    Event(n).recon = 0;      % no reconstruction.
    Event(n).process = 0;    % no processing
    Event(n).seqControl = [4 nsc];
    SeqControl(nsc).command = 'sync';
    nsc = nsc + 1;
    n = n+1;
end
    
    
end

Event(n).info = 'Return to Matlab';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 4;

n = n+1;

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
save('MatFiles/Matrix1024-3MHz-WideBeam_SRAC.mat');
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
    ReconRegion=evalin('base','ReconRegion');
    NumRndApe=evalin('base','NumRndApe');
    P=evalin('base','P');

    if UIState % 1 means pressed
        set(UI(6).handle,'String','3DRecon on');
        ReconRegion=1;
        Resource.Parameters.simulateMode = 2;
    else
        set(UI(6).handle,'String','3DRecon off');
        handle3D=evalin('base','handle3D');
        ReconRegion=5;
        close(handle3D)
        evalin('base','clear handle3D')
        Resource.Parameters.simulateMode = 0;

    end
    
    for w=1:NumRndApe
        k= P.numSyntheticRcvs*(w-1);
        for i = 1:P.numSyntheticRcvs
            if ReconRegion ==1 %ReconRegion = 1 reconstruct the full volume (every region, w), ReconRegion = 5, reconstruct the 3 orthogonal planes
                ReconInfo(i+k).regionnum =  w;
            else
                ReconInfo(i+k).regionnum =  NumRndApe+5;
            end
            
        end
        
    end
    

    assignin('base','ReconInfo',ReconInfo);
    assignin('base','Resource',Resource);
    assingin('base','ReconRegion',ReconRegion);

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


