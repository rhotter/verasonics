% SetupMatrix2048FlashAngles2DRndApe_2048Mux  This is a replica of
% the existing example script "SetupMatrix1024_3MHz_FlashAngles2DRndApe",
% but has been modified for use with the Vantage 512 system with SHI
% 2048-MUX, and a purely hypothetical 2048 element 3 MHz Matrix array
% having a rectangular aperture of 64 rows by 32 columns of elements.
%
% This basic functionality of this script is the same as the original, with
% the same conical volume for image data acquisition, over a range of
% acquisitions each with transmit steered at several different angles from
% the Z-axis and in a range of directions around the circular cross-section
% of the conical volume.  For each acquisition a randomly generated sparse
% aperture of 512 elements is used, with four synthetic receive aperture
% acquisition events each using 512 of the 2048 elements.
%
% Search for the string **V-512** for comments explaining the changes from
% the original V-256 script, to make it compatible with the V-512 sysem
% with SHI 2048-MUX.  Comments with the string **Matrix2048** cover changes
% to the Trans structure to redefine it as the 64 X 32 2048 element array.

% **V-512** The following comments are from the original V-256 script:
%
%  This script does Flash multiAngles using Randome Apertures using 256 elements
%  of the full apertures for each transmit. The receive is done with four
%  subapertues (1-256, 257-512, 513-768, 769-1024). Therefore every transmit
%  needs to be repeated four times so that it can be received with each subapertures.
%  A set of 100 Random Apertures taking into account the constrains of the MUX
%  switches have been created for this script which are loaded from a .mat
%  file "Saved100RndApod". However, the user is free to create new sets
%  of Random Apertures using the function "Rnd1024Aperture" for as many apertures
%  as needed.
%  For this script the user defines a set of Angles defined by the using a
%  vector name TxAngles. The user needs to specify the + and - angles
%  which should be symmetric. For each of this angles then a series of
%  transmits around 2*pi are calculated using different Random Apertures.
%  The number of transmits around the z axis is defined by na2. Where a na2 = 9
%  is a transmit every 45 degrees. If the  number of angles in TxAngles is odd,
%  the script will assume that the middle angle is zero for which it will create
% a diverging wave with no steering.

%  The user should be careful when choosing the TimeToNextAcq to determine
%  the frame rate and the number of frames acquired. Since having a high
%  frame rate with a few number of frames could lead to overwritting the
%  buffer while the reconstruction is taking place for one frame. In the
%  case of this scrip the average time for the reconstruction time for the
%  slices (Region = 4) was estimated around 600 ms given the host computer
%  that was used at the time of the testing. Taking into accound that the
%  TTNA1 is 0.25 ms, which is used for the 17 angles and the 4
%  synthetic receives, the Tx/Rcv for 1 frame takes 17 ms. Then the TTNA2
%  (pause between frames) is set to 20 ms (this allows the time to do the
%  transfer to host). The total period for 1 acquisition is 37 ms.
%  Therefore the number of frames acquired x Period needs to be larger than
%  the time it takes to do 1 reconstruction (600 ms). In this particular
%  case for this scrip the number of angles should be larger than 18
%  frames. If the script is only going to be run one time (no jump at the
%  end of it), the number of frames can be as low as 1 since there is no
%  risk of overwritting the buffer.
%  In order for the customer to determine the reconstruction time two
%  external functions (before and after the recon) have been added to time
%  the reconstruction (vstic and vstoc). The use should keep in mind that
%  the first %  reconstruction usually takes a significant amount of time
%  due to the initialization of the system.

% **V-512** Revision history
% Copyright 2020 Verasonics, Inc. Verasonics Registered U.S. Patent and Trademark Office.
% Created April 19, 2020 for internal V-512 development

clear all;clear functions;

P.startDepth = 0;   % start of Acquisition depth range in wavelengths
P.endDepth = 128; %   use longer lines for noise floor testing
P.viewAngle = 20 * pi/180;                    % angle between z-axis and surface of cone.


TxAngles =[-20 -10 0 10 20]; %Angles defined by user in degrees
na = length(TxAngles);      % Set na = number of angles.
na2= 9; %this parameter determines the number transmits around around 2*pi


P.numFrames = 18;
% **Matrix2048** four synthetic receive acquisitions of 512 elements each,
% to acquire data from all 2048 elements
P.numSyntheticRcvs = 4;

RcvProfile.AntiAliasCutoff = 15; % 10, 15, 20, 30
% RcvProfile.LnaGain = 18;
% RcvProfile.LnaZinSel = 31;

%% Specify system parameters.
% **V-512** Set numTransmit to 512, and note the field numRcvChannels is no longer used at all.
Resource.Parameters.numTransmit = 512;      % number of transmit channels.
% Resource.Parameters.numRcvChannels = 256;   % number of receive channels.
Resource.Parameters.verbose = 2; % 1 enables and 0 disables
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.speedCorrectionFactor = 1.0;
Resource.Parameters.fakeScanhead = 0; % allow system HW operation with nothing connected
Resource.Parameters.simulateMode = 0;
Resource.Parameters.Connector = 1;
Resource.VDAS.dmaTimeout=1000;


%% Specify Trans structure array.
% **Matrix2048** First call computeTrans to create Trans structure for the
% Vermon 3 MHz 1024 element 2D array.  Then we will redefine Trans fields
% as needed to define our hypothetical 2048 element 64 X 32 2D array; the
% transducer center frequency, bandwidth, element pitch, and other
% characteristics will all be inherited from the 'Matrix1024-3'.
Trans.name = 'Matrix1024-3';
Trans.units = 'mm';
Trans = computeTrans(Trans);
Trans.name = 'Matrix2048-3MHz-64x32'; % **Matrix2048** Probe name for our hypothetical array
% **Matrix2048** Note Trans.id for direct-connect probes is derived from
% the UTA or SHI to which they are connected; this has not been defined yet
% for SHI 2048-MUX so we will leave it unchanged for now
Trans.connType = 13; % **V-512** Redefine connType for probe connected to SHI 2048-MUX
Trans.numelements = 2048; % **Matrix2048** element count for our hypothetical array
Trans.ConnectorES = (1:2048)'; % **Matrix2048** assume one-to-one wiring for our hypothetical array
Trans.ElementPos = zeros(2048, 5); % **Matrix2048** Z-coordinate and element angular orientation are all zero
Trans.ElementPos(:, 1) = repmat((-31.5:31.5)' * Trans.spacingMm, 32, 1); % **Matrix2048** X-coordinate for each element in mm
for rownum = 1:32
    % **Matrix2048** for each row of 64 elements, assign the Y-coordinate for that row in mm
    rowOffset = (rownum-1)*64; % element number at the end of the previous row
    Trans.ElementPos(rowOffset + (1:64), 2) = (rownum-16.5) * Trans.spacingMm; % Y-coordinae for all elements in that row
end
Trans.HVMux = computeSHIMux2048; % **V-512** add the HV Mux programming for the SHI 2048-MUX

% **V-512** Instead of loading the precomputed random apertures used by
% the V-256 script we will create our own 512 element random apertures when
% creating the TX structures.
% load('Saved100RndApod'); % this files contains 100 Rnd Apertures (100x1024) that were calculated using an external function taking into account the limits set by the MUX switches

% Create Papod arrays to be pasted into the TX and Receive structures, with
% the fifth entry enabling all 1024 elements for a single flash transmit
% **V-512** Redefine all four apertures as 512 of the 2048 for receive
Papod =  [   ones(1,512) zeros(1,512) zeros(1,512) zeros(1,512); ...
            zeros(1,512)  ones(1,512) zeros(1,512) zeros(1,512); ...
            zeros(1,512) zeros(1,512)  ones(1,512) zeros(1,512); ...
            zeros(1,512) zeros(1,512) zeros(1,512)  ones(1,512)]; % **V-512**
%     RndApod];

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
% **Matrix2048** Limit extent to the 32 element dimension of the array
if ~exist('extent','var'), extent = max(Trans.ElementPosWL(:,2)); end
zApex = -extent/tan(P.viewAngle);

PData(1).Size(1) = ceil(2.0*(P.endDepth-zApex)*tan(P.viewAngle)/PData(1).PDelta(2));  if mod(PData(1).Size(1),2)==0, PData(1).Size(1) =  PData(1).Size(1)+1; end
PData(1).Size(2) = PData(1).Size(1);
PData(1).Size(3) = ceil((P.endDepth)/PData(1).PDelta(3));
PData(1).Origin = [-((PData(1).Size(2)-1)/2)*PData(1).PDelta(1), ((PData(1).Size(1)-1)/2)*PData(1).PDelta(1), 0];

PData(1).Region(1) = struct('Shape',struct('Name','Pyramid','Position',[0,0,zApex],'angle',P.viewAngle,'z1',P.startDepth,'z2',P.endDepth));

PData(1).Region(2) = struct('Shape',struct('Name','Slice','Orientation','yz','oPAIntersect',PData.Origin(1)+PData.PDelta(2)*(PData.Size(2)-1)/2,'andWithPrev',1));
PData(1).Region(3) = struct('Shape',struct('Name','Slice','Orientation','xz','oPAIntersect',PData.Origin(2)-PData.PDelta(1)*(PData.Size(1)-1)/2,'andWithPrev',1));
% PData(1).Region(4) = struct('Shape',struct('Name','Slice','Orientation','xy','oPAIntersect',90,'andWithPrev',1));

[PData(1).Region] = computeRegions(PData(1));

PData.Region(4).PixelsLA = unique([PData.Region(2).PixelsLA;PData.Region(3).PixelsLA]); %;PData.Region(+4).PixelsLA]);
PData.Region(4).Shape.Name = 'Custom';
PData.Region(4).numPixels = length(PData.Region(4).PixelsLA);



%% Specify Media.  Use point targets in middle of PData.
% for i=1:61
%     Media.MP(i,:)= [0,-30+i,50,1];
%     Media.MP(61+i,:)= [-30+i,0,72,1];
% end
% Media.MP(size(Media.MP,1)+1,:)=[18,32,82,1];


Media.MP(1,:) = [30,30,90,1.0];      % single point.
Media.MP(2,:) = [0,10,60,1.0];
Media.MP(3,:) = [-15,10,30,1.0];
Media.MP(4,:) = [30,30,60,1.0];
Media.MP(5,:) = [-30,-30,50,1.0];
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
Resource.RcvBuffer(1).rowsPerFrame = 4096*na*na2*P.numSyntheticRcvs;
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numTransmit;
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

%Calculating the AZ and EL steering angles for the TxAngles and the na2
%angles
rotang=linspace(0,2*pi,na2); % na2=9 creates a plane every 45 degrees (0, 45, 90, 135, 180, 225, 270, 315 and 360).
rotang=rotang(1:end-1); % this eliminates the last angle since 0 = 360

AZ=[];
EL=[];
txAngles=TxAngles.*pi./180;

for i=1:fix(na/2)
    AZ=[AZ atan(tan(txAngles(i)).*cos(rotang))];
    EL=[EL atan(tan(txAngles(i)).*sin(rotang))];
end

if fix(na/2) ~= na/2       % if na is odd, it assumes the script calls for a wave with no steering (plane wave or diverging wave)
    AZ=[AZ 0];
    EL=[EL 0];

end

nang=length(AZ);% number of angles transmitted combining TxAngles and na2 angles (and 0 if the TxAngles was odd)


TX = repmat(struct('waveform', 1, ...
    'Origin', [0.0,0.0,0.0], ...
    'focus', zApex, ...
    'Steer', [0.0,0.0], ...
    'Apod', zeros(1,1024), ...
    'Delay', zeros(1,1024),...
    'peakCutOff', 2,...
    'peakBLMax', 20,...
    'aperture',5,...
    'TXPD',[]), 1,nang);


for n=1:nang
    TX(n).Steer = [AZ(n),EL(n) ];
    % **V-512** **Matrix2048** create a new random aperture for each TX,using the
    % "Rnd512of2048Aperture" subfunction added at the end of this script
    TX(n).Apod = Rnd512of2048Aperture();
    % **V-512** call the "computeMuxAperture" function to create the HVMux
    % aperture for this Apod, and assign an aperture index to it
    TX(n).aperture = computeMuxAperture(TX(n).Apod, Trans);
    TX(n).Delay = computeTXDelays(TX(n));
end

%% allows to visualize the TX in 3D (similar to the EventAnalysisTool
% figure,for i=1:size(TX,2),plot3(Trans.ElementPos(find(TX(i).Apod==1),1),Trans.ElementPos(find(TX(i).Apod==1),2),TX(i).Delay(find(TX(i).Apod==1)),'k*'),zlim([0 10]),view(gca,-12,49),grid on,pause(.5),end

%%
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
    'callMediaFunc', 0), 1, P.numFrames*P.numSyntheticRcvs*nang);


for j = 1:P.numFrames
    tt=1;
    k=P.numSyntheticRcvs*nang*(j-1);
    for i = 1:P.numSyntheticRcvs*nang
        Receive(i+k).acqNum = i;
        Receive(i+k).framenum = j;
        Receive(i+k).aperture = tt;
        Receive(i+k).Apod = Papod(tt,:);
        if mod(i,P.numSyntheticRcvs)==0
            tt=1;
        else
            tt=tt+1;

        end
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
    'RINums', 1:nang*P.numSyntheticRcvs), 1, 1);


% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 'accumIQ', ...
    'txnum', 1, ...
    'rcvnum', 1, ...
    'scaleFactor', 1, ...
    'regionnum', 1), 1, nang*P.numSyntheticRcvs);

% - Set specific ReconInfo attributes.
% ReconInfo(1).Pre = 'clearImageBuf';
ReconInfo(1).Pre = 'clearInterBuf';

for w=1:nang
    k= P.numSyntheticRcvs*(w-1);
    for i = 1:P.numSyntheticRcvs
        ReconInfo(i+k).txnum = w;
        ReconInfo(i+k).rcvnum = i+k;
        ReconInfo(i+k).regionnum = 4; %1 for the whole volume, 4 for the slices
    end

end

ReconInfo(i+k).Post = 'IQ2IntensityImageBuf';

%%
pgain=5; % **V-512** changing from 20 to 5 since TX apertures now have 512 active elements
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
Process(4).method = 'volumetricPlot';
Process(4).Parameters = {'srcbuffer','image',...
    'srcbufnum',1,...
    'srcframenum',-1,...
    'dstbuffer','none'};

Process(5).classname = 'External';
Process(5).method = 'SlicePlot';
Process(5).Parameters = {'srcbuffer','image',...
    'srcbufnum',1,...
    'srcframenum',-1,...
    'dstbuffer','none'};

Process(6).classname = 'External';
Process(6).method = 'ModifVSXgui';
Process(6).Parameters = {'srcbuffer','none'};

Process(7).classname = 'External';
Process(7).method = 'vstic';
Process(7).Parameters = {'srcbuffer','none'};

Process(8).classname = 'External';
Process(8).method = 'vstoc';
Process(8).Parameters = {'srcbuffer','none'};


% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; % jump back to start
SeqControl(1).argument = 2;
SeqControl(2).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(2).argument = 250;
SeqControl(3).command = 'timeToNextAcq';  % time between frames
SeqControl(3).argument = 20000;%650000;  % 20 msec
% SeqControl(3).condition = 'ignore'; %Recon processing time
SeqControl(4).command = 'returnToMatlab';
nsc = 5; % nsc is count of SeqControl objects

n=1;
Event(n).info = 'Modif VSX gui';
Event(n).tx = 0;   % No TX structure.
Event(n).rcv = 0;  % NoRcv structure of frame.
Event(n).recon = 0;      % no reconstruction.
Event(n).process = 6;    % Processing
Event(n).seqControl = 0; %
n = n+1;


for j = 1:Resource.RcvBuffer(1).numFrames

    v=P.numSyntheticRcvs*nang*(j-1);
    for w=1:nang
        k= P.numSyntheticRcvs*(w-1);
        for i = 1:P.numSyntheticRcvs
            Event(n).info = 'Transmit/Receive.';
            Event(n).tx = w;   % Loop over the TX structure.
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
    %     SeqControl(nsc).command = 'waitForTransferComplete';
    %     SeqControl(nsc).argument = nsc-1;
    %     nsc = nsc + 1;
    n = n+1;

    if j==1 || j==2
        Event(n).info = 'tic ';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 7;
        Event(n).seqControl = 0;
        n = n+1;
    end

    Event(n).info = ['3D Reconstruct Frame ' num2str(j)];
    Event(n).tx = 0;    % No TX structure.
    Event(n).rcv = 0;   % No Rcv structure of frame.
    Event(n).recon = [1]; % Reconstruction of the slices in real time or the whole volume in "3DRecon" = "on"
    Event(n).process = 0;   % no processing
    Event(n).seqControl = 0;% [nsc];
    %           SeqControl(nsc).command = 'markTransferProcessed';  %
    %       SeqControl(nsc).argument = SeqControl(nsc-1).argument;
    %       nsc = nsc + 1;
    n = n+1;

    if j==1 || j==2
        Event(n).info = 'toc ';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 8;
        Event(n).seqControl = 0;
        n = n+1;
    end

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
    Event(n).process = 3;
    Event(n).seqControl = 0;
    n = n+1;


    Event(n).info = ['Slice Process  - Frame' num2str(j)];
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 5;
    Event(n).seqControl = 0;
    if floor(j/3) == j/3     % Exit to Matlab every 3rd frame reconstructed
        Event(n).seqControl = 4;
    end

    n = n+1;



    %             Event(n).info = 'HW/SW Sync';
    %         Event(n).tx = 0;         % no TX structure.
    %         Event(n).rcv = 0;        % no Rcv structure.
    %         Event(n).recon = 0;      % no reconstruction.
    %         Event(n).process = 0;    % no processing
    %         Event(n).seqControl = nsc;
    %         SeqControl(nsc).command = 'sync';
    %         nsc = nsc + 1;
    %         n = n+1;
    %

end
%     Event(n).info = 'Return to Matlab';
%     Event(n).tx = 0;
%     Event(n).rcv = 0;
%     Event(n).recon = 0;
%     Event(n).process = 0;
%     Event(n).seqControl = 4;
%     n = n+1;

Event(n).info = 'Jump';
Event(n).tx = 0;         % no TX structure.
Event(n).rcv = 0;        % no Rcv structure.
Event(n).recon = 0;      % no reconstruction.
Event(n).process = 0;    % call processing function
Event(n).seqControl = 1; %


% User specified UI Control Elements
% - Sensitivity Cutoff
UI(1).Control =  {'UserB7','Style','VsSlider','Label','Sens. Cutoff',...
    'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
    'SliderStep',[0.025,0.1],'ValueFormat','%1.3f'};
UI(1).Callback = text2cell('%SensCutoff');
% - Section number change
UI(2).Control = {'UserB1','Style','VsSlider','Label','Z-Section',...
    'SliderMinMaxVal',[1,P.endDepth,Resource.DisplayWindow(3).ReferencePt(3)],...
    'SliderStep',[1/(P.endDepth-1) 10/(P.endDepth-1)],'ValueFormat','%3.0f'};
UI(2).Callback = text2cell('%ZSectionChange');
% - Section number change
UI(3).Control = {'UserB2','Style','VsSlider','Label','Y-Section',...
    'SliderMinMaxVal',[PData(1).Origin(2)-(PData(1).Size(1)-1)*PData(1).PDelta(1), PData(1).Origin(2), Resource.DisplayWindow(1).ReferencePt(2)],...
    'SliderStep',[1/((PData(1).Size(2)-1)*PData(1).PDelta(2)-1) 10/((PData(1).Size(2)-1)*PData(1).PDelta(2)-1)],'ValueFormat','%3.0f'};
UI(3).Callback = text2cell('%YSectionChange');
% - Section number change
UI(4).Control = {'UserB3','Style','VsSlider','Label','X-Section',...
    'SliderMinMaxVal',[PData(1).Origin(1), PData(1).Origin(1)+(PData(1).Size(2)-1)*PData(1).PDelta(2), Resource.DisplayWindow(2).ReferencePt(1)],...
    'SliderStep',[1/((PData(1).Size(1)-1)*PData(1).PDelta(1)-1) 10/((PData(1).Size(1)-1)*PData(1).PDelta(1)-1)],'ValueFormat','%3.0f'};
UI(4).Callback = text2cell('%XSectionChange');
CompressionFactor=1.3;
% - compression
UI(5).Control =  {'UserA2','Style','VsSlider','Label','CompressionP',...
    'SliderMinMaxVal',[1,11,CompressionFactor],...
    'SliderStep',[1/40,4/40],'ValueFormat','%1.1f'};
UI(5).Callback = text2cell('%CompressionP');


% Recon3D switch button to change to the reconstruction of the whole volume of the previously
% acquired frames
UI(6).Control =  {'UserB4','Style','VsToggleButton','Label','3DRecon off'};
UI(6).Callback = text2cell('%SwitchRecon3D');

% EF(1).Function = text2cell('%-EF#1'); %Function for the Volumetric Plot
EF(1).Function = vsv.seq.function.ExFunctionDef('volumetricPlot', @volumetricPlot);

% EF(2).Function = text2cell('%-EF#2'); %Function for the slice plot)
EF(2).Function = vsv.seq.function.ExFunctionDef('SlicePlot', @SlicePlot);

EF(3).Function = text2cell('%-EF#3'); %Function to modify the VSX-GUI

EF(4).Function = text2cell('%-EF#4'); %Function to modify CLSlider

EF(5).Function = text2cell('%-EF#5'); %Function to modify CLValue

EF(6).Function = text2cell('%-EF#6'); %Function to modify CLValue

EF(7).Function = text2cell('%-EF#7'); %Function to modify CLValue


% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 3;


% Save all the structures to a .mat file.
save('MatFiles/Matrix2048FlashAngles2DRndApe_2048Mux.mat');
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

PData.Region(4).PixelsLA = unique([PData.Region(2).PixelsLA;PData.Region(3).PixelsLA]); %;PData.Region(+4).PixelsLA]);
PData.Region(4).Shape.Name = 'Custom';
PData.Region(4).numPixels = length(PData.Region(4).PixelsLA);

assignin('base','PData',PData);

Control = repmat(struct('Command',[],'Parameters',[]),1,2);
Control(1).Command = 'set&Run';
Control(1).Parameters = {'DisplayWindow',1,'ReferencePt',RefPt};
Control(2).Command = 'update&Run';
Control(2).Parameters = {'PData'};
assignin('base','Control', Control);
% keyboard
%YSectionChange

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
PData.Region(4).PixelsLA = unique([PData.Region(2).PixelsLA;PData.Region(3).PixelsLA]); %;PData.Region(+4).PixelsLA]);
PData.Region(4).Shape.Name = 'Custom';
PData.Region(4).numPixels = length(PData.Region(4).PixelsLA);

assignin('base','PData',PData);

Control = repmat(struct('Command',[],'Parameters',[]),1,2);
Control(1).Command = 'set&Run';
Control(1).Parameters = {'DisplayWindow',2,'ReferencePt',RefPt};
Control(2).Command = 'update&Run';
Control(2).Parameters = {'PData'};
assignin('base','Control', Control);
%XSectionChange

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
%ZSectionChange

%CompressionP
CompressionFactor = UIValue;
assignin('base','CompressionFactor',CompressionFactor);
%CompressionP



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

%SwitchRecon3D

%-EF#3
ModifVSXgui();
h4=findobj('Tag','CLSlider');
h4.Callback= @CLSliderCallback;

h4=findobj('Tag','CLValue');
h4.Callback= @CLValueCallback;


assignin('base','h4',h4);

%-EF#3

%-EF#4
CLSliderCallback(varargin);

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
%-EF#4

%-EF#5
CLValueCallback(varargin);

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

%-EF#5

%-EF#6
vstic(varargin)

tic
%-EF#6

%-EF#7
vstoc(varargin)

toc
%-EF#7


%% Subfunction for creating random apertures
function [Ap2048]=Rnd512of2048Aperture()
% This function creates a Random Aperture taking into account the
% constrains of the MUX switches for the V-512 SHI 2048MUX.
% No input is required for this function. The output
% of the function is a vector containing 2048 elements with "ones" for the
% 512 elements activated for that particular aperture. The function
% verifies that no two elements connected to the same MUX switch can be
% activated at the same time.

el512=0;
el512rep=0;

el1024=0;
el1024rep=0;

el1536=0;
el1536rep=0;

el2048=0;
el2048rep=0;

for i=1:512
    de=0;

    elem(i)=randi([1 2048],1,1);
    while de ~= 1
        if elem(i) <=512
            if sum(elem == elem(i))>1 || sum(elem == elem(i)+512)~=0  || sum(elem == elem(i)+1024)~=0 || sum(elem == elem(i)+1536)~=0
                elem(i)=randi([1 512],1,1);
                el512rep=el512rep+1;
            else
                de=1;
                el512=el512+1;
            end
        end

        if elem(i)>512 && elem(i) <=1024
            if sum(elem == elem(i))>1 || sum(elem == elem(i)-512)~=0   || sum(elem == elem(i)+512)~=0  || sum(elem == elem(i)+1024)~=0
                elem(i)=randi([513 1024],1,1);
                el1024rep=el1024rep+1;
            else
                de=1;
                el1024=el1024+1;
            end
        end

        if elem(i)>1024 && elem(i) <=1536
            if sum(elem == elem(i))>1 || sum(elem == elem(i)-1024)~=0  || sum(elem == elem(i)-512)~=0 || sum(elem == elem(i)+512)~=0
                elem(i)=randi([1025 1536],1,1);
                el1536rep=el1536rep+1;
            else
                de=1;
                el1536=el1536+1;
            end
        end

        if elem(i)>1536 && elem(i) <=2048
            if sum(elem == elem(i))>1 || sum(elem == elem(i)-1536)~=0 || sum(elem == elem(i)-1024)~=0  || sum(elem == elem(i)-512)~=0
                elem(i)=randi([1537 2048],1,1);
                el2048rep=el2048rep+1;
            else
                de=1;
                el2048=el2048+1;
            end
        end

    end
end

Ap2048=zeros(1,2048);
for i=1:512
Ap2048(elem(i))=1;
end
% ap1024=reshape(Ap1024,32,32);

end


%-EF#1
function volumetricPlot(ImgData)

persistent handle3D

PData = evalin('base','PData');
Trans = evalin('base','Trans');
P = evalin('base','P');
v_depth = 1;
CFactor = evalin('base','CompressionFactor');
ImgData=ImgData.^(1/2);
ImgData(:,:,1:v_depth)=0;
ImgData = flipdim(ImgData,3);

if isempty(handle3D) || ~ishandle(handle3D)
    figure('Position', [1085 116 414 391]);
    handle3D = axes;
end

set(handle3D,'NextPlot','replacechildren');
vol3d('cdata', ImgData,'texture', '3D','Parent', handle3D);
grid(handle3D, 'on'); colormap(handle3D,'gray');
set(handle3D,'NextPlot','add');
alphamap('rampup');
alphamap(5 .* alphamap);
alphamap('decrease')
alphamap('decrease')
alphamap('decrease')
alphamap('decrease')


xl=PData(1).Size(2);
yl=PData(1).Size(1);
zl=PData(1).Size(3);
dx=PData(1).PDelta(2);
dy=PData(1).PDelta(1);
dz=PData(1).PDelta(3);
plot3(handle3D,Trans.ElementPosWL(:,1)./dx+(xl-1)/2,Trans.ElementPosWL(:,2)./dy+(yl-1)/2,(Trans.ElementPosWL(:,3)+P.endDepth)./dz,'k.');
view(handle3D,-45,30);
end
%-EF#1

%-EF#2
function SlicePlot(ImgData)

persistent slicehandle3D

PData = evalin('base','PData');
Trans = evalin('base','Trans');
P = evalin('base','P');
RS=evalin('base','Resource');
UI=evalin('base','UI');

zApex=evalin('base','zApex');
waveLength=evalin('base','waveLength');
v_depth = 1;
CFactor = evalin('base','CompressionFactor');
tImgData=20*log10(ImgData(:,:,P.startDepth+1:P.endDepth)/max(max(max(ImgData(:,:,10:end)))));

if isempty(slicehandle3D) || ~ishandle(slicehandle3D)
    figure('Position', [640 659 560 420]);
    slicehandle3D = axes;
end

% keyboard
x=linspace(-(round(size(ImgData,1)/2)-1),round(size(ImgData,1)/2)-1,size(ImgData,1));
y=linspace(round(size(ImgData,2)/2)-1,-(round(size(ImgData,2)/2)-1),size(ImgData,2));
z=P.startDepth+1:P.endDepth;

ysect= RS.DisplayWindow(1).ReferencePt(2);
xsect= RS.DisplayWindow(2).ReferencePt(1);
set(slicehandle3D,'NextPlot','replacechildren');
slice(slicehandle3D,x,y,z,tImgData,xsect,ysect,[]),shading(slicehandle3D,'interp'),caxis(slicehandle3D,[-60 0]),set(slicehandle3D,'Zdir','reverse')

grid(slicehandle3D, 'on'); colormap(slicehandle3D,'gray');
set(slicehandle3D,'NextPlot','add');

xl=PData(1).Size(2);
yl=PData(1).Size(1);
zl=PData(1).Size(3);
dx=PData(1).PDelta(2);
dy=PData(1).PDelta(1);
dz=PData(1).PDelta(3);
plot3(slicehandle3D,Trans.ElementPosMm(:,1),Trans.ElementPosMm(:,2),Trans.ElementPosMm(:,3),'r.');
view(slicehandle3D,-45,30);

if ~isvalid(UI(6).handle) && ~ishandle(UI(6).handle)
    return
end

UIState=get(UI(6).handle,'Value');
if UIState==1 % 1 means pressed
    if evalin('base',[ 'exist(''handle3D'',''var'')'])
        handle3D=evalin('base','handle3D');
        handle3Daxes=handle3D.CurrentAxes;
    else
        handle3D=figure;
        set(handle3D,'Position', [1085 116 414 391]);
        handle3Daxes = axes;
        view(handle3Daxes,-45,30);
    end

    pl3D = findall(handle3Daxes, 'Tag', 'Transducer');
    set(handle3Daxes,'Zdir','reverse','NextPlot','replacechildren');
    vol3d('cdata', ImgData.^.5,'texture', '3D','Parent', handle3Daxes);
    grid(handle3Daxes, 'on'); colormap(handle3Daxes,'gray');
    alphamap(5 .* alphamap);
    alphamap(5 .* alphamap);
    alphamap('decrease')
    alphamap('decrease')
    alphamap('decrease')

    if isempty(pl3D)
        xl=PData(1).Size(2);
        yl=PData(1).Size(1);
        zl=PData(1).Size(3);
        dx=PData(1).PDelta(2);
        dy=PData(1).PDelta(1);
        dz=PData(1).PDelta(3);

        hold(handle3Daxes, 'on');
        plot3(handle3Daxes, Trans.ElementPosWL(:,1)./dx+(xl-1)/2,Trans.ElementPosWL(:,2)./dy+(yl-1)/2,(Trans.ElementPosWL(:,3))./dz,'k.', 'Tag', 'Transducer');
        hold(handle3Daxes, 'off');
    end
    assignin('base','handle3D',handle3D);
end
end

%-EF#2