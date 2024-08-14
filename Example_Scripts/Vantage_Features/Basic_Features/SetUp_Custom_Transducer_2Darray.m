% File name SetUp3DMtrixFlash.m:
%
% Set up file for Conical scan format, flash transmit.  Transducer has 256
%   elements arranged in a 2D grid (16x16).  The angle of the scan is
%   alpha, which is nominally around 20 degrees.
%
%   The 3D volume region is a cube with the transducer centered in one face.
%   The region imaged is an expanding cone, from a virtual apex behind
%      the array, as shown below.
%
%           +y                  .
%            |  +x         .  .   .
%            |/./     .      /     \
%            - |  .         .       .
%         . /|/.            .       .
% -z < + - . /--------------------------------- +z
%         ../|/             .       .
%          | -              .       .
%         /./|    .          \     /
%       -x   |          .     .   .
%            |                  .
%           -y
%
% Last update:
% 05/06/2020 - Update to SW 4.3 format for new user UIControls and External function definitions (VTS 1691).
%   More info:(?/Example_Scripts/Vantage_Features/New UI Scheme/SetUpL11_5vFlash_NewUI)

clear all

P.startDepth = 0;
P.endDepth = 96;
P.angle = 20 * pi/180; % P.angle between z axis and surface of cone.

% Define system parameters.
Resource.Parameters.numTransmit = 128;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;    % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Specify Trans structure array.
Trans.name = 'Matrix';
Trans.units = 'wavelengths';
Trans.frequency = 2.5;
Trans.type = 2;      % Array geometry is flat 2D array (z=0 for all elements).
Trans.numelements = 128;
Trans.elementWidth = 1.0;  % This is the x element width in wavelengths.
Trans.spacing = Trans.elementWidth;
Trans.impedance = 50;
Trans.maxHighVoltage = 50;
eleX = 1.0;  % element spacing in x direction
eleY = 1.0;  % element spacing in y direction
% Define a square array EM (Element Map) that shows the location of
% each element (numbered 1 to 128) on the physical grid of transducer
% elements.  Locations that represent either an element that does not exist
% or an element that is not connected are identified by setting the element
% number to zero.
GS = 16;  % GS (Grid Size) sets the size of the square array of elements, GS X GS
EM = zeros(GS,GS);  % define the array and fill it with zeros

% Now manually specify each element position in the following lines
EM(1,:)  = [  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0];
EM(2,:)  = [  0   0   0   0   0   0   0   1   2   0   0   0   0   0   0   0];
EM(3,:)  = [  0   0   0   0   0   3   4   5   6   7   8   0   0   0   0   0];
EM(4,:)  = [  0   0   0   0   9  10  11  12  13  14  15  16   0   0   0   0];
EM(5,:)  = [  0   0   0  17  18  19  20  21  22  23  24  25  26   0   0   0];
EM(6,:)  = [  0   0  27  28  29  30  31  32  33  34  35  36  37  38   0   0];
EM(7,:)  = [  0   0  39  40  41  42  43  44  45  46  47  48  49  50   0   0];
EM(8,:)  = [  0  51  52  53  54  55  56  57  58  59  60  61  62  63  64   0];
EM(9,:)  = [  0  65  66  67  68  69  70  71  72  73  74  75  76  77  78   0];
EM(10,:) = [  0   0  79  80  81  82  83  84  85  86  87  88  89  90   0   0];
EM(11,:) = [  0   0  91  92  93  94  95  96  97  98  99 100 101 102   0   0];
EM(12,:) = [  0   0   0 103 104 105 106 107 108 109 110 111 112   0   0   0];
EM(13,:) = [  0   0   0   0 113 114 115 116 117 118 119 120   0   0   0   0];
EM(14,:) = [  0   0   0   0   0 121 122 123 124 125 126   0   0   0   0   0];
EM(15,:) = [  0   0   0   0   0   0   0 127 128   0   0   0   0   0   0   0];
EM(16,:) = [  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0];

% Now test for obvious data entry errors
elcount = 0;
for x=1:GS
    for y=1:GS
        if(EM(y,x)<0 || EM(y,x)>128)
            error('data error: illegal entry in row %d column %d\n', y, x);
        end
        if(EM(y,x)>0)
            elcount = elcount + 1;
        end
    end
end
if(elcount~=128)
    error('data error: total number of active element positions is not 128');
end
clear elcount

% Next use EM array to assign xy position for each element in Trans.ElementPos array
Trans.ElementPos = zeros(128,5,'double');  % define array and fill with zeros
Trans.ElementPos(:,1) = -1;  % fill x,y positions with -1, to allow test for
Trans.ElementPos(:,2) = -1;  % uninitialized entries later
for i=1:GS
    for j=1:GS
        if(EM(i,j)>0)   % note positions are relative to center of array
            Trans.ElementPos(EM(i,j),2) = ((GS-1)*eleY)/2 - (i-1)*eleY;
            Trans.ElementPos(EM(i,j),1) = (j-1)*eleX - ((GS-1)*eleX)/2;
        end
    end
end
clear GS EM

% Test for uninitialized entries
for i=1:128
    if(Trans.ElementPos(i,1)==-1 || Trans.ElementPos(i,2)==-1)
        fprintf('data error: element %d not listed in Element Map array', i);
        return;
    end
end

% Set element sensitivity function (100 weighting values from -pi/2 to pi/2).
Theta = (-pi/2:pi/100:pi/2)';
Theta(51) = 0.0000001;
Trans.ElementSens = sin(Trans.elementWidth*pi*sin(Theta))./(Trans.elementWidth*pi*sin(Theta));
clear Theta

% Specify Media.  Use point targets in middle of PData.
%Media.MP(1,:) = [0,0,30,1.0];      % single point.
%Media.program = 'PointTargets3D';   % prog. to execute for creating media pts.
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

% Specify PData structure array.
PData.PDelta = [0.7,0.7,0.7];
% PData.Size computed so that PData x,y section contains max depth cross section.
radius = 7 * Trans.spacing;  % radius of transducer elements.
zApex = -radius/tan(P.angle);  % z value of cone apex.
PData.Size(1) = ceil(2.0*(P.endDepth-zApex)*tan(P.angle)/PData.PDelta(2));
PData.Size(2) = PData.Size(1);
PData.Size(3) = ceil((P.endDepth-P.startDepth)/PData.PDelta(3));
PData.Origin = [(-(PData.Size(2))/2)*PData.PDelta(1), ((PData.Size(1))/2)*PData.PDelta(1), P.startDepth];
PData.Region(1) = struct('Shape',struct('Name','Cone','Position',[0,0,zApex],'angle',P.angle,'z1',0,'z2',P.endDepth));
PData.Region(2) = struct('Shape',struct('Name','Slice','Orientation','xz','oPAIntersect',PData.Origin(2)-PData.PDelta(1)*PData.Size(1)/2,'andWithPrev',1));
PData.Region(3) = struct('Shape',struct('Name','Slice','Orientation','yz','oPAIntersect',PData.Origin(1)+PData.PDelta(2)*PData.Size(2)/2,'andWithPrev',1));
PData.Region(4) = struct('Shape',struct('Name','Slice','Orientation','xy','oPAIntersect',60,'andWithPrev',1));
% Calculate the Region data using 'computeRegions'.
PData.Region = computeRegions(PData);

% Set up Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 4096;
Resource.RcvBuffer(1).colsPerFrame = 128;  % this should match number of receive channels
Resource.RcvBuffer(1).numFrames = 4;
Resource.ImageBuffer(1).numFrames = 4;
Resource.DisplayWindow(1).Title = '3D Flash Image - XZ plane';
Resource.DisplayWindow(1).pdelta = 0.3;
Resource.DisplayWindow(1).Colormap = gray(256);
Resource.DisplayWindow(1).Position = [250,658, ...
    ceil(PData.Size(2)*PData.PDelta(2)/Resource.DisplayWindow(1).pdelta), ... % width
    ceil(PData.Size(3)*PData.PDelta(3)/Resource.DisplayWindow(1).pdelta)];    % height
Resource.DisplayWindow(1).Orientation = 'xz';
Resource.DisplayWindow(1).ReferencePt = [PData.Origin(1),0.0,0.0];
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(2).Title = '3D Flash Image - YZ plane';
Resource.DisplayWindow(2).pdelta = 0.3;
Resource.DisplayWindow(2).Colormap = gray(256);
Resource.DisplayWindow(2).Position = [640,658, ...
    ceil(PData.Size(1)*PData.PDelta(1)/Resource.DisplayWindow(2).pdelta), ... % width
    ceil(PData.Size(3)*PData.PDelta(3)/Resource.DisplayWindow(2).pdelta)];    % height
Resource.DisplayWindow(2).Orientation = 'yz';
Resource.DisplayWindow(2).ReferencePt = [0.0,-PData.Origin(2),0.0];
Resource.DisplayWindow(2).Type = 'Verasonics';
Resource.DisplayWindow(3).Title = '3D Flash Image - XY plane';
Resource.DisplayWindow(3).pdelta = 0.3;
Resource.DisplayWindow(3).Colormap = gray(256);
Resource.DisplayWindow(3).Position = [250,150, ...
    ceil(PData.Size(2)*PData.PDelta(2)/Resource.DisplayWindow(3).pdelta), ... % width
    ceil(PData.Size(1)*PData.PDelta(1)/Resource.DisplayWindow(3).pdelta)];    % height
Resource.DisplayWindow(3).Orientation = 'xy';
Resource.DisplayWindow(3).ReferencePt = [PData.Origin(1),-PData.Origin(2),60];
Resource.DisplayWindow(3).Type = 'Verasonics';


% Specify TW structure array.
TW.type = 'parametric';
TW.Parameters = [Trans.frequency,0.67,2,1];

% Specify TX structure array for flash transmit.
TX.waveform = 1;     % use 1st TW structure.
TX.Origin = [0.0,0.0,0.0];  % flash transmit origin at (0,0,0).
TX.focus = zApex;  % negative focus => concave TX.Delay profile.
TX.Steer = [0.0,0.0];       % theta, alpha = 0.
TX.Apod = ones(1,Resource.Parameters.numTransmit);
TX.Delay = computeTXDelays(TX);

% Specify Receive structure arrays -
%   P.endDepth - add additional acquisition depth to account for some channels
%              having longer path lengths.
maxAcqLength = ceil(sqrt(P.endDepth^2 + radius^2) - P.startDepth);
wlsPer128 = 128/(4*2); % wavelengths in 128 samples for 4 samplesPerWave
Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', P.startDepth + wlsPer128*ceil(maxAcqLength/wlsPer128), ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'mode', 0, ...
                        'callMediaFunc', 0),1,Resource.RcvBuffer(1).numFrames);
% Specify event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    Receive(i).framenum = i;
end

% Specify TGC Waveform structure.
TGC.CntrlPts = ([500,590,650,710,770,830,890,950]);
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Recon structure arrays.
Recon = struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame', -1, ...
               'ImgBufDest', [1,-1], ...  % auto-increment ImageBuffer each recon
               'RINums', 1:3);

% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 'replaceIntensity', ...  % replace Intensity.
                          'txnum', 1, ...
                          'rcvnum', 1, ...
                          'regionnum', 2),1,3);
% - Set specific ReconInfo attributes.
ReconInfo(2).regionnum = 3;
ReconInfo(3).regionnum = 4;

% Specify Process structure arrays.
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,... % number of buffer to process.
                         'framenum',-1,...
                         'pdatanum',1,...
                         'srcData','intensity3D',...
                         'pgain',1.0,...
                         'persistMethod','none',...
                         'persistLevel',30,...
                         'interpMethod','4pt',...    % method of interpolation
                         'compressMethod','power',...
                         'compressFactor',40,...
                         'mappingMethod','full',...
                         'display',1,...   % display image after processing
                         'displayWindow',1};
Process(2).classname = 'Image';
Process(2).method = 'imageDisplay';
Process(2).Parameters = {'imgbufnum',1,... % number of buffer to process.
                         'framenum',-1,...
                         'pdatanum',1,...
                         'srcData','intensity3D',...
                         'pgain',1.0,...
                         'persistMethod','none',...
                         'persistLevel',30,...
                         'interpMethod','4pt',...    % method of interpolation
                         'compressMethod','power',...
                         'compressFactor',40,...
                         'mappingMethod','full',...
                         'display',1,...   % display image after processing
                         'displayWindow',2};
Process(3).classname = 'Image';
Process(3).method = 'imageDisplay';
Process(3).Parameters = {'imgbufnum',1,... % number of buffer to process.
                         'framenum',-1,...
                         'pdatanum',1,...
                         'srcData','intensity3D',...
                         'pgain',1.0,...
                         'persistMethod','none',...
                         'persistLevel',30,...
                         'interpMethod','4pt',...    % method of interpolation
                         'compressMethod','power',...
                         'compressFactor',40,...
                         'mappingMethod','full',...
                         'display',1,...   % display image after processing
                         'displayWindow',3};

% Specify SeqControl structure arrays.
%  - Jump back to start.
SeqControl(1).command = 'jump';
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between frames
SeqControl(2).argument = 20000;  % 20ms
SeqControl(3).command = 'returnToMatlab';

nsc = 4; % nsc is count of SeqControl objects

n = 1; % n is count of Events

% Acquire all frames defined in RcvBuffer
for i = 1:Resource.RcvBuffer(1).numFrames
    Event(n).info = 'acquire frame';
    Event(n).tx = 1;
    Event(n).rcv = i;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = [2,nsc];
       SeqControl(nsc).command = 'transferToHost';
       nsc = nsc + 1;
    n = n+1;

    Event(n).info = 'Reconstruct';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 1;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = 'ProcessXY';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 2;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = 'ProcessXZ';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 3;
    Event(n).seqControl = 3;
    n = n+1;
end

Event(n).info = 'Jump back to first event';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 1;

% User specified UI Control Elements
% - Sensitivity cutoff
UI(1).Control =  vsv.seq.uicontrol.VsSliderControl(...
    'LocationCode','UserB7','Label','Sens. Cutoff',...
    'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
    'SliderStep',[0.025,0.1],'ValueFormat','%1.3f',...
    'Callback',@SensCutoff);

% - Section number change
UI(2).Control = vsv.seq.uicontrol.VsSliderControl(...
    'LocationCode','UserB2','Label','Section',...
    'SliderMinMaxVal',[1,P.endDepth,Resource.DisplayWindow(3).ReferencePt(3)],...
    'SliderStep',[1/(P.endDepth-1) 10/(P.endDepth-1)],...
    'ValueFormat','%3.0f','Callback',@SectionChange);

% Save all the structures to a .mat file.
save('MatFiles/3DMatrixFlash.mat');

%% Automatic VSX Execution:
% Uncomment the following line to automatically run VSX every time you run
% this SetUp script (note that if VSX finds the variable 'filename' in the
% Matlab workspace, it will load and run that file instead of prompting the
% user for the file to be used):

% filename = '3DMatrixFlash';  VSX;


%% **** Callback routines used by UIControls (UI) ****

function SensCutoff(~,~,UIValue)
% Sensitivity cutoff change
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

function SectionChange(~,~,UIValue)
    PData = evalin('base','PData');
    PData.Region(4).Shape.oPAIntersect = UIValue;
    PData.Region = computeRegions(PData);
    assignin('base','PData',PData);

    RS = evalin('base','Resource');
    RefPt = RS.DisplayWindow(3).ReferencePt;
    RefPt(3) = UIValue;
    RS.DisplayWindow(3).ReferencePt = RefPt;
    assignin('base','Resource',RS);

    Control = repmat(struct('Command',[],'Parameters',[]),1,2);
    Control(1).Command = 'set&Run';
    Control(1).Parameters = {'DisplayWindow',3,'ReferencePt',RefPt};
    Control(2).Command = 'update&Run';
    Control(2).Parameters = {'PData'};
    assignin('base','Control', Control);
end
