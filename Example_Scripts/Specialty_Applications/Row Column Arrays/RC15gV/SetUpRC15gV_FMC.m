% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: SetUpRC15gV_FMC.m - Example of full matrix capture imaging with row/col
%                                   array.
% Description:
%   Sequence programming file for 80x80 element Row/Column array, 80 transmits
%   on single elements in both the x and y directions.  160 elements are used on
%   receive for each acquisition.  Two slices are reconstructed - one in the XZ
%   plane and the other in the YZ plane.
%
%
% Last update:
%   03/27/23: VTS-5945 script adapted from RC6gV

clear all
P.startDepth = 0;   % Acquisition depth in wavelengths
P.endDepth = 192;   % EndDepth in wavelengths

ReconRegion=5;
xyplane=100; % XY plane for reconstruction (in waveleghts)

% Define system parameters.
Resource.Parameters.numTransmit = 256;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 256;    % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1; %forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2; %stops sequence and processes RcvData continuously.
Resource.Parameters.waitForProcessing=1; %Set waitForProcessing for synchronous script

% Specify Trans structure array.
Trans.name = 'RC15gV'; % Vermon 15 MHz Row Column Array
Trans.units = 'wavelengths';
Trans = computeTrans(Trans);

% Specify PData structure array.
PData.PDelta = [Trans.spacing, Trans.spacing, 0.5];
PData.Size(1) = 80;
PData.Size(2) = 80;
PData.Size(3) = ceil((P.endDepth-P.startDepth)/PData.PDelta(3)); % startDepth, endDepth and pdelta set PData.Size(3).
PData.Origin = [-Trans.spacing*39.5,Trans.spacing*39.5,P.startDepth]; % x,y,z of upper lft crnr of page 1.
PData.Region(1) = struct('Shape',struct('Name','PData'));
PData(1).Region(2) = struct('Shape',struct('Name','Slice','Orientation','yz','oPAIntersect',PData.Origin(1)+PData.PDelta(2)*(PData.Size(2)-1)/2,'andWithPrev',1));
PData(1).Region(3) = struct('Shape',struct('Name','Slice','Orientation','xz','oPAIntersect',PData.Origin(2)-PData.PDelta(1)*(PData.Size(1)-1)/2,'andWithPrev',1));
PData(1).Region(4) = struct('Shape',struct('Name','Slice','Orientation','xy','oPAIntersect',xyplane,'andWithPrev',1));

PData.Region = computeRegions(PData);

PData.Region(5).PixelsLA = unique([PData.Region(2).PixelsLA; PData.Region(3).PixelsLA; PData.Region(4).PixelsLA]);
PData.Region(5).Shape.Name = 'Custom';
PData.Region(5).numPixels = length(PData.Region(5).PixelsLA);

% Specify Media object.
Media.MP(1,:) = [-45.5,0,30,1.0];
Media.MP(2,:) = [-15.5,0,30,1.0];
Media.MP(3,:) = [15.5,0,30,1.0];
Media.MP(4,:) = [45.5,0,30,1.0];
Media.MP(5,:) = [-15.5,0,60,1.0];
Media.MP(6,:) = [-15.5,0,90,1.0];
Media.MP(7,:) = [-15.5,0,120,1.0];
Media.MP(8,:) = [-15.5,0,150,1.0];
Media.MP(9,:) = [-45.5,0,120,1.0];
Media.MP(10,:) = [15.5,0,120,1.0];
Media.MP(11,:) = [45.5,0,120,1.0];
Media.MP(12,:) = [-10.5,0,69,1.0];
Media.MP(13,:) = [-5.5,0,75,1.0];
Media.MP(14,:) = [0.5,0,78,1.0];
Media.MP(15,:) = [5.5,0,80,1.0];
Media.MP(16,:) = [10.5,0,81,1.0];
Media.MP(17,:) = [-75.5,0,120,1.0];
Media.MP(18,:) = [75.5,0,120,1.0];
Media.MP(19,:) = [-15.5,0,180,1.0];
k = 0; % -19
Media.MP(k+20,:) = [0,-45.5,30,1.0];
Media.MP(k+21,:) = [0,-15.5,30,1.0];
Media.MP(k+22,:) = [0,15.5,30,1.0];
Media.MP(k+23,:) = [0,45.5,30,1.0];
Media.MP(k+24,:) = [0,-15.5,60,1.0];
Media.MP(k+25,:) = [0,-15.5,90,1.0];
Media.MP(k+26,:) = [0,-15.5,120,1.0];
Media.MP(k+27,:) = [0,-15.5,150,1.0];
Media.MP(k+28,:) = [0,-45.5,120,1.0];
Media.MP(k+29,:) = [0,15.5,120,1.0];
Media.MP(k+30,:) = [0,45.5,120,1.0];
Media.MP(k+31,:) = [0,-10.5,69,1.0];
Media.MP(k+32,:) = [0,-5.5,75,1.0];
Media.MP(k+33,:) = [0,0.5,78,1.0];
Media.MP(k+34,:) = [0,5.5,80,1.0];
Media.MP(k+35,:) = [0,10.5,81,1.0];
Media.MP(k+36,:) = [0,-75.5,120,1.0];
Media.MP(k+37,:) = [0,75.5,120,1.0];
Media.MP(k+38,:) = [0,-15.5,180,1.0];
Media.MP(1,:) = [0,0,140,1.0];

Media.numPoints = size(Media.MP,1);
Media.attenuation = -0.5;
Media.function = 'movePoints';

Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 256*4096; % this size allows for maximum range
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 10;    % 10 frames stored in RcvBuffer.
Resource.InterBuffer(1).numFrames = 1;   % one intermediate buffer needed.
Resource.ImageBuffer(1).numFrames = 1;

Resource.DisplayWindow(1).Type = 'Matlab';
Resource.DisplayWindow(1).Title = '3D RC15gV FMC - XZ plane';
Resource.DisplayWindow(1).pdelta = 0.4;
Resource.DisplayWindow(1).Position = [0,40, ...
    ceil(PData(1).Size(2)*PData(1).PDelta(2)/Resource.DisplayWindow(1).pdelta), ... % width
    ceil(PData(1).Size(3)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta)];    % height
Resource.DisplayWindow(1).Orientation = 'xz';
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,0.0];
Resource.DisplayWindow(1).Colormap = gray(256);
Resource.DisplayWindow(1).AxesUnits = 'wavelengths';
Resource.DisplayWindow(1).numFrames = 20;

Resource.DisplayWindow(2).Type = 'Matlab';
Resource.DisplayWindow(2).Title = '3D RC15gV FMC - YZ plane';
Resource.DisplayWindow(2).pdelta = Resource.DisplayWindow(1).pdelta;
Resource.DisplayWindow(2).Position = [430,40, ...
    ceil(PData(1).Size(1)*PData(1).PDelta(1)/Resource.DisplayWindow(2).pdelta), ... % width
    ceil(PData(1).Size(3)*PData(1).PDelta(3)/Resource.DisplayWindow(2).pdelta)];    % height
Resource.DisplayWindow(2).Orientation = 'yz';
Resource.DisplayWindow(2).ReferencePt = [0,-PData(1).Origin(2),0];
Resource.DisplayWindow(2).Colormap = gray(256);
Resource.DisplayWindow(2).AxesUnits = 'wavelengths';
Resource.DisplayWindow(2).numFrames = 20;

Resource.DisplayWindow(3).Type = 'Matlab';
Resource.DisplayWindow(3).Title = '3D RC15gV FMC - XY plane';
Resource.DisplayWindow(3).pdelta = Resource.DisplayWindow(1).pdelta;
Resource.DisplayWindow(3).Position = [0,40, ...
    ceil(PData(1).Size(2)*PData(1).PDelta(2)/Resource.DisplayWindow(3).pdelta), ... % width
    ceil(PData(1).Size(1)*PData(1).PDelta(1)/Resource.DisplayWindow(3).pdelta)];    % height
Resource.DisplayWindow(3).Orientation = 'xy';
Resource.DisplayWindow(3).ReferencePt = [PData(1).Origin(1),-PData(1).Origin(2),xyplane];%PData.Region(end).Shape.oPAIntersect];
Resource.DisplayWindow(3).Colormap = gray(256);
Resource.DisplayWindow(3).AxesUnits = 'wavelengths';
Resource.DisplayWindow(3).numFrames = 20;

% Specify TW structure array.
TW.type = 'parametric';
TW.Parameters = [Trans.frequency,0.67,2,1];

% Specify TX structure array. Define steered plane wave transmits for the elements
% along the x and y axes.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'Apod', zeros(1,160), ...
                   'focus', 0.0, ...
                   'Steer', [0.0,0.0], ...
                   'Delay', zeros(1,Trans.numelements)), 1, 256);
% - Set event specific TX attributes.
for n = 1:160   % for 256 transmit events
    TX(n).Apod(n) = 1;
    TX(n).Delay = computeTXDelays(TX(n));
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [0,245,424,515,627,661,717,744];
TGC.CntrlPts = [0 785.2216 1023 1023 1023 1023 1023 1023];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Receive structure arrays.
% - We need na Receives for every frame.
maxAcqLength = ceil(sqrt(P.endDepth^2 + (80*Trans.spacing)^2));
Receive = repmat(struct('Apod', zeros(1,160), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength,...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, 160*Resource.RcvBuffer(1).numFrames);

% - Set event specific Receive attributes for each frame.
for i = 1:Resource.RcvBuffer(1).numFrames
%     Receive(na*(i-1)+1).callMediaFunc = 1;
    k = 160*(i-1);
    for j = 1:160
        if j<81
            Receive(k+j).Apod(81:160) = 1;
        else
            Receive(k+j).Apod(1:80) = 1;
        end 
        Receive(k+j).framenum = i;
        Receive(k+j).acqNum = j;
    end
end

% Specify Recon structure arrays.
% - We need one Recon structure which will be used for each frame.
Recon = repmat(struct('senscutoff', 0.6, ...
                      'pdatanum', 1, ...
                      'rcvBufFrame',-1, ...
                      'IntBufDest', [1,1], ...
                      'ImgBufDest', [1,-1], ...
                      'RINums', 1:160), 1, 1);
% Define ReconInfo structures.
% We need 256 ReconInfo structures.
ReconInfo = repmat(struct('mode', 'accumIQ', ...  % default is to accumulate IQ data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'scaleFactor', 1.0, ...
                   'regionnum', 1), 1, 160);
% - Set specific ReconInfo attributes.

na=80;
ReconInfo(1).mode = 'replaceIQ'; % replace IQ data
for j = 1:na  % For each transmit along x elements
    ReconInfo(j).txnum = j;
    ReconInfo(j).rcvnum = j;
    ReconInfo(j).regionnum = ReconRegion; %1 for the whole volume, 4 for the slices
end
ReconInfo(na).Post = 'IQ2IntensityImageBuf';
ReconInfo(na+1).Pre = 'clearInterBuf';
for j = 1:na  % For each transmit along y elements
    ReconInfo(na+j).txnum = na+j;
    ReconInfo(na+j).rcvnum = na+j;
    ReconInfo(na+j).regionnum = ReconRegion; %1 for the whole volume, 4 for the slices
end
ReconInfo(2*na).Post = 'IQ2IntensityImageBufAdd' ;

% Specify Process structure array.
pers = 20;
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'srcData','intensity3D',...
                         'pgain',8,...     % pgain is image processing gain
                         'reject',2,...      % reject level
                         'persistMethod','simple',...
                         'persistLevel',pers,...
                         'interpMethod','4pt',...
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',40,...
                         'mappingMethod','full',...
                         'display',1,...      % display image after processing
                         'displayWindow',1};
Process(2).classname = 'Image';
Process(2).method = 'imageDisplay';
Process(2).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'srcData','intensity3D',...
                         'pgain',8,...     % pgain is image processing gain
                         'reject',2,...      % reject level
                         'persistMethod','simple',...
                         'persistLevel',pers,...
                         'interpMethod','4pt',...
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',40,...
                         'mappingMethod','full',...
                         'display',1,...      % display image after processing
                         'displayWindow',2};
Process(3).classname = 'Image';
Process(3).method = 'imageDisplay';
Process(3).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'srcData','intensity3D',...
                         'pgain',8,...     % pgain is image processing gain
                         'reject',2,...      % reject level
                         'persistMethod','simple',...
                         'persistLevel',pers,...
                         'interpMethod','4pt',...
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',40,...
                         'mappingMethod','full',...
                         'display',1,...      % display image after processing
                         'displayWindow',3};

% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; % jump back to start
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(2).argument = 160;  % 160 usec
SeqControl(3).command = 'timeToNextAcq';  % time between frames
SeqControl(3).argument = 40000 - (256-1)*150;  % 50 msec (50 Hz)
SeqControl(4).command = 'returnToMatlab';
nsc = 5; % nsc is count of SeqControl objects

% Specify Event structure arrays.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    k = 160*(i-1);
    for j = 1:160                      % Acquire frame
        Event(n).info = 'xmit on each element';
        Event(n).tx = j;
        Event(n).rcv = k+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;
    end
    Event(n-1).seqControl = [3,nsc]; % modify last acquisition Event's seqControl
      SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
      nsc = nsc+1;

    Event(n).info = 'recon';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 0;
    Event(n).seqControl = 0;
    n = n+1;
%
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
    if floor(i/3) == i/3     % Exit to Matlab every 3rd frame reconstructed
        Event(n).seqControl = 4;
    end
    n = n+1;

    if ReconRegion==1
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

end

Event(n).info = 'Jump back';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 1;


%% User specified UI Control Elements
import vsv.seq.uicontrol.VsSliderControl;

% - Sensitivity Cutoff
UI(1).Control = VsSliderControl('LocationCode','UserB7','Label','Sens. Cutoff',...
                  'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                  'SliderStep',[0.025,0.1],'ValueFormat','%1.3f', ...
                  'Callback', @SensCutoffCallback );

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 3;

% Save all the structures to a .mat file.
save('MatFiles/RC15gV_FMC');
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