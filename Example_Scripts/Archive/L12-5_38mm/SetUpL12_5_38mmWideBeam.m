% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name SetUpL12_5_38mmWideBeam.m:
%
% Description:
%   Create wide beam sequence for L12-5_38mm 38mm Linear array on Vantage system.
%
% In the following aperture examples, each space represents 2 elements.
% Aperture 1 (P.numTx=48):
%   tttttttttttt---------------------------------------------------\--------------------------------
%   rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr--------------------------------
% Aperture for almost middle wide beam 24 (P.numTx=48, P.numRays=48):
%   ---------------/-------------------tttttttttttttttttttttttt-------------------\-----------------
%   ---------------rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr-----------------
% Aperture for last wide beam (P.numTx=48, P.numRays=48):
%   --------------------------------/---------------------------------------------------tttttttttttt
%   --------------------------------rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
%
% Last update:
% 12/13/2015 - modified for SW 3.0

clear all

P.startDepth = 2;
P.endDepth = 192;
P.numTx = 48;       % Number of active transmit elements in the aperture.
P.numRays = 48;     % Number of ray lines in frame.
P.txFocus = 600;    % transmit focal pt in wavelengths

% Define system parameters.
Resource.Parameters.numTransmit = 128;  % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;  % number of receive channels.
Resource.Parameters.speedOfSound = 1540;  % speed of sound in m/s
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Specify Trans structure array.
Trans.name = 'L12-5 38mm';
Trans.units = 'wavelengths';
Trans = computeTrans(Trans);  % L12-5_38mm transducer is 'known' transducer so we can use computeTrans.

% Specify PData structure array.
PData.PDelta = [1.0, 0, 0.5];
PData.Size(1) = ceil((P.endDepth-P.startDepth)/PData.PDelta(3)); % startDepth, endDepth and pdelta set PData.Size.
PData.Size(2) = ceil((Trans.numelements*Trans.spacing)/PData.PDelta(1));
PData.Size(3) = 1;      % single image page
PData.Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of upper lft crnr.
% Define P.numRays rectangular regions centered on TX beam origins. Define the width to be 3 times the
%   spacing between rays.
PData(1).Region = repmat(struct('Shape',struct('Name','Rectangle',...
                                               'Position',[0,0,P.startDepth],...
                                               'width',3*(192*Trans.spacing)/P.numRays,...
                                               'height',P.endDepth-P.startDepth)),1,P.numRays);
% -- Determine x coords. of TX aperture centers.
TxOrgX = (-95.5*Trans.spacing):(192*Trans.spacing/P.numRays):(95.5*Trans.spacing);
for n = 1:P.numRays, PData(1).Region(n).Shape.Position(1) = TxOrgX(n); end
PData(1).Region = computeRegions(PData(1));

% Specify Media object. 'pt1.m' script defines array of point targets.
pt1;
Media.attenuation = -0.5;
Media.function = 'movePoints';

% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = P.numRays*4096; % this size allows for max range
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 10;        % 10 frames used for RF cineloop.
Resource.InterBuffer(1).datatype = 'complex double';
Resource.InterBuffer(1).numFrames = 1;
Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).numFrames = 10;
Resource.DisplayWindow(1).Title = 'L12-5_38mmWideBeam';
Resource.DisplayWindow(1).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).numFrames = 20;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = gray(256);

% Specify Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,0.67,2,1];

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'aperture', 1, ...
                   'Apod', zeros(1,Resource.Parameters.numTransmit), ...
                   'focus', P.txFocus, ...
                   'Steer', [0.0,0.0], ...
                   'Delay', zeros(1,Resource.Parameters.numTransmit), ...
                   'TXPD', [], ...
                   'peakCutOff', 1.0, ...
                   'peakBLMax', 4.0), 1, P.numRays);

scaleToWvl = 1;
if strcmp(Trans.units, 'mm')
    scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);
end

% - Set event specific TX attributes.
for n = 1:P.numRays
    TX(n).Origin(1) = TxOrgX(n);
    % Compute transmit aperture
    % - find closest element to center of aperture.
    [Dummy,ce] = min(abs(scaleToWvl*Trans.ElementPos(:,1)-TxOrgX(n)));
    lft = ce - 64;
    if lft < 1, lft = 1; end
    if lft > 65, lft = 65; end
    TX(n).aperture = lft;
    % Compute TX.Apod within mux aperture.
    lft = ce - P.numTx/2;
    if lft < 1, lft = 1; end
    rt = ce + P.numTx/2;
    if rt > 192, rt = 192; end
    TX(n).Apod((lft-(TX(n).aperture-1)):(rt-(TX(n).aperture-1))) = 1;
    % Apply apodization function.
    [RIndices,CIndices,V] = find(TX(n).Apod);
    V = kaiser(size(V,2),1);
    TX(n).Apod(CIndices) = V;
    % Compute transmit delays
    TX(n).Delay = computeTXDelays(TX(n));
    % Compute transmit pixel data
    TX(n).TXPD = computeTXPD(TX(n),PData(1));
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [0,487,550,611,722,832,842,1001];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Receive structure arrays
%   endDepth - add additional acquisition depth to account for some channels
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
Receive = repmat(struct('aperture', 1, ...
                        'Apod', ones(1,128), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 0),1,P.numRays*Resource.RcvBuffer(1).numFrames);

% - Set event specific Receive attributes for each frame.
for i = 1:Resource.RcvBuffer(1).numFrames
    Receive(P.numRays*(i-1)+1).callMediaFunc = 1;
    for j = 1:P.numRays
        Receive(P.numRays*(i-1)+j).aperture = TX(j).aperture; % mux aperture same as transmit
        Receive(P.numRays*(i-1)+j).framenum = i;
        Receive(P.numRays*(i-1)+j).acqNum = j;
    end
end

% Specify Recon structure arrays.
% - We need one Recon structure which will be used for each frame.
Recon = struct('senscutoff', 0.5, ...
               'pdatanum', 1, ...
               'rcvBufFrame',-1, ...
               'IntBufDest',[1,1], ...
               'ImgBufDest',[1,-1], ...
               'RINums',1:P.numRays);

% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 'accumIQ', ...
                   'Pre',[],...
                   'Post',[],...
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'scaleFactor', 0.3, ...
                   'regionnum', 0, ...
                   'threadSync', 1), 1, P.numRays);
% - Set specific ReconInfo attributes.
ReconInfo(1).Pre = 'clearInterBuf';
for j = 1:P.numRays
    ReconInfo(j).txnum = j;
    ReconInfo(j).rcvnum = j;
    ReconInfo(j).regionnum = j;
end
ReconInfo(P.numRays).Post = 'IQ2IntensityImageBuf';

% Specify Process structure array.
pers = 20;
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',1.0,...            % pgain is image processing gain
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

% Specify SeqControl structure arrays.
%  - Jump back to start.
SeqControl(1).command = 'jump';
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between rayline acquisitions
SeqControl(2).argument = 160;  % 160 usec
SeqControl(3).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(3).argument = 1000;  % 1000 usec = 1 msec time between frames
SeqControl(4).command = 'returnToMatlab';
nsc = 5; % nsc is count of SeqControl objects

% Specify Event structure arrays.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    for j = 1:P.numRays        % Acquire frame
        Event(n).info = 'acquire rayline';
        Event(n).tx = j;
        Event(n).rcv = P.numRays*(i-1)+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;
    end
    Event(n-1).seqControl = [3,nsc]; % modify last acquisition Event's seqControl
      SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
      nsc = nsc+1;

    Event(n).info = 'recon and process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 1;
    Event(n).seqControl = 0;
    if (floor(i/5) == i/5)&&(i ~= Resource.RcvBuffer(1).numFrames)  % Exit to Matlab every 5th frame
        Event(n).seqControl = 4;
    end
    n = n+1;
end

Event(n).info = 'Jump back';
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
UI(1).Callback = text2cell('%SensCutOffCallback');

% - Range Change
MinMaxVal = [64,300,P.endDepth]; % default unit is wavelength
AxesUnit = 'wls';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
        AxesUnit = 'mm';
        MinMaxVal = MinMaxVal * (Resource.Parameters.speedOfSound/1000/Trans.frequency);
    end
end
UI(2).Control = {'UserA1','Style','VsSlider','Label',['Range (',AxesUnit,')'],...
                 'SliderMinMaxVal',MinMaxVal,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f'};
UI(2).Callback = text2cell('%RangeChangeCallback');

% - Peak CutOff
UI(3).Control = {'UserB2','Style','VsSlider','Label','Peak Cutoff',...
                  'SliderMinMaxVal',[0,20.0,TX(1).peakCutOff],...
                  'SliderStep',[0.005,0.020],'ValueFormat','%1.3f'};
UI(3).Callback = text2cell('%PeakCutOffCallback');

% - Max. Burst Length
UI(4).Control = {'UserB1','Style','VsSlider','Label','Max. BL',...
                  'SliderMinMaxVal',[0,20.0,TX(1).peakBLMax],...
                  'SliderStep',[0.005,0.020],'ValueFormat','%1.3f'};
UI(4).Callback = text2cell('%MaxBLCallback');

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 5;
% Save all the structures to a .mat file.
save('MatFiles/L12-5_38mmWideBeam');

return

% **** Callback routines to be converted by text2cell function. ****
%SensCutOffCallback - Sensitivity cutoff change
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
%SensCutOffCallback

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

PData = evalin('base','PData');
PData.Size(1) = ceil((P.endDepth-P.startDepth)/PData.PDelta(3));
PData.Region = repmat(struct('Shape',struct( ...
                    'Name','Rectangle',...
                    'Position',[0,0,P.startDepth],...
                    'width',3*(192*Trans.spacing)/P.numRays,...
                    'height',P.endDepth-P.startDepth)),1,P.numRays);
% - set position of regions to correspond to beam spacing and aperture centers.
TxOrgX = (-95.5*Trans.spacing):(192*Trans.spacing/P.numRays):(95.5*Trans.spacing);
for n = 1:P.numRays, PData(1).Region(n).Shape.Position(1) = TxOrgX(n); end
assignin('base','PData',PData);
evalin('base','PData.Region = computeRegions(PData);');
% Update TXPD data of TX structures.
TX = evalin('base','TX');
for i = 1:size(TX,2)
    TX(i).TXPD = computeTXPD(TX(i),PData);
end
assignin('base','TX',TX);
% Update Receive structures.
Receive = evalin('base', 'Receive');
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
for i = 1:size(Receive,2)
    Receive(i).endDepth = maxAcqLength;
end
assignin('base','Receive',Receive);
evalin('base','TGC.rangeMax = P.endDepth;');
evalin('base','TGC.Waveform = computeTGCWaveform(TGC);');
evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData.Size(1)*PData.PDelta(3)/Resource.DisplayWindow(1).pdelta);');
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'PData','Receive','TGC','TX','Recon','DisplayWindow','InterBuffer','ImageBuffer'};
assignin('base','Control', Control);
assignin('base', 'action', 'displayChange');
return
%RangeChangeCallback

%PeakCutOffCallback
TX = evalin('base', 'TX');
for i=1:size(TX,2)
    TX(i).peakCutOff = UIValue;
end
assignin('base','TX',TX);
% Set Control.Command to set TX
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'TX'};
assignin('base','Control', Control);
return
%PeakCutOffCallback

%MaxBLCallback
TX = evalin('base', 'TX');
for i=1:size(TX,2)
    TX(i).peakBLMax = UIValue;
end
assignin('base','TX',TX);
% Set Control.Command to set TX
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'TX'};
assignin('base','Control', Control);
return
%MaxBLCallback
