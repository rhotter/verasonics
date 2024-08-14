% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use.
%
% File name: SetUpL7_4SingleBeam_LE.m - Example of imaging with single focused
%                                      beam transmit
% Description:
%   Sequence programming file for L7-4 Linear array, using single beam
%   focused transmit with 2-1 synthetic aperture acquisation. The central
%   96 transmit and all 64 receive channels are active for each
%   synthetic acquisition. Processing are synchronous with respect to acquisition.
%
% Last update:
% 12/13/2015 - modified for SW 3.0

clear all
P.startDepth = 5;   % Acquisition depth in wavelengths
P.endDepth = 192;   % This should preferrably be a multiple of 128 samples.

% Define system parameters.
Resource.Parameters.numTransmit = 128;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 64;    % number of receive channels.
Resource.Parameters.speedOfSound = 1540;  % speed of sound in m/s
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Specify Trans structure array.
Trans.name = 'L7-4';
Trans.units = 'wavelengths'; % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);
Trans.maxHighVoltage = 50;  % set maximum high voltage limit for pulser supply.

% Specify PData structure array.
PData(1).PDelta = [Trans.spacing, 0, 0.5];
PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3)); % startDepth, endDepth and pdelta set PData(1).Size.
PData(1).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(3) = 1;      % single image page
PData(1).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of upper lft crnr.
% No PData.Region specified, so a default Region for the entire PData array will be created by computeRegions.

% Specify Media object. 'pt1.m' script defines array of point targets.
pt1;
Media.attenuation = -0.5;
Media.function = 'movePoints';

% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 4096;   % this size allows for maximum range
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numTransmit;
Resource.RcvBuffer(1).numFrames = 100;       % 100 frames used for RF cineloop.
Resource.InterBuffer(1).numFrames = 1;  % one intermediate buffer needed.
Resource.ImageBuffer(1).numFrames = 10;
Resource.DisplayWindow(1).Title = 'L7-4SingleBeam_LE';
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
TW(1).Parameters = [Trans.frequency,.67,2,1];

% Specify TX structure array.
TX(1).waveform = 1;            % use 1st TW structure.
TX(1).Origin = [0.0,0.0,0.0];  % flash transmit origin at (0,0,0).
TX(1).focus = 100;
TX(1).Steer = [0.0,0.0];       % theta, alpha = 0.
TX(1).Apod = zeros(1,Trans.numelements);
TX(1).Apod(17:112) = 1.0;
TX(1).Delay = computeTXDelays(TX(1));
TX(1).TXPD = computeTXPD(TX(1),PData);
TX(1).peakCutOff = 1.0;
TX(1).peakBLMax = 2.5;

% Specify TGC Waveform structure.
TGC.CntrlPts = [0,107,342,408,562,718,756,796];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Receive structure arrays -
%   endDepth - add additional acquisition depth to account for some channels
%              having longer path lengths.
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
Receive = repmat(struct('Apod', zeros(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength,...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 1),1,2*Resource.RcvBuffer(1).numFrames);

% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    % -- 1st synthetic aperture acquisition for full frame.
    Receive(2*i-1).Apod(1:Resource.Parameters.numRcvChannels) = 1.0;
    Receive(2*i-1).framenum = i;
    Receive(2*i-1).acqNum = 1;
    Receive(2*i-1).callMediaFunc = 1;
    % -- 2nd synthetic aperture acquisition for full frame.
    Receive(2*i).Apod((Resource.Parameters.numRcvChannels+1):Trans.numelements) = 1.0;
    Receive(2*i).framenum = i;
    Receive(2*i).acqNum = 2;   % two acquisitions per frame
end

% Specify Recon structure arrays.
Recon = struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame', -1, ...     % use most recently transferred frame
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...  % auto-increment ImageBuffer each recon
               'RINums', 1:2);

% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 'replaceIQ', ...
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'scaleFactor',0.2,...
                   'regionnum', 1),1,2);
ReconInfo(2).mode = 'accumIQ_replaceIntensity'; % accumulate IQ and detect
ReconInfo(2).rcvnum = 2;

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
SeqControl(1).command = 'jump'; % Jump back to start.
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between frames
SeqControl(2).argument = 20000;  % 20000usec = 20msec (~ 50 fps)
SeqControl(3).command = 'returnToMatlab';
nsc = 4; % nsc is count of SeqControl objects

n = 1; % n is count of Events

% Acquire all frames defined in RcvBuffer
for i = 1:Resource.RcvBuffer(1).numFrames

    Event(n).info = 'First half aperture.';
    Event(n).tx = 1;
    Event(n).rcv = 2*i-1;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = 'Second half aperture.';
    Event(n).tx = 1;
    Event(n).rcv = 2*i;
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
    if floor(i/5) == i/5     % Exit to Matlab every 5th frame
        Event(n).seqControl = 3;
    end
    n = n+1;
end

Event(n).info = 'Jump back to first event';
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
save('MatFiles/L7-4SingleBeam_LE');
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

evalin('base','PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));');
evalin('base','PData(1).Region = computeRegions(PData(1));');
evalin('base','TX(1).TXPD = computeTXPD(TX(1),PData);');
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
Control.Parameters = {'PData','InterBuffer','ImageBuffer','DisplayWindow','TX','Receive','TGC','Recon'};
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
%MaxBLCallback
