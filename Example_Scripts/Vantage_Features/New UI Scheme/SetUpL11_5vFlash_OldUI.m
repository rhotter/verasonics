% File name: SetUpL11_5vFlash_OldUI.m - Modified copy of L11_5vFlash to
% demonstrate the new approach added in 4.3.0 (VTS-1480) for defining GUI
% controls, callbacks and external function handles.  This script
% demonstrates the original "text2cell" approach; refer to companion script
% "SetUpL11_5vFlash_newUI.m" for an example of the same script using the
% new approach.
% 
% To demonstrate the differences between the new and old approaches, a new
% GUI control has been added, that simply changes the transmit burst
% duration. An external processing function has also been added, that
% records and reports how long the script has been running.
%
% This version of the script demonstrates the original "text2cell" scheme
% for defining GUI callbacks and external prcessing functions.  This
% original scheme is fully compatible with the 4.3.0 software, but since
% the external processing function defines a persistent variable this
% script would trigger a Matlab error when used with Matlab 2019b, but
% would function properly with Matlab 2018a, 2018b, and 2019a. Versonics
% reccommends redefining the external processing and GUI callback functions
% as genuine functions rather than preserving the text2cell approach.  The
% companion script "SetUpL11_5vFlash_NewUI.m" demonstrates the changes
% needed to go to the new scheme.  This will provide a much more robust
% solution as well as eliminating the Matlab compatibility problems.
%
% Version 1.0 | 2020-04-01 
% created for 4.3.0 release
% $Author: Ken Linkhart, Dr. Daniel Rohrbach
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent
% and Trademark Office.

clear all
P.startDepth = 5;   % Acquisition depth in wavelengths
P.endDepth = 192;   % This should preferrably be a multiple of 128 samples.

% Define system parameters.
Resource.Parameters.numTransmit = 128;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;    % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.fakeScanhead = 1;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Specify Trans structure array.
Trans.name = 'L11-5v';
Trans.units = 'mm'; % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);
Trans.maxHighVoltage = 50;  % much restricted than specified on data sheet


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
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 100;       % 100 frames used for RF cineloop.
Resource.InterBuffer(1).numFrames = 1;  % one intermediate buffer needed.
Resource.ImageBuffer(1).numFrames = 10;
Resource.DisplayWindow(1).Title = 'L11-5vFlash Old UI';
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
TX(1).focus = 0;
TX(1).Steer = [0.0,0.0];       % theta, alpha = 0.
TX(1).Apod = ones(1,Trans.numelements);
TX(1).Delay = computeTXDelays(TX(1));

% Specify TGC Waveform structure.
TGC.CntrlPts = [0,300,444,552,606,747,870,920];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Receive structure arrays -
%   endDepth - add additional acquisition depth to account for some channels
%              having longer path lengths.
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength,...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 1),1,Resource.RcvBuffer(1).numFrames);
% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    % -- Acquisition for full frame.
    Receive(i).framenum = i;
end

% Specify Recon structure arrays.
Recon = struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame', -1, ...     % use most recently transferred frame
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...  % auto-increment ImageBuffer each recon
               'RINums', 1);

% Define ReconInfo structures.
ReconInfo = struct('mode', 'replaceIntensity', ...
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 1);

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

% Process structure for calling the external processing function
Process(2).classname = 'External';
Process(2).method = 'runTimeMon';
Process(2).Parameters = {'srcbuffer','none',...  % name of buffer to process.
                         'dstbuffer','none'};

% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; % jump back to start.
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between frames
SeqControl(2).argument = 3000;  % 3 msec
SeqControl(3).command = 'returnToMatlab';
nsc = 4; % nsc is count of SeqControl objects

n = 1; % n is count of Events

% Acquire all frames defined in RcvBuffer
for i = 1:Resource.RcvBuffer(1).numFrames

    Event(n).info = 'Full aperture.';
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
    if floor(i/5) == i/5     % Exit to Matlab every 5th frame
        Event(n).seqControl = 3;
    else
        Event(n).seqControl = 0;
    end
    n = n+1;
    
    Event(n).info = 'call external processing function';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 2;
    Event(n).seqControl = 0;
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

% New UI Control for GUI Demo; modify TX burst duration
% In this script, the new Control is defined using the original UI
% structure and text2cell approach from releases prior to 4.3.0
UI(3).Control = {'UserB3', ...
                'Style','VsSlider', ...
                'Label','TX Duration cycles', ...
                'SliderMinMaxVal',[0 20 1], ...
                'SliderStep',[0.1,0.25], ...
                'ValueFormat','%3.0f'};
UI(3).Callback = text2cell('%TXdurChangeCallback');


% Create the External Processing function, using original text2cell approach.
EF(1).Function = text2cell('%EF#1');

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 5;

% Save all the structures to a .mat file.
save('MatFiles/L11-5vFlashOldUI');

%% Automatic VSX Execution:
% Uncomment the following line to automatically run VSX every time you run
% this SetUp script (note that if VSX finds the variable 'filename' in the
% Matlab workspace, it will load and run that file instead of prompting the
% user for the file to be used):

% filename = 'L11-5vFlashOldUI';  VSX;

return


% **** Callback function definitions. ****

% - Sensitivity cutoff change callback
%SensCutoffCallback
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

% - Range change callback
%RangeChangeCallback
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

% - TXdur change callback
%TXdurChangeCallback
TW = evalin('base', 'TW');
TW(1).Parameters(3) = 2*round(UIValue);
assignin('base','TW',TW);
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'TW'};
assignin('base','Control', Control);
return
%TXdurChangeCallback


%EF#1 External Processing function using original text2cell format
runTimeMon(varargin)
persistent startTime frameCount dispTic
% Initialize startTime if it does not exist.
if isempty(startTime)
    startTime = tic;
    frameCount = 1;
    dispTic = 1;
else
    % measure elapsed time so far
    % runTime is the total time in seconds the script has been running
    runTime = toc(startTime);
    assignin('base', 'runTime', runTime);
    % frameCount is total number of flash image frames that have been acquired
    frameCount = frameCount + 1;
    assignin('base', 'frameCount', frameCount);
    % notify at command prompt every 10 seconds
    if runTime > 10 * dispTic
        disp(['Script run time ', num2str(runTime), ' seconds.']);
        dispTic = dispTic + 1;
    end
end
return
%EF#1