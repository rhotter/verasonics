% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: SetUpL11_5vFlashAngles_Import_SequenceStorage_frames.m
% Example of plane wave imaging with steering angle transmits and saving RF data in realtime.
% This script shows how to use external function to save the RF data of all frames in a receive buffer
% in realtime. Here, only non-zeros channel data will be saved to reduce the
% saving time. The default frame limit is set to 100 RcvBuffers. Please change
% to a bigger number for more buffers.
%
% Description:
%   Sequence programming file for L11-5v Linear array, using plane wave
%   transmits with multiple steering angles. All 128 transmit and receive
%   channels are active for each acquisition. Processing is asynchronous
%   with respect to acquisition.
%
% In order to save each frame in a correct order, a "sync"
% command is required for the hardware to wait for the software to finish
% saving RF data, image reconstruction and image display.
% Therefore, a warning message "timeToNextAcq duration too short" might occur
% if the interval (SeqControl(3).argument in this script) between two frames is
% not long enough.
%
% Last update:
% 05/06/2020 - Update to SW 4.3 format for new user UIControls and External function definitions (VTS 1691).
%   More info:(?/Example_Scripts/Vantage_Features/New UI Scheme/SetUpL11_5vFlash_NewUI)
% 05/23/2016 - test with SW 3.0.7

clear all

P.startDepth = 5;   % Acquisition depth in wavelengths
P.endDepth = 192;   % This should preferably be a multiple of 128 samples.

na = 7;      % Set na = number of angles.
if (na > 1)
    dtheta = (36*pi/180)/(na-1);
    P.startAngle = -36*pi/180/2;
else
    dtheta = 0;
    P.startAngle=0;
end % set dtheta to range over +/- 18 degrees.

% Define system parameters.
Resource.Parameters.numTransmit    = 128;   % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;   % number of receive channels.
Resource.Parameters.speedOfSound   = 1540;  % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose        = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode   = 1;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Specify Trans structure array.
Trans.name = 'L11-5v';
Trans.units = 'wavelengths'; % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);  % L11-5v transducer is 'known' transducer so we can use computeTrans.
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

Resource.RcvBuffer(1).datatype      = 'int16';
Resource.RcvBuffer(1).rowsPerFrame  = na*4096; % this size allows for maximum range
Resource.RcvBuffer(1).colsPerFrame  = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames     = 30;    % 30 frames stored in RcvBuffer.

Resource.InterBuffer(1).numFrames   = 1;   % one intermediate buffer needed.
Resource.ImageBuffer(1).numFrames   = 10;
Resource.DisplayWindow(1).Title     = 'L11-5vFlashAngles';
Resource.DisplayWindow(1).pdelta    = 0.35;
ScrnSize = get(0,'ScreenSize');
DwWidth  = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt   = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type          = 'Verasonics';
Resource.DisplayWindow(1).AxesUnits     = 'mm';
Resource.DisplayWindow(1).Colormap      = gray(256);

% Specify Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,.67,2,1];

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'Apod', kaiser(Resource.Parameters.numTransmit,1)', ...
                   'focus', 0.0, ...
                   'Steer', [0.0,0.0], ...
                   'Delay', zeros(1,Trans.numelements)), 1, na);
% - Set event specific TX attributes.
if fix(na/2) == na/2       % if na even
    P.startAngle = (-(fix(na/2) - 1) - 0.5)*dtheta;
else
    P.startAngle = -fix(na/2)*dtheta;
end
for n = 1:na   % na transmit events
    TX(n).Steer = [(P.startAngle+(n-1)*dtheta),0.0];
    TX(n).Delay = computeTXDelays(TX(n));
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [0,141,275,404,510,603,702,782];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Receive structure arrays.
% - We need na Receives for every frame.
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
Receive      = repmat(struct(   'Apod',         ones(1,Trans.numelements), ...
                                'startDepth',   P.startDepth, ...
                                'endDepth',     maxAcqLength,...
                                'TGC',          1, ...
                                'bufnum',       1, ...
                                'framenum',     1, ...
                                'acqNum',       1, ...
                                'sampleMode',   'NS200BW', ...
                                'mode',          0, ...
                                'callMediaFunc', 0), 1, (na*Resource.RcvBuffer(1).numFrames) );

% - Set event specific Receive attributes for each frame.


for i = 1:Resource.RcvBuffer.numFrames
    Receive(na*(i-1)+1).callMediaFunc = 1;
    for j = 1:na
        Receive(na*(i-1)+j).bufnum   = 1;
        Receive(na*(i-1)+j).framenum = i;
        Receive(na*(i-1)+j).acqNum   = j;
    end
end


% Specify Recon structure arrays.
% - We need one Recon structures which will be used for each frame.
Recon = repmat(struct('senscutoff', 0.6, ...
               'newFrameTimeout', 10000, ...
               'pdatanum',      1, ...
               'rcvBufFrame',  -1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...
               'RINums', 1:na),1, Resource.RcvBuffer(1).numFrames);

% we create a recon for each frame
for N = 1:Resource.RcvBuffer(1).numFrames
    Recon(N).rcvBufFrame = N;
    Recon(N).RINums      = (1:na)+(N-1)*na;
end

% Define ReconInfo structures.
% We need na ReconInfo structures for na steering angles.
ReconInfo = repmat(struct('mode', 'accumIQ', ...  % default is to accumulate IQ data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 1), 1, na*Resource.RcvBuffer(1).numFrames);

% - Set specific ReconInfo attributes.
for N = 1:Resource.RcvBuffer(1).numFrames
    if na>1
        ReconInfo(1+na*(N-1)).mode = 'replaceIQ'; % replace IQ data
        for j = 1:na  % For each row in the column
            ReconInfo(na*(N-1)+j).txnum = j;
            ReconInfo(na*(N-1)+j).rcvnum = na*(N-1)+j;
        end
        ReconInfo(na*N).mode = 'accumIQ_replaceIntensity'; % accum and detect
    else
        ReconInfo(N).mode = 'replaceIntensity';
    end
end


% Specify Process structure array.
pers = 20;
Process(1).classname  = 'Image';
Process(1).method     = 'imageDisplay';
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

% obtain a local file path
str = which('VSX');
% save files in VSX folder under UserData
[fp, ~, ~,]   = fileparts(str);
storageFolder = fullfile( fullfile(fp, 'UserData'), 'StorageTest');
% make sure UserData folder exist
if( ~exist(storageFolder, 'dir'))
    error('No Storage Files available');
end
% The following will only work if start count was 0 for export
ls = vsv.file.FileTools.listFilesRegexp(storageFolder, '[0-9]{3}\.vrs');
numMaxFile = length(ls)-1;

% for N = 2:Resource.RcvBuffer.numFrames+1
    % Save Data Process
Process(2).classname  = 'Storage';
Process(2).method     = 'load'; % specify import for importing data
Process(2).Parameters = { ...
    'dstbuffer',    'receive',... % name of buffer to process.
    'dstframenum',  -1,       ... % process the last frame in a RcvBuffer
    'dstbufnum',    1,        ... % destination buffer
    'filename',    'TestFile',... % base file name needs to match th files in the folder
    'folder',       storageFolder, ... % folder where to look for frames
    'numdigits',    3, ...        % specify the number of digits to build the file name number e.g. if = 3; TestFile001.vrs, TestFile002.vrs, ...TestFile999.vrs
    'enabled',      1, ...        % disable the saving will add a button that allows to start the saving at a particular point, if set to one, saving will start with script first acquisitions
    'maxfiles',     numMaxFile, ...       % maximum number of files to import
    'repeat',       1, ...        % specify whether to repeat the file import when maximum files are reached
    'startcount',   0, ...        % specify the first file count to import
};

% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; % jump back to start
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(2).argument = 200;  % 200 usec
SeqControl(3).command = 'timeToNextAcq';  % time between frames
SeqControl(3).argument = 20000 - (na-1)*200;  % 25 msec
SeqControl(4).command = 'returnToMatlab';
SeqControl(5).command = 'sync';
SeqControl(5).argument = 1e6;  % 1s

nsc = 6; % nsc is count of SeqControl objects

% Specify Event structure arrays.
n = 1;
for N = 1:Resource.RcvBuffer.numFrames

    % Instead of simulating try to import RF data from disk
    Event(n).info = 'Import RF';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 2;
    Event(n).seqControl = 0;
    n=n+1;

    Event(n).info = 'Sync';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = 5;
    n = n+1;

    Event(n).info = 'recon and process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 1;
    Event(n).seqControl = 0;
    if floor(i/5) == i/5     % Exit to Matlab every 5th frame
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

import vsv.seq.uicontrol.VsSliderControl
import vsv.seq.uicontrol.VsButtonControl

% - Sensitivity Cutoff
UI(1).Control =  VsSliderControl('LocationCode','UserB7',...
    'Label','Sens. Cutoff',...
    'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
    'SliderStep',[0.025,0.1],'ValueFormat','%1.3f',...
    'Callback',@SensCutoffCallback);

% - Range Change
MinMaxVal = [64,300,P.endDepth]; % default unit is wavelength
AxesUnit = 'wls';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
        AxesUnit = 'mm';
        MinMaxVal = MinMaxVal * (Resource.Parameters.speedOfSound/1000/Trans.frequency);
    end
end
UI(2).Control = VsSliderControl('LocationCode','UserA1',...
    'Label',['Range (',AxesUnit,')'],...
    'SliderMinMaxVal',MinMaxVal,'SliderStep',[0.1,0.2],...
    'ValueFormat','%3.0f','Callback',@RangeChangeCallback);

UI(3).Control = VsButtonControl('LocationCode','UserB2',...
    'Label','Start','Callback',@StartCallback);


% External function definitions

import vsv.seq.function.ExFunctionDef

EF(1).Function = vsv.seq.function.ExFunctionDef('saveRFallFrames', @saveRFallFrames);


% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 5;

% Save all the structures to a .mat file.
save('MatFiles/L11-5vFlashAngles_saveRFallFrame');

%% Automatic VSX Execution:
% Uncomment the following line to automatically run VSX every time you run
% this SetUp script (note that if VSX finds the variable 'filename' in the
% Matlab workspace, it will load and run that file instead of prompting the
% user for the file to be used):

% filename = 'L11-5vFlashAngles_saveRFallFrame';  VSX;


%% **** Callback routines used by External function definition (EF) ****

function saveRFallFrames(RData)

    persistent RcvBufferNum RFfilename

    % file size can be reduced by elimating all zeros
    TXApod = evalin('base','TX.Apod');
    endSample = evalin('base','Receive(end).endSample');

    RcvNumLimit = 100;
    numLength = max(ceil(log10(abs(RcvNumLimit))),1)+1;

    if isempty(RcvBufferNum)
        RcvBufferNum = 1;
        RFfilename = ['RFdata_',datestr(now,'dd-mmmm-yyyy_HH-MM-SS')];
    else
        RcvBufferNum = RcvBufferNum + 1;
        if RcvBufferNum > RcvNumLimit
            RcvBufferNum = 1;
        end % set a limit for testing
    end

    fname = ['RcvBuffer_',num2str(RcvBufferNum,['%0',num2str(numLength),'.0f'])];
    newRData = RData(1:endSample, find(TXApod),:);
    eval([fname,' = newRData;']);

    tic
    if isequal(RcvBufferNum,1)
        save(RFfilename,fname,'-v6');
    else
        save(RFfilename,fname,'-v6','-append');
    end
    fprintf('saving time for buffer %g = %g s \n',RcvBufferNum,toc);
end


%% **** Callback routines used by UIControls (UI) ****

function StartCallback(~, ~, ~)
    if ~evalin('base','exist(''Process'',''var'')')
        return;
    end

    if evalin('base','exist(''Control'',''var'')')
        Cntrl = evalin('base','Control');
    else
        Cntrl = struct('Command',[],'Parameters',[]);
    end

    if isempty(Cntrl(1).Command)
        n=1;
    else
        n=length(Cntrl)+1;
    end

    Cntrl(n).Command    = 'set&Run';
    Cntrl(n).Parameters = {'Process', 2,...
                                 'enabled', 1};

    assignin('base','Control',Cntrl);
end

function SensCutoffCallback(~, ~, UIValue)
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

function RangeChangeCallback(hObject, ~, UIValue)
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
end