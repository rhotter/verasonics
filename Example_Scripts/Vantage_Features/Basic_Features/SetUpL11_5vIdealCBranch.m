% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use.
%
% File name: SetUpL11_5vIdealCBranch.m - Example of use of the "cBranch"
% sequence control command.
%
% Description:
%   This is the L11-5videal example script, modified to add a cBranch jump
%   to an alternate acquisition sequence which acquires an L11-5vFlash
%   image.  A gui button is added to trigger the cBranch, one Flash image
%   is acquired, processed, and displayed and then the script returns to
%   the "ideal" imaging sequence.
%
% Last update:
% 05/06/2020 - Update to SW 4.3 format for new user UIControls and External function definitions (VTS 1691).
%   More info:(?/Example_Scripts/Vantage_Features/New UI Scheme/SetUpL11_5vFlash_NewUI)
% Feb. 15, 2017 created for use with 3.2.1 release as cBranch example

clear all

P.startDepth = 5;   % Acquisition depth in wavelengths
P.endDepth = 192;   % This should preferrably be a multiple of 128 samples.

m = 128;  % number of synthetic aperture acquisitions per frame.

% Specify system parameters.
Resource.Parameters.numTransmit = 128;  % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;  % number of receive channels.
Resource.Parameters.speedOfSound = 1540;
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.


% Specify Trans structure array.
Trans.name = 'L11-5v';
Trans.units = 'mm'; % required in Gen3 to prevent default to mm units
Trans = computeTrans(Trans);  % L11-5v transducer is 'known' transducer so we can use computeTrans.
Trans.maxHighVoltage = 50;  % set maximum high voltage limit for pulser supply.

% Specify PData structure array.
PData(1).PDelta = [Trans.spacing, 0, 0.5];  % x, y, z pdeltas
PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));
PData(1).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(3) = 1;      % single image page
PData(1).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of upper lft crnr
% No PData.Region specified, so a default Region for the entire PData array will be created by computeRegions.

% Specify Media object.
pt1;
Media.attenuation = -0.5;
% Media.function = 'movePoints';

% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 2560*m;
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 10;     % 10 frames for RF data.
% Buffer 2 will be used by the Flash imaging sequence, when the cBranch is taken
Resource.RcvBuffer(2).datatype = 'int16';
Resource.RcvBuffer(2).rowsPerFrame = 2560; % single flash acquisition
Resource.RcvBuffer(2).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(2).numFrames = 1;     % single frame for cBranch flash image.
Resource.InterBuffer(1).numFrames = 1;  % one intermediate buffer needed.
Resource.InterBuffer(2).numFrames = 1;  % one intermediate buffer needed.
Resource.ImageBuffer(1).numFrames = 10;
Resource.ImageBuffer(2).numFrames = 10;
Resource.DisplayWindow(1).Title = 'L11-5vIdeal';
Resource.DisplayWindow(1).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [50,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).numFrames = 20;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow.Colormap = gray(256);

% display window 2 is for the Flash image, will only update when the
% cBranch is triggered by the GUI pushbutton
Resource.DisplayWindow(2).Title = 'L11-5vFlash';
Resource.DisplayWindow(2).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(2).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(2).pdelta);
Resource.DisplayWindow(2).Position = [500,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(2).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(2).Type = 'Verasonics';
Resource.DisplayWindow(2).numFrames = 1;
Resource.DisplayWindow(2).AxesUnits = 'mm';
Resource.DisplayWindow(2).Colormap = gray(256);

% Specify Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,.67,2,1];

% Specify m TX structure arrays. Transmit on element n in the array for event n.
TX = repmat(struct('waveform', 1, ...
                   'Origin', zeros(1,3), ...
                   'focus', 0, ...
                   'Steer', [0.0,0.0], ...
                   'Apod', zeros(1,Trans.numelements), ...
                   'Delay', zeros(1,Trans.numelements)), 1, m);

scaleToWvl = 1;
if strcmp(Trans.units, 'mm')
    scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);
end
% - Set event specific TX attributes.
for n = 1:m   % m transmit events
    % Set transmit Origins to positions of elements.
    TX(n).Origin = scaleToWvl*Trans.ElementPos(n,1:3);
    % Set transmit Apodization so that only one element is active.
    TX(n).Apod(n) = 1.0;    % Only one active transmitter for each TX.
end



% Specify TX structure array for cBranch Flash image.
cBtx = m+1;
TX(cBtx).waveform = 1;            % use 1st TW structure.
TX(cBtx).Origin = [0.0,0.0,0.0];  % flash transmit origin at (0,0,0).
TX(cBtx).focus = 0;
TX(cBtx).Steer = [0.0,0.0];       % theta, alpha = 0.
TX(cBtx).Apod = ones(1,Trans.numelements);
TX(cBtx).Delay = computeTXDelays(TX(1));

% Specify TGC Waveform structure.
TGC.CntrlPts = [189,314,457,698,770,911,948,976];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Receive structure arrays.
% - We need m Receive structures for each frame.
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode','BS100BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, m*Resource.RcvBuffer(1).numFrames+1);
% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    Receive(m*(i-1)+1).callMediaFunc = 1;
    for j = 1:m
        Receive(m*(i-1)+j).framenum = i;
        Receive(m*(i-1)+j).acqNum = j;
    end
end

% create one more Receive for the Flash image acquisition on the cBranch,
% using receive buffer 2 frame 1
cbRcv = m*Resource.RcvBuffer(1).numFrames + 1;
Receive(cbRcv).bufnum = 2;

% Specify Recon structure arrays.
% - We need one Recon structure.  Each frame will use
%   m ReconInfo structures, since we are using m
%   synthetic aperture acquisitions.
Recon = struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame', -1, ...     % use most recently transferred frame
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...  % auto-increment ImageBuffer each recon
               'RINums', 1:m);

% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 'accumIQ', ...  % accumulate IQ data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 1), 1, m);
% - Set specific ReconInfo attributes.
if m>1
    ReconInfo(1).mode = 'replaceIQ';  % replace IQ data
    for j = 1:m  % For each row in the column
        ReconInfo(j).txnum = j;
        ReconInfo(j).rcvnum = j;
    end
    ReconInfo(m).mode = 'accumIQ_replaceIntensity'; % accum & detect
else
    ReconInfo(1).mode = 'replaceIntensity';
end



% Specify cBranch processing Recon structure, using buffer #2 and
% explicitly processing only frame 1 in the receive buffer.  Note that the
% Recon processing will automatically wait for DMA transfer of receive
% buffer 2 frame 1 to complete.
Recon(2) = struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame', 1, ...     % use frame 1 (only one frame exists in buffer #2)
               'IntBufDest', [2,1], ...
               'ImgBufDest', [2,1], ...  % auto-increment ImageBuffer each recon
               'RINums', m+1);

% Define cBranch processing ReconInfo structure (only one Flash acquisition).
ReconInfo(m+1).mode = 'replaceIntensity';
ReconInfo(m+1).txnum = cBtx;
ReconInfo(m+1).rcvnum = cbRcv;
ReconInfo(m+1).regionnum = 1;

% Specify Process structure array for image display.
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


% Specify second Process structure array for the cBranch Flash image.
Process(2) = Process(1);
Process(2).Parameters = {'imgbufnum',2,...   % number of buffer to process.
                         'framenum',1,...   % (-1 => lastFrame)
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
                         'displayWindow',2};

% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; % jump back to start.
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(2).argument = 200;  % 200 usec
SeqControl(3).command = 'returnToMatlab';

SeqControl(4).command = 'cBranch';
SeqControl(4).condition = 'bFlag';

SeqControl(5).command = 'sync';

nsc = 6; % nsc is count of SeqControl objects


% Specify Event structure array.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    for j = 1:m                 % Acquire frames
        Event(n).info = 'Acquisition.';
        Event(n).tx = j;
        Event(n).rcv = m*(i-1)+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;
    end
    Event(n-1).seqControl = [2, nsc]; % modify last acquisition Event's seqControl
      SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
      nsc = nsc+1;

    Event(n).info = 'recon and process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 1;
    Event(n).seqControl = 3;
    n = n+1;

    Event(n).info = 'cBranch conditionaljump';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = 4;
    n = n+1;
end


Event(n).info = 'Jump back';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 1;
n = n+1;

% This is the start of the cBranch Flash image acquisition and processing
% events

SeqControl(4).argument = n;  % branch Event for the cBranch command

% The hardware and software sequencers are running asynchronously, so we
% don't know exactly when each of them will respond to the cBranch command
% and jump to this point.  The first thing we do here is a sync command
% so whoever gets here first will have to wait for the other and then they
% will both proceed together to the events following the sync
Event(n).info = 'sync after branch';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 5;
n = n+1;

Event(n).info = 'Acquisition event, single Flash acquisition with transfer to host';
Event(n).tx = cBtx;
Event(n).rcv = cbRcv;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = nsc;
      SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
      nsc = nsc+1;
n = n+1;

Event(n).info = 'Recon and display for flash image';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 2;
Event(n).process = 2;
Event(n).seqControl = 0;
n = n+1;

% after taking the cBranch to do the acquisition and processing of the
% Flash image, we jump back to event #1 to resume the 'ideal' imaging
% sequence
Event(n).info = 'Jump back';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 1;
n = n+1;


% User specified UI Control Elements

import vsv.seq.uicontrol.VsSliderControl
import vsv.seq.uicontrol.VsToggleButtonControl

% - Sensitivity Cutoff
UI(1).Control = VsSliderControl('LocationCode','UserB7',...
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
    'Label',['Range (',AxesUnit,')'],'SliderMinMaxVal',MinMaxVal,...
    'SliderStep',[0.1,0.2],'ValueFormat','%3.0f',...
    'Callback',@RangeChangeCallback);

% - cBranch Gui Button
UI(3).Control = VsToggleButtonControl('LocationCode','UserB4',...
    'Label','cBranch','Callback',@CBranchButtonCallback);

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 1;

% Save all the structures to a .mat file.
save('MatFiles/L11-5vIdealCBranch');

%% Automatic VSX Execution:
% Uncomment the following line to automatically run VSX every time you run
% this SetUp script (note that if VSX finds the variable 'filename' in the
% Matlab workspace, it will load and run that file instead of prompting the
% user for the file to be used):

% filename = 'L11-5vIdealCBranch';  VSX;


%% **** Callback routines used by UIControls (UI) ****

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

function CBranchButtonCallback(~, ~, ~)
    Control = evalin('base','Control');
    Control.Command = 'setBFlag';
    assignin('base','Control', Control);
end
