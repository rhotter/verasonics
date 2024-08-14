% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: SetUpL11_5vFlashPushP5Mon.m - Example of inserting a Push
% transmit burst within the L11-5vFlash imaging script.  GUI commands allow
% setting the Profile 5 Monitor limits for testing that feature.  Setting
% the limits to the same value disables the monitoring of the Profile 5
% voltage level.
%
% Description:
% Generate .mat Sequence Object file for L11-5v Linear array flash imaging,
% with a transition to TPC profile 5 for two Push transmit-only events
% using the Extended Transmit Option.
% 128 transmit channels and 128 receive channels are used.
%
% Push % Comment lines with this prefix identify the changes added to allow
% the script to exercise a TPC Profile 5 Push transmit-only event.  Other
% than these changes, this is a typical flash imaging example script.
%
% NOTE: This script should not be run with an actual L11-5v transducer
% connected; the Push transmit burst could be destructive to the
% transducer!!  It should be run with no transducer connected at all.  The
% purpose of this script is only to demonstrate how to define a script for Push
% transmit, not to actually do anything useful with an actual transducer.
%
% Revision history:
% 05/06/2020 - Update to SW 4.3 format for new user UIControls and External function definitions (VTS 1691).
%   More info:(?/Example_Scripts/Vantage_Features/New UI Scheme/SetUpL11_5vFlash_NewUI)
% Jan 2020 4.2.0 issue VTS-1583 P5 monitor limits control bugs
% July 9 2019 4.1.0 verification test support

clear all

P.startDepth = 5;   % Acquisition depth in wavelengths
P.endDepth = 192;   % This should preferrably be a multiple of 128 samples.

% Define system parameters.
Resource.Parameters.numTransmit = 128;      % number of transmit channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.fakeScanhead = 1;  % allows automatic operation with no probe connected

Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Resource.Parameters.SystemLED = {'activeTxAndOrRx', 'activeTxAndOrRxProfile5', 'off', 'starting'};
% Resource.Parameters.ProbeConnectorLED = [1 0 0 0];

% Specify Trans structure array.
Trans.name = 'L11-5v';
Trans.units = 'wavelengths'; % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);  % L11-5v transducer is 'known' transducer so we can use computeTrans.
Trans.maxHighVoltage = 40;  % set maximum high voltage limit for pulser supply.

% Push % Must specify maxHighVoltage explicitly whenever TPC profile 5 is
% to be used
TPC(5).maxHighVoltage = 40;

% Initialize the P5 monitor limit values
P5HVmin = 0;
P5HVmax = 0;


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
Resource.DisplayWindow(1).Title = 'L11-5vFlash';
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


% Push % A separate transmit waveform, TW(2), is defined for the Push
% transmit. This waveform will be used with TPC profile 5 and the Extended
% Burst transmit power supply, allowing high-power bursts of up to several
% milliseconds duration.
bdur = 0.1; % initial burst duration in msec
pushFreq = 4; % profile 5 push burst frequency MHz
numBursts = 1; % number of P5 transmit bursts- set to 1, 2, or 3
TW(2).type = 'parametric';
TW(2).Parameters = [pushFreq,.67,2000*bdur*pushFreq,1];
% to get number of half-cycles multiply bdur times 1000 to convert to usec,
% times pushFreq to get total number of cycles, times 2 to get total number
% of half-cycles

% Specify TX structure array.
TX(1).waveform = 1;            % use 1st TW structure.
TX(1).Origin = [0.0,0.0,0.0];  % flash transmit origin at (0,0,0).
TX(1).focus = 0;
TX(1).Steer = [0.0,0.0];       % theta, alpha = 0.
TX(1).Apod = ones(1,Trans.numelements);
TX(1).Delay = computeTXDelays(TX(1));

% Push % TX(2) is defined for the Push transmit, using the Push transmit
% waveform defined in TW(2).  Transmit focus control and TX.Apod aperture
% control both function the same for a Push burst as for any imaging burst.
TX(2) = TX(1);
TX(2).focus = 200;
TX(2).Apod(1:16) = 0;
TX(2).Delay = computeTXDelays(TX(2));
TX(2).waveform = 2;


% Specify TGC Waveform structure.
TGC.CntrlPts = [0,141,275,404,510,603,702,782];
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

% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; % jump back to start.
SeqControl(1).argument = 3; % Push % jump back to third event since we don't need to re-initialize TPC profile selection.
SeqControl(2).command = 'timeToNextAcq';  % time between frames
SeqControl(2).argument = 3000;  % 3 msec
SeqControl(3).command = 'returnToMatlab';

SeqControl(4).command = 'setTPCProfile';
SeqControl(4).argument = 5;
SeqControl(4).condition = 'next';

SeqControl(5).command = 'setTPCProfile';
SeqControl(5).argument = 1;
SeqControl(5).condition = 'next';

SeqControl(6).command = 'triggerOut';

SeqControl(7).command = 'timeToNextAcq';  % time between frames
SeqControl(7).argument = 40000;  % 40 msec

nsc = 8; % nsc is count of SeqControl objects

n = 1; % n is count of Events

% Push % For any script using multiple TPC profiles, and especially any
% script using Profile 5 for Push transmit, the script must explicitly
% specify an initial profile at the beginning of the script, prior to any
% transmit events.  This is to prevent the script from 'inheriting'
% whatever TPC Profile was in effect when some previous script was
% terminated.
Event(n).info = 'select TPC profile';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = nsc;
n = n+1;
   SeqControl(nsc).command = 'setTPCProfile';
   SeqControl(nsc).argument = 1;
   SeqControl(nsc).condition = 'immediate';
   nsc = nsc + 1;

% Push % After a set TPC profile command has been issued, enough time must
% be allowed for the TPC profile transition to complete.  10 msec is
% usually sufficient for this transition; the actual time required depends
% on how far apart the two voltage settings are.
Event(n).info = 'noop delay';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = nsc;
n = n+1;
   SeqControl(nsc).command = 'noop';
   SeqControl(nsc).argument = 5e5; % 500 msec delay
   nsc = nsc + 1;

% Acquire all frames defined in RcvBuffer
for i = 1:Resource.RcvBuffer(1).numFrames

    Event(n).info = 'Full aperture.';
    Event(n).tx = 1;
    Event(n).rcv = i;
    Event(n).recon = 0;
    Event(n).process = 0;
    if i < Resource.RcvBuffer(1).numFrames
        Event(n).seqControl = [2,nsc];
    else
        Event(n).seqControl = [7,nsc,4];
    end
       SeqControl(nsc).command = 'transferToHost';
       nsc = nsc + 1;
    n = n+1;

    % Push % In the last imaging acquisition event of the frame sequence
    % above, a set TPC Profile command to switch into profile 5 has been
    % added.  (See SeqControl(4) assignment above.)  THe 'next' condition
    % means this transition will not start until the imaging event is
    % completed.  In this case, the 40 msec timeToNextAcq delay
    % will allow plenty of time for the TPC profile transition to take
    % place before the Push event

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
end

% Push % Here we insert the Push transmit-only events, with a trigger out
% command to make it convenient to trigger an oscilloscope for examining
% the actual Push transmit output.  At the end of these events the transition
% back to TPC profile 1 will be started, for the next imaging acquisition
% sequence.  Another 40 msec timeToNextAcq delay will allow time for the
% transition to be completed.
Event(n).info = 'First Push event';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 0;
firstPushEvent = n; % Event index for use by numbursts GUI control
n = n+1;

Event(n).info = 'Second Push event';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 0;
n = n+1;

Event(n).info = 'Push event with trig out, and return to profile 1';
Event(n).tx = 2;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = [6, 5, 7];
n = n+1;

Event(n).info = 'Jump back to first event';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 1;

% User specified UI Control Elements

import vsv.seq.uicontrol.VsSliderControl

% User specified UI Control Elements
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

% P5 HV monitor max voltage limit
UI(3).Control = VsSliderControl('LocationCode','UserB5',...
    'Label','P5 Max Limit',...
    'SliderMinMaxVal',[0, Trans.maxHighVoltage, P5HVmax],...
    'SliderStep',[0.01,0.05],'ValueFormat','%3.1f',...
    'Callback',@P5MaxCallback);

% P5 HV monitor min voltage limit
UI(4).Control = VsSliderControl('LocationCode','UserB4',...
    'Label','P5 Min Limit',...
    'SliderMinMaxVal',[0, Trans.maxHighVoltage, P5HVmin],...
    'SliderStep',[0.01,0.05],'ValueFormat','%3.1f',...
    'Callback',@P5MinCallback);

% burst duration control
UI(5).Control = VsSliderControl('LocationCode','UserC3',...
    'Label','P5 Bdur msec','SliderMinMaxVal',[0, 1, bdur],...
    'SliderStep',[0.1,0.2],'ValueFormat','%3.2f',...
    'Callback',@BdurCallback);

% Num P5 bursts control
UI(6).Control = VsSliderControl('LocationCode','UserC2',...
    'Label','num P5 Bursts','SliderMinMaxVal',[1, 3, numBursts],...
    'SliderStep',[0.5,0.5],'ValueFormat','%3.0f',...
    'Callback',@numBurstCallback);


% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 5;

% Save all the structures to a .mat file.
save('MatFiles/L11-5vFlashPushP5Mon');

%% Automatic VSX Execution:
% Uncomment the following line to automatically run VSX every time you run
% this SetUp script (note that if VSX finds the variable 'filename' in the
% Matlab workspace, it will load and run that file instead of prompting the
% user for the file to be used):

% filename = 'L11-5vFlashPushP5Mon';  VSX;


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

function P5MaxCallback(~, ~, UIValue)
    P5HVmax = UIValue;
    assignin('base', 'P5HVmax', P5HVmax);
    P5HVmin = evalin('base', 'P5HVmin');
    % call the function to set new limit values
    rc = setProfile5VoltageLimits(P5HVmin, P5HVmax);
    if ~strcmp(rc, 'Success')
        fprintf(2, ['Error from setProfile5VoltageLimits: ', rc, ' \n'])
    end
end

function P5MinCallback(~, ~, UIValue)
    P5HVmin = UIValue;
    assignin('base', 'P5HVmin', P5HVmin);
    P5HVmax = evalin('base', 'P5HVmax');
    % call the function to set new limit values
    rc = setProfile5VoltageLimits(P5HVmin, P5HVmax);
    if ~strcmp(rc, 'Success')
        fprintf(2, ['Error from setProfile5VoltageLimits: ', rc, ' \n'])
    end
end

function BdurCallback(~, ~, UIValue)
    bdur = UIValue;
    assignin('base', 'bdur', bdur);  % Burst duration in milliseconds
    TW = evalin('base', 'TW');
    TW(2).Parameters(3) = 2000*bdur*TW(2).Parameters(1); % scale bdur * freq to number of halfcycles
    assignin('base', 'TW', TW);
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'TW'};
    assignin('base','Control', Control);
end

function numBurstCallback(~, ~, UIValue)
    numBursts = round(UIValue);
    assignin('base', 'numBursts', numBursts);
    Event = evalin('base', 'Event');
    firstPushEvent = evalin('base', 'firstPushEvent');
    switch numBursts
        case 1
            Event(firstPushEvent).tx = 0; % disable first push
            Event(firstPushEvent + 1).tx = 0; % disable second push
        case 2
            Event(firstPushEvent).tx = 0; % enable first push
            Event(firstPushEvent + 1).tx = 2; % disable second push
        case 3
            Event(firstPushEvent).tx = 2; % enable first push
            Event(firstPushEvent + 1).tx = 2; % enable second push
    end
    assignin('base', 'Event', Event);
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'Event'};
    assignin('base','Control', Control);
end