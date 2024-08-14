% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: SetUpL11-5vThermistorDemo.m - L11-5vFlash script modified to
% demonstrate thermal monitoring using the probe thermistor temperature
% sense feature of the system.
% Description:
%   An external processing function has been added to the script, to read
%   the thermistor values from the hardware system and display them.  Note
%   the L11-5v does not have thermistors, so the monitor will report an open
%   circuit if this script is run with the L11-5v.  It can be run in
%   "fakeScanhead" mode with nothing connected, and then you can connect a
%   resistor or short circuit to the thermistor pins at the connector, and
%   see the results from the monitor function.
%
%%%%%% THERMISTOR DEMO %%%%%%
% search for this comment string within the file, to find the code used to
% program the thermistor monitor function, and also to actually read the
% temperature sense thermistor A/D data from the HW system and process it:
%%%%%% THERMISTOR DEMO %%%%%%
%
% Last update:
% 05/06/2020 - Update to SW 4.3 format for new user UIControls and External function definitions (VTS 1691).
%   More info:(?/Example_Scripts/Vantage_Features/New UI Scheme/SetUpL11_5vFlash_NewUI)
% 4/14/2017 created as thermistor demo for 3.2.2 release

clear all

P.startDepth = 5;   % Acquisition depth in wavelengths
P.endDepth = 192;   % This should preferrably be a multiple of 128 samples.

% Define system parameters.
Resource.Parameters.numTransmit = 128;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;    % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

%%%%%% THERMISTOR DEMO %%%%%%
Resource.Parameters.fakeScanhead = 1; % allow system to run automatically, with nothing connected

% Uncomment following line to select connector 2; otherwise system will default to connector 1
% Resource.Parameters.Connector = 2;


% Set ProbeThermistor to all zeros (this duplicates what VSX would have
% done by default, but it allows us to set individual values as desired in
% the lines that follow)
Resource.Parameters.ProbeThermistor = ...
    repmat(struct('enable', 0, ...
                  'threshold', 0, ...
                  'reportOverThreshold', 0), 1, 2);

% the following three lines will enable monitoring on thermistor 1, with an
% automatic shutdown if the A/D reading is ever less than 40 (equivalent to
% approximately 20 Ohms), so if you jumper pin cc4 to thermistor return at
% cc8 the system will exit with an over-temperature error
Resource.Parameters.ProbeThermistor(1).enable = 1;
Resource.Parameters.ProbeThermistor(1).threshold = 40;
Resource.Parameters.ProbeThermistor(1).reportOverThreshold = 0;
% 0 directs the system to shut down if the thermistor value goes under the
% threshold value; to reverse this and shut down for values over the threshold,
% use reportOverThreshold = 1
%%%%%% THERMISTOR DEMO %%%%%%



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
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 4096;   % this size allows for maximum range
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 100;       % 100 frames used for RF cineloop.
Resource.InterBuffer(1).numFrames = 1;  % one intermediate buffer needed.
Resource.ImageBuffer(1).numFrames = 10;
Resource.DisplayWindow(1).Title = 'L11-5vThermistorDemo';
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


%%%%%% THERMISTOR DEMO %%%%%%
% External processing function.
Process(2).classname = 'External';
Process(2).method = 'ThermDemo';
Process(2).Parameters = {'srcbuffer','none',... % name of buffer to process.
                         'dstbuffer','none'};
%%%%%% THERMISTOR DEMO %%%%%%

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
end

%%%%%% THERMISTOR DEMO %%%%%%
Event(n).info = 'call external processing function';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 2;
Event(n).seqControl = 0;
n = n+1;
%%%%%% THERMISTOR DEMO %%%%%%

Event(n).info = 'Jump back to first event';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 1;

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 5;

%% **** Callback routines used by External function definition (EF) ****

import vsv.seq.function.ExFunctionDef

%%%%%% THERMISTOR DEMO %%%%%%
EF(1).Function = vsv.seq.function.ExFunctionDef('ThermDemo', @ThermDemo);
%%%%%% THERMISTOR DEMO %%%%%%


% Save all the structures to a .mat file.
save('MatFiles/L11-5vThermistorDemo');

%% Automatic VSX Execution:
% Uncomment the following line to automatically run VSX every time you run
% this SetUp script (note that if VSX finds the variable 'filename' in the
% Matlab workspace, it will load and run that file instead of prompting the
% user for the file to be used):

% filename = 'L11-5vThermistorDemo';  VSX;


%% **** Callback routines used by External function definition (EF) ****

function ThermDemo()
    persistent thermDemoFig labelHandles

    % initialize figures and parameters on the first call at startup
    
    VDAS     = vsv.seq.getBaseComp('VDAS'); % VDAS==0 for simulator, VDAS==1 on real hardware
    Resource = vsv.seq.getBaseComp( 'Resource' );
    
    if ~isempty(Resource) && ~isempty(VDAS)
        
        connSel        = Resource.Parameters.Connector(1);
        numThermistors = numel(Resource.Parameters.ProbeThermistor);

        % create figure with lables
        if isempty(thermDemoFig) || ~isvalid(thermDemoFig)
            thermDemoFig = uifigure('Position',[80 500 800 100],'Name','Thermistor Readings');

            % add a blank label for each thermistor (position is Left Bottom Width Height)
            labelHandlesList = cell(numThermistors, 1);
            for i=1:numThermistors
                labelHandlesList{i} = uilabel(thermDemoFig,'Text','','Position',[20 60-(i*20) 800 20]);
            end
            labelHandles = [labelHandlesList{:}];
        end

        % sample thermistors and convert to resistance values
        thermistorValues = vsv.hal.shi.thermo.readProbeThermistorValues(connSel, numThermistors,VDAS);
        resistances      = vsv.hal.shi.thermo.convertThermistorValues(thermistorValues);
        
        vsv.hal.shi.thermo.displayReadings(labelHandles, connSel, ...
                                           thermistorValues, resistances);

        drawnow('limitrate');
    end
    
end

