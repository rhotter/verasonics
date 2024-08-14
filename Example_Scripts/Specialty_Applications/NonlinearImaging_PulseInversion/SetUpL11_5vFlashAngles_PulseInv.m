% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: SetUpL11_5vFlashAngles_PulseInv.m - Example of harmonic imaging using pulse inversion with multiple plane waves.
%
% Description:
% This script provides an example of Pulse Inversion (PI) imaging using multiple angle unfocused waves,
%  for contrast agent imaging or tissue harmonic imaging. A broadband probe is required; this script is for the L11-5v.
%  Tissue nonlinearity is detectable but often too weak to be practical with unfocused transmits;
%  beam focusing greatly improves PI results in tissue.
%
% The 64-channel synthetic aperture example script is modified to do two transmits per pulse inversion pair,
%  using all 128 channels. Each acquisition is transferred to the host for examination, and this approach easily permits
%  adjusting the relative contribution of each pulse to explore the relative weighting of the fundamental and the second
%  harmonic images. Reconstruction processing would be faster if the summation was done in hardware accumulation before DMA;
%  in that implementation, the Receive gain could be used to adjust the relative weighting between individual and summed results.
%
% An external function provides plotting of 3 RF data channels, for the first acquisition and the sum of the 2 inverse pulses.
%  This slows down the frame rate due to time required for Matlab plotting, but it is easy to comment out the external
%  processing event for faster frame rates. See the function "plotThreeChannels(RData)" for which channels to plot.
%


%  This example also uses the nonlinear acoustic simulation capability added in V2.9.0, which permits
%  use of the simulator to help develop a nonlinear imaging script. This example uses cells (with the %% delimiter) and runs
%  VSX automatically after defining the variable "filename" as a string with the name of the matfile created by running
%  this script.
%
% Testing: Tested with software release 2.9.0 on Vantage 128,
%
%
% last update: 12/4/2014    - updated the header
%                           - updated the value of the NLCoefficient to work with 2.9.0 scaling of simulated RF data
% 12/15/2015 - Update to SW 3.0 format
% 05/06/2020 - Update to SW 4.3 format for new UIControls and function definitions (VTS 1691). 
%   More info:(.../Example_Scripts/Vantage_Features/New UI Scheme/SetUpL11_5vFlash_NewUI)

%
clear all
filename = 'L11-5vFlashAngles_PulseInv';

%% --- Common User Settings ------------------------------------------------------------------------------------------------

P.startDepth = 0;   % Acquisition depth in wavelengths
P.endDepth = 160;   % This should preferrably be a multiple of 128 samples.
na = 7;             % Set na = number of angles.
angleRange = 36;    % total angular range for plane wave angles (% set dtheta to range over +/- 18 degrees.)
Polarity =  -1;     % +1 or -1: Polarity adjusts the sign of the second pulse with respect to the first one
scaleFactor = 1.5;    % weight of the second acquisition wrt first one: allows adjusting ratio of fundamental to second harmonic in image
                    % 0=first acq only  1=equal weighting  2=second acq only (assuming the two cancel completely when =1)
NLCoefficient =  0.8;  % the nonlinearity coefficient adds harmonic content to the simulation waveform (0=linear,  0< NLCoefficient <= 1.0 is reasonable)
TWfreq = 5;       % transmit frequency (set lower than Trans.frequency for harmonic imaging)
% -------------------------------------------------------------------------------------------------------------------------

if (na > 1)
    dtheta = (angleRange*pi/180)/(na-1);
else
    dtheta=0;
end

% Define system parameters.
Resource.Parameters.numTransmit = 128;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;    % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

%% Specify Media object.
pt1;
Media.function = 'movePoints';

%% Specify Trans structure array.
Trans.name = 'L11-5v';
Trans.frequency = 7.812;
Trans.units = 'wavelengths';  % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);  % L11-5v transducer is 'known' transducer so we can use computeTrans.
Trans.maxHighVoltage = 50;  % set maximum high voltage limit for pulser supply.

% Specify PData structure array.
PData(1).PDelta = [0.5, 0, 0.25];
PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3)); % startDepth, endDepth and pdelta set PData(1).Size.
PData(1).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(3) = 1;      % single image page
PData(1).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of upper lft crnr.
% No PData.Region specified, so a default Region for the entire PData array will be created by computeRegions.

%% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 2*na*3072; % this size allows for maximum range
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 30;    % 30 frames stored in RcvBuffer.
Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).numFrames = 1;   % one intermediate buffer needed.
Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).numFrames = 10;
Resource.DisplayWindow(1).Title = filename;
Resource.DisplayWindow(1).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = gray(256);

%% Specify Transmit waveform structure.
% use a relatively low frequency (4.5) to be better able to detect the second harmonic (9 MHz)
TW(1).type = 'parametric';
TW(1).Parameters = [TWfreq,1,2,1];   % positive polarity
TW(1).nonLinear = NLCoefficient;   % can comment out to obtain usual linear result in simulation (has no effect when using hardware)

TW(2) = TW(1);
TW(2).Parameters(4) = Polarity*TW(1).Parameters(4);   % negative polarity

%% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'Apod', kaiser(Resource.Parameters.numTransmit,2)', ...
                   'focus', 0.0, ...
                   'Steer', [0.0,0.0], ...
                   'Delay', zeros(1,Trans.numelements)), 1, 2*na);
% - Set event specific TX attributes.
if fix(na/2) == na/2       % if na even
    startAngle = (-(fix(na/2) - 1) - 0.5)*dtheta;
else
    startAngle = -fix(na/2)*dtheta;
end
for n = 2:2:2*na   % 2*na transmit events
    TX(n-1).Steer = [(startAngle+(n/2-1)*dtheta),0.0]; % TW(1)
    TX(n-1).Delay = computeTXDelays(TX(n-1));

    TX(n).Steer = TX(n-1).Steer;                       % TW(2)
    TX(n).Delay = TX(n-1).Delay;
    TX(n).waveform = 2;
end

%% Specify Receive structure arrays.
% - We need 2*na Receives for every frame.
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
wlsPer128 = 128/(4*2); % wavelengths in 128 samples for 4 samplesPerWave
Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', P.startDepth + wlsPer128*ceil(maxAcqLength/wlsPer128), ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, 2*na*Resource.RcvBuffer(1).numFrames);
% - Set event specific Receive attributes for each frame.
for i = 1:Resource.RcvBuffer(1).numFrames
    k = 2*na*(i-1);
    Receive(k+1).callMediaFunc = 1;    % only move media points once per frame, not between angles
    for j = 1:2:2*na
        Receive(k+j).Apod(1:Resource.Parameters.numRcvChannels) = 1.0; % TW(1)
        Receive(k+j).framenum = i;
        Receive(k+j).acqNum = j;
        Receive(k+j+1).Apod(1:Resource.Parameters.numRcvChannels) = 1.0; % TW(2)
        Receive(k+j+1).framenum = i;
        Receive(k+j+1).acqNum = j+1;
    end
end

%% Specify TGC Waveform structure.
TGC.CntrlPts = [0,297,424,515,627,764,871,1000];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

%% Specify Recon structure arrays.
% - We need one Recon structure which will be used for each frame.
Recon = struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame',-1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...
               'RINums',1:2*na);

% Define ReconInfo structures.
% We need 2*na ReconInfo structures for 2-1 synthetic aperture and na steering angles.
ReconInfo = repmat(struct('mode', 'accumIQ', ...  % default is to accumulate IQ data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 1), 1, 2*na);
% - Set specific ReconInfo attributes.
ReconInfo(1).mode = 'replaceIQ'; % replace IQ data
for j = 1:2:2*na  % For each row in the column
    ReconInfo(j).txnum = j;     % first acq for each angle
    ReconInfo(j).rcvnum = j;

    ReconInfo(j+1).txnum = j+1; % second acq for each angle
    ReconInfo(j+1).rcvnum = j+1;
    ReconInfo(j+1).scaleFactor = scaleFactor;
end
ReconInfo(2*na).mode = 'accumIQ_replaceIntensity'; % accum and detect

%% Specify Process structure array.
pers = 0;
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

Process(2).classname = 'External';
Process(2).method = 'plotThreeChannels';
Process(2).Parameters = {'srcbuffer','receive',... % name of buffer to process.
                         'srcbufnum',1,...
                         'srcframenum',-1,...
                         'dstbuffer','none'};

%% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; % jump back to start
SeqControl(1).argument = 1;

SeqControl(2).command = 'timeToNextAcq';  % time between PI acquisitions
twoWayTravel = 2.5 * Receive(1).endDepth  / Trans.frequency; % microsecs
SeqControl(2).argument = twoWayTravel;  % formerly 160 usec, but now as fast as possible

SeqControl(3).command = 'timeToNextAcq';  % time between frames
SeqControl(3).argument = 30000;  % 30 msec
SeqControl(4).command = 'returnToMatlab';
nsc = 5; % nsc is next count of SeqControl objects

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 1;

%% Specify Event structure arrays.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    for j = 1:2:2*na                      % Acquire frame
        Event(n).info = 'First Waveform';
        Event(n).tx = j;
        Event(n).rcv = 2*na*(i-1)+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;

        Event(n).info = 'Second Waveform';
        Event(n).tx = j+1;
        Event(n).rcv = 2*na*(i-1)+j+1;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;
    end
    Event(n-1).seqControl = [3,nsc]; % modify last acquisition Event's seqControl
      SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
      nsc = nsc+1;

    Event(n).info = 'display RF line';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 2;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = 'recon and process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 1;
    Event(n).seqControl = 4;
    n = n+1;
end

Event(n).info = 'Jump back';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 1;


%% User specified UI Control Elements

import vsv.seq.uicontrol.VsSliderControl

% - Sensitivity Cutoff
UI(1).Control = VsSliderControl('LocationCode','UserB7','Label','Sens. Cutoff',...
                   'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                   'SliderStep',[0.025,0.1],'ValueFormat','%1.3f',...
                   'Callback', @SensCutoffCallback);


% - Range Change
MinMaxVal = [64,300,P.endDepth]; % default unit is wavelength
AxesUnit = 'wls';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
        AxesUnit = 'mm';
        MinMaxVal = MinMaxVal * (Resource.Parameters.speedOfSound/1000/Trans.frequency);
    end
end
UI(2).Control = VsSliderControl('LocationCode','UserA1','Label',['Range (',AxesUnit,')'],...
                    'SliderMinMaxVal',MinMaxVal,...
                    'SliderStep',[0.1,0.2],'ValueFormat','%3.0f',...
                    'Callback', @RangeChangeCallback);


% - Relative Weighting of second acquisition in Image
UI(3).Control = VsSliderControl('LocationCode','UserB1','Label','Acq 2 Weight',...
                    'SliderMinMaxVal',[0.05,2.0,scaleFactor],...
                    'SliderStep',[0.01,0.1],'ValueFormat','%1.3f',...
                    'Callback', @Weigh2PulseCallback);


% External function delimiter definition
EF(1).Function = vsv.seq.function.ExFunctionDef('plotThreeChannels',@plotThreeChannels);

% Save all the structures to a .mat file.
save(['MatFiles/',filename]);

VSX     % run VSX automatically!

% when VSX exits, the following plot will be created

figure ('Position', [680 560 860 535])  % plot the TW waveform for both individual pulses and their sum.
    plot(TW(1).Wvfm2Wy); str1 = ([ 'TWfreq = ' num2str(TWfreq) '   Trans.frequency = ' num2str(Trans.frequency) ...
        '   Trans.Bandwidth = [' num2str(Trans.Bandwidth(1)) ', ' num2str(Trans.Bandwidth(2)) ']  TW.peak = ' num2str(TW.peak)]);
    str2 = 'Simulated transmit waveforms'; title ({str2, str1})
    hold on, plot(TW(2).Wvfm2Wy, 'k'), plot(TW(1).Wvfm2Wy+TW(2).Wvfm2Wy, 'r'), legend('First transmission', 'Second transmission', 'Pulse Inversion sum'), hold off
    ylabel('Simulation Waveforms')

return


%% **** Callback routines used by UIControls (UI) ****

function SensCutoffCallback(~,~,UIValue)
%Sensitivity cutoff change
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

function RangeChangeCallback(hObject,~,UIValue)
%Range change
    UI = evalin('base','UI');
    set(UI(2).handle,'Interruptible','off');

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
    twoWayTravel = 2.5 * Receive(1).endDepth  / Trans.frequency; % microsecs
    SeqControl = evalin('base', 'SeqControl');
    SeqControl(2).argument = twoWayTravel;  % formerly 160 usec
    assignin('base','SeqControl',SeqControl);
    assignin('base','Receive',Receive);
    evalin('base','TGC.rangeMax = P.endDepth;');
    evalin('base','TGC.Waveform = computeTGCWaveform(TGC);');
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'PData','InterBuffer','ImageBuffer','DisplayWindow','Receive','TGC','Recon','SeqControl'};
    assignin('base','Control', Control);
    assignin('base', 'action', 'displayChange');
end

function Weigh2PulseCallback(~,~,UIValue)
%Weight of second pulse
    scaleFactor = UIValue;
    assignin ('base', 'scaleFactor', scaleFactor) % update the parameter
    ReconInfo = evalin('base', 'ReconInfo');
    numRI = size(ReconInfo,2); % IMPORTANT: the number of ReconInfo structures is changed by VSX (multiplied by the number of regions)
    for j = 1:2:numRI  % For every second ReconInfo
        ReconInfo(j+1).scaleFactor = scaleFactor;
    end
    assignin('base','ReconInfo',ReconInfo);
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'Recon'}; % to update the Reconstruction parameters, update Recon even though only ReconInfo has changed.
    assignin('base','Control', Control);
end

%% **** Callback routines used by External function definition (EF) ****

function plotThreeChannels(RData)
    % This plotting function is used to display the RF data from 3 channels.
    % The first acquisitions are in Blue, Red, and Black colored traces.
    % The sum of two consecutive pulses (pulse inversion summation) is drawn in Cyan for each channel.

    persistent myHandle10 endSample

    channel = 40;  % first Channel no. to plot
    xlims = [200,1600];
    xrange1 = xlims(1):xlims(2); % sample range for first acquisition
    ylims = 1000*[-1 1];
    dely = ylims(2);

    endSample = evalin ('base', 'Receive(1).endSample');
    xrange2 = xrange1 + endSample; % range of second acquisition

    if isempty(myHandle10) || ~ishandle(myHandle10)
        figure ('Position', [550 750 650 300]); %, 'MenuBar', 'none');
        myHandle10 = axes('XLim',xlims,'YLim',2*ylims);
        title(['Channel Numbers = ' num2str(channel)])
    end

    set(myHandle10,'NextPlot','replacechildren');
    plot(myHandle10, xrange1, RData(xrange1,channel)+dely, 'b');      %  first acq
    set(myHandle10,'NextPlot','add');
    plot(myHandle10, xrange1, RData(xrange2,channel)+dely, 'b--');           %  second acq
    plot(myHandle10, xrange1, RData(xrange1,channel)+RData(xrange2,channel)+dely, 'c'); % sum of first and second acqs

    % plot(myHandle10, xrange1, RData(xrange1,c2),      'r');
    % plot(myHandle10, xrange1, RData(xrange1,c2)+RData(xrange2,c2),      'c');
    % plot(myHandle10, xrange1, RData(xrange1,c3)-dely, 'k');
    % plot(myHandle10, xrange1, RData(xrange1,c3)+RData(xrange2,c3)-dely, 'c');

    % zoom on
    drawnow nocallbacks
end

