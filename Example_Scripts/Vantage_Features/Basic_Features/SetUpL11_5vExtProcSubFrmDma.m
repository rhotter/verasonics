% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: SetUpL11_5vExtProcSubFrmDma.m - Example of linear array imaging with
%                                     focused transmits modified to use DMA
%                                     subframes
% Description:
%   This script is a version of SetUpL7_4_128RyLns that has been updated to
%   use the L11-5v transducer and modified to serve as an example of using DMA
%   subframes while calling an external processing function instead of
%   using Recon/
%
%
% Last update:
% 05/06/2020 - Update to SW 4.3 format for new user UIControls and External function definitions (VTS 1691).
%   More info:(?/Example_Scripts/Vantage_Features/New UI Scheme/SetUpL11_5vFlash_NewUI)
% 8/27/18 - modified for 4.0 release to eliminate range change etc. not
% intended to be used in this demonstration script (VTS-837)

clear all

P.startDepth = 0;
P.endDepth = 192;
P.txFocus = 0;  % Initial transmit focus.
P.numTx = 32;  % number of transmit elements in TX aperture (where possible).

nf = 4; % number of receive data frames

% Define system parameters.
Resource.Parameters.numTransmit = 128;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;    % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.
Resource.Parameters.fakeScanhead = 1;
% Specify Trans structure array.
Trans.name = 'L11-5v';
Trans.units = 'wavelengths'; % Explicit declaration avoids warning message when selected by default
Trans.frequency = 5.2083;
Trans = computeTrans(Trans);  % L11-5v transducer is 'known' transducer so we can use computeTrans.
% note nominal center frequency from computeTrans is 6.25 MHz
Trans.maxHighVoltage = 50;  % set maximum high voltage limit for pulser supply.
Trans.ElementPos(:, 1) = 0;

% Specify PData structure array.
PData(1).PDelta = [Trans.spacing, 0, 0.5];  % x, y, z pdeltas
PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));
PData(1).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(3) = 1;      % single image page
PData(1).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of upper lft crnr
% - specify 128 Region structures.
PData(1).Region = repmat(struct('Shape',struct( ...
                    'Name','Rectangle',...
                    'Position',[0,0,P.startDepth],...
                    'width',Trans.spacing,...
                    'height',P.endDepth-P.startDepth)),1,128);
% - set position of regions to correspond to beam spacing.
for i = 1:128
    PData(1).Region(i).Shape.Position(1) = (-63.5 + (i-1))*Trans.spacing;
end
PData(1).Region = computeRegions(PData(1));

% Specify Media object.
Media.MP = [0 0 1 1];
Media.numPoints = 1;
% Media.attenuation = -0.5;
Media.function = 'movePoints';

% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 4096*128;%P.endDepth*2*128;   % this size allows for maximum range
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = nf;       % nf frames used for RF cineloop.
Resource.InterBuffer(1).numFrames = 1;  % one intermediate buffer needed.
Resource.ImageBuffer(1).numFrames = 1;
Resource.DisplayWindow(1).Title = 'L11-5v_128RyLns';
Resource.DisplayWindow(1).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).numFrames = 20;
Resource.DisplayWindow(1).AxesUnits = 'wavelengths';
Resource.DisplayWindow.Colormap = gray(256);

% Specify TW structure array.
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,0.67,2,1];

% Specify nf*128 TX structure arrays. Transmit centered on element n in the array for event n.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'focus', P.txFocus, ...
                   'Steer', [0.0,0.0], ...
                   'Apod', zeros(1,Trans.numelements), ...
                   'Delay', zeros(1,Trans.numelements)), 1, 128*nf);
% - Set event specific TX attributes.
for fr = 1:nf
    for n = 1:128   % 128 transmit events
        TX(n+128*(fr-1)).Apod(n) = 1.0;
        dly = 8*fr+n/2;
        if dly>72
            dly = dly-64;
        end
        TX(n+128*(fr-1)).Delay(n) = dly;
    end
end

% Specify Receive structure arrays.
% - We need 128 Receives for every frame.
% maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', P.endDepth/4, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW',...
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, 128*Resource.RcvBuffer(1).numFrames);
% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    k = 128*(i-1);
%     Receive(k+1).callMediaFunc = 1;
    for j = 1:128
        Receive(k+j).framenum = i;
        Receive(k+j).acqNum = j;
    end
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [0,138,260,287,385,593,674,810];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% % Specify Recon structure array.
% Recon = struct('senscutoff', 0.5, ...
%                'pdatanum', 1, ...
%                'rcvBufFrame',-1, ...
%                'IntBufDest', [1,1], ...
%                'ImgBufDest', [1,1], ...  % auto-increment ImageBuffer each recon
%                'RINums', 1:128);
%
% % Define ReconInfo structures.
% ReconInfo = repmat(struct('mode', 'replaceIntensity', ...
%                    'txnum', 1, ...
%                    'rcvnum', 1, ...
%                    'regionnum', 0), 1, 128);
% % - Set specific ReconInfo attributes.
% for j = 1:128
%     ReconInfo(j).txnum = j;
%     ReconInfo(j).rcvnum = j;
%     ReconInfo(j).regionnum = j;
% end

% Specify Process structure array.
pers = 20;
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',1.0,...            % pgain is image processing gain
                         'reject',2,...      % reject level
                         'persistMethod','simple',...
                         'persistLevel',pers,...
                         'interpMethod','4pt',...  %method of interp. (1=4pt)
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',40,...
                         'mappingMethod','full',...
                         'display',1,...      % display image after processing
                         'displayWindow',1};

Process(2).classname = 'External';
Process(2).method = 'SFDmaImg';
Process(2).Parameters = {'srcbuffer','receive',...  % name of buffer to process.
                         'srcbufnum',1,...
                         'srcframenum',-1,...
                         'dstbuffer','image',...
                         'dstbufnum',1,...
                         'dstframenum',1};

% Specify SeqControl structure arrays.
%  - Time between acquisitions in usec
t1 = round(2*768*(1/Trans.frequency)); % acq. time in usec for max depth
SeqControl(1).command = 'timeToNextAcq';
SeqControl(1).argument = t1;
%  - Time between frames at 20 fps at max endDepth.
SeqControl(2).command = 'timeToNextAcq';
SeqControl(2).argument = 500000;
%  - Return to Matlab
SeqControl(3).command = 'returnToMatlab';
%  - Jump back to start.
SeqControl(4).command = 'jump';
SeqControl(4).argument = 1;
%  - sync.
SeqControl(5).command = 'sync';
nsc = 6; % next SeqControl number

% Specify Event structure arrays.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    for j = [96:-1:37, 97:128, 36:-1:1] % Acquire frame jumbled order in three blocks
%     for j = 128:-1:1                    % Acquire frame reverse order 128 to 1
%     for j = 1:128                       % Acquire frame normal 'forward' order 1 to 128
        Event(n).info = 'Acquisition.';
        Event(n).tx =  j+(i-1)*128;   % use next TX structure.
        Event(n).rcv = j+(i-1)*128;
        Event(n).recon = 0;      % no reconstruction.
        Event(n).process = 0;    % no processing
        Event(n).seqControl = 1; % seqCntrl
        n = n+1;
    end
    % Replace last events SeqControl for inter-frame timeToNextAcq.
   Event(n-1).seqControl = 2;

   % insert a couple subframe DMA commands within the previous frame's
   % events
    Event(n-70).seqControl = [1, nsc]; % acqnum 59  (129-70)
       SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
       nsc = nsc+1;

    Event(n-30).seqControl = [1, nsc]; % acqnum 99  (129-30)
       SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
       nsc = nsc+1;

    Event(n).info = 'Transfer frame to host.';
    Event(n).tx = 0;        % no TX
    Event(n).rcv = 0;       % no Rcv
    Event(n).recon = 0;     % no Recon
    Event(n).process = 0;
    Event(n).seqControl = nsc;
       SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
       nsc = nsc+1;
    n = n+1;

    Event(n).info = 'external process funct call';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 0;      % reconstruction
    Event(n).process = 2;    % process
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = 'recon and process';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 0;      % reconstruction
    Event(n).process = 1;    % process
    Event(n).seqControl = 0;
    if (floor(i/5) == i/5)&&(i ~= Resource.RcvBuffer(1).numFrames)  % Exit to Matlab every 5th frame
        Event(n).seqControl = 3; % return to Matlab
    end
    n = n+1;

    Event(n).info = 'sync';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 0;      % reconstruction
    Event(n).process = 0;    % process
    Event(n).seqControl = 5;
    n = n+1;
end

Event(n).info = 'Jump back';
Event(n).tx = 0;        % no TX
Event(n).rcv = 0;       % no Rcv
Event(n).recon = 0;     % no Recon
Event(n).process = 0;
Event(n).seqControl = 4;


% User specified UI Control Elements

import vsv.seq.uicontrol.VsSliderControl

% % - Range Change
% wls2mm = 1;
% AxesUnit = 'wls';
% if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
%     if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
%         AxesUnit = 'mm';
%         wls2mm = Resource.Parameters.speedOfSound/1000/Trans.frequency;
%     end
% end
% UI(2).Control = VsSliderControl('LocationCode','UserA1',...
%     'Label',['Range (',AxesUnit,')'],...
%     'SliderMinMaxVal',[64,300,P.endDepth]*wls2mm,...
%     'SliderStep',[0.1,0.2],'ValueFormat','%3.0f',...
%     'Callback',@RangeChangeCallback);
%
% % - Transmit focus change
% UI(3).Control = VsSliderControl('LocationCode','UserB4',...
%     'Label',['TX Focus (',AxesUnit,')'],...
%     'SliderMinMaxVal',[20,320,P.txFocus]*wls2mm,...
%     'SliderStep',[0.1,0.2],'ValueFormat','%3.0f',...
%     'Callback',@TxFocusCallback);
%
% % - F number change
% UI(4).Control = VsSliderControl('LocationCode','UserB3',...
%     'Label','F Number',...
%     'SliderMinMaxVal',[1,20,round(P.txFocus/(P.numTx*Trans.spacing))],...
%     'SliderStep',[0.05,0.1],'ValueFormat','%2.0f',...
%     'Callback',@FNumCallback);


%% **** Callback routines used by External function definition (EF) ****

EF(1).Function = vsv.seq.function.ExFunctionDef('SFDmaImg',@SFDmaImg);

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 4;

% Save all the structures to a .mat file.
save('MatFiles/L11-5vExtProcSubFrmDma');

%% Automatic VSX Execution:
% Uncomment the following line to automatically run VSX every time you run
% this SetUp script (note that if VSX finds the variable 'filename' in the
% Matlab workspace, it will load and run that file instead of prompting the
% user for the file to be used):

% filename = 'L11-5vExtProcSubFrmDma';  VSX;


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

    PData = evalin('base','PData');
    PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));
    PData(1).Region = repmat(struct('Shape',struct( ...
                        'Name','Rectangle',...
                        'Position',[0,0,P.startDepth],...
                        'width',Trans.spacing,...
                        'height',P.endDepth-P.startDepth)),1,128);
    % - set position of regions to correspond to beam spacing.
    for i = 1:128
        PData(1).Region(i).Shape.Position(1) = (-63.5 + (i-1))*Trans.spacing;
    end
    assignin('base','PData',PData);
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

function TxFocusCallback(hObject, ~, UIValue)
    simMode = evalin('base','Resource.Parameters.simulateMode');
    % No focus change if in simulate mode 2.
    if simMode == 2
        set(hObject,'Value',evalin('base','P.txFocus'));
        return
    end
    Trans = evalin('base','Trans');
    Resource = evalin('base','Resource');
    scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);

    P = evalin('base','P');
    P.txFocus = UIValue;
    if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
        if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
            P.txFocus = UIValue*scaleToWvl;
        end
    end
    assignin('base','P',P);

    TX = evalin('base', 'TX');
    for n = 1:128   % 128 transmit events
        TX(n).focus = P.txFocus;
        TX(n).Delay = computeTXDelays(TX(n));
    end
    assignin('base','TX', TX);
    % Update Fnumber based on new P.txFocus
    evalin('base','set(UI(4).handle(2),''Value'',round(P.txFocus/(P.numTx*Trans.spacing)));');
    evalin('base','set(UI(4).handle(3),''String'',num2str(round(P.txFocus/(P.numTx*Trans.spacing))));');
    % Set Control command to update TX
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'TX'};
    assignin('base','Control', Control);
end

function FNumCallback(hObject,~,UIValue)
    simMode = evalin('base','Resource.Parameters.simulateMode');
    P = evalin('base','P');
    Trans = evalin('base','Trans');
    % No F number change if in simulate mode 2.
    if simMode == 2
        set(hObject,'Value',round(P.txFocus/(P.numTx*Trans.spacing)));
        return
    end
    P.txFNum = UIValue;
    P.numTx = round(P.txFocus/(P.txFNum*Trans.spacing));
    assignin('base','P',P);
    % - Redefine event specific TX attributes for the new P.numTx.
    TX = evalin('base', 'TX');
    for n = 1:128   % 128 transmit events
        % Set transmit Apodization.
        lft = n - floor(P.numTx/2);
        if lft < 1, lft = 1; end
        rt = n + floor(P.numTx/2);
        if rt > Trans.numelements, rt = Trans.numelements; end
        TX(n).Apod = zeros(1,Trans.numelements);
        TX(n).Apod(lft:rt) = 1.0;
        TX(n).Delay = computeTXDelays(TX(n));
    end
    assignin('base','TX', TX);
    % Set Control command to update TX
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'TX'};
    assignin('base','Control', Control);
end

function ImgDat = SFDmaImg(RData)
    ImgDat = zeros(384, 128);
    % ImgDat(100:110, :) = 1e6;
    RDataP = 1e3*abs(double(RData));
    % pks = max(RDataP(:,1))
    % pks = find(
    for i = 1:128
    %     ImgDat(i:(i+10), i) = 1e6;
        ImgDat(:, i) = RDataP((i-1)*384 + (1:384), i);
    end
    persistent myHandle
    Receive = evalin('base','Receive');
    % If myPlotChnl exists, read it for the channel to plot.
    if evalin('base','exist(''myPlotChnl'',''var'')')
        channel = evalin('base','myPlotChnl');
    else
        channel = 64;  % Channel no. to plot
    end

    % Create the figure if it does not exist.
    if isempty(myHandle)||~ishandle(myHandle)
        figure('name','Receive Signal','NumberTitle','off');
        myHandle = axes('XLim',[0,Receive(1).endSample],'YLim',[-2048 2048], ...
                        'NextPlot','replacechildren');
    end

    % Plot the element's RF data.
    plot(myHandle,RData(1:Receive(1).endSample,channel));
    drawnow limitrate
end
