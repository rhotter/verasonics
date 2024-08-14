% File name: SetUpL11_5vFlash_newUI.m - Modified copy of L11_5vFlash to
% demonstrate the new approach added in 4.3.0 (VTS-1480) for defining GUI
% controls, callbacks and external function handles.  This script
% demonstrates the new approach; refer to companion script
% "SetUpL11_5vFlash_OldUI.m" for an example of the same script using the
% original "text2cell" approach.
% 
% To demonstrate the differences between the new and old approaches, a new
% GUI control has been added, that simply changes the transmit burst
% duration. An external processing function has also been added, that
% records and reports how long the script has been running.
%
% The new approach uses a more advanced object oriented approach for
% defining GUI controls, which brings several benefits for the user. For
% example, it is harder now to make typing mistakes or entering values for
% parameters that would be invalid. In these cases the error will be thrown
% directly when script is executed. Furthermore, it is possible now to
% define real functions in your script and to pass real function handles to
% the callbacks of the GUI controls or the external functions. This has the
% benefit that the user can use Matlabs' tools to debug script functions
% directly. 
%
% To modify your script in order to take advantage of the new approach,
% only a few changes are required. To help our users to identify these
% changes we marked them in this script with **NewUI approach** and
% extended documentation.
%
% This script only demonstrates examples with sliders. However, we also
% support all other UIControls plus, with the new approach, a drop down
% menu. You can get easy help with the new approach: e.g., type
% >> help vsv.seq.uicontrol
% and you will get a list of possible UI controls. The following UIcontrols
% are currently supported:
%
%     VsButtonControl       - describes a use ui button control
%     VsButtonGroupControl  - a group of selectable options using radio buttons
%     VsDropDownControl     - defines a drop down menu that can be added to a GUI
%     VsSliderControl       - defines a slider control handle that can be added to a GUI
%     VsToggleButtonControl - uicontrol defining a toggle button with a state
% 
% You can also get help by using the following command
% >> help vsv.seq.uicontrol.VsSliderControl
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
Resource.DisplayWindow(1).Title = 'L11-5vFlash new UI';
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

% **NewUI approach**
% For this script, all of the older "text2cell" format GUI control
% definitions have been converted to the new package scheme
%
% You can use the import command to import controls. This will make it
% easier in the following ui control defintions because you won't need to
% type 'vsv.seq.uicontrol.' anymore. Instead you can use VsSliderControl
% directly. 
% 
% >> import vsv.seq.uicontrol.VsSliderControl
% >> UI(2).Control = VsSliderControl(...)
%
% However, you can also use the full name to define the UIcontrol:
% >> UI(2).Control = vsv.seq.uicontrol.VsSliderControl(...)
%
% This has the advantage that you can use the TAB key. When pressing the
% TAB key after the '.' (e.g., vsv.scrip.uicontro.<hit TAB>) Matlab will
% prompt some suggestion you can choose from. This has the advantage that
% you do not need to remember every command.
import vsv.seq.uicontrol.VsSliderControl; 

% **NewUI approach**
% define the slider and add the callback; note 'Style' field is no
% longer needed.
%
% Note, the callback is directly defined as part of the control object; the
% separate definition of the callback using the text2cell function has been
% eliminated. The callback is now a real function handle (as indicated by
% the @ sign). For support about function handles use help function_handle
% or see the Matlab documentation

% - Sensitivity Cutoff
UI(1).Control = VsSliderControl('LocationCode','UserB7','Label','Sens. Cutoff',...
                  'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                  'SliderStep',[0.025,0.1],'ValueFormat','%1.3f', ...
                  'Callback', @SensCutoffCallback );  

% - Range Change
MinMaxVal = [64,300,P.endDepth]; % default unit is wavelength
AxesUnit = 'wls';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
        AxesUnit = 'mm';
        MinMaxVal = MinMaxVal * (Resource.Parameters.speedOfSound/1000/Trans.frequency);
    end
end

% **NewUI approach**
% define the slider and add the callback; note 'Style' field is no longer
% needed.

UI(2).Control = VsSliderControl('LocationCode','UserA1','Label',['Range (',AxesUnit,')'],...
                 'SliderMinMaxVal',MinMaxVal,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f', ...
                 'Callback', @RangeChangeCallback);  

% New UI Control for GUI Demo; modify TX burst duration
% In this script, the new control is defined using the new scheme
% added in 4.3.0

% **NewUI approach**
% define the slider and add the callback; note 'Style' field is no
% longer needed.
UI(3).Control = VsSliderControl('LocationCode', 'UserB3', ...
                            'Label','TX Duration cycles', ...
                            'SliderMinMaxVal',[0 20 1], ...
                            'SliderStep',[0.1,0.25], ...
                            'ValueFormat','%3.0f', ...
                            'Callback', @TXdurChangeCallback ); 


% External function definitions

% Create the External Processing function, using the new function handle
% approach. Note that in this example no import was used. The user can
% directly call the new approach functions using the full package path
% You can use the TAB-key after you place a '.' to get suggestions for
% available functions. 
%
% Note the callback function is a real function handle as indicated by the
% @ sign. You can now mark this function handle, and right click and then
% select 'Open runTimeMon' and the editor will jump to the function
% definition in your script
EF(1).Function = vsv.seq.function.ExFunctionDef('runTimeMon', @()runTimeMon(obj));

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 5;

% Save all the structures to a .mat file.
save('MatFiles/L11-5vFlashNewUI');

%% Automatic VSX Execution:
% Uncomment the following line to automatically run VSX every time you run
% this SetUp script (note that if VSX finds the variable 'filename' in the
% Matlab workspace, it will load and run that file instead of prompting the
% user for the file to be used):

% filename = 'L11-5vFlashNewUI';  VSX;


%% **** Callback routines used by UIControls (UI) ****

% **NewUI approach**
%  Note the text2cell delimiter comments are no longer needed.
% Note also that no changes are needed to the actual code in the callback
% function at all; the only changes are to add the Matlab function
% definition line at the start of the function, and the required "end"
% statement at the end.  In these exaples the body of the funtion has been
% indented for readability.
%
% In the new callback function definition, three input arguments are
% required:  "(hObject, evt, UIValue)".
% However, the "hObject" and "evt" arguments are not always used in these
% examples; "evt" is a placeholder for additional features that will be
% added in the future.  In the examples given below, the unused input
% arguments have been replaced with "~".
% 
% SensCutoffCallback is a real nested function and is used with 
% @SensCutoffCallback as a function_handle the UI(1) callback

% - Sensitivity cutoff change callback
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

% **NewUI approach**
% RangeChangeCallback is a real nested function and is used with 
% @RangeChangeCallback as a function_handle the UI(2) callback

% - Range change callback
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

% **NewUI approach**
% TXdurChangeCallback is a real nested function and is used with 
% @TXdurChangeCallback as a function_handle the UI(3) callback

% - TXdur change callback
function TXdurChangeCallback(~, ~, UIValue)
    TW = evalin('base', 'TW');
    TW(1).Parameters(3) = 2*round(UIValue);
    assignin('base','TW',TW);
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'TW'};
    assignin('base','Control', Control);
end

% **NewUI approach**
% runTimeMon External Processing function using new nested functions

function runTimeMon(varargin)
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
end