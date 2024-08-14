% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: SetUpL11_5v_128RyLns_GPU_RF.m
%
% Description:
%   - Example of external FIR filtering of the raw RF data on the GPU
%      before being passed back to the CPU for reconstruction and display
%
%   - 3 different external functions performing FIR filtering on RF data
%      are called for comparison:
%       1. FIR filtering using a matlab script
%       2. FIR filtering using parallel computing toolbox on GPU
%       3. FIR filtering using GPU with cuda code
%
%   - When the script runs synchronously, the user can note the hardware
%      framerate change to match the software framerate.
%
%   - The gui button group allows the user to update&run to switch quickly
%      between methods
%
%  NOTE:  GPU/synchronization specific code are commented with @gpu for
%  easy searching


clear all

Resource.Parameters.waitForProcessing = 0;  % Toggle this to switch between async/sync acquisitions

Resource.VDAS.dmaTimeout = 1000;

P.startDepth = 5;
P.endDepth = 128;
P.txFocus = 100;  % Initial transmit focus.
P.numTx = 16;  % number of transmit elements in TX aperture (where possible).

% Define system parameters.
Resource.Parameters.numTransmit = 128;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;   % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Specify Trans structure array.
Trans.name = 'L11-5v';
Trans.units = 'wavelengths'; % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);  % L11-5v transducer is 'known' transducer so we can use computeTrans.
% note nominal center frequency from computeTrans is 6.25 MHz
Trans.maxHighVoltage = 50;  % set maximum high voltage limit for pulser supply.

% Specify PData structure array.
PData(1).PDelta = [Trans.spacing, 0, 1];  % x, y, z pdeltas
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
pt1;
Media.attenuation = -0.5;
Media.function = 'movePoints';

% Specify Resources in external input/output pairs
%--- Normal Host Memory ---%
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 15360*P.numTx;   % this size allows for maximum range
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 10;

%--- CUDA Host Memory ---%
Resource.RcvBuffer(2) = Resource.RcvBuffer(1);
Resource.RcvBuffer(2).gpuMemoryType = 'CPU';

%--- CUDA Unified Memory ---%
Resource.RcvBuffer(3)=Resource.RcvBuffer(1);
Resource.RcvBuffer(3).gpuMemoryType = 'Unified';
 
%@gpu we use a second RcvBuffer to pass the filtered RF data to for reconstruction
Resource.RcvBuffer(4)=Resource.RcvBuffer(1);
Resource.RcvBuffer(4).numFrames = 1;
Resource.RcvBuffer(4).gpuMemoryType = 'CPU';

Resource.InterBuffer(1).numFrames = 1;
Resource.ImageBuffer(1).numFrames = 10;
Resource.DisplayWindow(1).Title = 'L11-5v_128RyLns_GPU_RF';
Resource.DisplayWindow(1).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
    DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Matlab';
Resource.DisplayWindow(1).numFrames = 20;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow.Colormap = gray(256);

% Specify TW structure array.
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,0.67,2,1];

% Specify nr TX structure arrays. Transmit centered on element n in the array for event n.
TX = repmat(struct('waveform', 1, ...
    'Origin', [0.0,0.0,0.0], ...
    'focus', P.txFocus, ...
    'Steer', [0.0,0.0], ...
    'Apod', zeros(1,Trans.numelements), ...
    'Delay', zeros(1,Trans.numelements)), 1, 128);
% - Set event specific TX attributes.
for n = 1:128   % 128 transmit events
    % Set transmit Origins to positions of elements.
    TX(n).Origin = [(-63.5 + (n-1))*Trans.spacing, 0.0, 0.0];
    % Set transmit Apodization.
    lft = n - floor(P.numTx/2);
    if lft < 1, lft = 1; end
    rt = n + floor(P.numTx/2);
    if rt > Trans.numelements, rt = Trans.numelements; end
    TX(n).Apod(lft:rt) = 1.0;
    TX(n).Delay = computeTXDelays(TX(n));
end

% Specify Receive structure arrays.
% - We need 128 Receives for every frame.
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
    'startDepth', P.startDepth, ...
    'endDepth', maxAcqLength, ...
    'TGC', 1, ...
    'bufnum', 1, ...
    'framenum', 1, ...
    'acqNum', 1, ...
    'sampleMode', 'NS200BW',...
    'mode', 0, ...
    'callMediaFunc', 0), 1, 128*(Resource.RcvBuffer(1).numFrames+1));
% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    k = 128*(i-1);
    Receive(k+1).callMediaFunc = 1;
    for j = 1:128
        Receive(k+j).framenum = i;
        Receive(k+j).acqNum = j;
    end
end
% @gpu: A second set of receive structures for reconstructing the
% processed RF data is needed to point to the second RcvData Buffer
for i = 1:128
    Receive(128*Resource.RcvBuffer(1).numFrames + i) = Receive(i); %@gpu: make it identical to the 1st set
    Receive(128*Resource.RcvBuffer(1).numFrames + i).framenum = 1;
    Receive(128*Resource.RcvBuffer(1).numFrames + i).bufnum = 4;   %@gpu: Change the buffer number to 4 (the single frame RcvBuffer)
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [200,375,490,560,630,700,765,830];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Recon structure array.
%@gpu: Recon for frame 1
Recon = struct('senscutoff', 0.5, ...
    'pdatanum', 1, ...
    'rcvBufFrame', 1, ...
    'IntBufDest', [1,1], ...
    'ImgBufDest', [1,-1], ...  %@gpu: auto-increment ImageBuffer each recon
    'RINums', [1:128]);

% Define ReconInfo structures.
% - Set specific ReconInfo attributes.
for j = 1:128
    % @gpu: ReconInfo for frame k
    ReconInfo(j).mode = 'replaceIntensity';
    ReconInfo(j).txnum = j;
    ReconInfo(j).rcvnum = 128*Resource.RcvBuffer(1).numFrames + j;
    ReconInfo(j).regionnum = j;
end

% Specify Process structure array.
pers = 20;

Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
    'framenum',-1,...   % @gpu: (-1 => lastFrame)
    'pdatanum',1,...    % number of PData structure to use
    'pgain',1.0,...     % pgain is image processing gain
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

%@gpu: this is the external function that is used to FIR filter the data
Process(2).classname = 'External';
Process(2).method = 'externalFIR_CPU_M'; %@gpu: start with this external function
Process(2).Parameters = {'srcbuffer','receive',...
    'srcbufnum', 1,...
    'srcframenum', -1,...  %@gpu: auto increment the source frame
    'dstbuffer', 'receive',...
    'dstbufnum', 4,...    %@gpu: use rcvbuffer 4
    'dstframenum', 1, ... %@gpu: use frame 1 only
    };

Process(3).classname = 'External';
Process(3).method = 'externalFIR_GPU_PCT'; %@gpu: start with this external function
Process(3).Parameters = {'srcbuffer','receive',...
    'srcbufnum', 1,...
    'srcframenum', -1,...   %@gpu: use last frame transferred to buffer
    'dstbuffer', 'receive',...
    'dstbufnum', 4,...      %@gpu: use rcvbuffer 4
    'dstframenum', 1, ...   %@gpu: use frame 1 only
    };

Process(4).classname = 'External';
Process(4).method = 'externalFIR_GPU_MEXCUDA'; %@gpu: start with this external function
Process(4).Parameters = {'srcbuffer','receive',...
    'srcbufnum', 1,...
    'srcframenum', -1,...   %@gpu: use last frame transferred to buffer
    'dstbuffer', 'receive',...
    'dstbufnum', 4,...      %@gpu: use rcvbuffer 4
    'dstframenum', 1, ...   %@gpu: use frame 1 only
    };

% Specify SeqControl structure arrays.
%  - Time between acquisitions in usec
t1 = round(2*384*(1/Trans.frequency)); % acq. time in usec for max depth
SeqControl(1).command = 'timeToNextAcq';
SeqControl(1).argument = t1;
%  - Time between frames at 14 fps at max endDepth.
SeqControl(2).command = 'timeToNextAcq';
SeqControl(2).argument = 70000 - t1*127;
%  - Return to Matlab
SeqControl(3).command = 'returnToMatlab';
%  - Jump back to start.
SeqControl(4).command = 'jump';
SeqControl(4).argument = 1;
nsc = 5; % next SeqControl number

% Specify Event structure arrays.
externalProcessEvents=[];
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    for j = 1:128                      % Acquire frame
        Event(n).info = 'Acquisition.';
        Event(n).tx = j;
        Event(n).rcv = 128*(i-1)+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 1;
        n = n+1;
    end
    % Replace last events SeqControl for inter-frame timeToNextAcq.
    Event(n-1).seqControl = [2, nsc];
    SeqControl(nsc).command = 'transferToHost'; % @gpu: transfer frame to host buffer
    nsc = nsc+1;
    
    externalProcessEvents(end+1) = n;  %@gpu: keep track of which events contain the external process for switching from GUI
    Event(n).info = 'External Processing on RF';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 2;
    Event(n).seqControl = 0;
    n = n + 1;
    
    Event(n).info = 'recon and process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 1;
    Event(n).seqControl = 0;
    if floor(i/3) == i/3     % Exit to Matlab every 3rd frame
        Event(n).seqControl = nsc;
        SeqControl(nsc).command = 'returnToMatlab';
        nsc = nsc+1;
    end
    n = n+1;
end

Event(n).info = 'Jump back';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 4;

% User specified UI Control Elements
% @gpu: External Process Selection UI
import vsv.seq.uicontrol.VsButtonGroupControl;
UI(1).Control = VsButtonGroupControl('LocationCode','UserB5',...
    'Title','External Processing', ...
    'NumButtons', 3, ...
    'PossibleCases', {'CPU .m','GPU PCT', 'GPU .mexcuda'}, ...
    'Callback', @ExternalProcessMethodCallback);

UI(2).Control = VsButtonGroupControl('LocationCode','UserB7',...
    'Title','Rcv Source Type', ...
    'NumButtons', 3, ...
    'PossibleCases', {'Normal', 'CPU', 'Unified'}, ...
    'Callback', @RcvBufferNumCallback);


% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 3;

% Save all the structures to a .mat file.
save('MatFiles/L11-5v_128RyLns_GPU_RF');

return

%-ChangeExternalMethodCallback
function ExternalProcessMethodCallback(~, ~, UIState)
%@GPU: This is the callback for the UI to quickly switch between different
%external functions to compare their framerates without having to restart
%the sequence each time.
externalIndex = [];
switch UIState
    case 1
        disp('Switch to externalFIR_CPU_M')
        externalIndex = 2;
        % if it was in gpu, switch to 'normal' (1)
        RcvBufferNum = 1;
        set(findobj('tag','UserB7RadioButton1'), 'Value', 1);
        RcvBufferNumCallback([], [], RcvBufferNum);
    case 2
        v = ver;
        if (license('test','Distrib_Computing_Toolbox')  && any(strcmp('Parallel Computing Toolbox', {v.Name})) && (gpuDeviceCount()>0))
            disp('Switch to externalFIR_GPU_PCT')
            externalIndex = 3;
        else
            warning('Matlab Parallel Computing Toolbox needs to be installed to be able to test this function')
        end
    case 3
        if ~isempty(which('externalFIR_GPU_MEXCUDA'))
            disp('Switch to externalFIR_GPU_MEXCUDA')
            externalIndex = 4;
        else
            warning('externalFIR_CPU_MEX.mex does not exist or is not in the path.  Please type ''compileExampleFIRCode'' from the command line to generate the mex file')
        end
end

if ~isempty(externalIndex)
    Event=evalin('base','Event');
    externalProcessEvents = evalin('base','externalProcessEvents');
    for i = 1:length(externalProcessEvents)
        Event(externalProcessEvents(i)).process = externalIndex;  %replace external events with new process index
    end
    assignin('base','Event', Event);
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'Event'};
    assignin('base','Control', Control);
end

end


%-ChangeExternalMethodCallback
function RcvBufferNumCallback(~, ~, UIState)
%@GPU: This is the callback for the UI to quickly switch between different
%
switch UIState
    case 1
        disp('Switch Receive Buffer to Normal')
        RcvBufferNum = 1;
        
    case 2
        disp('Switch Receive Buffer to cuda-CPU')
        RcvBufferNum = 2;
        
    case 3
        disp('Switch Receive Buffer to cuda-Unified')
        RcvBufferNum = 3;

end

Process = evalin('base','Process');
for i = 2:4
    Process(i).Parameters{4} = RcvBufferNum;
end
assignin('base','Process', Process);

Receive = evalin('base','Receive');
numFrames = evalin('base','Resource.RcvBuffer(1).numFrames');
for i = 1:numFrames
    for j = 1:128
        Receive(128*(i-1) + j).bufnum = RcvBufferNum; %@gpu
    end
end
assignin('base','Receive', Receive);

Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'Process','Receive'};
assignin('base','Control', Control);

end
