% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use.
%
% File name SetUpL7_4PushTest.m:
% Description:
% Generate .mat Sequence Object file for L7-4 Linear array flash imaging,
% with a transition to TPC profile 5 for a single Push transmit-only event
% using the Extended Transmit Option.
% 128 transmit channels and 128 receive channels are used.
%
% Push % Comment lines with this prefix identify the changes added to allow
% the script to exercise a TPC Profile 5 Push transmit-only event.  Other
% than these changes, this is a typical flash imaging example script.
%
% NOTE: The Push transmit burst could be destructive to the % transducer!!
% It should be run with L7-4 or no transducer connected at all.
% The purpose of this script is only to demonstrate how to define a script for Push transmit
%
% Last update:
% April 2019 VTS-1142 IQBuffer replaced with separate I and Q Buffers
% 05/06/2020 - Update to SW 4.3 format for new UIControls and function definitions (VTS 1691). 
%   More info:(.../Example_Scripts/Vantage_Features/New UI Scheme/SetUpL11_5vFlash_NewUI)

clear all

P.startDepth = 5;   % Acquisition depth in wavelengths
P.endDepth = 192;   % This should preferrably be a multiple of 128 samples.

filename = ('L7-4PushTest');
na = 50;      % Set na = number of detect acquisitions.

% Define system parameters.
Resource.Parameters.numTransmit = 128;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;    % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Specify Trans structure array.
Trans.name = 'L7-4';
Trans.units = 'wavelengths';
Trans.frequency = 5.208; % nominal frequency in megahertz
Trans = computeTrans(Trans);  % L7-4 transducer is 'known' transducer so we can use computeTrans.
Trans.maxHighVoltage = 70;  % set maximum high voltage limit for pulser supply.

% Specify PData structure array.
PData(1).PDelta = [Trans.spacing, 0, 0.5];
PData.Size(1) = ceil((P.endDepth-P.startDepth)/PData.PDelta(3)); % startDepth, endDepth and pdelta set PData.Size.
PData.Size(2) = ceil((64*Trans.spacing)/PData.PDelta(1));
PData.Size(3) = 1;      % single image page
PData.Origin = [-Trans.spacing*31.5,0,P.startDepth]; % x,y,z of upper lft crnr.

% Specify Media object. 'pt1.m' script defines array of point targets.
pt1;
Media.function = 'movePoints';

% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = na*4096;
Resource.RcvBuffer(1).colsPerFrame = 64;  % RcvBuffer is 64 cols using syn aper.
Resource.RcvBuffer(1).numFrames = 10;     %
Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).numFrames = 1;  % one intermediate buffer needed.
Resource.InterBuffer(1).pagesPerFrame = na;
Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).numFrames = 10;
Resource.DisplayWindow(1).Title = filename;
Resource.DisplayWindow(1).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData.Size(2)*PData.PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData.Size(1)*PData.PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData.Origin(1),0,PData.Origin(3)]; % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).Colormap = gray(256);

% Specify Transmit waveform structure.
% - detect waveform
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,0.67,2,1];

% Push % A separate transmit waveform, TW(2), is defined for the Push
% transmit. This waveform will be used with TPC profile 5 and the Extended
% Burst transmit power supply, allowing high-power bursts of up to several
% milliseconds duration.
TW(2).type = 'parametric';
TW(2).Parameters = [Trans.frequency,1,2e4,1];

% Set TPC profile 5 high voltage limit.
TPC(5).maxHighVoltage = 70;

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'focus', 0.0, ...
                   'Steer', [0.0,0.0], ...
                   'Apod', ones(1,Trans.numelements), ...
                   'Delay', zeros(1,Trans.numelements)), 1, 2);

% Push % TX(2) is defined for the Push transmit, using the Push transmit
% waveform defined in TW(2).  Transmit focus control and TX.Apod aperture
% control both function the same for a Push burst as for any imaging burst.
TX(2).waveform = 2;
TX(2).focus = 50.0;
TX(2).Apod = [zeros(1,40),ones(1,48),zeros(1,40)];
TX(2).Delay = computeTXDelays(TX(2));

% Specify Receive structure arrays.
% -- Compute the maximum receive path length, using the law of cosines.
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
wlsPer128 = 128/(4*2); % wavelengths in 128 samples for 4 samplesPerWave
Receive = repmat(struct('Apod', zeros(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', P.startDepth + wlsPer128*ceil(maxAcqLength/wlsPer128), ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'mode', 0, ...
                        'InputFilter', repmat([0.0036,0.0127,0.0066,-0.0881,-0.2595,0.6494],Resource.Parameters.numRcvChannels,1), ...
                        'callMediaFunc', 0), 1, na*Resource.RcvBuffer(1).numFrames);
% - Set event specific Receive attributes for each frame.
for i = 1:Resource.RcvBuffer(1).numFrames
    Receive(na*(i-1)+1).callMediaFunc = 1;
    for j = 1:na  % na acquisitions per frame
        Receive(na*(i-1)+j).Apod(33:96) = 1.0;
        Receive(na*(i-1)+j).framenum = i;
        Receive(na*(i-1)+j).acqNum = j;
    end
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [500,590,650,710,770,800,850,950];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Recon structure arrays.
% - We need a Recon structure for the 2D image which will be used for each frame.
Recon(1) = struct(...
               'senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame',-1, ...
               'IntBufDest', [0,0], ...
               'ImgBufDest', [1,-1], ...
               'RINums',1);
Recon(2) = struct(...
               'senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame',-1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [0,0], ...
               'RINums',2:(na+1));

% Define ReconInfo structures.
% - ReconInfo for 2D frame.
ReconInfo = struct('mode','replaceIntensity', ...  % intensity output.
                   'txnum',1, ...
                   'rcvnum',1, ...
                   'regionnum',1);
k = 1; % k keeps track of index of last ReconInfo defined
% We need na ReconInfo structures for IQ reconstructions.
ReconInfo((k+1):(k+na)) = repmat(struct('mode', 'replaceIQ', ... % IQ output
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 1), 1, na);
% - Set specific ReconInfo attributes.
for j = 1:na  % For each row in the column
    ReconInfo(k+j).txnum = 1;
    ReconInfo(k+j).rcvnum = j;
    ReconInfo(k+j).pagenum = j;
end

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

Process(2).classname = 'External';
Process(2).method = 'processIQFunction';
Process(2).Parameters = {'srcbuffer','inter',... % name of buffer to process.
                         'srcbufnum',1,...
                         'srcframenum',1,...
                         'dstbuffer','none'};

% Specify SeqControl structure arrays.
% - Change to Profile 5 (high power)
SeqControl(1).command = 'setTPCProfile';
SeqControl(1).condition = 'immediate';
SeqControl(1).argument = 5;
% - Noop to allow time for charging external cap.
SeqControl(2).command = 'noop';
SeqControl(2).argument = 500000; % wait 100 msec.
% - time between push and detect acquisitions
SeqControl(3).command = 'timeToNextAcq';
SeqControl(3).argument = 500;               % 500usec
afterpush = 3;
% - time between detect acquisitions
SeqControl(4).command = 'timeToNextAcq';
SeqControl(4).argument = 100;               % 100usec (10KHz)
PRF=4;
% - time between frames
SeqControl(5).command = 'timeToNextEB';    % set time between extended bursts
% SeqControl(5).command = 'timeToNextAcq';    % set time between extended bursts
SeqControl(5).argument = 200000;            % 200000usec = 200msec (~ 5 fps)
TTNEB=5;
% - Return to Matlab
SeqControl(6).command = 'returnToMatlab';
% - Jump back to start.
SeqControl(7).command = 'jump';
SeqControl(7).argument = 1;
% - Trigger out
SeqControl(8).command = 'triggerOut';

nsc = 9;

% Specify Event structure arrays.
n = 1;

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
Event(n).info = 'noop for charging ext. cap.';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 2;
n = n+1;

for i = 1:Resource.RcvBuffer(1).numFrames
    % Push % Here we insert the Push transmit-only event, with a trigger out
    % command to make it convenient to trigger an oscilloscope for examining
    % the actual Push transmit output.  At the end of this event the transition
    % back to TPC profile 1 will be started, for the next imaging acquisition
    % sequence.  Another 40 msec timeToNextAcq delay will allow time for the
    % transition to be completed.
    Event(n).info = 'Switch to profile 5.';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = 1;
    n = n+1;

    % Push transmit
    Event(n).info = 'Push transmit';
    Event(n).tx = 2;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = [8,afterpush, TTNEB];
    n = n+1;

    for j = 1:na                      % Acquire frame
        Event(n).info = 'Acquire data';
        Event(n).tx = 1;
        Event(n).rcv = na*(i-1)+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = PRF;
        n = n+1;
    end

    Event(n-1).seqControl = 0; % do not want a TTNA here, since the next transmit is a EB

    Event(n).info = 'transfer data to Host';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = nsc;
      SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
      nsc = nsc+1;
    n = n+1;

    Event(n).info = 'recon and process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = [1,2];
    Event(n).process = 1;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = 'ext func to process IQ';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 2;
    Event(n).seqControl = 6;
    n = n+1;

end

% Comment out the Event below to pause at the end of 10 pushes.
Event(n).info = 'Jump back';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 7;


% User specified UI Control Elements

% - Sensitivity Cutoff
UI(1).Control = vsv.seq.uicontrol.VsSliderControl('LocationCode','UserB7','Label','Sens. Cutoff',...
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
UI(2).Control = vsv.seq.uicontrol.VsSliderControl('LocationCode','UserA1','Label',['Range (',AxesUnit,')'],...
                    'SliderMinMaxVal',MinMaxVal,...
                    'SliderStep',[0.1,0.2],'ValueFormat','%3.0f',...
                    'Callback', @RangeChangeCallback);
                
EF(1).Function = vsv.seq.function.ExFunctionDef('processIQFunction',@processIQFunction);

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 1;

% Save all the structures to a .mat file.
save(['MatFiles/',filename]);

% Uncomment the following line to automatically run VSX every time you run
% this SetUp script (note that if VSX finds the variable 'filename' in the
% Matlab workspace, it will load and run that file instead of prompting the
% user for the file to be used):

% VSX

return

%% **** Callback routines used by UIControls (UI) ****

function SensCutoffCallback(~,~,UIValue)
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

    
%% **** Callback routines used by External function definition (EF) ****

function processIQFunction(IBuffer, QBuffer)
    IQData = complex(IBuffer, QBuffer); % VTS-1142 separate I and Q buffers
    %processIQFunction: Computes power estimates from IQData
    %		Im = I(k) * Q(k+1) - I(k+1) * Q(k)
    %		Re = I(k) * I(k+1) + Q(k) * Q(k+1)
    %		Power = sqrt(Im*Im + Re*Re);

    persistent myHandle nf

    if isempty(myHandle)||~ishandle(myHandle)
        ImMean = (imag(IQData(1:150,:,1,3)) + imag(IQData(1:150,:,1,4)))/2;
        ReMean = (real(IQData(1:150,:,1,3)) + real(IQData(1:150,:,1,4)))/2;
        Im = (imag(IQData(1:150,:,1,3))-ImMean) .* (real(IQData(1:150,:,1,4))-ReMean) - ...
             (imag(IQData(1:150,:,1,4))-ImMean) .* (real(IQData(1:150,:,1,3))-ReMean);
        Re = (imag(IQData(1:150,:,1,3))-ImMean) .* (imag(IQData(1:150,:,1,4))-ImMean) + ...
             (real(IQData(1:150,:,1,3))-ReMean) .* (real(IQData(1:150,:,1,4))-ReMean);
        Power = (Im .* Im + Re .* Re).^0.125;
        figure;
        myHandle = imagesc(Power);
        colormap('gray');
        nf = 1;
    else
        for i = 3:49 % for all combinations of 2 pages
            ImMean = (imag(IQData(1:150,:,1,i)) + imag(IQData(1:150,:,1,i+1)))/2;
            ReMean = (real(IQData(1:150,:,1,i)) + real(IQData(1:150,:,1,i+1)))/2;
            Im = (imag(IQData(1:150,:,1,i))-ImMean) .* (real(IQData(1:150,:,1,i+1))-ReMean) - ...
                 (imag(IQData(1:150,:,1,i+1))-ImMean) .* (real(IQData(1:150,:,1,i))-ReMean);
            Re = (imag(IQData(1:150,:,1,i))-ImMean) .* (imag(IQData(1:150,:,1,i+1))-ImMean) + ...
                 (real(IQData(1:150,:,1,i))-ReMean) .* (real(IQData(1:150,:,1,i+1))-ReMean);
            Power = (Im .* Im + Re .* Re).^0.125;
            pause(0.05); % slow down the visulization of shearwave propagation
            set(myHandle,'CData',Power);
            drawnow
        end
        Resource = evalin('base','Resource');
        if (Resource.Parameters.simulateMode == 2)&&(nf<11)
            filename = ['IQData',num2str(nf)];
            IQDataS = squeeze(IQData);
            save(filename,'IQDataS');
            nf = nf + 1;
        end
    end
end

