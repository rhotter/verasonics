% File name: SetUpL2048_hFlash_SimOnly.m - Hypothetical 2048-element linear array
% Flash imaging script for development and testing of Vantage 512 system
% software using one acquisition per frame of all 2048 elements.  This is a
% simulate-only script with no HW support and no HVMux switching, intended
% only to verify software processing functionalit with a 2048 elment array.
% 
% Last update
%    3/15/20 Script created

clear all

numEl = 2048; % number of elements and channels to use throughout the script

P.startDepth = 2;   % Acquisition depth in wavelengths
P.endDepth = 192;   % This should preferrably be a multiple of 128 samples.

% Define system parameters.
Resource.Parameters.numTransmit = numEl;  % number of system channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

Resource.System.Product = 'SimulateOnly';

%% Trans Structure definition
    Trans.name = 'L512tst';
    Trans.id = -1;
    Trans.units = 'wavelengths'; % Explicit declaration avoids warning message when selected by default
    Trans.frequency = 3; % nominal frequency in MHz
    Trans.Bandwidth = [2 4];
    Trans.type = 0;     % Array geometry is linear (x values only).
    Trans.connType = 0; % Simulate Only
    Trans.numelements = numEl;
    Trans.spacing = 0.2*512/numEl; % spacing in wvl
    Trans.elementWidth = 0.9*Trans.spacing; 
%     Trans.elementWidth = .1703/2; % width in mm 512 elements at half the pitch of L12-5 50mm
%     Trans.spacingMm = .1953/2;   % Spacing between elements in mm.
        Trans.elevationApertureMm = 7.5; % active elevation aperture in mm (estimate)
        Trans.elevationFocusMm = 20; % nominal elevation focus depth from lens on face of transducer (estimate)
    Trans.ElementPos = zeros(Trans.numelements,5);
%     Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
    Trans.ElementPos(:,1) = Trans.spacing*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
    Trans.lensCorrection = 0;
    Trans.impedance = 51; % using default value for MUX probe
    Trans.maxHighVoltage = 50;
    % define mapping from transducer elements to element signals at
    % connector:
    Trans.ConnectorES = [1:numEl]';
    % Fill in the wavelengths-based fields
%     scaleToWvl = Trans.frequency/1.540; % conversion factor from mm to wavelengths
%     Trans.spacing = Trans.spacingMm * scaleToWvl;   % Spacing between elements in wavelengths
    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
    Theta = (-pi/2:pi/100:pi/2);
    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
    % note at this point elementWidth is in mm, so we have to
    % convert to wavelengths for the ElementSens calculation
    eleWidthWl = Trans.elementWidth;
    if eleWidthWl < 0.01
        % avoid the divide by zero for very small values (in this
        % case the sinc function will be extremely close to 1.0 for
        % all Theta, so we only need the cos term)
        Trans.ElementSens = abs(cos(Theta));
    else
        Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
    end
%     if strcmp(Trans.units, 'wavelengths')
        % convert all mm unit variables to wavelengths.  Note columns 4
        % and 5 of the ElementPos array are angles in radians, and do
        % not require units conversion
%         Trans.elementWidth = Trans.elementWidth * scaleToWvl;
%         Trans.ElementPos(:,1) = Trans.ElementPos(:,1) * scaleToWvl;
%         Trans.ElementPos(:,2) = Trans.ElementPos(:,2) * scaleToWvl;
%         Trans.ElementPos(:,3) = Trans.ElementPos(:,3) * scaleToWvl;
%         Trans.lensCorrection = Trans.lensCorrection * scaleToWvl;
%     end
% End of Trans definition





% Specify PData structure array.
PData.PDelta = [0.5, 0, 0.5];
PData.Size(1) = ceil((P.endDepth-P.startDepth)/PData.PDelta(3)); % startDepth, endDepth and pdelta set PData.Size.
PData.Size(2) = 512;
PData.Size(3) = 1;      % single image page
PData.Origin = [-255.5/2,0,P.startDepth]; % x,y,z of upper lft crnr.
% No PData.Region specified, so a default Region for the entire PData array will be created by computeRegions.

% Specify Media object. 'pt1.m' script defines array of point targets.
pt1;
Media.attenuation = -0.5;
Media.function = 'movePoints';

% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 4096; % this size allows for maximum range
Resource.RcvBuffer(1).colsPerFrame = numEl;
Resource.RcvBuffer(1).numFrames = 40;        % 40 frames used for RF cineloop.
Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).numFrames = 10;
Resource.DisplayWindow(1).Title = 'L2048-hFlash_SimOnly';
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
Resource.DisplayWindow(1).Colormap = gray(256);

% Specify Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,.67,1,1];

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'Apod', ones(1,numEl), ...
                   'focus', 0.0, ...
                   'Steer', [0.0,0.0], ...
                   'Delay', zeros(1,numEl)), 1, 1);

% Specify TGC Waveform structure.
TGC.CntrlPts = [234,368,514,609,750,856,1023,1023];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Receive structure arrays -
%   endDepth - add additional acquisition depth to account for some channels
%              having longer path lengths.
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
Receive = repmat(struct('Apod', ones(1,numEl), ...
                        'startDepth', P.startDepth, ...
                        'endDepth',maxAcqLength, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 0),1,Resource.RcvBuffer(1).numFrames);
% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames  % 3 acquisitions per frame
    Receive(i).callMediaFunc = 1;
    % acquisition for frame i
    Receive(i).framenum = i;
end

% Specify Recon structure arrays.
Recon = struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame', -1, ...     % use most recently transferred frame
               'ImgBufDest', [1,-1], ...  % auto-increment ImageBuffer each recon
               'RINums', 1);

% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 'replaceIntensity', ...  % replace IQ data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'scaleFactor', 1, ...
                   'regionnum', 1), 1, 1);

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
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(2).argument = 200;  % 200 usec
SeqControl(3).command = 'timeToNextAcq';  % time between frames
SeqControl(3).argument = 20000 - 200;  % 20000 usec = 20 msec
SeqControl(4).command = 'returnToMatlab';
nsc = 5; % nsc is count of SeqControl objects

n = 1; % n is count of Events

% Acquire all frames defined in RcvBuffer
for i = 1:Resource.RcvBuffer(1).numFrames
    Event(n).info = '1st aperture.';
    Event(n).tx = 1;
    Event(n).rcv = i;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = [3,nsc];
       SeqControl(nsc).command = 'transferToHost';
       nsc = nsc + 1;
    n = n+1;

    Event(n).info = 'Reconstruct';
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
MinMaxVal = [64,200,P.endDepth]; % default unit is wavelength
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

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 2;

% Save all the structures to a .mat file.
save('MatFiles/L2048-hFlash_SimOnly');
return

% **** Callback routines to be converted by text2cell function. ****
%SensCutoffCallback - Sensitivity cutoff change
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

%RangeChangeCallback - Range change
simMode = evalin('base','Resource.Parameters.simulateMode');
% No range change if in simulate mode 2.
if simMode == 2
    set(hObject,'Value',evalin('base','P.endDepth'));
    return
end
Trans = evalin('base','Trans');
Resource = evalin('base','Resource');
% scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);

P = evalin('base','P');
P.endDepth = UIValue;
% if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
%     if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
%         P.endDepth = UIValue*scaleToWvl;
%     end
% end
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
