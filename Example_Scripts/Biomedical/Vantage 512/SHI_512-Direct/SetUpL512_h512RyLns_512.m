% File name: SetUpL512_h512RyLns_512.m - Hypothetical 512-element linear array
% Ray line imaging script for development and testing of Vantage 512 system
% using 512 ray line acquisitions per frame
% This initial version demonstrates the programming techniques needed to
% run on the Vantage 512 hardware system but at this point can only be run
% in "software only" simulation mode.  When the Vantage 512 hardware
% configuration is actually available, some minor adjustments to this
% script may be needed for full compatiblility with the actual hardware
% system.

% Last update
%    10/19/2020 Initial version

clear all

numEl = 512; % number of elements and channels to use throughout the script

numRyLns = numEl; % number of ray line acquisitions per frame

P.startDepth = 2;   % Acquisition depth in wavelengths
P.endDepth = 192;   % This should preferrably be a multiple of 128 samples.
P.txFocus = 100;  % Initial transmit focus.
P.numTx = 64;  % number of transmit elements in TX aperture (where possible).

% Define system parameters.
Resource.Parameters.numTransmit = numEl;  % number of system channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

%% Trans Structure definition
    Trans.name = 'L512tst';
    Trans.id = -1; % Minus 1 means the probe has no ID code
    Trans.units = 'wavelengths'; % Explicit declaration avoids warning message when selected by default
    Trans.frequency = 7.813; % nominal frequency in MHz
    Trans.Bandwidth = [5, 11];
    Trans.type = 0;     % Array geometry is linear (x values only).
    Trans.connType = 13; % V-512 Captive SHI connector type
    Trans.numelements = numEl;
    Trans.elementWidth = .1703/2; % width in mm 512 elements at half the pitch of L12-5 50mm
    Trans.spacingMm = .1953/2;   % Spacing between elements in mm.
        Trans.elevationApertureMm = 7.5; % active elevation aperture in mm (estimate)
        Trans.elevationFocusMm = 20; % nominal elevation focus depth from lens on face of transducer (estimate)
    Trans.ElementPos = zeros(Trans.numelements,5);
    Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
    Trans.lensCorrection = 2.365; % in mm units; was 12 wavelengths;
    Trans.impedance = 51; % using default value for MUX probe
    Trans.maxHighVoltage = 50;
    % define mapping from transducer elements to element signals at
    % connector:
    Trans.ConnectorES = [1:numEl]';
    % Fill in the wavelengths-based fields
    scaleToWvl = Trans.frequency/1.540; % conversion factor from mm to wavelengths
    Trans.spacing = Trans.spacingMm * scaleToWvl;   % Spacing between elements in wavelengths
    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
    Theta = (-pi/2:pi/100:pi/2);
    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
    % note at this point elementWidth is in mm, so we have to
    % convert to wavelengths for the ElementSens calculation
    eleWidthWl = Trans.elementWidth * scaleToWvl;
    if eleWidthWl < 0.01
        % avoid the divide by zero for very small values (in this
        % case the sinc function will be extremely close to 1.0 for
        % all Theta, so we only need the cos term)
        Trans.ElementSens = abs(cos(Theta));
    else
        Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
    end
    if strcmp(Trans.units, 'wavelengths')
        % convert all mm unit variables to wavelengths.  Note columns 4
        % and 5 of the ElementPos array are angles in radians, and do
        % not require units conversion
        Trans.elementWidth = Trans.elementWidth * scaleToWvl;
        Trans.ElementPos(:,1) = Trans.ElementPos(:,1) * scaleToWvl;
        Trans.ElementPos(:,2) = Trans.ElementPos(:,2) * scaleToWvl;
        Trans.ElementPos(:,3) = Trans.ElementPos(:,3) * scaleToWvl;
        Trans.lensCorrection = Trans.lensCorrection * scaleToWvl;
    end
% End of Trans definition





% Specify PData structure array.
PData.PDelta = [Trans.spacing, 0, 0.5];
PData.Size(1) = ceil((P.endDepth-P.startDepth)/PData.PDelta(3)); % startDepth, endDepth and pdelta set PData.Size.
PData.Size(2) = ceil((numEl*Trans.spacing)/PData.PDelta(1));
PData.Size(3) = 1;      % single image page
PData.Origin = [-Trans.spacing*(numEl-1)/2,0,P.startDepth]; % x,y,z of upper lft crnr.
% - specify numRyLns Region structures.
PData(1).Region = repmat(struct('Shape',struct( ...
                    'Name','Rectangle',...
                    'Position',[0,0,P.startDepth],...
                    'width',Trans.spacing,...
                    'height',P.endDepth-P.startDepth)),1,numRyLns);
% - set position of regions to correspond to beam spacing.
for i = 1:numRyLns
    PData(1).Region(i).Shape.Position(1) = (-(numEl-1)/2 + (i-1))*Trans.spacing;
end
PData(1).Region = computeRegions(PData(1));

% Specify Media object. 'pt1.m' script defines array of point targets.
pt1;
Media.attenuation = -0.5;
Media.function = 'movePoints';

% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = numRyLns*4096; % this size allows for maximum range
Resource.RcvBuffer(1).colsPerFrame = numEl;
Resource.RcvBuffer(1).numFrames = 4;        % 4 frames used for RF cineloop.
Resource.InterBuffer(1).numFrames = 1;  % one intermediate buffer needed.
Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).numFrames = 10;
Resource.DisplayWindow(1).Title = 'L512-h512RyLns_512';
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
TW(1).Parameters = [Trans.frequency,.67,1,1];

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'Apod', zeros(1,numEl), ...
                   'focus', P.txFocus, ...
                   'Steer', [0.0,0.0], ...
                   'Delay', zeros(1,numEl)), 1, numRyLns);
% - Set event specific TX attributes.
for n = 1:numRyLns   % numRyLns transmit events
    % Set transmit Origins to positions of elements.
    TX(n).Origin = [(-(numEl-1)/2 + (n-1))*Trans.spacing, 0.0, 0.0];
    % Set transmit Apodization.
    lft = n - floor(P.numTx/2);
    if lft < 1, lft = 1; end
    rt = n + floor(P.numTx/2);
    if rt > numEl, rt = numEl; end
    TX(n).Apod(lft:rt) = 1.0;
    TX(n).Delay = computeTXDelays(TX(n));
end

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
                        'callMediaFunc', 0),1,numRyLns*Resource.RcvBuffer(1).numFrames);
% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    k = numRyLns*(i-1);
    Receive(k+1).callMediaFunc = 1;
    for j = 1:numRyLns
        Receive(k+j).framenum = i;
        Receive(k+j).acqNum = j;
    end
end

% Specify Recon structure arrays.
Recon = struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame', -1, ...     % use most recently transferred frame
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...  % auto-increment ImageBuffer each recon
               'RINums', 1:numRyLns);

% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 'replaceIntensity', ...
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 0), 1, numRyLns);
% - Set specific ReconInfo attributes.
for j = 1:numRyLns
    ReconInfo(j).txnum = j;
    ReconInfo(j).rcvnum = j;
    ReconInfo(j).regionnum = j;
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

% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; % jump back to start.
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between ray line acquisitions
SeqControl(2).argument = 200;  % 200 usec
SeqControl(3).command = 'timeToNextAcq';  % time between frames
SeqControl(3).argument = 20000 - 200;  % 20000 usec = 20 msec
SeqControl(4).command = 'returnToMatlab';
nsc = 5; % nsc is count of SeqControl objects

n = 1; % n is count of Events

% Acquire all frames defined in RcvBuffer
for i = 1:Resource.RcvBuffer(1).numFrames
    for j = 1:numRyLns                      % Acquire frame
        Event(n).info = 'Acquisition.';
        Event(n).tx = j;
        Event(n).rcv = numRyLns*(i-1)+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;
    end
    % Replace last events SeqControl for inter-frame timeToNextAcq.
    Event(n-1).seqControl = 3;

    Event(n).info = 'Transfer frame to host.';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = nsc;
       SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
       nsc = nsc+1;
    n = n+1;

    Event(n).info = 'Reconstruct';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 1;
    Event(n).seqControl = 4; % Exit to Matlab every 5th frame
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
MinMaxVal = [64,300,P.endDepth]; % default unit is wavelength
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
save('MatFiles/L512-h512RyLns_512');
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
return
%RangeChangeCallback
