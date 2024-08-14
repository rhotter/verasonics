% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: SetUpL38-22v_128RyLnsI.m - Example of single focal zone line
%                                      line mode imaging w interleaved sampling.
% Description:
%    Generate .mat file for L38-22v linear array for Vantage system.
%    Transmit/Receive is performed with a 128 element mux aperture that is
%    moved across the 256 element aperture. The Tx aperture is moved with each
%    ray line, but the mux aperture is moved only for ray lines in the
%    central portion of the transducer aperture where the full 128 channel
%    aperture can be centered around the Tx aperture.
%
% In the following aperture examples, each space represents 4 elements.
% Aperture 1 (P.numTx=24):
%   ttt-------------------------------------------------------------
%   rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr--------------------------------
% Aperture for middle beam 128 (TX.aperture=65, P.numTx=24, numRay=128):
%   -----------------------------tttttt-----------------------------
%   ----------------rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr----------------
% Aperture for last wide beam (TX.aperture=129, P.numTx=24, numRay=256):
%   -------------------------------------------------------------ttt
%   --------------------------------rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
%
%   This version does asynchronous acquisition and processing. Interleaved
%   sampling is used to effectively sample RF data at 125MHz.
%
% Last update
%    11/01/2019

clear all

P.startDepth = 5;  % Define startDepth and endDepth at top for use in defining other parameters.
P.endDepth = 200;
P.numRays = 128;  % number of rays in a single frame.
P.numTx = 24;    % number of active transmitters in TX aperture.
P.txFocus = 100;   % transmit focal pt in wavelengths

% Define system parameters.
Resource.Parameters.numTransmit = 128;  % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;  % number of receive channels.
Resource.Parameters.speedOfSound = 1540.;
Resource.Parameters.verbose = 2;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Specify Trans structure array.
Trans.name = 'L38-22v';
Trans.units = 'wavelengths'; % required to prevent default to mm units
Trans = computeTrans (Trans);
Trans.name = 'Custom';  % only needed for prototype without personality eprom
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);

% Specify PData structure array.
PData.PDelta = [Trans.spacing,0,0.5];
PData.Size(1) = ceil((P.endDepth-P.startDepth)/PData.PDelta(3)); % P.startDepth, P.endDepth and pdelta set PData.Size.
PData.Size(2) = ceil((Trans.numelements*Trans.spacing)/PData.PDelta(1));
PData.Size(3) = 1;      % single image page
PData.Origin = [-((Trans.numelements-1)/2)*Trans.spacing,0,P.startDepth]; % x,y,z of upper lft crnr.
% Define P.numRays rectangular regions, positioned around the origin of the wide beams.
txdx = 255*Trans.spacing/(P.numRays-1);
TxOrgX = (-127.5*Trans.spacing):txdx:(127.5*Trans.spacing); % x coords of beam centers
for n = 1:P.numRays
    PData.Region(n) = struct('Shape',struct('Name','Rectangle','Position',[TxOrgX(n),0,0],'width',txdx,'height',P.endDepth));
end
PData.Region = computeRegions(PData);

% Specify Media object. 'pt1.m' script defines array of point targets.
pt1;
Media.function = 'movePoints';
Media.attenuation = -0.5;

% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 2*P.numRays*400*8; % size allows for interleaved rays, with range of up to 400 wvlngths
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 20;        % 20 frames in RcvBuffer
Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).numFrames = 1;  % one intermediate buffer needed.
Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).numFrames = 10;
Resource.DisplayWindow(1).Title = 'L38-22v_128RyLnsI';
Resource.DisplayWindow(1).pdelta = 0.4;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData.Origin(1),0,PData.Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).numFrames = 200;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = gray(256);

% Specify Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,0.67,2,1];

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'aperture', 1, ...
                   'Apod', zeros(1,Resource.Parameters.numTransmit), ...
                   'focus', P.txFocus, ...
                   'Steer', [0.0,0.0], ...
                   'Delay', zeros(1,Resource.Parameters.numTransmit)), 1, 2*P.numRays);

% - Set event specific TX attributes.
for n = 1:P.numRays  % specify P.numRays transmit events
    % Set transmit Origins.
    TX(n).Origin(1) = TxOrgX(n);
    % Compute available transmit mux aperture
    [Dummy,ce] = min(abs(Trans.ElementPos(:,1)-TxOrgX(n))); % ce is closest ele to cntr of aper.
    lft = round(ce - 64);
    if lft < 1, lft = 1; end
    if lft > 129, lft = 129; end
    TX(n).aperture = lft;
    % Compute TX.Apod within mux aperture.
    lft = round(ce - P.numTx/2);
    if lft < 1, lft = 1; end
    rt = round(ce + P.numTx/2);
    if rt > 256, rt = 256; end
    TX(n).Apod((lft-(TX(n).aperture-1)):(rt-(TX(n).aperture-1))) = 1;
    % Apply apodization function.
    [RIndices,CIndices,V] = find(TX(n).Apod);
    V = kaiser(size(V,2),1);
    TX(n).Apod(CIndices) = V;
    % Define 2nd TXs for interleave.
    TX(n+P.numRays) = TX(n);
    % Compute transmit delays for normal and interleave
    TX(n+P.numRays).Delay = computeTXDelays(TX(n+P.numRays));
    % Add the delay offset representing one half of A/D sample period (1/62.5MHz = 16ns) to first acq's delay.
    TX(n).Delay = TX(n+P.numRays).Delay + 0.5*Trans.frequency/62.5; % 0.24 wls of Trans.frequency=30 equivalent to 8 nsec.
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [400,622,798,891,925,965,1000,1023];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Receive structure arrays
RcvProfile.LnaZinSel = 0;
%  - P.endDepth - add additional acquisition depth to account for some channels
%              having longer path lengths.
% Coefficients for high pass filter with 23MHz -3dB cutoff. We need a high pass
%   filter for interleaved sampling since the non-interleaved sample rate is half the
%   interleaved rate and the Nyquist frequency is near the transducer center frequency.
HPFilter = [-0.00128 +0.00189 -0.00061 -0.00256 +0.00537 -0.00443 -0.00165 ...
            +0.00961 -0.01251 +0.00488 +0.01114 -0.02441 +0.02124 +0.00348 ...
            -0.03760 +0.05469 -0.02856 -0.04800 +0.15350 -0.24527 +0.28162];
maxAcqLength = ceil(sqrt(P.endDepth^2 + (Trans.numelements*Trans.spacing)^2));
Receive = repmat(struct('aperture', 1, ...
                        'Apod', ones(1,128), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BWI', ...
                        'InputFilter', HPFilter, ...
                        'mode', 0, ...
                        'callMediaFunc', 0),1,2*P.numRays*Resource.RcvBuffer(1).numFrames);

% - Set event specific Receive attributes for each frame.
for i = 1:Resource.RcvBuffer(1).numFrames
    k = 2*P.numRays*(i-1);
    Receive(k+1).callMediaFunc = 1;
    for j = 1:P.numRays
        % first Receive of interleaved pair.
        Receive(k+2*j-1).aperture = TX(j).aperture; % mux aperture same as transmit
        Receive(k+2*j-1).framenum = i;
        Receive(k+2*j-1).acqNum = 2*j-1;
        % 2nd Receive of interleaved pair
        Receive(k+2*j).aperture = TX(j).aperture; % mux aperture same as transmit
        Receive(k+2*j).framenum = i;
        Receive(k+2*j).acqNum = 2*j;
    end
end

% Specify Recon structure arrays.
Recon = repmat(struct('senscutoff', 0.5, ...
               'pdatanum', 1, ...
               'rcvBufFrame',-1, ...
               'ImgBufDest', [1,-1], ...
               'RINums',1:P.numRays), 1, 1);

% Define ReconInfo structures.
ReconInfo = repmat(struct('mode','replaceIntensity', ...
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'scaleFactor', 1.0, ...
                   'regionnum', 0), 1, P.numRays);
% - Set specific ReconInfo attributes.
for j = 1:P.numRays
    ReconInfo(j).txnum = j;
    ReconInfo(j).rcvnum = 2*j-1;
    ReconInfo(j).regionnum = j;
end

% Specify Process structure array.
pers = 20;
cmpFactor = 45;
rej = 2;
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',1,...            % pgain is image processing gain
                         'reject',rej,...
                         'persistMethod','simple',...
                         'persistLevel',pers,...
                         'interp','4pt',...      % method of interpolation (1=4pt interp)
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',cmpFactor,...
                         'mappingMethod','full',...
                         'display',1,...      % display image after processing
                         'displayWindow',1};

% Specify SeqControl structure arrays.
%  - Jump back to start.
SeqControl(1).command = 'jump';
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(2).argument = 50;  % 50 usec
SeqControl(3).command = 'timeToNextAcq';  % time to next wide beam acquisition
SeqControl(3).argument = 100;  % 100 usec
SeqControl(4).command = 'timeToNextAcq';
SeqControl(4).argument = 20000; %20000 usec = 20msec time between frames
SeqControl(5).command = 'returnToMatlab';
nsc = 6; % nsc is count of SeqControl objects

% Specify Event structure arrays.
n = 1;
for i = 1:Resource.RcvBuffer.numFrames
    % Acquire frame
    k = 2*P.numRays*(i-1);
    for j = 1:P.numRays        % Acquire frame
        Event(n).info = 'acquire 1st interleave pair';
        Event(n).tx = j;   % use 1st set of TX structures.
        Event(n).rcv = k+2*j-1; % 1st set of Receives
        Event(n).recon = 0;      % no reconstruction.
        Event(n).process = 0;    % no processing
        Event(n).seqControl = 2;
        n = n+1;
        Event(n).info = 'acquire 2nd interleave pair';
        Event(n).tx = P.numRays+j; % use 2nd set of TX structures.
        Event(n).rcv = k+2*j; % 1st set of Receives
        Event(n).recon = 0;      % no reconstruction.
        Event(n).process = 0;    % no processing
        Event(n).seqControl = 3;
        n = n+1;
    end
    Event(n-1).seqControl = [4,nsc]; % modify last acquisition Event's seqControl
      SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
      nsc = nsc+1;

    Event(n).info = 'recon and process';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 1;      % reconstruction
    Event(n).process = 1;    % process
    Event(n).seqControl = 0;
    if (floor(i/4) == i/4)&&(i ~= Resource.RcvBuffer(1).numFrames) % Exit to Matlab every 3rd frame
        if (i~=Resource.RcvBuffer.numFrames)
            Event(n).seqControl = 5; % return to Matlab
        end
    end
    n = n+1;
end

Event(n).info = 'Jump back';
Event(n).tx = 0;        % no TX
Event(n).rcv = 0;       % no Rcv
Event(n).recon = 0;     % no Recon
Event(n).process = 0;
Event(n).seqControl = 1;


% User specified UI Control Elements
% - Sensitivity Cutoff
UI(1).Control =  {'UserB7','Style','VsSlider','Label','Sens. Cutoff',...
                  'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                  'SliderStep',[0.025,0.1],'ValueFormat','%1.3f'};
UI(1).Callback = text2cell('%SensCutOffCallback');

% - Range Change
wls2mm = 1;
AxesUnit = 'wls';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
        AxesUnit = 'mm';
        wls2mm = Resource.Parameters.speedOfSound/1000/Trans.frequency;
    end
end
UI(2).Control = {'UserA1','Style','VsSlider','Label',['Range (',AxesUnit,')'],...
                 'SliderMinMaxVal',[64,300,P.endDepth]*wls2mm,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f'};
UI(2).Callback = text2cell('%RangeChangeCallback');

% - Transmit focus change
UI(3).Control = {'UserB4','Style','VsSlider','Label',['TX Focus (',AxesUnit,')'],...
                 'SliderMinMaxVal',[20,320,P.txFocus]*wls2mm,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f'};

UI(3).Callback = text2cell('%TxFocusCallback');
% - TX Aperture change
UI(4).Control = {'UserB3','Style','VsSlider','Label','TX Aper',...
                 'SliderMinMaxVal',[1,128,P.numTx],'SliderStep',[1/128,1/64],'ValueFormat','%3.0f'};
UI(4).Callback = text2cell('%ApertureCallback');

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 4;
% Save all the structures to a .mat file.
save('MatFiles/L38-22v_128RyLnsI');

return
% =========================================================================
% **** Callback routines to be converted by text2cell function. ****
%SensCutOffCallback - Sensitivity cutoff change
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
%SensCutOffCallback

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
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
        P.endDepth = UIValue*scaleToWvl;
    end
end
assignin('base','P',P);

PData = evalin('base','PData');
PData.Size(1) = ceil((P.endDepth-P.startDepth)/PData.PDelta(3)); % P.startDepth, P.endDepth and pdelta set PData.Size(1).
for i = 1:size(PData.Region,2)
    PData.Region(i).Shape.height = P.endDepth;
end
PData.Region = computeRegions(PData);
assignin('base','PData',PData);
evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData.Size(1)*PData.PDelta(3)/Resource.DisplayWindow(1).pdelta);');
Receive = evalin('base', 'Receive');
maxAcqLength = sqrt(P.endDepth^2 + (Trans.numelements*Trans.spacing)^2);
for i = 1:size(Receive,2)
    Receive(i).endDepth = maxAcqLength;
end
assignin('base','Receive',Receive);
evalin('base','TGC.rangeMax = P.endDepth;');
evalin('base','TGC.Waveform = computeTGCWaveform(TGC);');
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'PData','TX','Receive','Recon','InterBuffer','ImageBuffer','DisplayWindow'};
assignin('base','Control', Control);
assignin('base', 'action', 'displayChange');
return
%RangeChangeCallback

%TxFocusCallback
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
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
        P.txFocus = UIValue*scaleToWvl;
    end
end
assignin('base','P',P);
% - Redefine event specific TX attributes for the new focus.
TX = evalin('base', 'TX');
PData = evalin('base','PData');
m = P.numRays;
for n = 1:m
    % write new focus value to TX
    TX(n+m).focus = P.txFocus;
    TX(n+m).Delay = computeTXDelays(TX(n+m));
    % Compute transmit pixel data
    TX(n).focus = P.txFocus;
    TX(n).Delay = TX(n+m).Delay + 0.5*Trans.frequency/62.5;
    % Compute transmit pixel data
end
assignin('base','TX', TX);
% Set Control command to update TX
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'TX'};
assignin('base','Control', Control);
return
%TxFocusCallback

%ApertureCallback
simMode = evalin('base','Resource.Parameters.simulateMode');
% No focus change if in simulate mode 2.
if simMode == 2
    set(hObject,'Value',evalin('base','P.numTx'));
    return
end
Trans = evalin('base', 'Trans');
P = evalin('base','P');
P.numTx = UIValue;
assignin('base','P',P);

TX = evalin('base', 'TX');
TxOrgX = evalin('base','TxOrgX');
PData = evalin('base','PData');
scaleToWvl = evalin('base','scaleToWvl');

% - Redefine event specific TX attributes for the new aperture.
% - Set event specific TX attributes.

for n = 1:P.numRays  % specify P.numRays transmit events
    TX(n).Apod = zeros(1,128);
    % Set transmit Origins.
    TX(n).Origin(1) = TxOrgX(n);
    % Compute available transmit mux aperture
    [Dummy,ce] = min(abs(Trans.ElementPos(:,1)-TxOrgX(n))); % ce is closest ele to cntr of aper.
    lft = round(ce - 64);
    if lft < 1, lft = 1; end
    if lft > 129, lft = 129; end
    TX(n).aperture = lft;
    % Compute TX.Apod within mux aperture.
    lft = round(ce - P.numTx/2);
    if lft < 1, lft = 1; end
    rt = round(ce + P.numTx/2);
    if rt > 256, rt = 256; end
    TX(n).Apod((lft-(TX(n).aperture-1)):(rt-(TX(n).aperture-1))) = 1;
    % Apply apodization function.
    [RIndices,CIndices,V] = find(TX(n).Apod);
    V = kaiser(size(V,2),1);
    TX(n).Apod(CIndices) = V;
    % Define 2nd TXs for interleave.
    TX(n+P.numRays) = TX(n);
    % Compute transmit delays for normal and interleave
    TX(n+P.numRays).Delay = computeTXDelays(TX(n+P.numRays));
    % Add the delay offset representing one half of A/D sample period (1/62.5MHz = 16ns) to first acq's delay.
    TX(n).Delay = TX(n+P.numRays).Delay + 0.5*Trans.frequency/62.5; % 0.24 wls of Trans.frequency=30 equivalent to 8 nsec.
end
assignin('base','TX', TX);
% Set Control command to update TX
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'TX'};
assignin('base','Control', Control);
return
%ApertureCallback
