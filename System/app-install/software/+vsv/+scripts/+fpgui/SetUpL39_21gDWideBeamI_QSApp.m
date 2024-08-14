% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use.
%
% File name: SetUpL39_21gDWideBeamI.m - Example of wide beam imaging
%                                       with interleaved sampling.
% Description:
%    Generate .mat file for L39_21gD linear array for Vantage system.
%    Transmit/Receive is performed with a 128 element aperture. The Tx
%    aperture is moved with each ray line and is centered around the beam
%    origin. The TX beam is focused just beyond the endDepth and forms a
%    triangular region, which is matched by the PData region specification.
%
%   This version uses asynchronous acquisition and processing. Interleaved
%   sampling is used to effectively sample RF data at 125MHz.
%
% Last update
%    April 2021 VTS-2152 computeTrans entry for L39-21gD
%    03/11/2021 Mods for L39-21gD with new lens
%    12/31/2020 Working with GE408 UTA

%clear[ 	]+all

P.startDepth = 0;  % Define startDepth and endDepth at top for use in defining other parameters.
P.endDepth = 200;
P.numRays = 128;   % number of rays in a single frame.
P.numTx = 50;      % number of active transmitters in TX aperture.
P.txFocus = P.endDepth+10;   % transmit focal pt in wavelengths

% Define system parameters.
Resource.Parameters.numTransmit = 128;  % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;  % number of receive channels.
Resource.Parameters.speedOfSound = 1540.;
Resource.Parameters.verbose = 2;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Specify Trans structure array.
Trans.name = 'L39-21gD';
Trans.units = 'mm';
Trans = computeTrans(Trans);

mm2wl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);

% Specify PData structure array.
PData.PDelta = [0.75*Trans.spacing,0,0.5];
PData.Size(1) = ceil((P.endDepth-5)/PData.PDelta(3)); % 5, P.endDepth and pdelta set PData.Size.
PData.Size(2) = ceil((Trans.numelements*Trans.spacing)/PData.PDelta(1));
PData.Size(3) = 1;      % single image page
PData.Origin = [-((Trans.numelements-1)/2)*Trans.spacing,0,5]; % x,y,z of upper lft crnr.
% Define P.numRays trpezoidal regions, positioned around the origin of the wide beams.
txdx = 127*Trans.spacing/(P.numRays-1);
TxOrgX = (-63.5*Trans.spacing):txdx:(63.5*Trans.spacing); % x coords of beam centers
for n = 1:P.numRays
    PData.Region(n) = struct('Shape',struct('Name','Trapezoid','Position',[TxOrgX(n),0,0],'top',...
                             P.numTx*Trans.spacing,'bottom',2*txdx,'height',P.endDepth));
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
Resource.RcvBuffer(1).numFrames = numRcvFrames;        % 20 frames in RcvBuffer
Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).numFrames = 1;  % one intermediate buffer needed.
Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).numFrames = 5;
Resource.DisplayWindow(1).Title = 'L39-21gD WideBeamI';
Resource.DisplayWindow(1).pdelta = 0.25;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData.Origin(1),PData.Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = gray(256);
Resource.DisplayWindow(1).numFrames = numImageFrames;

% Specify Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [31.25,0.67,2,1];

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'Apod', zeros(1,Resource.Parameters.numTransmit), ...
                   'focus', P.txFocus, ...
                   'Steer', [0.0,0.0], ...
                   'Delay', zeros(1,Resource.Parameters.numTransmit), ...
                   'TXPD', [],...
                   'peakCutOff', 0.8, ...
                   'peakBLMax', 10.0),1,2*P.numRays);

% - Set event specific TX attributes.
for n = 1:P.numRays  % specify P.numRays transmit events
    % Set transmit Origins.
    TX(n).Origin(1) = TxOrgX(n);
    % Compute available aperture
    [Dummy,ce] = min(abs(mm2wl*Trans.ElementPos(:,1)-TxOrgX(n))); % ce is closest ele to cntr of aper.
    lft = round(ce - P.numTx/2);
    if lft < 1, lft = 1; end
    rt = round(ce + P.numTx/2);
    if rt > 128, rt = 128; end
    TX(n).Apod(lft:rt) = 1;
    % Apply apodization function.
    [RIndices,CIndices,V] = find(TX(n).Apod);
    V = kaiser(size(V,2),1);
    TX(n).Apod(CIndices) = V;
    % Define 2nd TXs for interleave.
    TX(n+P.numRays) = TX(n);
    % Compute transmit delays for normal and interleave
    TX(n+P.numRays).Delay = computeTXDelays(TX(n+P.numRays));
    % Add the delay offset representing one half of A/D sample period (1/62.5MHz = 16ns) to first acq's delay.
    TX(n).Delay = TX(n+P.numRays).Delay + 0.5*Trans.frequency/62.5; % 0.25 wls of Trans.frequency=31.25 equivalent to 8 nsec.
    % Compute transmit pixel data only for the first TX in the interleave pair.
    TX(n).TXPD = computeTXPD(TX(n),PData);
    TX(n+P.numRays).TXPD = TX(n).TXPD;
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [200,524,798,891,925,965,1000,1023];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Receive structure arrays
RcvProfile.LnaZinSel = 31; % put LNA in "high-z" input state for best sensitivity

%  Specify Receive structure
%  - compute maxAcqLength to account for some channels having longer path lengths.
maxAcqLength = ceil(sqrt(P.endDepth^2 + (Trans.numelements*Trans.spacing)^2));
Receive = repmat(struct('Apod', ones(1,128), ...
                        'startDepth', 5, ...
                        'endDepth', maxAcqLength, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BWI', ...
                        'mode', 0, ...
                        'callMediaFunc', 0),1,2*P.numRays*Resource.RcvBuffer(1).numFrames);

% - Set event specific Receive attributes for each frame.
for i = 1:Resource.RcvBuffer(1).numFrames
    k = 2*P.numRays*(i-1);
    Receive(k+1).callMediaFunc = 1;
    for j = 1:P.numRays
        % first Receive of interleaved pair.
        Receive(k+2*j-1).framenum = i;
        Receive(k+2*j-1).acqNum = 2*j-1;
        % 2nd Receive of interleaved pair
        Receive(k+2*j).framenum = i;
        Receive(k+2*j).acqNum = 2*j;
    end
end

% Specify Recon structure arrays.
Recon = repmat(struct('senscutoff', 0.5, ...
               'pdatanum', 1, ...
               'rcvBufFrame',-1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...
               'RINums',1:P.numRays), 1, 1);

% Define ReconInfo structures.
ReconInfo = repmat(struct('mode','accumIQ', ...  % default is to accumulate IQ data.
               'Pre', [], ...
               'Post', [], ...
               'txnum', 1, ...
               'rcvnum', 1, ...
               'scaleFactor', 1.0, ...
               'regionnum', 0), 1, P.numRays);
% - Set specific ReconInfo attributes.
ReconInfo(1).Pre = 'clearInterBuf';
for j = 1:P.numRays
    ReconInfo(j).txnum = j;
    ReconInfo(j).rcvnum = 2*j-1;
    ReconInfo(j).regionnum = j;
end
ReconInfo(P.numRays).Post = 'IQ2IntensityImageBuf';

% Specify Process structure array.
pers = 20;
cmpFactor = 48;
rej = 2;
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',3,...       % pgain is image processing gain
                         'reject',rej,...
                         'persistMethod','simple',...
                         'persistLevel',pers,...
                         'interp','4pt',...
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
SeqControl(3).argument = 80;  % 80 usec
SeqControl(4).command = 'timeToNextAcq';
SeqControl(4).argument = 10000; %10000 usec = 10msec time between frames
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
    if (floor(i/4) == i/4)&&(i ~= Resource.RcvBuffer(1).numFrames) % Exit to Matlab every 4th frame
        if i~=Resource.RcvBuffer.numFrames, Event(n).seqControl = 5; end % return to Matlab
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
MinMaxVal = [64,200,P.endDepth]; % default unit is wavelength
AxesUnit = 'wls';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
        AxesUnit = 'mm';
        MinMaxVal = MinMaxVal/mm2wl;
    end
end
UI(2).Control = {'UserA1','Style','VsSlider','Label',['Range (',AxesUnit,')'],...
                 'SliderMinMaxVal',MinMaxVal,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f'};
UI(2).Callback = text2cell('%RangeChangeCallback');

% - TX Aperture change
UI(3).Control = {'UserB3','Style','VsSlider','Label','TX Aper',...
                 'SliderMinMaxVal',[1,128,P.numTx],'SliderStep',[1/128,1/64],'ValueFormat','%3.0f'};
UI(3).Callback = text2cell('%ApertureCallback');

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 4;
% Save all the structures to a .mat file.
filename = fullfile( vsv.file.getVSXDir, 'MatFiles/SetUpL39_21gDWideBeamI_QSApp.mat');

save(filename);

return

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
mm2wl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);
P = evalin('base','P');
P.endDepth = UIValue;
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
        P.endDepth = UIValue*mm2wl;
    end
end
P.txFocus = P.endDepth;

% PData
PData = evalin('base','PData');
PData(1).Size(1) = ceil((P.endDepth-5)/PData(1).PDelta(3)); % rows
for n = 1:size(PData(1).Region,2)
    PData(1).Region(n).Shape.height = P.endDepth;
end
PData.Region = computeRegions(PData);
assignin('base','PData',PData);
evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData.Size(1)*PData.PDelta(3)/Resource.DisplayWindow(1).pdelta);');

% Update TXPD data of TX structures.
TX = evalin('base','TX');
h = waitbar(0,'Program TX parameters, please wait!');
for n = 1:P.numRays
    TX(n).focus = P.txFocus;
    TX(n+P.numRays).focus = P.txFocus;
    % Compute transmit delays for normal and interleave
    TX(n+P.numRays).Delay = computeTXDelays(TX(n+P.numRays));
    % Add the delay offset representing one half of A/D sample period (1/62.5MHz = 16ns) to first acq's delay.
    TX(n).Delay = TX(n+P.numRays).Delay + 0.5*Trans.frequency/62.5; % 0.24 wls of Trans.frequency=30 equivalent to 8 nsec.
    TX(n).TXPD = computeTXPD(TX(n),PData);
    waitbar(n/P.numRays)
end
close(h)
assignin('base','TX',TX);

% Receive
Receive = evalin('base', 'Receive');
maxAcqLength = sqrt(P.endDepth^2 + (Trans.numelements*Trans.spacing)^2);
for i = 1:size(Receive,2)
    Receive(i).endDepth = maxAcqLength;
end
assignin('base','Receive',Receive);
evalin('base','TGC.rangeMax = P.endDepth;');
evalin('base','TGC.Waveform = computeTGCWaveform(TGC);');
Control.Command = 'update&Run';
Control.Parameters = {'PData','InterBuffer','ImageBuffer','DisplayWindow','TX','Receive','TGC','Recon'};
assignin('base','Control', Control);
assignin('base', 'action', 'displayChange');

return
%RangeChangeCallback

%ApertureCallback
simMode = evalin('base','Resource.Parameters.simulateMode');
% No focus change if in simulate mode 2.
if simMode == 2
    set(hObject,'Value',evalin('base','P.numTx'));
    return
end
P = evalin('base','P');
P.numTx = UIValue;
assignin('base','P',P);

Trans = evalin('base', 'Trans');
TX = evalin('base', 'TX');
PData = evalin('base','PData');
for n = 1:size(PData.Region,2)
    PData.Region(n).Shape.top = P.numTx*Trans.spacing;
end
PData.Region = computeRegions(PData);
assignin('base','PData',PData);

TxOrgX = evalin('base','TxOrgX');
h = waitbar(0,'Program TX parameters, please wait!');
% - Redefine event specific TX attributes for the new aperture.
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
    % Compute transmit pixel data only for the first TX in the interleave pair.
    TX(n).TXPD = computeTXPD(TX(n),PData);
    waitbar(n/P.numRays)
end
close(h)
assignin('base','TX', TX);
% Set Control command to update TX
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'PData','TX','Recon'};
assignin('base','Control', Control);
return
%ApertureCallback
