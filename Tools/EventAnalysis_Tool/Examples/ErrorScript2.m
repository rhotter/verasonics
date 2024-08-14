% ErrorScript2 has a bug about unused TX and Receive
% EventAnalysisTool shows a warning message if any TX, Rcv, and Seq is not
% listed in the event sequence. In some cases, if you have dummy TX,
% Receive, or SeqControl on purpose, please ignore the warning message

clear all

P.startDepth = 2;   % Acquisition depth in wavelengths
P.endDepth = 192;   % This should preferrably be a multiple of 128 samples.

na = 7;      % Set na = number of angles.
if (na > 1), dtheta = (36*pi/180)/(na-1); startAngle = -36*pi/180/2; else dtheta = 0; startAngle=0; end % set dtheta to range over +/- 18 degrees.

% Define system parameters.
Resource.Parameters.numTransmit = 128;  % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;  % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Specify Trans structure array.
Trans.name = 'L12-3v';
Trans.units = 'mm'; % Explicit declaration avoids warning message when selected by default
% note nominal center frequency in computeTrans is 7.813 MHz
Trans = computeTrans(Trans);  % L12-3v transducer is 'known' transducer so we can use computeTrans.

% Specify PData structure array.
PData.PDelta = [Trans.spacing, 0, 0.5];
PData.Size(1) = ceil((P.endDepth-P.startDepth)/PData.PDelta(3)); % startDepth, endDepth and pdelta set PData.Size.
PData.Size(2) = ceil((Trans.numelements*Trans.spacing)/PData.PDelta(1));
PData.Size(3) = 1;      % single image page
PData.Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of upper lft crnr.
% No PData.Region specified, so a default Region for the entire PData array will be created by computeRegions.

% Specify Media object. 'pt1.m' script defines array of point targets.
pt1;

% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 5120*2*na; % this size allows for 3 acqs, maximum range
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 10;        % 10 frames used for RF cineloop.
Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).numFrames = 1;  % one intermediate buffer needed.
Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).numFrames = 10;
Resource.DisplayWindow(1).Title = 'Error2';
Resource.DisplayWindow(1).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).numFrames = 20;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = gray(256);

% Specify Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [4.5,.67,2,1];   % A, B, C, D - transmit at 4.5MHz

TW(2).type = 'parametric';
TW(2).Parameters = [4.5,.67,2,-1];   % A, B, C, D

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'aperture', 1, ...
                   'Apod', ones(1,Resource.Parameters.numTransmit), ...
                   'focus', 0.0, ...
                   'Steer', [0.0,0.0], ...
                   'Delay', zeros(1,Resource.Parameters.numTransmit)), 1, 4*na);
% - Set event specific TX attributes.
for n = 1:2:2*na   % 2*na transmit events for NORMAL POLARITY
    angle = startAngle+((n+1)/2-1)*dtheta;
    %Normal Polarity for aperture 1
    TX(n).waveform = 1;
    TX(n).Steer = [angle,0.0];
    TX(n).aperture = 1; % Use the tx aperture that starts at element 1.
    TX(n).Delay = computeTXDelays(TX(n),'TOAE'); % use 'TransmitOnAllElements' flag
    %Normal Polarity for aperture 2
    TX(n+1).waveform = 1;
    TX(n+1).Steer = [angle,0.0];
    TX(n+1).aperture = 65; % Use the tx aperture that starts at element 65.
    TX(n+1).Delay = computeTXDelays(TX(n+1),'TOAE');
end
m=2*na;
for n = 1:2:2*na   % 2*na transmit events for INVERTED POLARITY
    %Inverted Polarity for aperture 1
    TX(m+n).waveform = 2;
    TX(m+n).Steer = TX(n).Steer;
    TX(m+n).aperture = TX(n).aperture; % Use the tx aperture that starts at element 1.
    TX(m+n).Delay = computeTXDelays(TX(m+n),'TOAE'); % use 'TransmitOnAllElements' flag
    %Inverted Polarity for aperture 2
    TX(m+n+1).waveform = 2;
    TX(m+n+1).Steer = TX(n+1).Steer;
    TX(m+n+1).aperture = TX(n+1).aperture; % Use the tx aperture that starts at element 65.
    TX(m+n+1).Delay = computeTXDelays(TX(m+n+1),'TOAE');
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [442,599,728,795,863,930,997,1023];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Receive structure arrays -
%   endDepth - add additional acquisition depth to account for some channels
%              having longer path lengths.
%   InputFilter - The same coefficients are used for all channels. The
%              coefficients below give a broad bandwidth bandpass filter.
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
Receive = repmat(struct('Apod', zeros(1,128), ...
                        'aperture', 1, ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'demodFrequency', 9.0, ...
                        'mode', 0, ...
                        'callMediaFunc', 0),1,4*na*Resource.RcvBuffer(1).numFrames);
% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames  % 6*na acquisitions per frame
    k = 4*na*(i-1);
    m = k + 2*na;
    Receive(k+1).callMediaFunc = 1;
    for j = 1:2:2*na;
        % -- 1st synthetic aperture acquisition for aperture 1 Normal Polarity.
        Receive(k+j).Apod(1:96) = 1.0;
        Receive(k+j).aperture = 1;
        Receive(k+j).framenum = i;
        Receive(k+j).acqNum = j;
        Receive(k+j).mode = 0;
        % -- 2nd synthetic aperture acquisition for aperture 2 Normal Polarity.
        Receive(k+j+1).Apod(33:128) = 1.0;
        Receive(k+j+1).aperture = 65;
        Receive(k+j+1).framenum = i;
        Receive(k+j+1).acqNum =j+1;
        Receive(k+j+1).mode = 0;
    end

    %Accumulate Receive Structure
    for j = 1:2:2*na
        % -- 1st synthetic aperture acquisition for aperture 1 Inverted Polarity.
        Receive(m+j).Apod(1:96) = 1.0;
        Receive(m+j).aperture = 1;
        Receive(m+j).framenum = i;
        Receive(m+j).acqNum = j;
        Receive(m+j).mode = 1;
        % -- 2nd synthetic aperture acquisition for aperture 2 Inverted Polarity.
        Receive(m+j+1).Apod(33:128) = 1.0;
        Receive(m+j+1).aperture = 65;
        Receive(m+j+1).framenum = i;
        Receive(m+j+1).acqNum =j+1;
        Receive(m+j+1).mode = 1;
    end
end

% Specify Recon structure arrays.
Recon = struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame', -1, ...     % use most recently transferred frame
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...  % auto-increment ImageBuffer each recon
               'RINums', 1:2*na);

% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 'accumIQ', ...  % default is to accumulate IQ data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 1), 1, 2*na);
% - Set specific ReconInfo attributes.
ReconInfo(1).mode = 'replaceIQ';
for j = 1:2:2*na
    ReconInfo(j).txnum = j;
    ReconInfo(j).rcvnum = j;
    ReconInfo(j+1).txnum = j+1;
    ReconInfo(j+1).rcvnum = j+1;
end
ReconInfo(2*na).mode = 'accumIQ_replaceIntensity';

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
                         'interpMethod','4pt',...  %method of interp. (1=4pt)
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
SeqControl(2).argument = 220;  % 220 usec
SeqControl(3).command = 'timeToNextAcq';  % time between frames
SeqControl(3).argument = 40000 - (4*na-1)*220;  % 40000 usec = 40 msec
SeqControl(4).command = 'returnToMatlab';
nsc = 5; % nsc is count of SeqControl objects

n = 1; % n is count of Events

% Acquire all frames defined in RcvBuffer
for i = 1:Resource.RcvBuffer(1).numFrames
    k = 4*na*(i-1);
    for j = 1:na  % Acquire all angles for first half of aperture
        Event(n).info = '1st aperture NORMAL POLARITY';
        Event(n).tx = j;
        Event(n).rcv = k+j;
        Event(n).recon = 0;      % no reconstruction.
        Event(n).process = 0;    % no processing
        Event(n).seqControl = 2; % time between syn. aper. acqs.
        n = n+1;

        Event(n).info = '1st aperture INVERTED POLARITY';
        Event(n).tx = 2*na +j;
        Event(n).rcv = 2*na +k+j;
        Event(n).recon = 0;      % no reconstruction.
        Event(n).process = 0;    % no processing
        Event(n).seqControl = 2; % time between syn. aper. acqs.
        n = n+1;
    end
    for j = 1:na  % Acquire all angles for second third of aperture
        Event(n).info = '2nd aperture NORMAL POLARITY';
        Event(n).tx = j+1;
        Event(n).rcv = k+j+1;
        Event(n).recon = 0;      % no reconstruction.
        Event(n).process = 0;    % no processing
        Event(n).seqControl = 2; % time between syn. aper. acqs.
        n = n+1;

        Event(n).info = '2nd aperture INVERTED POLARITY';
        Event(n).tx = 2*na +j+1;
        Event(n).rcv = 2*na +k+j+1;
        Event(n).recon = 0;      % no reconstruction.
        Event(n).process = 0;    % no processing
        Event(n).seqControl = 2; % time between syn. aper. acqs.
        n = n+1;
    end
    Event(n-1).seqControl = [nsc,3]; % use SeqControl struct defined below.
       SeqControl(nsc).command = 'transferToHost';
       nsc = nsc + 1;

    Event(n).info = 'Reconstruct & process';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 1;      % reconstruction
    Event(n).process = 1;    % processing
    Event(n).seqControl = 0;
    if floor(i/2) == i/2     % Exit to Matlab every 2nd frame
        Event(n).seqControl = 4;
    end
    n = n+1;
end

Event(n).info = 'Jump back to first event';
Event(n).tx = 0;        % no TX
Event(n).rcv = 0;       % no Rcv
Event(n).recon = 0;     % no Recon
Event(n).process = 0;
Event(n).seqControl = 1; % jump command


% User specified UI Control Elements
% - Sensitivity Cutoff
UI(1).Control =  {'UserB7','Style','VsSlider','Label','Sens. Cutoff',...
                  'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                  'SliderStep',[0.025,0.1],'ValueFormat','%1.3f'};
UI(1).Callback = text2cell('%-UI#1Callback');

% - Range Change
MinMaxVal = [64,300,P.endDepth]; % default unit is wavelength
AxesUnit = 'wls';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
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
filename = 'MatFiles/Error2';
save(filename);
% VSX;

return


% **** Callback routines to be converted by text2cell function. ****
%-UI#1Callback - Sensitivity cutoff change
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
%-UI#1Callback

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
