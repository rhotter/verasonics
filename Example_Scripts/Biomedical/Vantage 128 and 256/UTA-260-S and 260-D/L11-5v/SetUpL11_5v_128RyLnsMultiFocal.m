% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: SetUpL11_5v_128RyLnsMultiFocal.m - Imaging with dynamic focus transmits
%
% Description:
%   Sequence programming file for L11-5v Linear array, using 128 ray lines
%   (focus transmits) and 128 receive acquisitions. Of the 128 transmit
%   channels with 3 different focal zones, the active transmit aperture is
%   limited based on user-entered transmit focus and f-number. All 128
%   receive channels are active for each acquisition. The receive
%   acquisitions use 100% bandwidth to improve DMA transfers. Processing is
%   asynchronous with respect to acquisition.
%
% Last update:
% 11/10/2015 - modified for SW 3.0

clear all

P.startDepth = 5;
P.endDepth = 192;
P.numRays = 128;

nr = P.numRays;  % nr is number of raylines.

% Number of Focal Zones
fz = 3;

% Transmit focus in wavelengths
TXFocus(1) = 160;
TXFocus(2) = 120;
TXFocus(3) = 90;

% Receive zone endDepth in wavelengths
RcvZone(1) = 192;
RcvZone(2) = 130;
RcvZone(3) = 90;

% Specify system parameters.
Resource.Parameters.connector = 1;
Resource.Parameters.numTransmit = 128;  % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;  % number of receive channels.
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Specify Trans structure array.
Trans.name = 'L11-5v';
Trans.units = 'wavelengths'; % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);  % L11-5v transducer is 'known' transducer so we can use computeTrans.
Trans.maxHighVoltage = 50;  % set maximum high voltage limit for pulser supply.

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
pt1;
Media.attenuation = -0.5;
Media.function = 'movePoints';

% Specify Resources.
Resource.RcvBuffer.datatype = 'int16';
Resource.RcvBuffer.rowsPerFrame = 196608; % this should be larger than 128*Receive.endDepth*2(round trip)*2(smpls/wave) for max depth
Resource.RcvBuffer.colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer.numFrames = 10;
Resource.InterBuffer(1).numFrames = 1;  % one intermediate buffer needed.
Resource.ImageBuffer(1).numFrames = 10;
Resource.DisplayWindow.Title = 'L11-5v_128RyLnsMultiFocal';
Resource.DisplayWindow.pdelta = 0.3;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).numFrames = 20;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow.Colormap = gray(256);

% Specify TW structure array.
TW.type = 'parametric';
TW.Parameters = [Trans.frequency,0.67,2,1];

% Specify nr TX structure arrays. Transmit on element n in the array for event n.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'focus', TXFocus(1), ...
                   'Steer', [0.0,0.0], ...
                   'Apod', zeros(1,Trans.numelements), ...
                   'Delay', zeros(1,Trans.numelements)), 1, fz*nr);

% Determine TX aperture based on focal point and desired f number.
txFNum = 3.5;  % set to desired f-number value for transmit (range: 1.0 - 20)
txNumEl = zeros(1,fz);
for j = 1:fz
    txNumEl(j)=round((TXFocus(j)/txFNum)/Trans.spacing/2); % no. of elements in 1/2 aperture.
    if txNumEl(j) > (Trans.numelements/2 - 1), txNumEl(j) = floor(Trans.numelements/2 - 1); end
end
% txNumEl is the number of elements to include on each side of the
% center element, for the specified focus and sensitivity cutoff.
% Thus the full transmit aperture will be 2*txNumEl + 1 elements.
%display('Number of elements in transmit aperture:');
%disp(2*txNumEl+1);

% - Set event specific TX attributes.
for j = 1:128   % 128 transmit events
    k = fz*(j-1);
    for n = 1:fz
        % Set transmit Origins to positions of elements.
        TX(n+k).Origin = [(-63.5 + (j-1))*Trans.spacing, 0.0, 0.0];
        % Set transmit Apodization so (1 + 2*TXnumel) transmitters are active.
        lft = j - txNumEl(n);
        if lft < 1, lft = 1; end;
        rt = j + txNumEl(n);
        if rt > Trans.numelements, rt = Trans.numelements; end;
        TX(n+k).Apod(lft:rt) = 1.0;
        TX(n+k).focus = TXFocus(n);
        TX(n+k).Delay = computeTXDelays(TX(n+k));
    end
end

% Specify Receive structure arrays.
% - We need nr Receives for every frame.
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'BS100BW',...
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, fz*nr*Resource.RcvBuffer(1).numFrames);
% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    k = nr*fz*(i-1);
    Receive(k+1).callMediaFunc = 1;
    for j = 1:P.numRays
        w = fz*(j-1);
        for z = 1:fz
            %First half of aperture - Set max acq length for first TX in each set of focal zones
            if z == 1
                Receive(k+w+z).endDepth = maxAcqLength;
            else
                Receive(k+w+z).endDepth = RcvZone(z);
            end
            Receive(k+w+z).framenum = i;
            Receive(k+w+z).acqNum = j;      % two acquisitions per frame
        end
    end
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [0,137,149,287,333,644,827,1023];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Recon structure array.
Recon = struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame',-1, ...
               'IntBufDest', [0,0], ...
               'ImgBufDest', [1,-1], ...
               'RINums', 1:nr);

% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 'replaceIntensity', ...  % replace intensity data
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 0), 1, nr);
% - Set specific ReconInfo attributes.
for j = 1:nr
    ReconInfo(j).txnum = 1+fz*(j-1); % only need the 1st tx out of each 3 tx (for different focal zones) as they have the same origin and characteristics (focus)
    ReconInfo(j).rcvnum = 1+fz*(j-1);
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
SeqControl(1).command = 'timeToNextAcq';
SeqControl(1).argument = 200;  % 200 usec between ray lines; should be > 2*Receive.enDepth*T(period)
SeqControl(2).command = 'timeToNextAcq';
SeqControl(2).argument = round(80000 - 127*3*200); % 12.5 frames per second, fz=3
SeqControl(3).command = 'returnToMatlab';
SeqControl(4).command = 'jump'; % Jump back to start.
SeqControl(4).argument = 1;
nsc = 5;

% Specify Event structure arrays.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    k = nr*fz*(i-1);
    for j = 1:nr                 % Acquire all ray lines for frame
        w = fz*(j-1);
        for z = 1:fz
            Event(n).info = 'Acquire ray line';
            Event(n).tx = z+w;
            Event(n).rcv = k+z+w;
            Event(n).recon = 0;
            Event(n).process = 0;
            Event(n).seqControl = 1;
            n = n+1;
        end
    end
    % Replace last events SeqControl with inter-frame timeToNextAcq and transfer to host.
    Event(n-1).seqControl = [2,nsc];
       SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
       nsc = nsc+1;

    Event(n).info = 'recon and process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 1;
    if floor(i/4) == i/4     % Exit to Matlab every 4th frame reconstructed
        Event(n).seqControl = 3;
    else
        Event(n).seqControl = 0;
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
% - Sensitivity Cutoff
UI(1).Control =  {'UserB7','Style','VsSlider','Label','Sens. Cutoff',...
                  'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                  'SliderStep',[0.025,0.1],'ValueFormat','%1.3f'};
UI(1).Callback = text2cell('%SensCutoffCallback');

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 4;

% Save all the structures to a .mat file.
save('MatFiles/L11-5v_128RyLnsMultiFocal');

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
