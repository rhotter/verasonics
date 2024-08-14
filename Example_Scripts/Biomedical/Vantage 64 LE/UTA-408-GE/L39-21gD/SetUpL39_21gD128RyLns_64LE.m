% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: SetUpL39_21gD128RyLns_64LE.m - Example of scanline imaging with 
%                                 focused transmits at a single focal depth.
%
% Description:
%   Sequence programming file for Dax L39-21gD Linear array, using 128 ray lines
%   (focus transmits) and 64 receive acquisitions. Of the 128 transmit
%   channels, the active transmit aperture is limited based on user-entered
%   transmit focus and f-number. All 64 receive channels are active for
%   each acquisition. 
%
% Last update:
%   April 2021 VTS-2152 computeTrans entry for L39-21gD
%   12/28/2020 - New transducer on GE408 UTA
%   07/07/2021 - created for 64LE 

clear all

P.startDepth = 0;
P.endDepth = 192;    % Acquisition depth in wavelengths
P.txFocus = 80;      % Initial transmit focus.
P.numRays = 128;     % no. of Rays

% Specify system parameters.
Resource.Parameters.numTransmit = 128;  % number of transmit channels.
Resource.Parameters.numRcvChannels = 64;  % number of receive channels.
Resource.Parameters.speedOfSound = 1540;
Resource.Parameters.speedCorrectionFactor = 1.0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Specify Trans structure array.
Trans.name = 'L39-21gD';
Trans.units = 'mm';
Trans = computeTrans(Trans);

RcvProfile.LnaZinSel = 31; % put LNA in "high-z" input state for best sensitivity

% Specify PData structure array.
PData.PDelta = [Trans.spacing/2, 0, 0.2];
PData.Size(1) = ceil((P.endDepth-P.startDepth)/PData.PDelta(3));
PData.Size(2) = ceil((Trans.numelements*Trans.spacing)/PData.PDelta(1));
PData.Size(3) = 1;      % single image page
PData.Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of upper lft crnr
% - specify 128 Region structures.
PData.Region = repmat(struct('Shape',struct( ...
                    'Name','Rectangle',...
                    'Position',[0,0,P.startDepth],...
                    'width',Trans.spacing,...
                    'height',P.endDepth-P.startDepth)),1,128);
% - set position of regions to correspond to beam spacing.
for i = 1:128
    PData.Region(i).Shape.Position(1) = (-63.5 + (i-1))*Trans.spacing;
end
PData.Region = computeRegions(PData(1));

% Specify Media object.
pt1;
Media.attenuation = -0.5;
Media.function = 'movePoints';

% Specify Resources.
Resource.RcvBuffer.datatype = 'int16';
Resource.RcvBuffer.rowsPerFrame = 4*184320; % this should be larger than 128*Receive.endDepth*4 for max depth (doubled for 4X sampling)
Resource.RcvBuffer.colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer.numFrames = 10;
Resource.InterBuffer.numFrames = 1;  % one intermediate buffer needed.
Resource.ImageBuffer.numFrames = 50;
Resource.DisplayWindow.Title = 'L39-21gD 128RyLns 64LE interleaved sampling';
Resource.DisplayWindow.pdelta = 0.25;
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
TW = struct('type','parametric',...
            'Parameters',[31.25,.67,2,1]);

% Specify P.numRays TX structure arrays. Transmit centered on element n in the array for event n.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'focus', P.txFocus, ...
                   'Steer', [0.0,0.0], ...
                   'Apod', zeros(1,Trans.numelements), ...
                   'Delay', zeros(1,Trans.numelements)), 1, 2*P.numRays);

% Determine TX aperture based on focal point and desired f number.
txFNum = 2;  % set to desired f-number value for transmit (range: 1.0 - 20)
P.numTx = round((P.txFocus/txFNum)/Trans.spacing); % no. of elements in 1/2 aperture.
txNumEl = floor(P.numTx/2);
if txNumEl > (Trans.numelements/2 - 1), txNumEl = floor(Trans.numelements/2 - 1); end
% txNumEl is the number of elements to include on each side of the
% center element, for the specified focus and sensitivity cutoff.
% Thus the full transmit aperture will be 2*txNumEl + 1 elements.
%display('Number of elements in transmit aperture:');
%disp(2*txNumEl+1);

% - Set event specific TX attributes.
for n = 1:2:2*P.numRays   % 2*P.numRays transmit events
    % Set transmit Origins to positions of elements.
    TX(n).Origin = [(-63.5 + (n/2 + 0.5 -1))*Trans.spacing, 0.0, 0.0];
    % Set transmit Apodization so (1 + 2*TXnumel) transmitters are active.
    lft = n/2 + 0.5 - txNumEl;
    if lft < 1, lft = 1; end
    rt = n/2 + 0.5 + txNumEl;
    if rt > Trans.numelements, rt = Trans.numelements; end
    TX(n).Apod(lft:rt) = 1.0;
    % interleave TXs
    TX(n+1).Origin = TX(n).Origin;
    TX(n+1).Apod = TX(n).Apod;
    TX(n+1).Delay = computeTXDelays(TX(n+1));
    TX(n).Delay = TX(n+1).Delay + 0.25*(Trans.frequency/31.25);    
end

% Specify Receive structure arrays.
% - We need 4*P.numRay Receives for every acquisition angle, since two acquisitions are needed to
%   produce one interleaved acquisition line of data and two acquisitons for 2-1 synthetic aperture. For the interleaved acquisition approach with
%   2X interleave, we must take into account the actual sample rate at which the digital filters in
%   the hardware will be operating: twice the center frequency set by Trans.frequency, not the
%   typical 4X factor.  This means that the Nyquist limit for the filters will be at
%   Trans.frequency; the higher half of the transducer signal frequency spectrum will be folded over
%   the lower half due to aliasing.  (The subsequent interleave combination of the two acquisition
%   events will unfold this aliasing).  Therefore the filter actually used for the Input Filter must
%   be defined as a high-pass filter. The net effect after interleave will be a symmetric bandpass
%   filter centered at 31.25 MHz.
% - A Highpass coefficient array has been defined below, representing a fractional bandwidth of
%   100%, relative to Fc at 31.25 MHz.  The coefficient array listed below was developed using the
%   matlab fdatool, with a Hamming window.
% - Note that for the L30 we would ideally use a bandpass filter centered near the transducer's
%   nominal center frequency of 30 MHz, not the 31.25 MHz forced on the CGD filters.
HighPassCoef100 = [-0.0000    0.0014   -0.0000   -0.0024   -0.0000    0.0046   -0.0000...
                   -0.0081   -0.0000    0.0136   -0.0000   -0.0217   -0.0000    0.0341...
                   -0.0000   -0.0551   -0.0000    0.1009   -0.0000   -0.3169    0.5006];

maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
Receive = repmat(struct('Apod', zeros(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BWI', ...
                        'demodFrequency', 31.25, ...
                        'InputFilter', HighPassCoef100, ...
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, 4*P.numRays*Resource.RcvBuffer(1).numFrames);
% - Set event specific Receive attributes for each frame.
for i = 1:Resource.RcvBuffer(1).numFrames
    k = 4*P.numRays*(i-1);
    Receive(k+1).callMediaFunc = 1;
    for j = 1:4:4*P.numRays
        % left aperture
        % pairs of acquisition receives to buffer 1 for interleaved acquisition at 1/2 sample rate:
        % First of each pair
        Receive(k+j).Apod(1:64) = 1;
        Receive(k+j).framenum = i;
        Receive(k+j).acqNum = j;
        % second of each pair
        Receive(k+j+1).Apod(1:64) = 1;
        Receive(k+j+1).framenum = i;
        Receive(k+j+1).acqNum = j+1;
        
        % right aperture
        % pairs of acquisition receives to buffer 1 for interleaved acquisition at 1/2 sample rate:
        % First of each pair
        Receive(k+j+2).Apod(65:128) = 1;
        Receive(k+j+2).framenum = i;
        Receive(k+j+2).acqNum = j+2;
        % second of each pair
        Receive(k+j+3).Apod(65:128) = 1;
        Receive(k+j+3).framenum = i;
        Receive(k+j+3).acqNum = j+3;
    end
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [250,460,665,765,818,869,900,910];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Recon structure array.
Recon = struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame',-1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...
               'RINums', 1:2*P.numRays);

% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 'replaceIQ', ...  % replace intensity data
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 0), 1, 2*P.numRays);
% - Set specific ReconInfo attributes.
for j = 1:2:2*P.numRays
    ReconInfo(j).txnum = j;
    ReconInfo(j).rcvnum = 2*j-1;
    ReconInfo(j).regionnum = (j+1)/2;
     
    ReconInfo(j+1).mode = 'accumIQ_replaceIntensity'; % accum and detect
    ReconInfo(j+1).txnum = j;
    ReconInfo(j+1).rcvnum = 2*j+1;
    ReconInfo(j+1).regionnum = (j+1)/2;
end

% Specify Process structure array.
pers = 20;
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',3.0,...            % pgain is image processing gain
                         'reject',2,...      % reject level
                         'persistMethod','simple',...
                         'persistLevel',pers,...
                         'interpMethod','4pt',...
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',45,...
                         'mappingMethod','full',...
                         'display',1,...      % display image after processing
                         'displayWindow',1};

% Specify SeqControl structure arrays.
SeqControl(1).command = 'timeToNextAcq';
SeqControl(1).argument = 150;  % 150 usec between ray lines
SeqControl(2).command = 'timeToNextAcq';
SeqControl(2).argument = round(44400 - 127*SeqControl(1).argument); % 40 frames per second
SeqControl(3).command = 'returnToMatlab';
SeqControl(4).command = 'jump'; % Jump back to start.
SeqControl(4).argument = 1;
nsc = 5;

% Specify Event structure arrays.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    k = 4*P.numRays*(i-1);
    for j = 1:2:2*P.numRays                 % Acquire all ray lines for frame
        Event(n).info = 'Acquire ray line, first interleave, left aperture';
        Event(n).tx = j;
        Event(n).rcv = k+2*j-1;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 1;
        n = n+1;
        
        Event(n).info = 'Acquire ray line, 2nd interleave samples,left aperture';
        Event(n).tx = j+1;   % use TX structure with interleave offset delay added.
        Event(n).rcv = k+2*j;
        Event(n).recon = 0;      % no reconstruction.
        Event(n).process = 0;    % no processing
        Event(n).seqControl = 1;
        n = n+1;
        
        Event(n).info = 'Acquire ray line, first interleave, right aperture';
        Event(n).tx = j;
        Event(n).rcv = k+2*j+1;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 1;
        n = n+1;
        
        Event(n).info = 'Acquire ray line, 2nd interleave samples,right aperture';
        Event(n).tx = j+1;   % use TX structure with interleave offset delay added.
        Event(n).rcv = k+2*j+2;
        Event(n).recon = 0;      % no reconstruction.
        Event(n).process = 0;    % no processing
        Event(n).seqControl = 1;
        n = n+1;
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
import vsv.seq.uicontrol.VsSliderControl; 

% - Sensitivity Cutoff
UI(1).Control = VsSliderControl('LocationCode','UserB7','Label','Sens. Cutoff',...
                  'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                  'SliderStep',[0.025,0.1],'ValueFormat','%1.3f', ...
                  'Callback', @SensCutoffCallback );

% - Range Change
wls2mm = 1;
AxesUnit = 'wls';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
        AxesUnit = 'mm';
        wls2mm = Resource.Parameters.speedOfSound/1000/Trans.frequency;
    end
end
UI(2).Control = VsSliderControl('LocationCode','UserA1','Label',['Range (',AxesUnit,')'],...
                 'SliderMinMaxVal',[64,300,P.endDepth]*wls2mm,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f', ...
                 'Callback', @RangeChangeCallback);  
             

% - Transmit focus change
UI(3).Control = VsSliderControl('LocationCode','UserB4','Label',['TX Focus (',AxesUnit,')'],...
                 'SliderMinMaxVal',[64,300,P.endDepth]*wls2mm,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f', ...
                 'Callback', @TxFocusCallback);  


% - F number change
UI(4).Control = VsSliderControl('LocationCode','UserB3','Label','F Number',...
                 'SliderMinMaxVal',[1,20,round(P.txFocus/(P.numTx*Trans.spacing))],'SliderStep',[0.05,0.1],'ValueFormat','%2.0f', ...
                 'Callback', @FNumCallback);

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 4;

% Save all the structures to a .mat file.
save('MatFiles/L39-21gD128RyLns_64LE');
% filename = ('L39-21gD128RyLns_64LE'); % VSX    % permits immediately running VSX without specifying the matfile name


%% **** Callback routines used by UIControls (UI) ****
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

% - RangeChangeCallback - Range change
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


% - TxFocusCallback - TX focus changel
function TxFocusCallback(hObject,~,UIValue)
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
    for n = 1:2:2*P.numRays   % 2*P.numRays transmit events
        TX(n).focus = P.txFocus;
        TX(n+1).focus = TX(n).focus;
        TX(n+1).Delay = computeTXDelays(TX(n+1));
        TX(n).Delay = TX(n+1).Delay + 0.25*(Trans.frequency/31.25); 
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


% - FNumCallback - F number change
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
    txNumEl = floor(P.numTx/2);
    for n = 1:2:2*P.numRays   % 2*P.numRays transmit events
        % Set transmit Apodization.
        lft = n/2 + 0.5 - txNumEl;
        if lft < 1, lft = 1; end
        rt = n/2 + 0.5 + txNumEl;
        if rt > Trans.numelements, rt = Trans.numelements; end
        TX(n).Apod = zeros(1,Trans.numelements);
        TX(n).Apod(lft:rt) = 1.0;
        TX(n+1).Apod = TX(n).Apod;
        TX(n+1).Delay = computeTXDelays(TX(n+1));
        TX(n).Delay = TX(n+1).Delay + 0.25*(Trans.frequency/31.25);
    end
    assignin('base','TX', TX);
    % Set Control command to update TX
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'TX'};
    assignin('base','Control', Control);
end

