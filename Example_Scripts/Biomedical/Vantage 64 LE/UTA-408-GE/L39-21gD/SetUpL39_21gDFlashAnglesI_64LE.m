% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use.
%
% File name SetUpL39_21gDFlashAnglesI_64LE.m 
%   This script is derived from SetUpL39_21gDFlashAnglesI.m.
%
% Description: Sequence programming file for L39-21gD Linear array, using plane wave
%   transmits with multiple steering angles and 2-1 synthetic aperture
%   (two acquisitions per line to acquire data from all 128 receive
%   elements). All 128 transmit and 64 receive channels are active for each
%   synthetic aperture acquisition. Processing is asynchronous with
%   respect to acquisition.
%
% Last Update
%   April 2021 VTS-2152 computeTrans entry for L39-21gD
%   07/13/2021 created for 64LE


clear all

% Specify parameters for settings and presets.
P.startDepth = 5;
P.endDepth = 192;
P.numAngles = 7;
% - set angle increment, dtheta, to space lines over +/- 18 degrees.
if (P.numAngles > 1)
    dtheta = (36*pi/180)/(P.numAngles-1);
    startAngle = -36*pi/180/2;
else
    dtheta = 0;
    startAngle=0;
end

% Specify system parameters.
Resource.Parameters.connector = 1;
Resource.Parameters.numTransmit = 128;  % number of transmit channels.
Resource.Parameters.numRcvChannels = 64;  % number of receive channels.
Resource.Parameters.speedOfSound = 1540;
Resource.Parameters.simulateMode = 0; % 0 means no simulation, if hardware is present.
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Specify Trans structure array.
Trans.name = 'L39-21gD';
Trans.units = 'mm';
Trans = computeTrans(Trans);

RcvProfile.LnaZinSel = 31; % put LNA in "high-z" input state for best sensitivity

% Specify PData structure array.
PData.PDelta = [Trans.spacing/2, 0, 0.5]; % [pdeltaX, pdeltaY, pdeltaZ]
PData.Size(1) = ceil((P.endDepth-P.startDepth)/PData.PDelta(3)); % range, origin and pdelta set PData.Size.
PData.Size(2) = ceil((Trans.numelements*Trans.spacing)/PData.PDelta(1));
PData.Size(3) = 1;      % single image page
PData.Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of uppr lft crnr.
% No PData.Region specified, so a default Region for the entire PData array will be created by computeRegions.

% Specify Media object.
pt1;
Media.function = 'movePoints';

% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 2*P.numAngles*6144; % this size allows for maximum range
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 20;    % frames stored in RcvBuffer.
Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).numFrames = 1;   % one intermediate buffer needed.
Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).numFrames = 10;
Resource.DisplayWindow(1).Title = 'L39-21gD 64LE Flash Angles, Interleaved Sampling';
Resource.DisplayWindow(1).pdelta = 0.3;
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
TW.type = 'parametric';
TW.Parameters = [31.25, 0.67, 4, 1];   % 31.25 MHz center frequency, two cycle burst

% Specify TX structure array.  We need two TX definitions for each interleave acquisition.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'focus', 0.0, ...
                   'Steer', [0.0,0.0], ...
                   'Apod', ones(1,Trans.numelements), ...
                   'Delay', zeros(1,Trans.numelements)), 1, 2*P.numAngles);
% - Set event specific TX attributes.
%   Define two transmit events at each of the desired steering angles.
%   The two transmits will be identical with the exception of the first
%   having an extra delay of 0.25 wavelengths. TX(1) added delay will be
%   1/4 of the wls, based on the demodulation frequency.
%   TX.delay will be converted to time based on the Trans.frequency

for n = 1:2:2*P.numAngles
    TX(n).Steer = [(startAngle+((n-1)/2)*dtheta),0.0];
    TX(n+1).Steer = TX(n).Steer;
    TX(n+1).Delay = computeTXDelays(TX(n+1));
    TX(n).Delay = TX(n+1).Delay + 0.25*(Trans.frequency/31.25);
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [0,288,490,613,716,818,869,920];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Receive structure arrays.
% - We need 4*P.numRay Receives for every acquisition angle, since two acquisitions are needed to
%   produce one interleaved acquisition line of data and two for 2-1 synthetic aperture. For the interleaved acquisition approach with
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
% - Note that for the L35-16vX we would ideally use a bandpass filter centered near the transducer's
%   nominal center frequency of 28 MHz, not the 31.25 MHz forced on the CGD filters.
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
                        'callMediaFunc', 0), 1, 4*P.numAngles*Resource.RcvBuffer(1).numFrames);
% - Set event specific Receive attributes for each frame.
for i = 1:Resource.RcvBuffer(1).numFrames
    k = 4*P.numAngles*(i-1);
    Receive(k+1).callMediaFunc = 1;
    for j = 1:4:4*P.numAngles
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

% Specify Recon structure arrays.
% - We need one Recon structure which will be reused for all frames.
Recon = struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame',-1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...
               'RINums',1:2*P.numAngles);

% Define ReconInfo structures.
% We need na ReconInfo structures for the na steering angles.
ReconInfo = repmat(struct('mode', 'accumIQ', ...  % default is to accumulate IQ data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 1), 1, 2*P.numAngles);
% - Set specific ReconInfo attributes.
ReconInfo(1).mode = 'replaceIQ';

for j = 1:2:2*P.numAngles
    
    ReconInfo(j).txnum = j;
    ReconInfo(j).rcvnum = 2*j-1;
     
    ReconInfo(j+1).txnum = j;
    ReconInfo(j+1).rcvnum = 2*j+1;
end
ReconInfo(2*P.numAngles).mode = 'accumIQ_replaceIntensity';

% Specify Process structure array.
pers = 20;
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,... % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',5.0,...     % pgain is image processing gain
                         'reject',2,...
                         'persistMethod','simple',... % none, simple, dynamic
                         'persistLevel',pers,...
                         'interpMethod','4pt',...  % method of interpolation
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',45,...
                         'mappingMode','full',...
                         'display',1,...      % display image after processing
                         'displayWindow',1};

% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump';
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';
SeqControl(2).argument = 80;  % 80usec
SeqControl(3).command = 'timeToNextAcq';
SeqControl(3).argument = 20000;  % 20000 = 20msec (~ 50 fps)
SeqControl(4).command = 'returnToMatlab';
nsc = 5; % nsc is count of SeqControl objects

% Specify Event structure arrays.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    k = 4*P.numAngles*(i-1);
    for j = 1:2:2*P.numAngles  % Acquire frame
        Event(n).info = 'Acquire ray line, first interleave, left aperture';
        Event(n).tx = j;
        Event(n).rcv = k+2*j-1;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;
        
        Event(n).info = 'Acquire ray line, 2nd interleave samples,left aperture';
        Event(n).tx = j+1;   % use TX structure with interleave offset delay added.
        Event(n).rcv = k+2*j;
        Event(n).recon = 0;      % no reconstruction.
        Event(n).process = 0;    % no processing
        Event(n).seqControl = 2;
        n = n+1;
        
        Event(n).info = 'Acquire ray line, first interleave, right aperture';
        Event(n).tx = j;
        Event(n).rcv = k+2*j+1;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;
        
        Event(n).info = 'Acquire ray line, 2nd interleave samples,right aperture';
        Event(n).tx = j+1;   % use TX structure with interleave offset delay added.
        Event(n).rcv = k+2*j+2;
        Event(n).recon = 0;      % no reconstruction.
        Event(n).process = 0;    % no processing
        Event(n).seqControl = 2;
        n = n+1;
    end
        Event(n-1).seqControl = 3; % set frame rate in last acquisition

    Event(n).info = 'Transfer acquisitions to host';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = nsc; % SeqControl struct defined below.
       SeqControl(nsc).command = 'transferToHost';
       nsc = nsc + 1;
    n = n+1;

    Event(n).info = 'recon and process';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 1;      % reconstruction
    Event(n).process = 1;    % process
    if floor(i/5) == i/5     % Exit to Matlab every 5th frame reconstructed
        Event(n).seqControl = 4;
    else
        Event(n).seqControl = 0;
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
import vsv.seq.uicontrol.VsSliderControl; 

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
UI(2).Control = VsSliderControl('LocationCode','UserA1','Label',['Range (',AxesUnit,')'],...
                 'SliderMinMaxVal',MinMaxVal,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f', ...
                 'Callback', @RangeChangeCallback);

clear i j n

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 5;

% Save all the structures to a .mat file.
save('MatFiles/L39-21gDFlashAnglesInterleave_64LE');
% filename=('L39-21gDFlashAnglesInterleave_64LE'); VSX


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
    mm2wl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);
    P = evalin('base','P');
    P.endDepth = UIValue;
    if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
        if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
            P.endDepth = UIValue*mm2wl;
        end
    end
    assignin('base','P',P);
    % Modify PData rows and adjust DisplayWindow height.
    evalin('base','PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));');
    evalin('base','PData(1).Region = computeRegions(PData(1));');
    evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);');
    % Modify Receive endDepth
    Receive = evalin('base', 'Receive');
    maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
    for i = 1:size(Receive,2)
        Receive(i).endDepth = maxAcqLength;
    end
    assignin('base','Receive',Receive);
    % Adjust TGC range
    evalin('base','TGC.rangeMax = P.endDepth;');
    evalin('base','TGC.Waveform = computeTGCWaveform(TGC);');
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'PData','InterBuffer','ImageBuffer','DisplayWindow','Receive','TGC','Recon'};
    assignin('base','Control', Control);
    assignin('base', 'action', 'displayChange');
end
