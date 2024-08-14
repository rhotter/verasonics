% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use.
%
% File name SetUpL39_21gDWideBeamISC.m:
%    Generate .mat Sequence Object file for L39-21gD Linear array wide beam
%    scan, with 3 steering angles for spatial compounding using 4X interleaved sampling.
%    Aquisitions for each steering angle are acquired into separate RcvBuffers, and
%    the most recent frame from each buffer is reconstructed. A sliding 3 frame
%    averaging window is used on the processed image frames to combine the images.
%
% Description: This script uses the "interleaved sampling" mechanism for RF data
%    acquisition.  The L39-21gD is operated at a nominal center frequency for
%    Recon processing of 31.25 MHz, with 4X sampling so the required RF data
%    acquisition sample rate is 125 MHz.  Since the HW system cannot acquire
%    data at sample rates above 62.5 MHz, we acquire two sets of receive data
%    for each line, at an A/D sample rate of 62.5 MHz.  These pairs of
%    acquisition data lines are then interleaved in a 'Pre' routine for the
%    reconstruction processing to produce a single composite line sampled at
%    125 MHz.  The data acquired in the second acquistion of each pair
%    must have its data shifted by half a sample period, so it can be
%    interleaved with the first line.  This is accomplished by adding an
%    offset to the transmit delays for the first acquisition equal to that
%    half-sample interval; this in effect shifts each receive sample's data
%    earlier by the desired one half sample period.
%
% Last update:
%   April 2021 VTS-2152 computeTrans entry for L39-21gD
%   3-11-2021 Initial version of script
%

clear all

% Specify parameters for settings and presets.
P.startDepth = 0;  % Interleave acquisitions should have startDepth=0
P.endDepth = 192;
P.numTx = 32;   % no. of elements in TX aperture.
P.numRays = 48; % no. of wide beam rays for each steering angle
P.txFocus = 2*P.endDepth;
P.dtheta = 10*(pi/180);  % angle delta for beams steered left or right

% Specify system parameters.
Resource.Parameters.connector = 1;
Resource.Parameters.numTransmit = 128;  % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;  % number of receive channels.
Resource.Parameters.speedOfSound = 1540;
Resource.Parameters.simulateMode = 0; % 0 means no simulation, if hardware is present.
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Specify Trans structure array.
Trans.name = 'L39-21gD';
Trans.units = 'mm';
Trans = computeTrans(Trans);

mm2wl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);

% Specify PData structure arrays.
% - PData contains rectangular regions for no beam steering.
PData.PDelta = [Trans.spacing,0,0.5];
PData.Size(1) = ceil((P.endDepth-5)/PData.PDelta(3)); % startDepth, endDepth and PDelta set PData.Size.
PData.Size(2) = ceil((Trans.numelements*Trans.spacing)/PData.PDelta(1));
PData.Size(3) = 1;      % single image page
PData.Origin = [-Trans.spacing*(Trans.numelements/2),0,5]; % x,y,z of upper lft crnr.
% - Predefine all the Region structures, then modify them below.
PData.Region = repmat(struct('Shape',struct('Name','Parallelogram',...
                                            'Position',[0,0,5],...
                                            'width',14*Trans.spacing,...
                                            'height',P.endDepth-5,...
                                            'angle', 0.0)),1,3*P.numRays); % default is no steering
TxOrgX = (-63.5*Trans.spacing):(127*Trans.spacing/(P.numRays-1)):(63.5*Trans.spacing); % x coords of beam centers
% Specify P.numRays rectangular regions centered on TX beam origins (use default angle of 0.0).
for n = 1:P.numRays, PData.Region(n).Shape.Position(1) = TxOrgX(n); end
m = P.numRays;
% Define numRays steered left parallelogram regions, centered on TX beam origins. Adjust the angle
%   so that the steering goes to zero over 8 beams at the left and right edge.
for n = 1:P.numRays
    if n<=8
        angle = -((n-1)/8)*P.dtheta;
    elseif n>(P.numRays-8)
        angle = -((P.numRays-n)/8)*P.dtheta;
    else
        angle = -P.dtheta;
    end
    PData.Region(n+m).Shape.Position(1) = TxOrgX(n);
    PData.Region(n+m).Shape.height = (P.endDepth-5)/cos(angle);
    PData.Region(n+m).Shape.angle = angle;
end
m = m + P.numRays;
% Define numRays steered right parallelogram regions, centered on TX beam origins. Adjust the angle
%   so that the steering goes to zero over 8 beams at the left and right edge.
for n = 1:P.numRays
    if n<=8
        angle = ((n-1)/8)*P.dtheta;
    elseif n>(P.numRays-8)
        angle = ((P.numRays-n)/8)*P.dtheta;
    else
        angle = P.dtheta;
    end
    PData.Region(n+m).Shape.Position(1) = TxOrgX(n);
    PData.Region(n+m).Shape.height = (P.endDepth-5)/cos(angle);
    PData.Region(n+m).Shape.angle = angle;
end
PData.Region = computeRegions(PData);

% Specify Media object. 'pt1.m' script defines array of point targets.
% Media.MP = rand(20000,4);
% Media.MP(:,4) = 0.03*Media.MP(:,3) + 0.015;  % Random amplitude
% Media.MP(:,1) = (Trans.spacing*Trans.numelements)*(Media.MP(:,1)-0.5);
% Media.MP(:,2) = 0;
% Media.MP(:,3) = 200*Media.MP(:,3);
pt1;
Media.function = 'movePoints';

% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = P.numRays*6144; % this size allows for maximum range
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 12;    % frames stored in RcvBuffer.
Resource.RcvBuffer(2).datatype = 'int16';
Resource.RcvBuffer(2).rowsPerFrame = P.numRays*6144; % this size allows for all rays, with range of up to 400 wvlngths
Resource.RcvBuffer(2).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(2).numFrames = 12;    % should be multiple of 3
Resource.RcvBuffer(3).datatype = 'int16';
Resource.RcvBuffer(3).rowsPerFrame = P.numRays*6144; % this size allows for all rays, with range of up to 400 wvlngths
Resource.RcvBuffer(3).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(3).numFrames = 12;    % should be multiple of 3
Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).numFrames = 1;   % one intermediate buffer needed.
Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).numFrames = 10;
Resource.DisplayWindow(1).Title = 'L39-21gD WideBeamISC, Interleaved Sampling';
Resource.DisplayWindow(1).pdelta = 0.25;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).numFrames = 200;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = gray(256);

% Specify Transmit waveform structure.
TW.type = 'parametric';
TW.Parameters = [31.25, 1, 2, 1];   % 31.25 MHz center frequency, one cycle burst

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'Apod', zeros(1,Resource.Parameters.numTransmit), ...
                   'focus', P.txFocus, ...
                   'Steer', [0.0,0.0], ...
                   'Delay', zeros(1,Resource.Parameters.numTransmit), ...
                   'TXPD', [], ...
                   'peakCutOff', 1.0, ...
                   'peakBLMax', 4.0), 1, 3*2*P.numRays);

% - Set event specific TX attributes.
%   Define two transmit events for all beams of the three steering angles.
%   The two transmits will be identical with the exception of the first
%   having an extra delay of 0.25 wavelengths. TX(1) added delay will be
%   1/4 of the wls, based on the demodulation frequency.
for n = 1:2:2*P.numRays  % specify P.numRays transmit events
    TX(n).Origin(1) = TxOrgX((n+1)/2);
    TX(n+1).Origin(1) = TX(n).Origin(1);
    % Compute transmit aperture apodization
    TX(n).Apod = +(((mm2wl*Trans.ElementPos(:,1))>(TxOrgX((n+1)/2)-Trans.spacing*P.numTx/2))& ...
                 ((mm2wl*Trans.ElementPos(:,1))<(TxOrgX((n+1)/2)+Trans.spacing*P.numTx/2)))';
    [RIndices,CIndices,V] = find(TX(n).Apod);
    V = kaiser(size(V,2),1);
    TX(n).Apod(CIndices) = V;
    TX(n+1).Apod = TX(n).Apod;
    % Compute transmit delays
    TX(n+1).Delay = computeTXDelays(TX(n));
    TX(n).Delay = TX(n+1).Delay + 0.25; % 0.25 wls of xmit frequency=31.25 equivalent to 8 nsec.
    % Compute transmit pixel data
    TX(n).TXPD = computeTXPD(TX(n),PData(1));
end
m = 2*P.numRays;
for n = 1:2:2*P.numRays
    TX(n+m).Origin(1) = TX(n).Origin(1);
    TX(n+m+1).Origin(1) = TX(n+m).Origin(1);
    TX(n+m).Apod = TX(n).Apod;
    TX(n+m+1).Apod = TX(n+m).Apod;
    if n<=16
        TX(n+m).Steer = [-(((n+1)/2-1)/8)*P.dtheta,0.0];
    elseif n>(2*P.numRays-16)
        TX(n+m).Steer = [-((P.numRays-(n+1)/2)/8)*P.dtheta,0.0];
    else
        TX(n+m).Steer = [-P.dtheta,0.0];
    end
    TX(n+m+1).Apod = TX(n+m).Apod;
    % Compute transmit delays
    TX(n+m+1).Delay = computeTXDelays(TX(n+m));
    TX(n+m).Delay = TX(n+m+1).Delay + 0.25*(Trans.frequency/31.25);
    % Compute transmit pixel data
    TX(n+m).TXPD = computeTXPD(TX(n+m),PData);
end
m = m + 2*P.numRays;
for n = 1:2:2*P.numRays
    TX(n+m).Origin(1) = TX(n).Origin(1);
    TX(n+m+1).Origin(1) = TX(n+m).Origin(1);
    TX(n+m).Apod = TX(n).Apod;
    TX(n+m+1).Apod = TX(n+m).Apod;
    if n<=16
        TX(n+m).Steer = [(((n+1)/2-1)/8)*P.dtheta,0.0];
    elseif n>(2*P.numRays-16)
        TX(n+m).Steer = [((P.numRays-(n+1)/2)/8)*P.dtheta,0.0];
    else
        TX(n+m).Steer = [P.dtheta,0.0];
    end
    TX(n+m+1).Apod = TX(n+m).Apod;
    % Compute transmit delays
    TX(n+m+1).Delay = computeTXDelays(TX(n+m));
    TX(n+m).Delay = TX(n+m+1).Delay + 0.25*(Trans.frequency/31.25);
    % Compute transmit pixel data
    TX(n+m).TXPD = computeTXPD(TX(n+m),PData);
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [0,335,570,712,832,951,995,1020];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

RcvProfile.LnaZinSel = 31; % put LNA in "high-z" input state for best sensitivity

% Specify Receive structure arrays.
% - We need 2*P.numRay Receives for every acquisition angle, since two acquisitions are needed to
%   produce one interleaved acquisition line of data. For the interleaved acquisition approach with
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
Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
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
                        'callMediaFunc', 0), 1, 3*2*P.numRays*Resource.RcvBuffer(1).numFrames);
% - Set event specific Receive attributes for each frame.
for i = 1:Resource.RcvBuffer(1).numFrames
    k = 3*2*P.numRays*(i-1);
    Receive(k+1).callMediaFunc = 1;
    for j = 1:2:2*P.numRays
        % Receives for straight ahead.
        % First of each pair
        Receive(k+j).framenum = i;
        Receive(k+j).acqNum = j;
        % second of each pair
        Receive(k+j+1).framenum = i;
        Receive(k+j+1).acqNum = j+1;
        % Receives for steered left.
        Receive(k+2*P.numRays+j).bufnum = 2;
        Receive(k+2*P.numRays+j).framenum = i;
        Receive(k+2*P.numRays+j).acqNum = j;
        Receive(k+2*P.numRays+j+1).bufnum = 2;
        Receive(k+2*P.numRays+j+1).framenum = i;
        Receive(k+2*P.numRays+j+1).acqNum = j+1;
        % Receives for steered right
        Receive(k+4*P.numRays+j).bufnum = 3;
        Receive(k+4*P.numRays+j).framenum = i;
        Receive(k+4*P.numRays+j).acqNum = j;
        Receive(k+4*P.numRays+j+1).bufnum = 3;
        Receive(k+4*P.numRays+j+1).framenum = i;
        Receive(k+4*P.numRays+j+1).acqNum = j+1;
    end
end

% Specify Recon structure arrays.
% - We need three Recon structures, one for each steering direction.
Recon = repmat(struct('senscutoff', 0.5, ...
               'pdatanum', 1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...
               'rcvBufFrame',-1, ...
               'RINums',1:P.numRays), 1, 3);
% - Set specific Recon attributes.
Recon(2).RINums = (P.numRays+1):(2*P.numRays);
Recon(3).RINums = (2*P.numRays+1):(3*P.numRays);

% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 'accumIQ', ...  % default is to accumulate IQ data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'scaleFactor', 1.0, ...
                   'regionnum', 0), 1, 3*P.numRays);
% - Set specific ReconInfo attributes.
ReconInfo(1).Pre = 'clearInterBuf';
for j = 1:P.numRays
    ReconInfo(j).txnum = 2*j-1;
    ReconInfo(j).rcvnum = 2*j-1;
    ReconInfo(j).regionnum = j;
end
ReconInfo(P.numRays).Post = 'IQ2IntensityImageBuf';
k = P.numRays;
ReconInfo(k+1).Pre = 'clearInterBuf';
for j = 1:P.numRays
    ReconInfo(j+k).txnum = 2*j-1+2*k;
    ReconInfo(j+k).rcvnum = 2*j-1+2*k;
    ReconInfo(j+k).regionnum = j+k;
end
ReconInfo(k+P.numRays).Post = 'IQ2IntensityImageBuf';
k = k + P.numRays;
ReconInfo(k+1).Pre = 'clearInterBuf';
for j = 1:P.numRays
    ReconInfo(j+k).txnum = 2*j-1+4*P.numRays;
    ReconInfo(j+k).rcvnum = 2*j-1+4*P.numRays;
    ReconInfo(j+k).regionnum = j+2*P.numRays;
end
ReconInfo(k+P.numRays).Post = 'IQ2IntensityImageBuf';

% Specify Process structure array.
pers = 20;
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',2.5,...     % pgain is image processing gain
                         'reject',2,...      % reject level
                         'persistMethod','simple',...
                         'persistLevel',pers,...
                         'interpMethod','4pt',...  %method of interp. (1=4pt)
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','runAverage3',...
                         'compressMethod','power',...
                         'compressFactor',48,...
                         'mappingMethod','full',...
                         'display',1,...      % display image after processing
                         'displayWindow',1};

% Specify SeqControl structure arrays.
%  - Jump back to start.
SeqControl(1).command = 'jump';
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between ray line acquisitions
SeqControl(2).argument = 200;  % 200 usec
SeqControl(3).command = 'timeToNextAcq';
SeqControl(3).argument = 20000; %20000 usec = 20msec time between frames
SeqControl(4).command = 'returnToMatlab';
nsc = 5; % nsc is count of SeqControl objects

% Specify Event structure arrays.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    % Acquire frame with straight ahead wide beams
    k = 3*2*P.numRays*(i-1);
    for j = 1:2:2*P.numRays  % Acquire straight ahead frame
        Event(n).info = 'Acquire angle, 1st half interleave samples';
        Event(n).tx = j;   % use steered TX structure.
        Event(n).rcv = k+j;
        Event(n).recon = 0;      % no reconstruction.
        Event(n).process = 0;    % no processing
        Event(n).seqControl = 2;
        n = n+1;

        Event(n).info = 'Acquire angle, 2nd half interleave samples';
        Event(n).tx = j+1;   % use TX structure with interleave offset delay added.
        Event(n).rcv = k+j+1;
        Event(n).recon = 0;      % no reconstruction.
        Event(n).process = 0;    % no processing
        Event(n).seqControl = 2;
        n = n+1;
    end
    Event(n-1).seqControl = [3,nsc]; % modify last acquisition Event's seqControl
      SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
      nsc = nsc+1;

    Event(n).info = 'recon and process';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 1;      % reconstruction
    Event(n).process = 1;    % process
    Event(n).seqControl = 0;
    n = n+1;

    % Acquire frame with steered left wide beams
    for j = 1:2:2*P.numRays        % Acquire frame
        Event(n).info = 'Acquire ray, 1st half interleave samples';
        Event(n).tx = j+2*P.numRays;   % 2nd set of TX
        Event(n).rcv = k+2*P.numRays+j;  % 2nd set of Receives.
        Event(n).recon = 0;      % no reconstruction.
        Event(n).process = 0;    % no processing
        Event(n).seqControl = 2;
        n = n+1;

        Event(n).info = 'Acquire ray, 2nd half interleave samples';
        Event(n).tx = j+1+2*P.numRays;   % 2nd set of TX
        Event(n).rcv = k+2*P.numRays+j+1;  % 2nd set of Receives.
        Event(n).recon = 0;      % no reconstruction.
        Event(n).process = 0;    % no processing
        Event(n).seqControl = 2;
        n = n+1;
    end
    Event(n-1).seqControl = [3,nsc]; % modify last acquisition Event's seqControl
      SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
      nsc = nsc+1;

    Event(n).info = 'recon and process';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 2;      % reconstruction
    Event(n).process = 1;    % process
    Event(n).seqControl = 0;
    n = n+1;

    % Acquire frame with steered right wide beams
    for j = 1:2:2*P.numRays        % Acquire frame
        Event(n).info = 'acquire aperture.';
        Event(n).tx = j+4*P.numRays;   % use 3rd set of TX structures.
        Event(n).rcv = k+4*P.numRays+j; % use 3rd set of Receives
        Event(n).recon = 0;      % no reconstruction.
        Event(n).process = 0;    % no processing
        Event(n).seqControl = 2;
        n = n+1;

        Event(n).info = 'acquire aperture.';
        Event(n).tx = j+1+4*P.numRays;   % use 3rd set of TX structures.
        Event(n).rcv = k+4*P.numRays+j+1; % use 3rd set of Receives
        Event(n).recon = 0;      % no reconstruction.
        Event(n).process = 0;    % no processing
        Event(n).seqControl = 2;
        n = n+1;
    end
    Event(n-1).seqControl = [3,nsc]; % modify last acquisition Event's seqControl
      SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
      nsc = nsc+1;

    Event(n).info = 'recon and process';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 3;      % reconstruction
    Event(n).process = 1;    % process
    if i~=Resource.RcvBuffer(1).numFrames, Event(n).seqControl = 4;
    else Event(n).seqControl = 0;
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
wls2mm = 1;
AxesUnit = 'wls';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
        AxesUnit = 'mm';
        wls2mm = Resource.Parameters.speedOfSound/1000/Trans.frequency;
    end
end

UI(2).Control = VsSliderControl('LocationCode','UserA1','Label',['Range (',AxesUnit,')'],...
                 'SliderMinMaxVal',[64,200,P.endDepth]*wls2mm,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f', ...
                 'Callback', @RangeChangeCallback);

% - Transmit focus change
UI(3).Control = VsSliderControl('LocationCode','UserB4','Label',['TX Focus (',AxesUnit,')'],...
                 'SliderMinMaxVal',[20,500,P.txFocus]*wls2mm,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f', ...
                 'Callback', @txFocusCallback);

% - TX Aperture change
UI(4).Control = VsSliderControl('LocationCode','UserB3','Label','TX Aper',...
                 'SliderMinMaxVal',[1,128,P.numTx],'SliderStep',[1/128,1/64],'ValueFormat','%3.0f', ...
                 'Callback', @ApertureCallback);

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 3;
% Save all the structures to a .mat file.
save('MatFiles/L39-21gDWideBeamISC');


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

% - Range change
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

    PData = evalin('base','PData');
    PData.Size(1) = ceil((P.endDepth-5)/PData.PDelta(3)); % P.startDepth, P.endDepth and pdelta set PData.Size(1).
    for i = 1:size(PData.Region,2)
        PData.Region(i).Shape.height = P.endDepth-5;
    end
    PData.Region = computeRegions(PData);
    assignin('base','PData',PData);
    % Update TXPD data of TX structures.
    TX = evalin('base','TX');
    h = waitbar(0,'Program TX parameters, please wait!');
    steps = 3*2*P.numRays;
    for ind = 1:2:steps
        TX(ind).TXPD = computeTXPD(TX(ind),PData);
        waitbar(ind/steps)
    end
    close(h)
    assignin('base','TX',TX);
    evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData.Size(1)*PData.PDelta(3)/Resource.DisplayWindow(1).pdelta);');
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
    Control.Parameters = {'PData','InterBuffer','ImageBuffer','DisplayWindow','TX','Receive','TGC','Recon'};
    assignin('base','Control', Control);
    assignin('base', 'action', 'displayChange');
end


function txFocusCallback(hObject,~,UIValue)
    simMode = evalin('base','Resource.Parameters.simulateMode');
    % No focus change if in simulate mode 2.
    if simMode == 2
        set(hObject,'Value',evalin('base','P.txFocus'));
        return
    end
    mm2wl = evalin('base','mm2wl');
    Resource = evalin('base','Resource');
    P = evalin('base','P');
    P.txFocus = UIValue;
    if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
        if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
            P.txFocus = UIValue*mm2wl;
        end
    end
    assignin('base','P',P);
    % - Redefine event specific TX attributes for the new focus.
    TX = evalin('base', 'TX');
    PData = evalin('base','PData');
    Trans = evalin('base','Trans');
    h = waitbar(0,'Program TX parameters, please wait!');
    steps = 3*2*P.numRays;
    for ind = 1:2:steps
        % write new focus value to TX
        TX(ind).focus = P.txFocus;
        TX(ind+1).focus = TX(ind).focus;
        TX(ind+1).Delay = computeTXDelays(TX(ind));
        TX(ind).Delay = TX(ind+1).Delay + 0.25*(Trans.frequency/31.25);
        % Compute transmit pixel data
        TX(ind).TXPD = computeTXPD(TX(ind),PData);
        waitbar(ind/steps)
    end
    close(h)
    assignin('base','TX',TX);
    % Set Control command to update TX
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'TX'};
    assignin('base','Control', Control);
end

function ApertureCallback(hObject,~,UIValue)
    simMode = evalin('base','Resource.Parameters.simulateMode');
    % No aperture change if in simulate mode 2.
    if simMode == 2
        set(hObject,'Value',evalin('base','P.numTx'));
        return
    end
    Trans = evalin('base', 'Trans');
    P = evalin('base','P');
    P.numTx = UIValue;
    assignin('base','P',P);

    TX = evalin('base', 'TX');
    mm2wl = evalin('base','mm2wl');
    TxOrgX = evalin('base','TxOrgX');
    PData = evalin('base','PData');
%     rayDelta = Trans.numelements*Trans.spacing/P.numRays;
%     firstRayLocX = -((Trans.numelements-1)/2)*Trans.spacing;
    % - Redefine event specific TX attributes for the new aperture.
    h = waitbar(0,'Program TX parameters, please wait!');
    steps = 3*2*P.numRays;
    for n = 1:2:2*P.numRays  % specify P.numRays transmit events
        TX(n).Origin(1) = TxOrgX((n+1)/2);
        TX(n+1).Origin(1) = TX(n).Origin(1);
        % Compute transmit aperture apodization
        TX(n).Apod = +(((mm2wl*Trans.ElementPos(:,1))>(TxOrgX((n+1)/2)-Trans.spacing*P.numTx/2))& ...
                     ((mm2wl*Trans.ElementPos(:,1))<(TxOrgX((n+1)/2)+Trans.spacing*P.numTx/2)))';
        [~,CIndices,V] = find(TX(n).Apod);
        V = kaiser(size(V,2),1);
        TX(n).Apod(CIndices) = V;
        TX(n+1).Apod = TX(n).Apod;
        % Compute transmit delays
        TX(n+1).Delay = computeTXDelays(TX(n));
        TX(n).Delay = TX(n+1).Delay + 0.25*(Trans.frequency/31.25);
        % Compute transmit pixel data
        TX(n).TXPD = computeTXPD(TX(n),PData(1));
        waitbar(n/steps);
    end
    m = 2*P.numRays;
    for n = 1:2:2*P.numRays
        TX(n+m).Origin(1) = TX(n).Origin(1);
        TX(n+m+1).Origin(1) = TX(n+m).Origin(1);
        TX(n+m).Apod = TX(n).Apod;
        TX(n+m+1).Apod = TX(n+m).Apod;
        if n<=16
            TX(n+m).Steer = [-(((n+1)/2-1)/8)*P.dtheta,0.0];
        elseif n>(2*P.numRays-16)
            TX(n+m).Steer = [-((P.numRays-(n+1)/2)/8)*P.dtheta,0.0];
        else
            TX(n+m).Steer = [-P.dtheta,0.0];
        end
        TX(n+m+1).Apod = TX(n+m).Apod;
        % Compute transmit delays
        TX(n+m+1).Delay = computeTXDelays(TX(n+m));
        TX(n+m).Delay = TX(n+m+1).Delay + 0.25*(Trans.frequency/31.25);
        % Compute transmit pixel data
        TX(n+m).TXPD = computeTXPD(TX(n+m),PData);
        waitbar((n+2*P.numRays)/steps);
    end
    m = m + 2*P.numRays;
    for n = 1:2:2*P.numRays
        TX(n+m).Origin(1) = TX(n).Origin(1);
        TX(n+m+1).Origin(1) = TX(n+m).Origin(1);
        TX(n+m).Apod = TX(n).Apod;
        TX(n+m+1).Apod = TX(n+m).Apod;
        if n<=16
            TX(n+m).Steer = [(((n+1)/2-1)/8)*P.dtheta,0.0];
        elseif n>(2*P.numRays-16)
            TX(n+m).Steer = [((P.numRays-(n+1)/2)/8)*P.dtheta,0.0];
        else
            TX(n+m).Steer = [P.dtheta,0.0];
        end
        TX(n+m+1).Apod = TX(n+m).Apod;
        % Compute transmit delays
        TX(n+m+1).Delay = computeTXDelays(TX(n+m));
        TX(n+m).Delay = TX(n+m+1).Delay + 0.25*(Trans.frequency/31.25);
        % Compute transmit pixel data
        TX(n+m).TXPD = computeTXPD(TX(n+m),PData);
        waitbar((n+4*P.numRays)/steps);
    end
    close(h)
    assignin('base','TX', TX);
    % Set Control command to update TX
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'TX'};
    assignin('base','Control', Control);
end
