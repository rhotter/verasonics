% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use.
%
% File name SetUpDualXdcr.m:
% Example script to illustrate the use of two separate transducers imaging
% simultaneously on Vantage-256 system.  For this example, an L11-5v
% transducer is connected to the left connector and a P4-2v to the right
% connector.  This script is a combination of the "L11-5vFlash" and
% "P4-2vFlash" example scripts, with both of them running in parallel but
% independently of each other in the same script.  To accomplish this, we
% define separate structures for PData, Recon, and ImageDisplay for each
% transducer.
%
% Note that this example is not intended to demonstrate any level of
% imaging or RF data acquistion that involves interaction between the two
% probes (such as transmitting on one and receiving on the other in the
% same event), or any defined physical relationship between the two probes
% and thus between their resulting images.
%
% What this example does illustrate is the tradeoff associated with
% combining two different probes of significantly different center
% frequencies in a single script.  A fundamental constraint of the current
% Vantage software is that only one Trans structure can be defined for the
% entire script. Wavelength units are used in user scripts to define PData
% structures, TX delays, etc.- and the wavelength unit is defined as the
% wavelength of Trans.frequency at the speed of sound specified in the
% script.  Since the L11-5v center frequency is 7.6 MHz while the P4-2v is
% at 2.7 MHz, we define a compromise shared center frequency at 4.5 MHz as
% the basis for wavelength units.  By assigning this value to
% Trans.frequency before calling computeTrans for the two individual
% probes, the shared compromise definition of a wavelength will be applied
% consistently to all wavelength units throughout the script.  The transmit
% waveform and receive processing center frequencies can still be set
% properly for both probes in the separate TW, TX, and Receive
% structures used for each probe.  For fields specified in the shared Trans
% structure we do not have this flexibility, however.  For these fields
% (such as Trans.Bandwidth and Trans.lensCorrection) a shared compromise
% value must be chosen, as illustrated in this example.
%
% The script requires a single Trans structure, defined as if it were a 192
% element transducer using both connectors- but within that shared
% structure elements 1:128 represent the L11-5v and elements 129:192
% represent the P4-2v.  The TX and Receive .Apod arrays for each transducer
% will have zeros in the element positions representing the other
% transducer.  Note that the GUI controls will generally modify settings
% for the L11-5v only, since it is indexed first in the relevant
% structures. Additional GUI control object(s) can be created to switch
% between probes.
%
% Two separate TPC profiles are used, to allow independent control of the
% transmit voltage for each transducer using P1 and P2 HV sliders.
%
% Throughout this script, comments starting with: % ***Dual Xdcr***
% have been inserted to illustrate the changes made for the dual-imaging
% example
%
% Last update:
% 05/06/2020 - Update to SW 4.3 format for new user UIControls and External function definitions (VTS 1691).
%   More info:(?/Example_Scripts/Vantage_Features/New UI Scheme/SetUpL11_5vFlash_NewUI)
% 4-30-2019 - Modified for use with 4.0 software, and modified to use
%   Verasonics L11-5v and P4-2v probes

clear all

% --- Commonly Changed User Parammeters -------------------------
numRFFrames = 10;        % RF Data frames
frameRateFactor = 5;     % Factor for converting sequenceRate to frameRate.
Fc = 4.5;                % ***Dual Xdcr*** This is the compromise shared
% frequency to be used by both L11-5v (7.6 MHz) and P4-2v (2.7 MHz), set to
% roughly the geometric mean of those two frequencies
% ---------------------------------------------------------------

% Specify system parameters.
Resource.Parameters.numTransmit = 256; % ***Dual Xdcr*** number of transmit channels.
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% ***Dual Xdcr*** We must explicitly specify a Connector array
% value of [1 2], to inform the system that both connectors will be used to
% support an overall "transducer" with 192 elements
Resource.Parameters.Connector = [1 2];
% ***Dual Xdcr*** Also set the 'fakeScanhead' parameter, to allow this
% script to run on the HW system with no transducer actually connected and
% avoid confusion when the two distinct transducers are present.
Resource.Parameters.fakeScanhead = 1;

% Specify Trans structure array.
% ***Dual Xdcr*** First we will use computeTrans to define two separate
% TransL11 and TransP4 structures, and then we will use parameters from
% both of those to create the shared Trans structure that will actually be
% used by the script.

% ***Dual Xdcr*** Create the TransL11 structure
TransL11.name = 'L11-5v';
TransL11.frequency = Fc; % ***Dual Xdcr*** This is the shared frequency to be used by both L11-5v and P4-2v
TransL11.units = 'wavelengths';
TransL11 = computeTrans(TransL11);    % L11-5v transducer is 'known' transducer so we can use computeTrans.

% ***Dual Xdcr*** Create the TransP4 structure
TransP4.name = 'P4-2v';
TransP4.frequency = Fc; % ***Dual Xdcr*** This is the shared frequency to be used by both L11-5v and P4-2v
TransP4.units = 'wavelengths';
TransP4 = computeTrans(TransP4);

% ***Dual Xdcr*** Now use TransL11 and TransP4 to create the shared Trans structure ***
Trans.name = 'custom';      % Must be 'custom' to prevent confusion from the two unique transducer ID's that will actually be connected
Trans.id = -1; % id of -1 means ignore the actual ID value
Trans.units = 'wavelengths';
Trans.frequency = Fc;      % ***Dual Xdcr*** This is the shared frequency to be used by both L11-5v and P4-2v
Trans.Bandwidth = [TransP4.Bandwidth(1), TransL11.Bandwidth(2)]; % set bandwidth to span both probes
Trans.type = 0;             % 1D straight array geometry applies to both L11-5v and P4-2v
Trans.numelements = 192;    % total over both connectors
% Concatenate the two element position and ConnectorES arrays
Trans.ElementPos = [TransL11.ElementPos; TransP4.ElementPos];
Trans.ConnectorES = [TransL11.ConnectorES; (128+TransP4.ConnectorES)]; % add 128 to get channel numbers for second connector
% For the following parameters just copy the L11-5v values for the first
% image processing
Trans.lensCorrection = TransL11.lensCorrection;
Trans.spacing = TransL11.spacing;
Trans.elementWidth = TransL11.elementWidth;
Trans.ElementSens = TransL11.ElementSens;
% For the following use an appropriate shared value
Trans.impedance = 50;
Trans.maxHighVoltage = 50;  % set maximum high voltage limit for pulser supply.
Trans.connType = 1;

% Specify PData(1) structure array for L11-5v.
P(1).startDepth = 5;   % Acquisition depth in wavelengths
P(1).endDepth = 192;   % This should preferrably be a multiple of 128 samples.
PData(1).PDelta = [Trans.spacing, 0, 0.5];
PData(1).Size(1) = ceil((P(1).endDepth-P(1).startDepth)/PData(1).PDelta(3)); % startDepth, endDepth and pdelta set PData(1).Size.
PData(1).Size(2) = ceil((TransL11.numelements*Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(3) = 1;      % single image page
PData(1).Origin = [-Trans.spacing*(TransL11.numelements-1)/2,0,P(1).startDepth]; % x,y,z of upper lft crnr.
% No PData.Region specified, so a default Region for the entire PData array will be created by computeRegions.

% Specify PData(2) structure array for P4-2v.
P(2).startDepth = 0;
P(2).endDepth = 160;
P(2).theta = -pi/4;
P(2).rayDelta = 2*(-P(2).theta);
P(2).aperture = TransP4.numelements*TransP4.spacing;  % P.aperture in wavelengths
P(2).radius = (P(2).aperture/2)/tan(-P(2).theta); % dist. to virt. apex

% Set up PData structure.
PData(2).PDelta = [0.875, 0, 0.5];
PData(2).Size(1) = 10 + ceil((P(2).endDepth-P(2).startDepth)/PData(2).PDelta(3));
PData(2).Size(2) = 10 + ceil(2*(P(2).endDepth + P(2).radius)*sin(-P(2).theta)/PData(2).PDelta(1));
PData(2).Size(3) = 1;
PData(2).Origin = [-(PData(2).Size(2)/2)*PData(2).PDelta(1),0,P(2).startDepth];
PData(2).Region = struct(...
            'Shape',struct('Name','SectorFT', ...
            'Position',[0,0,-P(2).radius], ...
            'z',P(2).startDepth, ...
            'r',P(2).radius+P(2).endDepth, ...
            'angle',P(2).rayDelta, ...
            'steer',0));
PData(2).Region = computeRegions(PData(2));

% Specify Media object.
pt1;
Media.attenuation = -0.5;
Media.function = 'movePoints';

% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 2*4096; % ***Dual Xdcr*** doubles since there will be two acquisitions per frame, one for each xdcr
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numTransmit;
Resource.RcvBuffer(1).numFrames = numRFFrames;  % e.g., 10 frames used for RF cineloop.
% ***Dual Xdcr***  define the first ImageBuffer and DisplayWindow as usual, for the
% L11-5v
Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).numFrames = 10;
Resource.DisplayWindow(1).Title = 'L11-5vFlash';
Resource.DisplayWindow(1).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [50,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).numFrames = 20;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = gray(256);

% ***Dual Xdcr*** Now define a second ImageBuffer and DisplayWindow for the
% P4-2v
Resource.ImageBuffer(2).datatype = 'double';
% Resource.ImageBuffer(2).rowsPerFrame = 1024; % this is for maximum depth
% Resource.ImageBuffer(2).colsPerFrame = PData(2).Size(2);
Resource.ImageBuffer(2).numFrames = 10;
Resource.DisplayWindow(2).Title = 'P4-2vFlash';
Resource.DisplayWindow(2).pdelta = 0.35;
DwWidth = ceil(PData(2).Size(2)*PData(2).PDelta(1)/Resource.DisplayWindow(2).pdelta);
DwHeight = ceil(PData(2).Size(1)*PData(2).PDelta(3)/Resource.DisplayWindow(2).pdelta);
Resource.DisplayWindow(2).Position = [600,200, DwWidth, DwHeight];
Resource.DisplayWindow(2).ReferencePt = [PData(2).Origin(1),0,PData(2).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(2).Type = 'Verasonics';
Resource.DisplayWindow(2).numFrames = 20;
Resource.DisplayWindow(2).AxesUnits = 'mm';
Resource.DisplayWindow(2).Colormap = gray(256);

% Specify Transmit waveform structure.
% ***Dual Xdcr***
TW(1).type = 'parametric';
TW(1).Parameters = [7.6,.67,2,1];
TW(2).type = 'parametric';
TW(2).Parameters = [2.7,.67,2,1];

% ***Dual Xdcr*** Specify TX structure array for the L11-5v.
TX(1).waveform = 1;            % use 1st TW structure.
TX(1).Origin = [0.0,0.0,0.0];  % flash transmit origin at (0,0,0).
TX(1).focus = 0;
TX(1).Steer = [0.0,0.0];       % theta, alpha = 0.
TX(1).Apod = [ones(1,128), zeros(1, 64)]; % ***Dual Xdcr*** L11-5v uses channels 1:128
TX(1).Delay = computeTXDelays(TX(1));

% ***Dual Xdcr*** Separate TX structure for the P4-2v
TX(2).waveform = 2;
TX(2).Origin = [0,0,0];             % set origin to 0,0,0 for flat focus.
TX(2).focus = -P(2).radius;  	% set focus to negative for concave TX.Delay profile.
TX(2).Steer = [0,0];
TX(2).Apod = [zeros(1,128), ones(1,64)]; % ***Dual Xdcr*** P4-2v uses channels 129:192
TX(2).Delay = computeTXDelays(TX(2));

% Specify TPC structures ... creates two TPC profiles and two HV control sliders.
TPC(1).name = 'L11-5v';
TPC(1).maxHighVoltage = 50;
TPC(2).name = 'P4-2v';
TPC(2).maxHighVoltage = 50;

% Specify TGC Waveform structure.
% ***Dual Xdcr*** Separate TGC controls for each transducer
% - L11-5v TGC
TGC(1).CntrlPts = [0,141,275,404,550,700,850,1000];
TGC(1).rangeMax = P(1).endDepth;
TGC(1).Waveform = computeTGCWaveform(TGC(1));
% - P4-2v TGC
TGC(2).CntrlPts = [87,580,639,698,750,844,929,1023];
TGC(2).rangeMax = P(2).endDepth;
TGC(2).Waveform = computeTGCWaveform(TGC(2));

% Specify Receive structure arrays -
%   endDepth - add additional acquisition depth to account for some channels
%              having longer path lengths.
% ***Dual Xdcr*** For our simultaneous dual transducer acquisition scheme,
% we define two interleaved sets of Receive structures.  Each acquisition
% frame will consist of one L11-5vFlash acquisition followed by one P4-2vFlash
% acquisition.  This same concept can be easily extended to
% multi-acquisition formats such as FlashAngles, Doppler ensembles, or ray
% line imaging.
maxAcqLengthL11 = ceil(sqrt(P(1).endDepth^2 + ((TransL11.numelements-1)*TransL11.spacing)^2));
maxAcqLengthP4 = sqrt(P(2).aperture^2 + P(2).endDepth^2 - 2*P(2).aperture*P(2).endDepth*cos(P(2).theta-pi/2)) - P(2).startDepth;
wlsPer128 = 128/(4*2); % wavelengths in 128 samples for 4 samplesPerWave
Receive = repmat(struct('Apod', zeros(1,Trans.numelements), ...
                        'startDepth', P(1).startDepth, ...
                        'endDepth', P(1).startDepth + wlsPer128*ceil(maxAcqLengthL11/wlsPer128), ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 1),1,2*Resource.RcvBuffer(1).numFrames); % ***Dual Xdcr*** two Receive sturctures per frame

% ***Dual Xdcr***  - Set event-specific and transducer-specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    % ***Dual Xdcr*** For each frame number i,
    % Receive(2*i-1) will be for the L11-5v, and
    % Receive(2*i) will be for the P4-2v.
    % Note that most of the L11-5v specific values were used in the initial
    % Receive structure definition so they don't need to be updated here

    % L11-5v Receive
    Receive(2*i-1).framenum = i;
    Receive(2*i-1).Apod = [ones(1,128), zeros(1,64)]; % L11-5v uses elements 1:128
    Receive(2*i-1).demodFrequency = TW(1).Parameters(1);
    % P4-2v Receive
    Receive(2*i).framenum = i;
    Receive(2*i).Apod = [zeros(1,128), ones(1,64)]; % P4-2v uses channels 129:192
    Receive(2*i).startDepth = P(2).startDepth;
    Receive(2*i).endDepth = P(2).startDepth + wlsPer128*ceil(maxAcqLengthP4/wlsPer128);
    Receive(2*i).TGC = 2;
    Receive(2*i).acqNum = 2; % P4-2v is second acquisition in each frame
    Receive(2*i).demodFrequency = TW(2).Parameters(1);
    Receive(2*i).callMediaFunc = 0; % only move the media points once per frame
end

% Specify Recon structure for L11-5v.
Recon(1) = struct('senscutoff', 0.5, ...
               'pdatanum', 1, ...
               'rcvBufFrame', -1, ...
               'ImgBufDest', [1,-1], ...
               'RINums', 1);

% Define ReconInfo structure for L11-5v.
ReconInfo(1) = struct('mode', 'replaceIntensity', ...          % replace amplitude.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum',1);

% Specify Recon structure for P4-2v.
% Copy Recon(1) and then modify values that are different
Recon(2) = Recon(1);
Recon(2).pdatanum = 2;
Recon(2).ImgBufDest = [2,-1];
Recon(2).RINums = 2;

% Define ReconInfo structure for P4-2v using same approach as for Recon.
ReconInfo(2) = ReconInfo(1);
ReconInfo(2).txnum = 2;
ReconInfo(2).rcvnum = 2;

% Specify Process structure array for the L11-5v.
pers = 20;
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',5.0,...            % pgain is image processing gain
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

% Specify separate Process structure array for the P4-2v.
Process(2).classname = 'Image';
Process(2).method = 'imageDisplay';
Process(2).Parameters = {'imgbufnum',2,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',2,...    % number of PData structure to use
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
                         'displayWindow',2};

% Specify SeqControl structure arrays.
% Set the frame interval using timeToNextAcq
SeqControl(1).command = 'timeToNextAcq';
SeqControl(1).argument = 50000;

% return to matlab for GUI updates every fifth frame
SeqControl(2).command = 'returnToMatlab';

% at the end, jump back and start over
SeqControl(3).command = 'jump';
SeqControl(3).argument = 3; % don't need to repeat the first to events that made initial TPC profile selection

% select TPC profile 2 for the P4-2v
SeqControl(4).command = 'setTPCProfile';
SeqControl(4).argument = 2;
SeqControl(4).condition = 'next';

% select TPC profile 1 for the L7
SeqControl(5).command = 'setTPCProfile';
SeqControl(5).argument = 1;
SeqControl(5).condition = 'next';

nsc = 6; % nsc is index of next SeqControl object to be defined

% Specify the Event sequence
n = 1; % n is count of Events

Event(n).info = 'select TPC profile 1 at startup';
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

Event(n).info = 'noop delay for initial profile selection';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = nsc;
n = n+1;
   SeqControl(nsc).command = 'noop';
   SeqControl(nsc).argument = 100e3/.2; % 100 msec delay
   nsc = nsc + 1;

% Acquire all frames defined in RcvBuffer
for i = 1:Resource.RcvBuffer(1).numFrames
    Event(n).info = 'acquisition for L11-5v';
    Event(n).tx = 1;
    Event(n).rcv = 2*i-1;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = [4, 1];
    n = n+1;

    Event(n).info = 'acquisition for P4-2v';
    Event(n).tx = 2;
    Event(n).rcv = 2*i;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = [5, 1, nsc];
       SeqControl(nsc).command = 'transferToHost';
       nsc = nsc + 1;
    n = n+1;

    Event(n).info = 'Reconstruct & display L11-5v';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 1;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = 'Reconstruct & display P4-2v';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 2;
    Event(n).process = 2;
    if floor(i/frameRateFactor) == i/frameRateFactor     % Exit to Matlab every 5th frame
        Event(n).seqControl = 2;
    else
        Event(n).seqControl = 0;
    end
    n = n+1;

end

Event(n).info = 'Jump back to third event to repeat the sequence';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 3;

% Save all the structures to a .mat file.
save('MatFiles/DualXdcr');

%% Automatic VSX Execution:
% Uncomment the following line to automatically run VSX every time you run
% this SetUp script (note that if VSX finds the variable 'filename' in the
% Matlab workspace, it will load and run that file instead of prompting the
% user for the file to be used):

% filename = 'DualXdcr';  VSX;