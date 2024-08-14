% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: SetUpL11_5vFlashAngles_saveRFallFrames.m Example of plane wave
% imaging with steering angle transmits and saving RF data in realtime. This
% script shows how to use external function to save the RF data of all frames in
% a receive buffer in realtime. Here, only non-zeros channel data will be saved
% to reduce the saving time. The hardware and software are synchronized with a
% sync command after data acquisition of each buffer to prevent the data from
% being overwritten. Therefore, please note that the interval between the last
% frame of the current buffer and the first frame of the next buffer will be
% longer than the timeToNextAcq defined in the SeqControl(2). In addition, the
% default frame limit is set to 100 RcvBuffers. Please change to a bigger number
% for more buffers.
%
% Last update:
% 05/06/2020 - Update to SW 4.3 format for new user UIControls and External function definitions (VTS 1691).
%   More info:(?/Example_Scripts/Vantage_Features/New UI Scheme/SetUpL11_5vFlash_NewUI)
% 12/28/2018 - test with SW 4.0

clear all

P.startDepth = 5;   % Acquisition depth in wavelengths
P.endDepth   = 192; % This should preferably be a multiple of 128 samples.
P.maxDepth   = 195;

NumBuffer = 4;

na = 7;      % Set na = number of angles.
if (na > 1)
    dtheta = (36*pi/180)/(na-1);
    P.startAngle = -36*pi/180/2;
else
    dtheta = 0;
    P.startAngle=0;
end % set dtheta to range over +/- 18 degrees.

% Define system parameters.
Resource.Parameters.numTransmit    = 128;   % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;   % number of receive channels.
Resource.Parameters.speedOfSound   = 1540;  % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose        = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode   = 1;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Specify Trans structure array.
Trans.name  = 'L11-5v';
Trans.units = 'wavelengths'; % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);  % L11-5v transducer is 'known' transducer so we can use computeTrans.
Trans.maxHighVoltage = 50;  % set maximum high voltage limit for pulser supply.

% Specify PData structure array.
PData(1).PDelta  = [Trans.spacing, 0, 0.5];
PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3)); % startDepth, endDepth and pdelta set PData(1).Size.
PData(1).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(3) = 1;      % single image page
PData(1).Origin  = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of upper lft crnr.
% No PData.Region specified, so a default Region for the entire PData array will be created by computeRegions.

% Specify Media object. 'pt1.m' script defines array of point targets.
pt1;
Media.attenuation = -0.5;
Media.function = 'movePoints';

% define the max buffer length to prevent writing too much data
maxBufLength     = ceil(sqrt(P.maxDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
maxBufSizePerAcq = 128*ceil(maxBufLength*8/128);

% Specify Resources.
for N = 1:NumBuffer
    Resource.RcvBuffer(N).datatype     = 'int16';
    Resource.RcvBuffer(N).rowsPerFrame = na*maxBufSizePerAcq; % this size allows for maximum range
    Resource.RcvBuffer(N).colsPerFrame = Resource.Parameters.numRcvChannels;
    Resource.RcvBuffer(N).numFrames    = 30;    % 30 frames stored in RcvBuffer.
end

% Specify Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,.67,2,1];

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'Apod', kaiser(Resource.Parameters.numTransmit,1)', ...
                   'focus', 0.0, ...
                   'Steer', [0.0,0.0], ...
                   'Delay', zeros(1,Trans.numelements)), 1, na);
% - Set event specific TX attributes.
if fix(na/2) == na/2       % if na even
    P.startAngle = (-(fix(na/2) - 1) - 0.5)*dtheta;
else
    P.startAngle = -fix(na/2)*dtheta;
end
for n = 1:na   % na transmit events
    TX(n).Steer = [(P.startAngle+(n-1)*dtheta),0.0];
    TX(n).Delay = computeTXDelays(TX(n));
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [0,141,275,404,510,603,702,782];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Receive structure arrays.
% - We need na Receives for every frame.
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength,...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, (na*Resource.RcvBuffer(1).numFrames)*NumBuffer);

% - Set event specific Receive attributes for each frame.
for N = 1:NumBuffer
    ind = Resource.RcvBuffer(N).numFrames*na*(N-1);
    for i = 1:Resource.RcvBuffer(N).numFrames
        Receive(ind+na*(i-1)+1).callMediaFunc = 1;
        for j = 1:na
            Receive(ind+na*(i-1)+j).bufnum = N;
            Receive(ind+na*(i-1)+j).framenum = i;
            Receive(ind+na*(i-1)+j).acqNum = j;
        end
    end
end

% obtain a local file path
str = which('VSX');
% save files in VSX folder under UserData
fp = fileparts(str);
storageFolder = fullfile( fullfile(fp, 'UserData'), 'StorageTest' );
% make sure UserData folder exist
if ~exist(storageFolder, 'dir')
    mkdir(storageFolder);
end

for N = 1:NumBuffer
% Define a storage class for each Buffer
    % Save Data Process
Process(N).classname  = 'Storage'; % define the storage class name
Process(N).method     = 'save';    % to export files specify 'save'
Process(N).Parameters = { ...
    'srcbuffer','receive',... % name of buffer to export.
    'srcbufnum',    N,    ... % buffer number of buffer to export
    'srcframenum',  0,    ... % process the last frame in a RcvBuffer
    'filename',         'TestFileAllFrames',   ... % the name of the file, do not provide the file tag such as .vrs or similar. It will be removed
    'folder',           storageFolder... % the folder where to save the storage data, this should be folder that exist on the file system
    'compression', 'BL',  ... % specify the type of compression, either NONE or BL, future releases may provide other compression algorithms
    'doPostCompression', 1, ... % set to 1 if the files should be compressed in a post processing process, 0 if not
    'timestamp', 1, ...  % indicate whether to store a time stamp in the header
    'numdigits', 3, ...  % specify the number of digits to build the file name number e.g. if = 3; TestFile001.vrs, TestFile002.vrs, ...TestFile999.vrs
    'enabled',   1, ...  % disable the saving will add a button that allows to start the saving at a particular point, if set to one, saving will start with script first acquisitions
    'maxfiles',  999, ... % Can be used to specify the maximum number of files
    'studyID',   'Test ID',      ... % will save 'Test ID' as a studyID in the file header of each stored file
    'comment',   'Some Notes',  ... % will save 'Some Notes' as a comment in the file header of each stored file
    'sampleID',  'A Sample ID', ... % will save 'A Sample ID' as a sampleID in the file header
    };
end

% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; % jump back to start
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(2).argument = 300;  % usec
SeqControl(3).command = 'sync';

nsc = 4; % nsc is count of SeqControl objects

% Specify Event structure arrays.
n = 1;
for N = 1:NumBuffer
    ind = Resource.RcvBuffer(N).numFrames*na*(N-1);
    for i = 1:Resource.RcvBuffer(N).numFrames
        for j = 1:na                      % Acquire frame
            Event(n).info = 'Full aperture.';
            Event(n).tx = j;
            Event(n).rcv = ind+na*(i-1)+j;
            Event(n).recon = 0;
            Event(n).process = 0;
            Event(n).seqControl = 2;
            n = n+1;
        end
        Event(n-1).seqControl = [2,nsc]; % modify last acquisition Event's seqControl
        SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
        nsc = nsc+1;
    end

    Event(n-1).seqControl = [3,nsc-1]; % modify last acquisition Event's seqControl
    SeqControl(nsc-1).command = 'transferToHost'; % transfer frame to host buffer

    Event(n).info = 'Save RF data';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = N;
    Event(n).seqControl = 0;
    n=n+1;
end

Event(n).info = 'Jump back';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 1;

% Save all the structures to a .mat file.
save('MatFiles/L11-5vFlashAngles_saveRFallFrame');

%% Automatic VSX Execution:
% Uncomment the following line to automatically run VSX every time you run
% this SetUp script (note that if VSX finds the variable 'filename' in the
% Matlab workspace, it will load and run that file instead of prompting the
% user for the file to be used):

% filename = 'L11-5vFlashAngles_saveRFallFrame';  VSX;