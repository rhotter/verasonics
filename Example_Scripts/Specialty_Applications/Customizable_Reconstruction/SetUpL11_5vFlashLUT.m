% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: SetUpL11_5vFlashLUT.m - Example of imaging with single plane wave
%                                 transmit
% Description:
%   This is a duplicate of the "SetUpL11_5vFlash" example script, modified
%   to demonstrate the use of the user-programmable Recon LUT feature.  The
%   Recon LUT tables are defined explicitly in this script (and redefined
%   in the callback function for the range change control), but are defined
%   to duplicate what Recon would have created automatically if the
%   user-programmable LUT feature was not being used.  Changes from the
%   L11-5vFlash script are highlighted with the comment line **Recon LUT**.
%
% Last update:
% 06/17/2020 - created for demonstration of user-defined Recon LUT

clear all
% **Recon LUT** Note that for scripts with a user-defined Recon LUT the
% Receive startDepth must be set to zero to ensure the resulting combined
% LUT will not have a negative receive sample data index, since that can
% cause a Matlab crash
P.startDepth = 0;   % Acquisition depth in wavelengths
P.endDepth = 192;   % This should preferrably be a multiple of 128 samples.

% Define system parameters.
Resource.Parameters.numTransmit = 128;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;    % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;

%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Specify Trans structure array.
Trans.name = 'L11-5v';
Trans.units = 'mm'; % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);
Trans.maxHighVoltage = 50;  % much restricted than specified on data sheet
% **Recon LUT**  mm to wavelength conversion factor for Recon LUT calculations
mm2wls = 1000*Trans.frequency/Resource.Parameters.speedOfSound;


% Specify PData structure array.
PData(1).PDelta = [Trans.spacing, 0, 0.5];
PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3)); % startDepth, endDepth and pdelta set PData(1).Size.
PData(1).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(3) = 1;      % single image page
PData(1).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of upper lft crnr.
% **Recon LUT** No PData.Region specified, so a default Region for the
% entire PData array will be created by computeRegions for use in Recon LUT
% calculations.
PData.Region = computeRegions(PData);

% Specify Media object. 'pt1.m' script defines array of point targets.
pt1;
Media.attenuation = -0.5;
Media.function = 'movePoints';

% Specify Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 4096;   % this size allows for maximum range
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 100;       % 100 frames used for RF cineloop.
Resource.InterBuffer(1).numFrames = 1;  % one intermediate buffer needed.
Resource.ImageBuffer(1).numFrames = 10;
Resource.DisplayWindow(1).Title = 'L11-5vFlashLUT';
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
TW(1).Parameters = [Trans.frequency,.67,2,1];

% Specify TX structure array.
TX(1).waveform = 1;            % use 1st TW structure.
TX(1).Origin = [0.0,0.0,0.0];  % flash transmit origin at (0,0,0).
TX(1).focus = 0;
TX(1).Steer = [0.0,0.0];       % theta, alpha = 0.
TX(1).Apod = ones(1,Trans.numelements);
TX(1).Delay = computeTXDelays(TX(1));

% Specify TGC Waveform structure.
TGC.CntrlPts = [0,300,444,552,606,747,870,920];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Receive structure arrays -
%   endDepth - add additional acquisition depth to account for some channels
%              having longer path lengths.
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
                        'callMediaFunc', 1),1,Resource.RcvBuffer(1).numFrames);
% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    % -- Acquisition for full frame.
    Receive(i).framenum = i;
end

% Specify Recon structure arrays.
Recon = struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame', -1, ...     % use most recently transferred frame
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...  % auto-increment ImageBuffer each recon
               'rcvLUT', zeros(128*PData.Size(1)*PData.Size(2),1,'uint16'),...
               'RINums', 1);
% **Recon LUT** Compute Recon.rcvLUT
if (isfield(Resource.Parameters,'speedCorrectionFactor'))&&(~isempty(Resource.Parameters.speedCorrectionFactor))
    spcf = Resource.Parameters.speedCorrectionFactor;
else
    spcf = 1.0;
end
angle = 0;
for j = 51:101
    if Trans.ElementSens(j)<Recon.senscutoff, break; end
    angle = angle + 0.0314;
end
cosangle = cos(angle);
Recon.rcvLUT = zeros(Trans.numelements*PData.Size(1)*PData.Size(2),1,'uint16');
for j = 1:(PData.Size(1)*PData.Size(2))
    x = PData.Origin(1) + floor((j-1)/PData.Size(1))*PData.PDelta(1);
    z = PData.Origin(3) + rem((j-1),PData.Size(1))*PData.PDelta(3);
    zsq = z*z;
    k = Trans.numelements*(j-1);
    for m = 1:Trans.numelements
        pl = sqrt((x - mm2wls*Trans.ElementPos(m,1))*(x - mm2wls*Trans.ElementPos(m,1)) + zsq);
        if (z/pl) > cosangle
            pl = spcf*pl + Trans.lensCorrection;
            if pl > 0.0
                Recon.rcvLUT(k+m) = uint16(16*pl);
            else
                Recon.rcvLUT(k+m) = uint16(0);
            end
        else
            Recon.rcvLUT(k+m) = uint16(0);
        end
    end
end

% Define ReconInfo structures.
ReconInfo = struct('mode', 'replaceIntensity', ...
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 1, ...
                   'LUT', zeros(3*PData.Size(1)*PData.Size(2),1,'int32'));
% **Recon LUT** Compute ReconInfo.LUT
        ReconInfo.LUT = zeros(PData.Region(ReconInfo.regionnum).numPixels,1,'int32');
        n = 1;
        for j = 1:(PData.Size(1)*PData.Size(2)) % reconstruct all PData pixels in region 1
            ReconInfo.LUT(n) = int32(j-1); % linear adr. from 0
            n = n+1;
            ReconInfo.LUT(n) = int32(256); % weight of 1.0 in 24.8 format.
            n = n+1;
            x = PData.Origin(1) + floor((j-1)/PData.Size(1))*PData.PDelta(1);
            z = PData.Origin(3) + rem((j-1),PData.Size(1))*PData.PDelta(3);
            zsq = z*z;
            % For unfocused beams, transmit path length is shortest path from pixel to transducer.
            plmin = 1000;
            m = 0;
            for k = 1:Trans.numelements % for all active aperture elements
                if TX(1).Apod(k) == 0.0; continue; end
                pl = spcf*sqrt((x-mm2wls*Trans.ElementPos(k,1))*(x-mm2wls*Trans.ElementPos(k,1)) + zsq) + TX(1).Delay(k);
                if pl < plmin
                    plmin = pl;
                    m = k;
                end
            end
            ReconInfo.LUT(n) = int32(round((plmin + 1.0 + Trans.lensCorrection)*16));
            n = n+1;
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
SeqControl(2).command = 'timeToNextAcq';  % time between frames
SeqControl(2).argument = 3000;  % 3 msec
SeqControl(3).command = 'returnToMatlab';
nsc = 4; % nsc is count of SeqControl objects

n = 1; % n is count of Events

% Acquire all frames defined in RcvBuffer
for i = 1:Resource.RcvBuffer(1).numFrames

    Event(n).info = 'Full aperture.';
    Event(n).tx = 1;
    Event(n).rcv = i;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = [2,nsc];
       SeqControl(nsc).command = 'transferToHost';
       nsc = nsc + 1;
    n = n+1;

    Event(n).info = 'Reconstruct';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 1;
    if floor(i/5) == i/5     % Exit to Matlab every 5th frame
        Event(n).seqControl = 3;
    else
        Event(n).seqControl = 0;
    end
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
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
        AxesUnit = 'mm';
        MinMaxVal = MinMaxVal * (Resource.Parameters.speedOfSound/1000/Trans.frequency);
    end
end
UI(2).Control = {'UserA1','Style','VsSlider','Label',['Range (',AxesUnit,')'],...
                 'SliderMinMaxVal',MinMaxVal,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f'};
UI(2).Callback = text2cell('%RangeChangeCallback');

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 5;

% Save all the structures to a .mat file.
save('MatFiles/L11-5vFlashLUT');
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
% **Recon LUT**  mm to wavelength conversion factor for Recon LUT calculations
mm2wls = Trans.frequency/(Resource.Parameters.speedOfSound/1000);

P = evalin('base','P');
P.endDepth = UIValue;
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
        P.endDepth = UIValue*mm2wls;
    end
end
assignin('base','P',P);

% **Recon LUT**  redefine PData.Region for use in Recon LUT calculations
PData = evalin('base','PData');
PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));
PData(1).Region = computeRegions(PData(1));
assignin('base','PData',PData);

evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);');
Receive = evalin('base', 'Receive');
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
for i = 1:size(Receive,2)
    Receive(i).endDepth = maxAcqLength;
end
assignin('base','Receive',Receive);
evalin('base','TGC.rangeMax = P.endDepth;');
evalin('base','TGC.Waveform = computeTGCWaveform(TGC);');
% **Recon LUT**  compute new Recon.rcvLUT
spcf = evalin('base','spcf');
Recon = evalin('base','Recon');
angle = 0;
for j = 51:101
    if Trans.ElementSens(j)<Recon.senscutoff, break; end
    angle = angle + 0.0314;
end
cosangle = cos(angle);
Recon.rcvLUT = zeros(Trans.numelements*PData.Size(1)*PData.Size(2),1,'uint16');
for j = 1:(PData.Size(1)*PData.Size(2))
    x = PData.Origin(1) + floor((j-1)/PData.Size(1))*PData.PDelta(1);
    z = PData.Origin(3) + rem((j-1),PData.Size(1))*PData.PDelta(3);
    zsq = z*z;
    k = Trans.numelements*(j-1);
    for m = 1:Trans.numelements
        pl = sqrt((x - mm2wls*Trans.ElementPos(m,1))*(x - mm2wls*Trans.ElementPos(m,1)) + zsq);
        if (z/pl) > cosangle
            pl = spcf*pl + Trans.lensCorrection;
            if pl > 0.0
                Recon.rcvLUT(k+m) = uint16(16*pl);
            else
                Recon.rcvLUT(k+m) = uint16(0);
            end
        else
            Recon.rcvLUT(k+m) = uint16(0);
        end
    end
end
assignin('base','Recon',Recon);
% **Recon LUT**  compute ReconInfo.LUT - assume only one ReconInfo
ReconInfo = evalin('base','ReconInfo');
TX = evalin('base','TX');
ReconInfo.LUT = zeros(PData.Region(ReconInfo.regionnum).numPixels,1,'int32');
n = 1;
for j = 1:(PData.Size(1)*PData.Size(2)) % reconstruct all PData pixels in region 1
    ReconInfo.LUT(n) = int32(j-1); % linear adr. from 0
    n = n+1;
    ReconInfo.LUT(n) = int32(256); % weight of 1.0 in 24.8 format.
    n = n+1;
    x = PData.Origin(1) + floor((j-1)/PData.Size(1))*PData.PDelta(1);
    z = PData.Origin(3) + rem((j-1),PData.Size(1))*PData.PDelta(3);
    zsq = z*z;
    % For unfocused beams, transmit path length is shortest path from pixel to transducer.
    plmin = 1000;
    m = 0;
    for k = 1:Trans.numelements % for all active aperture elements
        if TX(1).Apod(k) == 0.0; continue; end
        pl = spcf*sqrt((x-mm2wls*Trans.ElementPos(k,1))*(x-mm2wls*Trans.ElementPos(k,1)) + zsq) + TX(1).Delay(k);
        if pl < plmin
            plmin = pl;
            m = k;
        end
    end
    ReconInfo.LUT(n) = int32(round((plmin + 1.0 + Trans.lensCorrection)*16));
    n = n+1;
end
assignin('base','ReconInfo',ReconInfo);

Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'PData','InterBuffer','ImageBuffer','DisplayWindow','Receive','TGC','Recon'};
assignin('base','Control', Control);
assignin('base', 'action', 'displayChange');
return
%RangeChangeCallback
