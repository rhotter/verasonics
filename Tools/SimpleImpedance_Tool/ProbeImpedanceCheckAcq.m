function [ProbeZdata] = ProbeImpedanceCheckAcq(txFreqIndx, Trans)
% File name ProbeImpedanceCheckAcq.m: RF data Acquisition for making
% approximate measurements of probe impedance, using BIST driver.

% Copyright 2001-2021, Verasonics, Inc. All worldwide rights and remedies under
% all intellectual property laws and industrial property laws are reserved.

%% Set preferences

TestAcq = 0; % set to 1 to allow free-running with GUI controls; 0 for normal use

%% Get system configuration

SysConfig = hwConfigCheck(1);

VDAS = SysConfig.VDAS; % will be one if HW system is present, zero if not
if VDAS == 0
    fprintf(2, 'HW system not detected. This Utility cannot be used.\n')
    return
end

if nargin == 0
    % set default frequency index for 5.2 MHz
    txFreqIndx = 18;
end

% Set freqRng to match the HW configuration provided by hwConfigCheck
freqRng = SysConfig.Frequency;
if strcmpi(freqRng, 'LF')
    fprintf(2, 'Low Frequency configuration not supported; will be added in a future release.\n')
    return
end

if ~TestAcq
    Mcr_GuiHide = 1; % prevent display of vsx_gui
    showPlots = 0;
else
    showPlots = 1; % plot RF data etc
end

%% Specify initial operating state parameters

% Select frequency for testing 
Freqin = 2.^((-5:22)/5); % 30 steps from 0.5 to 31.25 MHz, approx 15% per step
rcvFs = zeros(1, 28);
AAF = zeros(1, 28);
ADRateDiv = zeros(1, 28);
AAFlist = [5 10 15 20 30 35];
FreqOut = zeros(1, 28);
% Convert to the actual TX frequencies we wil use
Fs = 125;
for freqnum = 28:-1:1
    div = round(Fs/Freqin(freqnum));
    if div > 9
        div = div/2;
        Fs = Fs/2;
    end
    FreqOut(freqnum) = Fs/div;
    rcvFs(freqnum) = Fs;
    ADRateDiv(freqnum) = div;
    aafIndx = find((FreqOut(freqnum) < 0.9*AAFlist), 1);
    AAF(freqnum) = AAFlist(aafIndx);
end
numFrames = 2; % number of RF frames in the data acquisition sequence

numAper = 4; % number of Receive subaperture acquisitions
% using four subapertures will allow use of V32LE with UTA 260-Mux, and 4:1
% HVMux aperture selection with UTA 1024-Mux
numBist = 2; % separate BIST acquistion driven from each end,
% to null out drop along length of the BIST network

PRF = 0.005; % PRF in KHz for the numAper acquisition events per frame (note will be slower if accumulate>1)
% Note that this script runs synchronously, so frame rate will be set by SW
% event sequence processing time unless PRF is set very low.

initTGC = 1023;%500; % initial TGC setting (range 0 to 1023)


rcvES = 30; % initial receive element signal to select & plot
    % a value from 1 to numES plots that specific channel

autoscale = 1; % initial setting for autoscale control.  set to 1 to enable autoscaling of RF data plots; 0 to disable

nullDC = 1; % initial state for removal of per channel DC offset

startDepthUsec = 5;%128;%512 / 15.625; % derived from original startdepth of 512 wavelengths at 15.625 MHz
acqRangeUsec = 128 / 15.625;

plotdepth = 0.27; % initial depth for RF data plot, as a fraction of linelength from 0 to 1
plotrange = 1.0;%0.1; % initial range for RF data plot, as a fraction of linelength from 0 to 1
bdurUsec = 12; % burst duration in usec

% RcvProfile.DCsubtract = 'off';

%% Initialize frequency range dependent values and channel count

switch freqRng
    case 'SF'
        txFreqList = FreqOut; % 14 choices for the transmit burst frequency
        samplesPerWave = 4;
    case 'HF'
        txFreqList = FreqOut; % 14 choices for the transmit burst frequency
        samplesPerWave = 4;
    case 'LF'
        RcvProfile.PgaHPF = 0; % disable 80 KHz highpass from PGA integrator
        RcvProfile.LnaGain = 18; % 15, 18, 24 (dB) are the only options
        RcvProfile.LnaHPF = 100; % 0, 50, 100, 150, 200 (kHz) with 0 disabling the filter
        RcvProfile.LnaZinSel = 20; % 0-31 (int) 0=lowest 30=high 31=highest(off; no feedback)
        txFreqList = FreqOut/10; % 14 choices for the transmit burst frequency
        samplesPerWave = 4;
end



%% Initialize and Set system parameters

% Receive Profile
RcvProfile.LnaGain = 24;
RcvProfile.AntiAliasCutoff = AAF(txFreqIndx);
RcvProfile.LnaZinSel = 31;

rcvSampleRate = 62.5;

%% Specify Trans structure array.
connectorSelect = 1;
if isempty(SysConfig.UTAtype)
    SysConfig.UTAtype = [1, 1, 1, 0];
end
UTA = computeUTA(SysConfig.UTAtype, connectorSelect);
% this call to computeUTA will tell us how many channels are usable on this
% system configuration, if we have to create a dummy Trans structure
numES = length(UTA.TransConnector); % number of element signals supported by this system configuration

% If a Trans structure has not been provided, create a dummy one here
if nargin ~= 2 || ~isstruct(Trans)
    Trans.name = 'custom';
    Trans.id = -1;
    Trans.numelements = numES;
    Trans.ElementPos = zeros(numES, 4); % for simulation, put all elements at same location for identical round trip path length
    Trans.maxHighVoltage = 90;
    Trans.impedance = 53; % value presented by the load box at low frequencies
    Trans.units = 'mm';
    Trans.connType = -1;
end
if ~isfield(Trans, 'numelements') || isempty(Trans.numelements)
    Trans.numelements = numES;
    Trans.ElementPos = zeros(numES, 4);
    Trans.id = -1;
end
Trans.frequency = rcvSampleRate / 4;  % override Trans.freqency with test freqency we will be using.
Trans.ElementSens = ones(1, 101); % dummy array required by VSX but will not be used since there is no recon function
Trans.type = 0; % another dummy value that VSX requires but irrelevant for this application

assignin('base', 'Trans', Trans); % will be accessed by computeUTA

startDepth = startDepthUsec * Trans.frequency;   % start of Acquisition depth range in wavelengths
acqRange = acqRangeUsec * Trans.frequency;
endDepth = startDepth + acqRange;
linelength = 2*samplesPerWave*acqRange; % acquisition duration in samples

% Specify Media object. Single point target.
Media.MP = [0 0 2 1];
Media.numPoints = 1;
Media.function = 'movePoints';

UTA = computeUTA(SysConfig.UTAtype, connectorSelect);
% this call to computeUTA will tell us how many channels are usable on this
% system configuration, based on the Trans structure.
numCH = UTA.numCh;
numES = Trans.numelements; % number of elements used by this system configuration; size of Apod arrays
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.Connector = connectorSelect;
Resource.Parameters.verbose = 0;

Resource.Parameters.numTransmit = numCH;      % number of transmit channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.speedCorrectionFactor = 1.0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.
Resource.VDAS.halDebugLevel=0; % 0 to disable, 5 to dump all HAL output

Resource.Parameters.fakeScanhead = 1; % allow system HW operation with nothing connected

% connectorSelect = 1:SysConfig.UTAtype(3); % select all connectors that exist on UTA module
%% Specify Resources.
% Receive buffer 1 is for the Live Acquisition and RF data plot
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 2*numAper * numBist * linelength; %  acq. per frame
Resource.RcvBuffer(1).colsPerFrame = numCH;
Resource.RcvBuffer(1).numFrames = numFrames;   % as set above

Resource.VDAS.dmaTimeout = 10000;

%% Specify Transmit TW and TX structures.

% TW using Parametric waveform
txFreq = txFreqList(txFreqIndx);
TW(1).type = 'parametric';
TWA = round(125/txFreq); % "TWA" value that will go into pulsecode
txFreq = 125/TWA; % initial transmit frequency value, with quantization applied
TWC = 6;%ceil(2 * bdurUsec * txFreq); % number of half-cycles in waveform bdurUsec long
TW(1).Parameters = [txFreq, 0.67, TWC, 1];   % A, B, C, D  for specified transmit frequency

Trans.Bandwidth = [0.6, 1.4] * txFreq;
% Specify TX structure array
delay = zeros(1,numES) + 1.5 * (2*startDepthUsec + 2) * Trans.frequency;


TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'Apod', ones(1,numES), ...
                   'focus', 0.0, ...
                   'Steer', [0.0,0.0], ...
                   'VDASBistDriveEnable', 1, ...
                   'Delay', delay), 1, numBist); % numBist TX structures


% enable Bist drivers, alternating polarity from board to board
for bistnum = 1:numBist
%     if VDAS && numES >= 128
%         TX(bistnum).VDASBistDriveSelect = bistnum; % Bist drive from one end only
%         TX(bistnum).Apod = repmat([ones(1, 64), -ones(1,64)], 1, numES/128);
%     else
        TX(bistnum).VDASBistDriveSelect = bistnum; % Bist drive from one end only
%     end
end




%% Specify TGC Waveform structure.
TGC.CntrlPts = initTGC*ones(1,8); % default to maximum gain flat TGC curve
TGC.rangeMax = endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

%% Specify Receive structure arrays -

Receive = struct('Apod', zeros(1, numES), ...
                'startDepth', startDepth, ...
                'endDepth', endDepth, ...
                'TGC', 1, ...
                'bufnum', 1, ...
                'framenum', 1, ...
                'acqNum', 1, ... %'samplesPerWave', samplesPerWave, ...
                'mode', 0, ...
                'callMediaFunc', 0);

Receive.InputFilter = [zeros(1,20), 1]; %all-pass

% populate all of the Receive structures
Receive = repmat(Receive,1,(numAper*numBist*numFrames));

% define subaperture Apod arrays
SubApod = zeros(numAper, numES);
activeES = ceil(numES/numAper); % number of active elements in each subaperture
% note fourth subaperture will be smaller if numES is not a multiple of 4
for i = 1:numAper
    aperst = (i-1)*activeES +1; % first ES in subaperture 
    aperend = min(i*activeES, numES); % last ES in subaperture
    SubApod(i, aperst:aperend) = 1;
end
% - Set event specific Receive attributes
for frnum = 1:numFrames
    for bistnum = 1:numBist
        rindx = (frnum-1)*numAper*numBist + (bistnum-1)*numAper;
        for apnum = 1:numAper
            % numAper acquisitions per frame
            Receive(rindx + apnum).framenum = frnum;
            Receive(rindx + apnum).acqNum = (bistnum-1)*numAper + apnum;
            Receive(rindx + apnum).Apod = SubApod(apnum, :);
        end
    end
end

clear endDepth LPC BPF samplesPerWave

%% Define Processing functions
% External RF plot function.
Process(1).classname = 'External';
Process(1).method = 'ZtstPlots';
Process(1).Parameters = {'srcbuffer','receive',... % name of buffer to process.
                         'srcbufnum',1,...
                         'srcframenum',-1, ...
                         'dstbuffer','none'};



%% Specify SeqControl structure arrays.
%  - Jump back to start.
SeqControl(1).command = 'jump';
SeqControl(1).argument = 1;

SeqControl(2).command = 'timeToNextAcq';  % Frame rate
SeqControl(2).argument = 15000;  % 40 msec or 25 Hz acquisition rate

SeqControl(3).command = 'triggerOut';

SeqControl(4).command = 'timeToNextAcq';  % time between numAper acquisitions
SeqControl(4).argument = 1000/PRF;  % convert PRF in KHz to usec

SeqControl(5).command = 'returnToMatlab';

SeqControl(6).command = 'setRcvProfile';  % select a receive profile
SeqControl(6).argument = 1;  % select ibitial profile

SeqControl(7).command = 'sync';
SeqControl(7).argument = 1e7; % 10 second timeout



nsc = 8; % nsc is count of SeqControl objects


%% Specify Event Sequence
n = 1; % n is count of Events

% Acquire all frames defined in RcvBuffer
for frnum = 1:numFrames
    for bistnum = 1:numBist
        rindx = (frnum-1)*numAper*numBist + ((bistnum-1))*numAper;
        for apnum = 1:numAper
            Event(n).info = 'sub-aperture accumulate acquisitions';
            Event(n).tx = bistnum;
            Event(n).rcv = rindx + apnum;    % unique Rcv each frame and subaperture.
            Event(n).recon = 0;      % no reconstruction.
            Event(n).process = 0;    % no processing
            Event(n).seqControl = 4; % time between the numAper acqs.
            n = n+1;
        end
    end
    Event(n-1).seqControl = [2, nsc]; % time between frames and transfer to host
       SeqControl(nsc).command = 'transferToHost';
       nsc = nsc + 1;

    Event(n).info = 'Dummy transmit';
    Event(n).tx = 1;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = 4; % ttna
    n = n+1;

    Event(n).info = 'sync';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = 7; % sync.
    n = n+1;

    Event(n).info = 'Process RF data';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 0;      % reconstruction
    Event(n).process = 1;    % processing
    Event(n).seqControl = 0;
    n = n+1;


    Event(n).info = 'return to matlab';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 0;      % reconstruction
    Event(n).process = 0;    % processing
    Event(n).seqControl = 5;
    n = n+1;


end

Event(n).info = 'Jump back to repeating sequence';
Event(n).tx = 0;        % no TX
Event(n).rcv = 0;       % no Rcv
Event(n).recon = 0;     % no Recon
Event(n).process = 0;
Event(n).seqControl = 1; % jump command
n = n+1;

clear  n nsc frnum

%% User specified UI Control Elements

if ~TestAcq
    % create GUI TX Frequency control only if TestAcq is true
    import vsv.seq.uicontrol.VsSliderControl;


    UI(1).Control = VsSliderControl('LocationCode','UserC7','Label','TX Frequency',...
                      'SliderMinMaxVal',[1, 28, txFreqIndx],...
                      'SliderStep',[1/27, 1/27],'ValueFormat','%1.2f', ...
                      'Callback', @TxFreqCallback );
end



%% Cleanup, save .mat file, run VSX

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 1;

% External function definitions.
EF(1).Function = vsv.seq.function.ExFunctionDef('ZtstPlots', @ZtstPlots);

% Save all the structures to a .mat file.
save('MatFiles/probeZcheck');

% The following two lines will automatically invoke VSX, after the setup
% script has finished creating the .mat file
evalin('base', 'filename = ''probeZcheck''; VSX');
ProbeZdata.medianRMSlvl = evalin('base', 'medianChRMS');
ProbeZdata.ESrmsTot = evalin('base', 'ESrmsTot');
ProbeZdata.Trans = evalin('base', 'Trans');
ProbeZdata.txFreq = txFreq;
delete(findobj('tag','UI'));
if ~TestAcq
    evalin('base', 'clear'); % clear base workspace from having run VSX
end
end

%% Callback routines for GUIcontrols.

% - Transmit frequency callback (used only when TestAcq is true)
function TxFreqCallback(~, ~, UIValue)
    txFreqList = evalin('base', 'txFreqList');
    bdurUsec = evalin('base', 'bdurUsec');
    txFreqIndx = round(UIValue);
    assignin('base','txFreqIndx', txFreqIndx);
    txFreq = txFreqList(txFreqIndx);
    TWA = round(125/txFreq);
    assignin('base','TWA', TWA);
    txFreq = 125/TWA;
    txFreqtxt = num2str(txFreq,'%2.3f');
    h = findobj('tag','UserC7Edit');
    set(h,'String',txFreqtxt);
    assignin('base','txFreq',txFreq);
    TW = evalin('base','TW');
    TW(1).Parameters(:,1) = txFreq;
    TW(1).Parameters(:,2) = .67;
    TW(1).Parameters(:,3) = 6;%ceil(2 * bdurUsec * txFreq); % number of half-cycles
    assignin('base','TW',TW);
    evalin('base', 'RcvProfile.AntiAliasCutoff = AAF(txFreqIndx);');
    evalin('base', 'RcvProfile = computeRcvProfile(RcvProfile);');
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'TW', 'RcvProfile'};
    assignin('base','Control', Control);
end


%% External Function ZtstPlots

function ZtstPlots(RDatain)
    persistent RFsampleplothndl RFplotrange probeZplothndl ...
         chPlotRMSavg probeZplotRange medianChRMS
    drawnow;

    numCH = evalin('base', 'numCH');
    numES = evalin('base', 'numES');
    SubApod = evalin('base', 'SubApod');
    numAper = size(SubApod, 1);

    ADfreq = evalin('base', 'Receive(1).decimSampleRate');
    numBist = evalin('base', 'numBist');
    numAcq = numBist;


    rcvES = evalin('base','rcvES');  % Channel no. to plot
    showPlots = evalin('base', 'showPlots');

    PinNames = evalin('base', 'UTA.ChPinNames');
    TransConnector = evalin('base', 'UTA.TransConnector');

    Trans = evalin('base', 'Trans');
    linelength = evalin('base','linelength'); % length of a line of rcv data

    RDataPerAcq = zeros(linelength, numES, numAcq); % Define RData buffer with columns for all receive channels
    EnvP = RDataPerAcq;
    EnvN = RDataPerAcq;
% %     RDataB2 = RDataB1;
    % data is converted from 16 bit integers to doubles, and scaled to A/D
    % sample units with full scale = +/- 8192 (14 bit A/D).  Since Receive data
    % comes out with the 14 bit A/D data mapped to the middle 14 of the 16
    % bits, we need to divide by two to get to an integer representation of the
    % A/D sample values.
    for acqnum = 1:numAcq
        for apernum = 1:numAper
            offset = ((acqnum-1)*numAper + apernum-1)*linelength;
            for es = 1:numES
                if SubApod(apernum, es)
                    ch = Trans.Connector(es);
                    RDataPerAcq(:, es, acqnum) = double(RDatain(offset+(1:linelength), ch))/2; % convert to double and scale
                end
            end
        end
    end
    
    
    



    X = 1:linelength; % X spans full line length; plot range is set by Xlim parameters

    plotdepth = round(linelength * evalin('base','plotdepth')); % depth and range to be plotted, converted to samples
    plotrange = round(max( 50, linelength * evalin('base','plotrange'))); % minimum range of 50 samples
    xlim1 = round(max(1, min(plotdepth-plotrange/2, linelength-plotrange))); % plotdepth sets center sample, but restricted to preserve full plotrange
    xlim2 = round(min(linelength, xlim1+plotrange-1)); %
    plotxlim = [xlim1, xlim2];

    assignin('base','plotxlim', plotxlim);
    
    if (isempty(medianChRMS)||~exist('medianChRMS', 'var'))
        chPlotRMSavg = 0;
%         chPlotPeakavg = 0;
        
        medianChRMS = 0;
    end

    % initialize figures and parameters on the first call at startup:
    if showPlots && (isempty(RFsampleplothndl)||~ishandle(RFsampleplothndl))
        RFplotrange = 13;
        figure('Position',[20 500 800 300],'Name','Receive Buffer Data Plot for selected channel');
        RFsampleplothndl = axes('XLim',plotxlim,'YLim',[-1.1*2^RFplotrange 1.1*2^RFplotrange],'NextPlot','replacechildren');
        figure('Position',[20 80 800 300],'Name','RMS signal Level from each Element Signal');
        probeZplotRange = 16000;
        probeZplothndl = axes('XLim',[1 numCH], 'YLim',[0 probeZplotRange], 'NextPlot','replacechildren');


    end

    txFreq = evalin('base', 'txFreq');


    if evalin('base','nullDC')
        for acqnum = 1:numAcq
            for es=1:numES
                DC = sum(RDataPerAcq(:, es, acqnum))/linelength;
                RDataPerAcq(:, es, acqnum) = RDataPerAcq(:, es, acqnum)-DC;
            end
        end
    end
    
    % Compute envelope signals
    for acqnum = 1:numAcq
        for es=1:numES
            [EnvP(:, es, acqnum), EnvN(:, es, acqnum)] = envelope(RDataPerAcq(:, es, acqnum),50,'analytic');
        end
    end

    chplot1 = RDataPerAcq(:, rcvES, 1);
% %     chplot1FltB1 = RDataFltB1(1:intrpLnLgth,rcvES);
% %     chplot1FltB2 = RDataFltB2(1:intrpLnLgth,rcvES);

    %  We want to evaluate RMS level over an integer number of cycles that
    %  represent about 4 usec total duration
    sampPerCyc = ADfreq / txFreq; % samples per cycle (should be an integer)
    numCyc = round(4 * txFreq); % number of cycles in 4 usec, to nearest integer
    Nmeas = round(numCyc * sampPerCyc); % total samples to measure

    % Burst starts about 2 usec into acquisition interval and is about 12
    % usec long.  Within the burst we want to position the 4 usec
    % measurement interval from about usec 6 to 10, so start at 8 usec into
    % the interval:
    Ms1 = round(8 * ADfreq); % delay in samples to first measurement sample, for 8 usec
    Ms2 = Ms1 + Nmeas - 1; % last measurement sample

    % compute RMS level over the Nmeas samples for each element signal
    ESrmsTot = zeros(1, numES);
    ESenvPk = zeros(1, numES);
%     ESrmsB2 = ESrmsB1;
% %     EstZtot = ESrmsB1;
    for acqnum = 1:numAcq
        for ESnum = 1:numES
            ESrmsTot(ESnum) = ESrmsTot(ESnum) + sqrt(sum(RDataPerAcq(Ms1:Ms2, ESnum, acqnum).*RDataPerAcq(Ms1:Ms2, ESnum, acqnum))/Nmeas);
            ESenvPk(ESnum) =  ESenvPk(ESnum) + max(EnvP(:, ESnum, acqnum)) - min(EnvN(:, ESnum, acqnum));
        end
    end
%     ESRMStot = ESrmsB1 + ESrmsB2;
    assignin('base', 'ESrmsTot', ESrmsTot);
    assignin('base', 'ESenvPk', ESenvPk);




    RFpeak=max(16, max(abs(chplot1(xlim1:xlim2)))); % for autoscaling the sample plot Y axis

    Txstr = 'Transmit through BIST network driver';

%     chPlotRMS = sqrt(sum(chplot1(xlim1:xlim2).*chplot1(xlim1:xlim2))/Nsamp);
%     chPlotPeak = (max(chplot1(xlim1:xlim2)) - min(chplot1(xlim1:xlim2)))/2; % avg of pos and neg peaks
    if chPlotRMSavg == 0
        chPlotRMSavg = ESrmsTot(rcvES);
    else
        chPlotRMSavg = chPlotRMSavg + (1/10)*(ESrmsTot(rcvES) - chPlotRMSavg);
    end

    if ~isempty(PinNames)
        rcvch = TransConnector(rcvES);
        rcvPinTxt = ['  Connector Pin ', PinNames{rcvch, 2}, '  ', PinNames{rcvch, 1}];
    else
        rcvPinTxt = ' ';
    end
    

    newmedianChRMS = median(ESrmsTot);
    
    if abs((newmedianChRMS - medianChRMS)/newmedianChRMS) > .02
        medianChRMS = 0.985 * newmedianChRMS;
    else
        medianChRMS = 0.80 * medianChRMS + 0.2 * newmedianChRMS;
    end
    assignin('base', 'medianChRMS', medianChRMS);

    filtDif = abs((newmedianChRMS - medianChRMS)/newmedianChRMS);

    if filtDif < .1 || newmedianChRMS < 1
        if showPlots
    
            plot(RFsampleplothndl,X, RDataPerAcq(:, rcvES, 1), X, EnvP(:, rcvES, 1),  [Ms1, Ms1], [-2^14, 2^14], [Ms2, Ms2], [-2^14, 2^14]);
        %     plot(RFsampleplothndl, X, chplot1, [4*txDly, 4*txDly], [-2^13, 2^13], [Ms2, Ms2], [-2^13, 2^13]);% % %, X, env); % plot channel data from pulse
            title(RFsampleplothndl,{['Receive data plot for Receive Element Signal ', num2str(rcvES,'%.0f'), rcvPinTxt];...
                Txstr; ['Plot data RMS level ', num2str(chPlotRMSavg,'%.1f'), ' LSB steps']});
            set(RFsampleplothndl,'XLim',plotxlim);

            if evalin('base', 'autoscale')
                if RFpeak > 1.1*2^RFplotrange
                    RFplotrange = min(13, RFplotrange+1);
                    set(RFsampleplothndl,'YLim',[-1.1*2^RFplotrange 1.1*2^RFplotrange]);
                elseif RFpeak < 0.4*2^RFplotrange
                    RFplotrange = RFplotrange-1;
                    set(RFsampleplothndl,'YLim',[-1.1*2^RFplotrange 1.1*2^RFplotrange]);
                end
            end





            plot(probeZplothndl, 1:numES, ESrmsTot,  1:numES, 0.125*ESenvPk);

                title(probeZplothndl, {['RMS level for each Element Signal, at ', num2str(txFreq,'%.2f'), ' MHz']; ...
                    ['Median level ', num2str(medianChRMS,'%.0f'), ' A/D steps.']});

            if evalin('base', 'autoscale')
                maxZ = max(ESrmsTot);
                if maxZ > 0.9 * probeZplotRange
                    probeZplotRange = 2 * probeZplotRange;
                    set(probeZplothndl,'YLim',[0, probeZplotRange]);
                elseif maxZ < 0.4 * probeZplotRange
                    probeZplotRange = 0.5 * probeZplotRange;
                    set(probeZplothndl,'YLim',[0, probeZplotRange]);
                end
            end
        end

        if evalin('base', '~TestAcq')
            evalin('base', 'vsExit = 1;');
        end
    end

end
