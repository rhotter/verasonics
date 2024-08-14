function varargout = spectralDoppler(varargin)
%% spectralDoppler  (main) Spectral Doppler processor for pulse Doppler data.
% Copyright 2001-2017 Verasonics, Inc.  All world-wide rights and remedies under all intellectual property laws and industrial property laws are reserved.  Verasonics Registered U.S. Patent and Trademark Office.
%
% DESCRIPTION
%     This routine processes pulse data to obtain classical spectral
%     Doppler output.
%       Usage:
%       1. Environment: Asynchronous processing runacq.c paradigm; called
%       with external process/m-file protocol (See Sequence Programming Manual).
%
%       2. >> spectralDoppler('cleanup');
%       Sets this program to uninitialized state and destroys audio driver
%       object.
%
%       3. >> p=spectralDoppler('spectrogramParameters',framePeriod,dopSweepTimeIndex);
%       Returns spectrogram parameters given frame period in seconds and
%       sweep time control number of 1,2,3, or 4 (fastest to slowest).
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%       Uses portaudio-based audio driver, mexaudiostream.c.
%
% Assumptions, conditions, and behavior:
%   1. This function is invoked at fixed intervals (usually 20 mSec).
%   2. This function draws lines to the spectral display MatLab window at
%   sweep intervals of 10ms, 20ms, 30ms etc. (multiples of invocation).
%   3. Operating Notes:
%   MACI - Leopard, MatLab R2009b(32bit): must run from shell, invoking matlab with -nojvm option
%       Otherwise "drawnow" causes excessive delays and glitching
%   MACI - SnowLeopard, 32bit, Matlab R2009b(32bit):  runs well from MatLab Gui launch
%   MACI - SnowLeopard, 32bit, Matlab R2010a(32bit):
%               drawnow - 20ms/invocation, occasionally 30ms
%               pause(.0001) - 13ms per invocation, occasionally 30 ms
%               no pause or drawnow - no glitches, but not very smooth.
%%    Vista 64 - runs well.
%
% REVISION HISTORY
%       April 2019 VTS-1142 IQ buffer replaced with separate I and Q buffers
%       V1.0:  12/18/2009 JAF stripped down from pwrtas_vs.m
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

startTime_clock = clock;
tic;
drawnowUpdateFrameTrigger = 5;

verbose = 0;
testDataGen = 0; %set for fake test data gen
captureData = 0; %set to one to capture IQ data columns in iqsave_xxx.mat file
%(is reset later in this file)
%.
audioDataSave = 0;  %save audio data path diagnostics
%can define other capture data sets in switch statement

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent HdWhp
persistent k_chunk sgp subFrameInd IQstate
persistent sounddata sounddataRS sounddataLR %make persistent so reallocation time is reduced
persistent rsParams channelizerState
persistent NoccStartThresh
persistent FsDAC dtDAC
persistent wfState rsState lpfState stateSpecPersist statePhasor
persistent bLPF aLPF startedAudioStream
persistent PWdopPRF
persistent baselineSave dopWFCutoffNormSave dopSweepTimeIndexSave %use to auto generate the update flag
persistent noiseFloorHistory noiseFloor
persistent spectralDoppler_state launchTime captureVarsString
persistent specWindow SdatAll
persistent computerType invokeCount
persistent sounddataRSCell sounddataCell invocationCount    %debug
persistent prevStartTime startTimeHistory

FrameCountMax = 700;

if isempty(invocationCount)
    invocationCount = 0;
    startTimeHistory = zeros(1, FrameCountMax);
end
invocationCount = invocationCount+1;


startTime = tic;
if ~isempty(prevStartTime)
    invocationInterval = toc(prevStartTime);
else
    invocationInterval = 0;
end
prevStartTime = startTime;

startTimeHistory(mod(invocationCount,FrameCountMax )+1 ) = invocationInterval;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%signature parsing
comm = '';
if ischar(varargin{1})
    %DUAL COMMANDS
    %first argument is a command:
    comm = varargin{1};
    if length(varargin)>1
        %second argument is a command:
        if ischar(varargin{2})
            comm = varargin{2};   %override
        end
    end
else
    %SINGLE COMMAND but in second argument
    if length(varargin)>1
        %second argument is a command:
        if ischar(varargin{2})
            comm = varargin{2};   %override
        end
    end
end

%execute subfunction specified in command string variable
if ~isempty(comm)

    %is comm a subfunction in this file?
    [~, fff]=fileparts(which(comm));
    issubfunc=isequal(fff, mfilename);
    if ~issubfunc
        error('unrecognized command')
    end

    [varargout{1:nargout}]=feval(comm,varargin{2:end});

    return %end of command processing mode

end %signature check



if isempty(spectralDoppler_state)
    spectralDoppler_state = 'not.initialized';
end
if isempty(computerType)
    computerType=computer;
    invokeCount = 0;
end
invokeCount=invokeCount+1;


%.
% .. . . . . . . . . .

wkSpace = 'caller';
scalePortAudio = 10.0; %nominal volume control.
ScaleSpec = 1.0; %adjust spectral
narrowbanding = 0;
NumNoiseHist = 30;%noise floor estimation median window size

switch nargin

    case 2 % - - - - - - - - -
        %Called by extern object in the script

        %get the IQ data: (VTS-1142 separates I and Q to two buffers, two separate input arguments)
        IQin = complex(varargin{1}, varargin{2});
        IQdata = squeeze(IQin);
%         IQdata = squeeze(varargin{1}); % original scheme with complex IQ input buffer
        sizeIQ = size(IQdata);
        %need to get because PRF change makes different sizes:
        nPulses = evalin(wkSpace,'nDopPRIs'); % get no. of PRIs in a 'frame'
        SDop.nPulses = nPulses;
        nRange = sizeIQ(1);
        if narrowbanding
            NumDepthSummation =nRange;
            IQdata = sum(IQdata(1:nRange,1:nPulses),1);
        else
            rangeSubSampleIndex = 1;
            NumDepthSummation =1;
            IQdata = IQdata(rangeSubSampleIndex,1:nPulses);
        end
        %resized now so get new size:
        sizeIQ = size(IQdata);

        if testDataGen
            IQdata = generateTestData('tones.clutter.slow',nPulses);
            if rand(1)< 0*.01
                disp([mfilename,' %%%%%%%%45 '])
                keyboard
            end
        end

        if strcmp(spectralDoppler_state,'initialized')
            SDop = evalin(wkSpace, 'SDop');
        else
            SDop = initSDop;
        end

        %parameter passing:
        %assuming that SDop is not being used by script callbacks to pass
        % parameters to this program (legacy interface):
        SDop.TFrame = evalin(wkSpace, 'framePeriod');
        SDop.PWdopPRF = evalin(wkSpace, 'dopPRF');
        dopWFCutoffNorm = evalin(wkSpace, 'dopWFCutoffNorm');
        sweepSpeedIndex = evalin(wkSpace, 'dopSweepTimeIndex');
        baseLineShiftNorm = evalin(wkSpace, 'baseLineShiftNorm');
        SDDynRangeDB = evalin(wkSpace, 'SDDynRangeDB');
        SDop.nIqfft = evalin(wkSpace, 'NWindow');
        SDop.despeckle = evalin(wkSpace, 'SDDespeckle');

        %set initialization mode flags:
        PWPRFnew = ~isequal(PWdopPRF,SDop.PWdopPRF);
        PWdopPRF = SDop.PWdopPRF;
        %.
        WFnew = ~isequal(dopWFCutoffNorm,dopWFCutoffNormSave);
        dopWFCutoffNormSave =  dopWFCutoffNorm;
        %.
        Sweepnew = ~isequal(sweepSpeedIndex,dopSweepTimeIndexSave);
        dopSweepTimeIndexSave =  sweepSpeedIndex;
        %.
        BLnew = ~isequal(baseLineShiftNorm,baselineSave);
        baselineSave = baseLineShiftNorm;

        nuIQRot = baseLineShiftNorm + 0.5;

    otherwise
        error('bad signature switch')
end %nargin switch

if ~strcmp(spectralDoppler_state,'initialized')
    disp([mfilename,': initializing state.'])
    PWPRFnew = 1; %override

    [sounddataRSCell, sounddataCell] = deal({}); %DEBUG


    clear mexaudiostream
    clear audiostream
    spectralDoppler_state = 'initialized';
    launchTime = startTime_clock;
end

enableAudio = SDop.audioOn;

%disable reads/writes from DAC fifo (for testing in VDAS simulation mode).
disableDACDriver = enableAudio~=1 ; %set for no audio
disableSpectralDisplay  = SDop.displayOn~=1; %set for no display

FrameRate = 1.0./SDop.TFrame;
displayTime = NaN;
Nwind = SDop.nIqfft;

%% --- Initialize the sweeping display window and Wall filter
if PWPRFnew || Sweepnew || BLnew || WFnew  % First time through, or any time PRF changes
    if verbose
        fprintf('\n')
    end

    k_chunk = 0;
    if PWPRFnew
        rsState = [];
    end

    %calc spectrogram parameters:
    sgp=spectrogramParameters(SDop.TFrame,sweepSpeedIndex);
    THist = sgp.THist;
    %partition the data chunk into subrames, one for each line:
    subFrameInd = subframeInd100(sgp,PWdopPRF); %get subframe indices within the iq chunk
    stateSpecPersist =[];
    if Sweepnew
        specWindow=[];
    end

    [noiseFloorHistory,noiseFloor]=deal([]);

    if verbose
        fprintf('I')
    end
    if enableAudio
        FsDAC = mexaudiostream(908); %get output sample rate for audio
        dtDAC = 1/FsDAC;

        DACFrameSize = mexaudiostream(920); %size of dac frame in stereo samples
        TDACFrame = DACFrameSize * dtDAC ;
        DACFramesPerInvoke =  SDop.TFrame/TDACFrame ; %non integer value
        %two floats per audio sample:
        NoccStartThresh = 1.0*ceil(DACFramesPerInvoke)*DACFrameSize; %how full we want the fifo.
        audioIsOpen = mexaudiostream(3);
        if ~audioIsOpen
            audiostream( 'open' )
        end
    end
    startedAudioStream = 0;

    [ channelizerState , wfState , rsState ] = deal([]);

    NumDataFrames =  sgp.NumDataFrames;
    nLines =  sgp.NumLines; %actual interpolated lines drawn

    %init spectrogram history if size has changed
    if ~isequal(size(SdatAll),[SDop.nfft , nLines])
        SdatAll = zeros( SDop.nfft , nLines);
    end

    %make a new figure
    hFig = findobj('Tag', 'sdop_figure');
    if ~isempty(hFig), close(hFig), end
    TXFreq_MHz = evalin(wkSpace, 'Trans.frequency');
    speedOfSound_mps = evalin(wkSpace, 'Resource.Parameters.speedOfSound');
    [hFig,hAxs,HImg] = sdop_gui_lite(THist, nLines, SDop.nfft, PWdopPRF, baseLineShiftNorm ,TXFreq_MHz,speedOfSound_mps );
    SDop.hFigure = hFig;
    SDop.hAxes = hAxs;
    SDop.HSpectrogram = HImg;

    %sliding FFT window buffer: make sure it's big enough
    [IQstate]=circbuff('initialize', SDop.nfft + 160 ); %corresponds to 8000 max prf

    SDop.isInit = 1;

    assignin(wkSpace, 'SDop', SDop); %pass SDop struct back to caller (legacy interface)

    switch captureData %initialize:
        case 0
            %do nothing
        case 1
            IQsave(zeros(nPulses,1),500);%init
        case 2
            captureVars = {...
                'invokeTime','runningTime','displayTime', ...
                'PWdopPRF' ...
                };
            captureVarsString = vars2VecEvalStr(captureVars);
            IQsave(zeros(length(captureVars),1),100,captureVars);%init

        case 3
            captureVars = {...
                'invokeTime','runningTime','displayTime', ...
                'PWdopPRF' ...
                };
            captureVarsString = vars2VecEvalStr(captureVars);
            IQsave(zeros(length(captureVars),1),256,captureVars); %init
        case 4
            captureVars = {...
                'invokeTime','runningTime','displayTime', ...
                'drawnowTime' , 'PWdopPRF' ...
                };
            captureVarsString = vars2VecEvalStr(captureVars);
            IQsave(zeros(length(captureVars),1),256,captureVars);%init

        otherwise
            disp([mfilename,':error: bad capture option'])
    end

    audMute = muter(FrameRate,1); %mute audio

else %initialize
    audMute = muter(FrameRate);
end %initialize
IQSlideBlock = zeros(Nwind,sgp.R);
% --- High pass WALL filter design
if WFnew
    [fproto]= getFilterPrototype('wallfilter');
    [HdWhp.Numerator,HdWhp.Denominator]=filterTranslate(...
        fproto.tf.b,fproto.tf.a,fproto.nuNominal,dopWFCutoffNorm,'high');
end
if ~isequal(size(specWindow),[Nwind,sgp.R])
    w = hammingvs(Nwind);
    w = w/sum(w);
    specWindow = w(:)*ones(1,sgp.R);
end

SDProcGainDB = 10*log10(Nwind*NumDepthSummation);
baseLineBin = min(SDop.nfft,max(1,round(nuIQRot*SDop.nfft)+1));

audCompressionFactor= SDop.audCompressionFactor ;%higher is more compression (low level gain).
SDCompressionMethod=SDop.SDCompressionMethod ;%spectral display compression method
if enableAudio
    %Update resampler params to PRF= PRFNominal*adjustment
    %until fifo accupancy is ok.
    %Then update resampler params to PRF = PRFEst.
    %%%% fifoState = audiostream('fifostate');

    %initialize resampler state:
    if isempty(rsState)

        rsParams = resamplingParams(PWdopPRF,FsDAC);
        if verbose
            disp([mfilename,':initializing resampler state.'])
        end

        analogResamplerDesign = 'rc.23'; %reasonable (but not optimal) analog resampler design.
        rsState = rs_irrat('initialize', ...
            PWdopPRF , rsParams.FsIntermediate ,nPulses,analogResamplerDesign);
        %second-stage resampler (integer upsample) reconstruction filter
        [fproto]= getFilterPrototype('audio_resampler');
        [bLPF,aLPF]=filterTranslate(...
            fproto.tf.b,fproto.tf.a,fproto.nuNominal,1/(2*rsParams.Qupsample),'low');
        %%%% [bLPF,aLPF] = ellip(5,1.5,60,2/(2*rsParams.Qupsample));

    end %init resampler state

end% if audio on


if exist('SDop','var')~=1
    SDop = evalin(wkSpace, 'SDop');
end

audScale = SDop.audioGain * audMute * scalePortAudio;

%finished initialization - - - - - -  - - - - - - - - - - - - -

%SECTION 1: WALL FILTERING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- High pass filter the entire frame of data (Wall Filter)
[IQ,wfState] = filter(HdWhp.Numerator, HdWhp.Denominator, IQdata,wfState);

%SECTION 2: AUDIO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a. channelize doppler in to pos and neg freqs.
% b. convert sample rate by arbitrary fractional amount to nearest rate
% related rationally to DAC sample rate.
% c. upsample by rational fraction rate to DAC sample rate.

% send sound samples for one frame to audio driver
if enableAudio
    %%%%%%%%%%%%%%%%%%
    % --- channelizer operations
    [sounddata,channelizerState]= doppler_channelize(IQ,channelizerState);

    % --- convert PRF to next higher freq. that divides the DAC
    % sample playback rate
    [sounddataRS ,rsState] = rs_irrat('resample',sounddata.',rsState);

    if audioDataSave
        sounddataRSCell{end+1}=sounddataRS; %debug
        sounddataCell{end+1}=sounddata;
        if invocationCount==FrameCountMax
            disp([mfilename,': debug save RS audio: audiors_debug.mat'])
            fifodata=mexaudiostream(901,15000);
            save('audiors_debug','sounddataRSCell','sounddataCell','fifodata','rsState')
            invocationCount = 0;
            [sounddataRSCell, sounddataCell] = deal({}); %DEBUG

        end
    end

    % --- upsample and LPF by the integer amount to give the DAC sample rate:
    sounddataRS = upsamplevs(sounddataRS,rsParams.Qupsample);

    [sounddataRS,lpfState] = filter(audScale*bLPF,aLPF, sounddataRS,lpfState);

    %audio compression
    sounddataLR = limiter('soft', ...
        [real(sounddataRS(:)),imag(sounddataRS(:))] , ...
        audCompressionFactor );

    %send data to DAC interface fifo:
    if disableDACDriver
        %for testing on non-realtime environment
        disp([mfilename,':DEBUG: disabled DAC fifo write.'])
    else
        audiostream('load' , sounddataLR);
    end

    if ~startedAudioStream
        if disableDACDriver
            %for testing on non-realtime environment
            disp([mfilename,': DAC start call disabled.'])
        else
            %check size to start:
            Nocc = 0.5*mexaudiostream(903); %in stereo samples
            if Nocc>=NoccStartThresh
                if verbose
                    %disp([mfilename,': starting audio stream.'])
                    fprintf('S')
                end
                Ncap = (1/2)*mexaudiostream(905); %capacity in stereo samples

                startFifoFraction=NoccStartThresh/Ncap ; %fraction of fifo occupancy needed to start audio
                audiostream('start', startFifoFraction );
                startedAudioStream = 1;
            end
        end
    end

end %if audioOn

% SECTION 3: SPECTRAL DISPLAY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~disableSpectralDisplay

    k_chunk = rem(k_chunk,sgp.NumDataFrames)+1;

    %rotate the frame of data
    [fazor,statePhasor]= phasor( nPulses, nuIQRot, statePhasor);
    IQrotated = IQ.*fazor;

    %spectral processing
    for r=1:sgp.R %collect enough data to draw lines for input chunk
        %R is the spectral slowtime interpolation factor
        % adjust the baseline, using a circular shift
        IQseg = IQrotated(subFrameInd{r});
        %write a subsegment to the buffer
        %(load buffer structure with newest iq data subchunk):
        IQstate = circbuff('write',IQstate,IQseg(:).');
        IQsegSlide = circbuff('read',IQstate,Nwind);%get most recent iq data from buffer
        IQSlideBlock(:,r)=IQsegSlide(:); %load multihistory matrix for FFT processing
    end
    %fft on all window slides:
    sfBlock = abs(fft(IQSlideBlock.*specWindow,SDop.nfft)).^2;

    % flip direction of flow
    if isequal(SDop.revFlowDir, 1)
        sfBlock = flipud(sfBlock);
    end

    %normalizing spec display data:
    %calc noise floor
    noiseFloorHistory = [median(sfBlock(:,end),1),noiseFloorHistory];
    noiseFloorHistory(NumNoiseHist:end)=[];
    nfhMid = median(noiseFloorHistory);
    if isempty(noiseFloor)
        noiseFloor = nfhMid;
    else
        noiseFloor = noiseFloor*SDop.noisePersist + ...
            nfhMid*(1-SDop.noisePersist);
    end

    sfBlock = sfBlock.'/noiseFloor; %normalize to noise floor

    %spectral persistence parameters:
    aSpecPers = [1, -SDop.specPersist];
    bSpecPers = 1-SDop.specPersist;
    %spectral persistence, each bin over slow time:
    [sfBlock,stateSpecPersist] = filter(bSpecPers,aSpecPers, sfBlock,stateSpecPersist) ;

    % compress the spectral power
    switch lower(SDCompressionMethod)
        case {'db'}
            noiseFloorColorIndex = -20;
            sfBlock = 10*log10(sfBlock);
            dynRangeComp = SDDynRangeDB+SDProcGainDB;
        case 'power'
            noiseFloorColorIndex = 1;
            sfBlock = sfBlock .^ SDop.compression;
            dynRangeComp = (10^((SDDynRangeDB+SDProcGainDB)/10)).^SDop.compression;
        otherwise
            error('bad switch ')
    end
    %hole-filling processing:
    sfBlock = holefill(sfBlock,SDop.despeckle^(Nwind/SDop.nfft));

    %transform to colormap space:
    Sxx = sfBlock.'/dynRangeComp*(SDop.cMapLen - noiseFloorColorIndex)+noiseFloorColorIndex;

    %clip to colormap:
    Sxx = min(SDop.cMapLen,max(1,Sxx));

    %write baseline:
    Sxx(baseLineBin,:) = SDop.cMapLen/2;

    %Sxx(maxFbgn:maxFend) = SDop.cMapLen; %trace of max Frequency in flow
    if ~disableSpectralDisplay
        imUpdate= paintImage( ...
            'get.line.indices' , SdatAll, Sxx , k_chunk );
        SdatAll(:,imUpdate.columnIndicesRaw) = Sxx;
        if ~isempty(imUpdate.cursorData)
            SdatAll(:,imUpdate.columnIndicesCursor) = imUpdate.cursorData;
        end
        if ishandle(SDop.HSpectrogram) % Check that figure wasn't closed
            set(SDop.HSpectrogram,'CData', SdatAll); % Paint Spectrogram
        end

        tempT = clock;
        switch computerType
            case 'MACI'

                if invokeCount==drawnowUpdateFrameTrigger
                    pause(.0001) %update gui
                    invokeCount=0;
                end
            case 'MACI64'
                if invokeCount==drawnowUpdateFrameTrigger
                    pause(.0001) %update gui
                    invokeCount=0;
                end
            otherwise
                drawnow %update gui
        end
        drawnowTime = etime(clock,tempT);

    end

end %  spectral display code

invokeTime = toc;
runningTime = etime(clock,launchTime);

if nargout>0
    varargout{1} = invokeTime;
end

%diagnostics:
switch captureData
    case 1
        iqdatasavename = IQsave(IQdata(:));
    case 2
        iqdatasavename = IQsave(eval(captureVarsString),100);
    case 3
        iqdatasavename = IQsave(eval(captureVarsString),256);
    case 4
        iqdatasavename = IQsave(eval(captureVarsString),256);
    otherwise
        iqdatasavename='';
end
if ~isempty(iqdatasavename)
    disp([mfilename,':saved diagnostics in file: ',iqdatasavename ])
end



return

end %main
%% End of spectralDoppler.m

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% functions % functions % functions % functions %
% functions % functions % functions % functions %
% functions % functions % functions % functions %
% functions % functions % functions % functions %
% functions % functions % functions % functions %
% functions % functions % functions % functions %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function a = muter(FrameRate,reset)

persistent a0
if isequal(FrameRate,'cleanup')
    [a0]=deal([]);
    return
end

%muter: mute audio
if nargin<2
    reset =[];
end
if isempty(reset)
    reset = 0;
end

if isempty(a0)
    a0 = 1.0;
end

alpha = exp(-5/FrameRate);

if reset
    a0 = 1.0;
    a = 0.0;
    return
end

a0 = a0*alpha;

a = 1.0 - a0;

end %muter


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s=vars2VecEvalStr(varsList)
%used by data capture utility
s=strcat(varsList,';');
s=[ '[',cat(2,s{:}),']' ];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xf=holefill(X,alpha)
%spectral hole filler - set to zero for no holefill
% X must be a stack of row vectors
Xf = zeros(size(X));
Nr = size(X,1);
%note: envdet requires "forgetting factor"
for k=1:Nr
    [y1,y2]=envdet(1-alpha,X(k,:));
    Xf(k,:) = min([y1;y2]);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function iiCell = subframeInd100(sgp,prf)
%subframeInd100
%partition the data chunk into subrames, one for each line:
%splits up data indices in a frame interval
%R = interpolation factor
R=sgp.R;
TFrame = sgp.TFrame;


ii = round(linspace(1/R,1,R)'*prf*TFrame);
II = [[1;ii(1:end-1)+1],ii];
for k=1:R
    iiCell{k} = II(k,1):II(k,2);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout=paintImage(method,varargin)

persistent Zn2

if isequal(method,'cleanup')
    [Zn2]=deal([]);
    return
end

switch lower(method)

    case {'get.line.indices'}%method %%%%%%%%%%%%%%%%%
        [ SdatAll,Sdat,lineIndex]=deal(varargin{:});

        Nfft = size(Sdat,1);

        if size(Zn2,1)~=Nfft
            Zn2 =     ones(Nfft,1)*[255 0];
        end
        %convert frame to line0
        R = size(Sdat,2); %num lines per data frame
        NumDataFrames = size(SdatAll,2)/R;
        klineBase = (lineIndex-1)*R;
        klineIndStart = klineBase+ (1:R);

        %draw cursor lines:
        cursorIndices = [];
        cursorData = [];
        if lineIndex<NumDataFrames-2
            cursorIndices = (klineIndStart(end) + (1:2));
            cursorData = repmat(Sdat(:,R)*2,1,2);
        end

        %collect image update info:
        imUpdate.columnIndicesRaw = klineIndStart ;
        imUpdate.columnIndicesCursor = cursorIndices ;
        imUpdate.cursorData = cursorData;

        varargout{1} = imUpdate;

    case {'set.image.data.local'}%method %%%%%%%%%%%%%%%%%%%%%%%%%%
        %note: only for illustration.  Do this by inline code:
        entireSpectrogramData = varargin{1};
        newSpecData = varargin{2};
        imUpdate = varargin{3};

        SdatAll = entireSpectrogramData ;
        SdatAll(:,imUpdate.columnIndicesRaw) = newSpecData;
        SdatAll(:,imUpdate.columnIndicesCursor) = imUpdate.cursorData;

        varargout{1} = SdatAll;

    case {'set.image'}%method %%%%%%%%%%%%%%%%%%%%%%%%%%

        imageHandle = varargin{2};

        set(imageHandle,'CData',varargin{1});

    otherwise %method %%%%%%%%%%%%%%%%%%%%%
        error('bad switch')
end %method switch %%%%%%%%%%%%%%%%%%%%%%%%
end %paintimage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hFig, hAxs, HImg] = sdop_gui_lite(varargin)
%sdop_gui_lite Spectral Doppler GUI
%.
% [hFig, hAxs, HImg] = sdop_gui_lite(tSweep, nLine, nFreq, pwPrf,<baselineshift>)

m2cm = 100;MHz2Hz = 1e6;


yUnits = 'frequency';
[tSweep, nLine, nFreq, pwPrf, baselineshift, TXFreq_MHz,speedOfSound_mps]=deal(varargin{:});

% Specified parameters
pX = (600/nLine)*1.25; %pixels per x-unit
pY = 1.0; %pixels per y-unit
nX = nLine;
nY = nFreq;

% Y-axis is scaled from -nY/2+1:nY/2 and labels are in units of velocity
freq2vel = m2cm*speedOfSound_mps/(2*TXFreq_MHz*MHz2Hz); %frequency (Hz) to velocity (cm/s) scale factor
YLim = [-1/2 1/2]*nY;
nYTick = 9; % must be odd

% Colormapping
cMapLen = 128;

CMap = bone(cMapLen );

hAxs = [];

%
if isempty(baselineshift)
    baselineshift = 0.0;
end
hFig = figure;
set(hFig,'Toolbar','none','MenuBar','none', ...
    'Position',[100 20 pX*nX+100 pY*nY+100],...
    'Tag','sdop_figure' );
switch yUnits
    case 'frequency'
        yvec = (((0:nY-1)-(nY/2))/nY - baselineshift)* pwPrf;
        yLabelStr = 'Doppler Frequency (Hz)';
    case 'velocity'
        yvec = (((0:nY-1)-(nY/2))/nY - baselineshift)* pwPrf*freq2vel;
        yLabelStr = 'Doppler Velocity (cm/s)';
    otherwise
        error('bad switch')
end
tvec = linspace(0,tSweep,nX);
imageHandle = image(tvec,yvec,zeros(nY,nX));

ylabel(yLabelStr)
xlabel('time (s)')
set(gca,'YDir','normal');
colormap(CMap);

HImg = imageHandle;

drawnow('expose');

end%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sgp=spectrogramParameters(TFrame,sweepTimeInd)
% spectrogramParameters: helper fcn for pwrtas_vs and spectral doppler script
% sgp=spectrogramParameters(TFrame,sweepTimeInd);
% Example:
% sgp=spectrogramParameters(.02,1);
% sgp
%
% sgp =
%
%            TFrame: 0.0200
%        THistAllow: [2 3 4 6]
%          NumLines: 600
%             THist: 2
%     NumDataFrames: 100
%                 R: 6

%jaf 11/1/2009

%round to nearest microsecond:

sgp.TFrame = TFrame;

TFrame_uSec =  round(TFrame*1e6);


%sweeptime parameters are adapted here:
% Want parameters to satisfy the following restriction:
% NumLines/NumDataFrames = NumLines/(THist/TFrame) = integer only

switch TFrame_uSec

    case 10000
        sgp.THistAllow = [2 3  6];
        sgp.NumLines = 600;

    case { 20000 }
        sgp.THistAllow = [2 3 4 6];
        sgp.NumLines = 600;

    case { 30000 }
        sgp.THistAllow = [2  3   6];
        sgp.NumLines = 600;

    case { 40000 }
        sgp.THistAllow = [ 2 3 4 6 ];
        sgp.NumLines = 600;

    case { 50000 }
        sgp.THistAllow = [ 2 3 5 6 ];
        sgp.NumLines = 600;

    otherwise
        THAnom_uSec=(2:6)*1e6;
        ndfnom=THAnom_uSec./TFrame_uSec;
        ndf = round(ndfnom);
        LPSnom = 600;
        R =  round(LPSnom./ndfnom) ;
        LPS = R.*ndf;
        THA_uSec = TFrame_uSec.*ndf;
        THA = THA_uSec*1e-6;

        sgp.THistAllow = round(THA);
        sgp.NumLines = LPS(sweepTimeInd);


end

THistAllow = sgp.THistAllow;
numAllowedTHist = length(THistAllow);
if  sweepTimeInd > numAllowedTHist
    error('sweep speed out of range.'),
end
THist = THistAllow(sweepTimeInd);
sgp.THist = THist; %seconds per sgram
NumDataFramesPerSGram = fix(THist/sgp.TFrame); %frames/sgram
sgp.NumDataFrames =  NumDataFramesPerSGram;
sgp.R = round(sgp.NumLines/NumDataFramesPerSGram);  %lines/frame

if ~any(THist==sgp.THistAllow)
    error(['Time history length must be one of: ',num2str(THistAllow)])
end

end %func
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test
%test audio and spectrogram with siren IQ data

specDoppFunc = @spectralDoppler; %simplified for delivery

load('L7-4SpectralDoppler') %get parameters

clear mexaudiostream

WavesPerFrame = 6;
NsampPF = framePeriod*dopPRF;
iqchunk = getAudioTestDataChunk('siren.1',NsampPF,[],500); %prime
iqchunk=[1 1i]*reshape(iqchunk,2,NsampPF);%convert frome stereo to analytic signal
starttime = clock;
schedFrame=0;
for k=1:800 %frame loop


    %generate test IQ signal:
    iqchunk = getAudioTestDataChunk('siren.1');
    iqchunk=[1 1i]*reshape(iqchunk,2,NsampPF);%convert frome stereo to analytic signal
    IQtest= ones(4,1)*iqchunk;
    %add noise and scale:
    IQtest = 8000*(3*IQtest + 1*(randn(size(IQtest)) + 1i*randn(size(IQtest)))/sqrt(2) );

    extime = feval(specDoppFunc,IQtest);

    %wait for next sched. frame time before generating next test frame:
    while (k-1)>(schedFrame-15)
        runningTime = etime(clock,starttime);
        schedFrame = runningTime/framePeriod;
        switch computer
            case {'PCWIN','PCWIN64'}
                if rand(1)<.0001
                    fprintf('.')
                end
            case {'MACI','MACI64'}
                fprintf('.')
            otherwise
                fprintf('.')
        end
    end
    if rem(k,50)==0,disp(['test frame: ',num2str(k)]),end

end

clear mexaudiostream


disp('...test done.')
%%%%%%%%%%%%%%%%%%%%

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function x=getAudioTestDataChunk(testSignalType,varargin)
%getAudioTestDataChunk: returns stereo audio data chunk as row vector
%usage:
% initial call:
% x = getAudioTestDataChunk(testSignalType,Nload,<Nblocks>);
% x ~ [1 x Nload ]
% default Nblocks is  500
% subsequent calls:
% x = getAudioTestDataChunk(testSignalType);
% methods:
% testSignalType = 'siren.1' - continuous
%

% john flynn 11dec2009

persistent PhzVecs kload Nload NFrameCyc vibratoAmount
if isequal(testSignalType,'cleanup')
    [PhzVecs ,kload ,Nload ,NFrameCyc ,vibratoAmount]=deal([]);
    return
end


%parse inputs:
Kvai = length(varargin);kvai = 1; %use kvai and template below:
if kvai<=Kvai, tempVar=varargin{kvai} ; %%%%%%%%%%%%%%
else, tempVar=[];
end;kvai=kvai+1; NloadIn = tempVar;
if kvai<=Kvai, tempVar=varargin{kvai} ; %%%%%%%%%%%%%%
else, tempVar=[];
end;kvai=kvai+1; NFrameCycIn = tempVar;
if kvai<=Kvai, tempVar=varargin{kvai} ; %%%%%%%%%%%%%%
else, tempVar=[];
end;kvai=kvai+1; vibratoAmountIn = tempVar;

if isempty(PhzVecs)
    disp([mfilename,': initializing audio test signal gen...'])

    if isempty(NFrameCycIn)
        NFrameCycIn = 500;
    end
    if isempty(vibratoAmountIn)
        vibratoAmountIn = 2;
    end
    if isempty(NloadIn)
        error('must specify load size')
    else
        Nload = NloadIn; % samples per DAC frame
    end
    NFrameCyc = NFrameCycIn;
    vibratoAmount = vibratoAmountIn;

    %siren phase history:
    PhzVecs = 2*pi*3*(0:NFrameCyc*Nload-1)/(Nload)+ ...
        (sin(2*pi*3*(0:NFrameCyc*Nload-1)/(NFrameCyc*Nload))*vibratoAmount)*2*pi;
    PhzVecs = reshape(PhzVecs,Nload,NFrameCyc)';

    kload = 0;

end

Nblock = size(PhzVecs,1);
switch lower(testSignalType)
    case {'siren.1'}%testSignalType %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        kload = mod(kload+1-1,Nblock)+1;
        phzvec=PhzVecs(kload,:);
        x = [cos(phzvec);sin(phzvec)];
        x=x(:).';
    otherwise
        error('bad test signal type name')
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y,stateOut]=doppler_channelize(x,state)
% doppler_channelize: split Doppler IQ into L/R by spectral halves.
% Output is coded into real/imag as left/right
% usage:
%  [y,stateOut]=doppler_channelize(x);
%  [y,stateOut]=doppler_channelize(x,state);
%
% Output    Spectrum
% - - - - - - - - -
% real      positive
% imag      negative
%
% ----------------------------------------------------
% John Flynn 16Feb2009
% (c) 2009 VeraSonics, Inc.

%initialize
ChannelizerBWFactor = 1.0; % range = (0.0 -  1.0)

if ischar(x)
    %special command processing
    command = x; clear x
    switch lower(command)
        case {'test'} %comm switch %%%%%%%%%%%%%%%%%%%%%%%

            x = dopplersignal(30000);
            X = reshape(x ,3,10000);

            X = conj(X); %switches audio to opposite channel

            [Y(1,:),dcs] = doppler_channelize(X(1,:));
            [Y(2,:),dcs] = doppler_channelize(X(2,:),dcs);
            [Y(3,:),dcs] = doppler_channelize(X(3,:),dcs);

            y = Y.';
            y = Y(:);

            y2 = upsamplevs(y,5);

            ylr = [real(y2(:)),imag(y2(:))];

            %should have sound in left or right channel only
            soundsc(ylr,44100);

        otherwise
            error('bad switch ')
    end %comm switch %%%%%%%%%%%%%%%%%%%%%%%%%

    return;

end%command check

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normal processing


if nargin<2
    state = [];
end

if isempty(state)
    % init state:
    state.fshifters = [] ; %note: all stages are same since all sizes are same
    state.lpfLeft = [];
    state.lpfRight = [];

    %    filter design:
    [fproto]= getFilterPrototype('doppler_channelize');
    [b,a]=filterTranslate(...
        fproto.tf.b,fproto.tf.a,fproto.nuNominal,  0.25*ChannelizerBWFactor ,'low');

    %%%%[b,a]=ellip(LPFOrd,LPFPassBAtten,LPFStopBAtten,2*(0.25*ChannelizerBWFactor));
    state.b = b;
    state.a = a;

end
b = state.b;
a = state.a;

stateOut = state;

[xhi,fshifters] = freqshift(x,-0.25,state.fshifters);
[xlo] = freqshift(x,0.25,state.fshifters);
stateOut.fshifters = fshifters;

[xhif , stateOut.lpfRight]= filter(b,a,xhi,state.lpfRight);
[xlof , stateOut.lpfLeft] = filter(b,a,xlo,state.lpfLeft);
[yr] = freqshift(xhif,0.25,state.fshifters);
[yi] = freqshift(xlof,-0.25,state.fshifters);

y = real(yr)+1i*imag(yi);

end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xc = limiter(method,x,varargin)
%limiter: saturation appropriate for audio compression/compressor
% compression/limiting is w.r.t output range of +/- 1.0
%usage:
% xc = limiter('soft',x,varargin);
% xc = limiter('hard',x,varargin);

%jaf 14apr2009

if isempty(method)
    method = 'hard';
end

Kvai = length(varargin);kvai = 1; %use kvai and template below:

switch lower(method)
    case {'hard'} %method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        xc = x;
        [posLim,negLim]=deal(1,-1);

        x(find(x>posLim)) = posLim;
        x(find(x<negLim)) = negLim;

    case {'soft'} %method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if kvai<=Kvai
            tempVar=varargin{kvai} ;
        else
            tempVar=[];
        end
        kvai=kvai+1; compressionFactor = tempVar;

        if isempty(compressionFactor)
            compressionFactor = 1.0;
        end

        xc = atand( x * compressionFactor) /(90 );

    otherwise %method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        error('bad switch ')

end %method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [y,stateOut,s] = freqshift(x,nu,state)
%freqshift: shift freq of input vector by nu.
%usage:
% [y,stateOut,<shifterVec>] = freqshift(x,nu,state);
% x is a vector of signal samples.
% nu is normalized frequency.
%

%jaf 16feb2009


[m,n]=size(x);
isRow = m==1;

if min(m,n)>1
    error('matrix input not supported')
end
if nargin<3
    state = [];
end

if isempty(state)
    state.nextStartPhasor = 1;
end

L = length(x);


if isRow
    s = exp((2*pi*1i*nu)*[0:L])*state.nextStartPhasor;
else
    s = exp((2*pi*1i*nu)*[0:L]')*state.nextStartPhasor;
end

y = x.*s(1:end-1);

stateOut.nextStartPhasor = s(end);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,stateOut]=phasor(N,nu,state)
%phasor: complex phasor with frequency nu, N samples
%usage: [x,state]=phasor(N,nu,state);
%N - number of samples
%nu normalized freq.
%state- init for recurrence

%john flynn 12/7/2009
if nargin<3
    state = [];
end

if ischar(N)
    comm = N;
    switch lower(comm)
        case 'test'
            [x1, s]=phasor(12, 0.134);
            [x2, ~]=phasor(6, 0.134, s);
            [x12]=phasor(18, 0.134);

            dx=x12-[x1,x2];
            err=norm(dx)/norm(x2);

            x = err;
            disp([mfilename,':',comm,': test errors: ',num2str(err)])

        otherwise
            error('bad switch')
    end

    return

end

if isempty(state)
    state = 0.0;
end
%   2*pi*i*nu*[0:N] + 2*pi*i*phi
% = 2*pi*i*(nu*[] + phi)
% = 2*pi*i*(nu*([]+phi/nu))
phz = (2*pi*nu*1i)*((0+state):(N-1+state)) ;
stateOut = N+state;
x=exp(phz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
function rsParams = resamplingParams(PRF,FsDAC)
% resamplingParams: resampling parameters
%.
% Given PRF and FsDAC, returns a structure with:
% FsIntermediat = Intermediate resample rate;
% Qupsample = upsampling factor
%

%jaf 17feb2009

if FsDAC<PRF
    error('PRF is higher than specified DAC sample rate')
end

MaxQ = 48;

Qset = 2:MaxQ;
FsDACsub = FsDAC./Qset;
ff = find(FsDACsub>PRF) ;
Index_Fs_leastGT_PRF = ff(end);

rsParams.Qupsample= Qset(Index_Fs_leastGT_PRF);
%upsampler rate from smallest subsample GT PRF

rsParams.FsIntermediate = FsDACsub(Index_Fs_leastGT_PRF);
rsParams.FsDAC = FsDAC;
rsParams.PRF = PRF; %input Fs

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [varargout]=circbuff(method,varargin)
%circbuff: manage circular buffer.
% read does not advance pointers.
%signatures -
% [state]=circbuff('initialize',Nbuff);
% [state]=circbuff('write',state,dataIn);
% [dataOut]=circbuff('read',state,Nread);

%john flynn 11/5/2009
Kvai = length(varargin);kvai = 1; %use kvai and template below:

switch lower(method)
    case {'initialize'} %--------------------------
        if kvai<=Kvai, tempVar=varargin{kvai} ; %%%%%%%%%%%%%%
        else tempVar=[];
        end;kvai=kvai+1; NBuff = tempVar;

        state.NBuff = NBuff;
        state.isInitialized = 1;
        state.nextWriteIndex = 1;
        state.data = zeros(1,state.NBuff);

        varargout{1}=state;

    case {'write'} %--------------------------
        if kvai<=Kvai, tempVar=varargin{kvai} ; %%%%%%%%%%%%%%
        else tempVar=[];
        end;kvai=kvai+1; state = tempVar;
        if kvai<=Kvai, tempVar=varargin{kvai} ; %%%%%%%%%%%%%%
        else tempVar=[];
        end;kvai=kvai+1; dataIn = tempVar;

        Nnew = length(dataIn);
        wrInd = [0:Nnew-1]+state.nextWriteIndex;
        wrInd = restrictIndex(state,wrInd);

        state.data(wrInd) = dataIn;
        state.nextWriteIndex = restrictIndex(state,state.nextWriteIndex+Nnew);
        state.isInitialized = 0;

        varargout{1}=state;

    case {'read'} %--------------------------
        if kvai<=Kvai, tempVar=varargin{kvai} ; %%%%%%%%%%%%%%
        else tempVar=[];
        end;kvai=kvai+1; state = tempVar;
        if kvai<=Kvai, tempVar=varargin{kvai} ; %%%%%%%%%%%%%%
        else tempVar=[];
        end; kvai=kvai+1; Nread = tempVar;

        if state.isInitialized,error('attempt to read from empty circular buffer '),end

        K=state.nextWriteIndex;
        ri = K-1-Nread+1:K-1;
        ri = restrictIndex(state,ri);

        varargout{1} = state.data(ri);


    otherwise %--------------------------
        error( 'bad switch' )
end %--------------------------


    function ir=restrictIndex(state,indVec)
        ir = mod(indVec-1,state.NBuff)+1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dataFileName]=IQsave(x,NframesIn,saveVarNames)
%IQsave: capture IQ data into buffer and save into file.
%.
%initialize:
% IQsave(dataColumnSnapshop,NumHist);
%.
% accumulate and auto save if history is full:
% saveFileName = IQsave(dataColumnSnapshop);
% saveFileName is empty unless save occurred.

%jaf 05apr2009
%

verbose = 1;

persistent X nframes Nframes captureVars

if nargin<3
    saveVarNames=[];
end

if isequal(x,'cleanup')
    [X ,nframes, Nframes]=deal([]);
    return
end

if isempty(captureVars)
    if isempty(saveVarNames)
        [MX,~]=size(X);
        saveVarNames = cell(1,MX);
        saveVarNames(:) = {''};
    end
    captureVars = saveVarNames;
end

doInit = length(x)~=size(X,1);
if doInit
    %initialization
    if nargin<2
        NframesIn = [];
    end
    if isempty(NframesIn)
        NframesIn = 300;
    end
    Nframes = NframesIn;
    X = single(zeros(length(x),Nframes));
    nframes = 0;

else

nframes = nframes+1;
nframes = rem(nframes-1,Nframes)+1;

X(:,nframes) = single(x(:));

end

saveBuffer = nframes==Nframes ;
if saveBuffer
    fmt=30; %datetime format
    timeAtSave = now;
    dataFileName = ['iqsave_',datestr(now,fmt)];
    creator = mfilename;
    callingStack = dbstack;
    savelist = {'X' , 'Nframes' , 'creator', 'dataFileName' , 'callingStack' ,'captureVars', 'timeAtSave'};
    save(dataFileName, savelist{:});
    if verbose
        disp(' ')
        disp([mfilename,':captured data in file: ', dataFileName,'.mat'])
    end
else
    dataFileName = '';
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SDop = initSDop
%defaults
SDop = struct('PWdopPRF', [], ...
    'TFrame' , [] , ...
    'nSampleVolume', 1, ...
    'nPulses', [], ...
    'nfft', 256 , ... %FFT size
    'nIqfft', 35 , ... %FFT time window
    'revFlowDir', 0 , ... %flips axis
    'specPersist', 0.3 , ... %spectral slowtime smoother, 0.0 to 1.0
    'despeckle',0.45 , ...  %hole filler, 0.0 to 1.0
    'noisePersist', 0.95 , ... %noise floor smoother, 0.0 to 1.0
    'displayOn', 1, ...
    'cMapLen', 128, ...
    'audioOn', 1 , ... %set to 0 to disable audio
    'audioGain', 3e-6, ...
    'audCompressionFactor' , 1.5, ... %higher is more compression (low level gain).
    'compression', 0.5, ... %spectral display compression level for power method,
    ... % if set to 1.0 ==> makes the colormap index units be watts into 1 Ohm load
    'SDCompressionMethod' , 'power' , ... %set to 'power' or 'db'
    'isInit', 0);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fproto]= getFilterPrototype(filterPurpose)
%getFilterPrototype: retrieve digital filter prototype for specific
% spectral doppler task.

persistent filterDesigns
if isempty(filterDesigns)
    load spectraldoppler_filtercoeff
end

fproto =filterDesigns.(filterPurpose);

end%
%%%%%%%%%%%%%%%%%%%%%%%%
function filterDesigns = filterDesignSpecs
%this is called by specDoppDesigner, which can be isolated from remainder
%of spectral doppler package for a compilation not dependent on the Signal Processing Toolbox

%call filter design functions with nominal critical freq,
% and desired attenuation params, in lowpass form.
% Will transform as needed to desired critical frequency at
% run-time.
filterDesigns = struct;
nuNominal = 0.25;
style = 'low'; %will be transformed to desired style in runtime.
% - - - -
%custom part:
filterPurpose = 'audio_resampler'; %at top level
callingSignature= { 5,1.5,60, nuNominal*2 ,style};
filterHandle = 'ellip';
%common:
fds.nuNominal = nuNominal;
fds.filterPurpose = filterPurpose;
fds.filterHandle = filterHandle;
fds.callingSignature = callingSignature;
filterDesigns.(filterPurpose) = fds;

% - - - -
%custom part:
filterPurpose = 'doppler_channelize'; %in a subroutine
callingSignature={5,1.5,30, nuNominal*2 ,style};
filterHandle = 'ellip';
%common:
fds.filterPurpose = filterPurpose;
fds.filterHandle = filterHandle;
fds.callingSignature = callingSignature;
filterDesigns.(filterPurpose) = fds;

% - - - -
%custom part:
filterPurpose = 'wallfilter';
callingSignature= { 7, 65  nuNominal*2 ,style };
filterHandle = 'cheby2';
%common:
fds.filterPurpose = filterPurpose;
fds.filterHandle = filterHandle;
fds.callingSignature = callingSignature;
filterDesigns.(filterPurpose) = fds;

varargout{1}= filterDesigns;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = generateTestData(method, ppf)
t = (0:(ppf-1))/ppf  ; %time vector to wrap every frame
iphz = 2*pi*1i*t ;
audioScale = 40000;
clutterLevel = 330;
switch method
    case 'tones.1'
        flow = exp(iphz*20) +  exp(-iphz*18)  ;
        y=audioScale*(  flow );
    case 'tones.clutter.slow'
        flow = exp(iphz*8) +  exp(-iphz*9)  ;
        clutter = clutterLevel*exp(iphz);
        y=audioScale*( clutter +  flow  );
    case 'tones.clutter.dc'
        flow = exp(iphz*20) +  exp(-iphz*18)  ;
        clutter = clutterLevel*1;
        y=audioScale*( clutterLevel*clutter +  flow  );
    otherwise
        error('bad switch')
end
end % generateTestData
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [varargout] = restrictAcqPRF( varargin )

%parse inputs:
Kvai = length(varargin);kvai = 1; %use kvai and template below:

method = 'prf-constraints.1';

switch method
    case 'legacy'
        %inputs: prfnom,tframenom
        [prfnom,tframenom] =deal(varargin{[2 1]} );
        % parameter in seconds or hz
        %produces even number of acq pulses per frame
        ttnanom = 1/prfnom ;
        ttnanom_microsec = ttnanom*1e6 ;
        ttna_microsec = round(ttnanom_microsec);
        ttna = ttna_microsec/1e6;
        prf = 1/ttna;
        if nargin>1
            ppffnom = prf*tframenom;
            ppf=round(ppffnom);
            %round to higher even
            ppf = ppf + double(rem(ppf,2)==1);
            tframe = ppf/prf;
        end
        varargout{:} = deal( prf,ttna_microsec,ttnanom_microsec,tframe,ppf );

    case 'prf-constraints.1' %-------------------
%signature:  FTNom_uSec, PRFDopNom_Hz
% where:
%  class(FTNom_uSec) == 'char'  ==> use exactly
%  class(FTNom_uSec) == numeric ==> will quantize

        if kvai<=Kvai
            tempVar=varargin{kvai} ; %%%%%%%%%%%%%%
        else
            tempVar=[];
        end
        kvai=kvai+1; FTNom_uSec = tempVar;
        if kvai<=Kvai
            tempVar=varargin{kvai} ; %%%%%%%%%%%%%%
        else
            tempVar=[];
        end
        kvai=kvai+1; PRFDopNom = tempVar;
        if kvai<=Kvai
            tempVar=varargin{kvai} ; %%%%%%%%%%%%%%
        else
            tempVar=[];
        end
        kvai=kvai+1; verbose = tempVar;

        nc=(190:600);  %candidate clocks/PRI

        if isempty(FTNom_uSec)
%default
            FTNom_uSec =  '2*2*2*2*3*3*5*7*11' ; % 55440
        end
        if isempty(verbose)
            verbose =  0 ;
        end

        %choose closest allowed Frametime
        %frame time in microseconds
        if ischar(FTNom_uSec)
            %use exact
            FT_uSec = eval(FTNom_uSec);
        else
            %will quantize to "best" nearest
            %best allowing choice of prfs without changing frame time
            FTvec_uSec = [ ...
                eval('2*2*2*2*2*2*2*3*3*5*7') ... 40320
                eval('2*2*2*3*3*5*11*13') ... 51480
                eval('2*2*2*2*3*3*5*7*11') ... 55440
                eval('2*2*2*2*3*3*5*7*12') ... 60480
                ];
            [~,indft]= min(abs(FTNom_uSec-FTvec_uSec));
            FT_uSec = FTvec_uSec(indft);
        end

        pphalffnom=(FT_uSec/2)./nc;
        pphalffc=round(pphalffnom);
        iii=find(abs(pphalffnom-pphalffc)./abs(pphalffc)<1e-6);
        pphalff = pphalffc(iii);
        pria=nc(iii)*1e-6;

        prfa=1.0./pria;
        prfd = (prfa/2);
        prfdnom=round(prfd);

        pulseSched.PRFAcq = prfa;
        pulseSched.PRFDop = prfd;
        pulseSched.maxPRFAcq = max(prfa);
        pulseSched.maxPRFDop = max(prfd);
        pulseSched.minPRFAcq = min(prfa);
        pulseSched.minPRFDop = min(prfd);
        pulseSched.pulsesPerFrameAcq  = round(pphalff*2);
        pulseSched.pulsesPerFrameAcqMax  = max(pulseSched.pulsesPerFrameAcq);
        pulseSched.pulsesPerFrameDop  = round(pphalff);
        pulseSched.pulsesPerFrameDopMax  = max(pulseSched.pulsesPerFrameDop);
        pulseSched.PRIAcq_uSec  = nc(iii);
        pulseSched.PRIDop_uSec  = nc(iii)*2;
        pulseSched.FrameTime_uSec  = round(FT_uSec);

        if ~isempty(PRFDopNom)
            [~,indPRF] = min(abs(pulseSched.PRFDop-PRFDopNom));
        else
            indPRF = [];
        end
        pulseSched.indexPRF = indPRF;

        if ~isempty( pulseSched.indexPRF)
             fnames=fieldnames(pulseSched);
             for k=1:length(fnames)
                 varname = fnames{k};
                 var = pulseSched.(varname);
                 L=length(var);
                 if L>1
                     pulseSched.(varname)=var(pulseSched.indexPRF);
                 end
             end
        end

        if verbose
            FT_uSec;
            plot(prfdnom,'o'),
            grid on
        end
        varargout{1}= pulseSched;

    otherwise
        error('bad switch')
end %method switch
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function unInitPersistentInCaller
% assign empty to all persistent vars in caller wkspc
w=evalin('caller','whos;');

wp={w.persistent};
wp=cat(2,wp{:});
wn = {w(wp==1).name};

for k=1:length(wn)
    evalin('caller',[wn{k},'= [];']);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cleanup
%close audio driver:
audiostream('cleanup')
%need to init persist in caller (not in this function)
evalin('caller', 'unInitPersistentInCaller;')

%uninit. subfunctions with persistent vars:
muter('cleanup')
paintImage('cleanup')
getAudioTestDataChunk('cleanup')
IQsave('cleanup')

end
% [EOF] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
