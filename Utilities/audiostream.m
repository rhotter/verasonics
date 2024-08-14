function varargout=audiostream(method,varargin)
%audiostream: wrapper/caller for audiostream PortAudio c mexfile
% Copyright 2001-2017 Verasonics, Inc.  All world-wide rights and remedies under all intellectual property laws and industrial property laws are reserved.  Verasonics Registered U.S. Patent and Trademark Office.
% Type ">> audiostream " at the command line to get help.

%John Flynn 08feb2009
% (c) 2009 VeraSonics, Inc.
% - - -
% revisions:
% 31mar2009 JAF - added more info/parameters to "fifostate" method.

if nargin<1
    disp([mfilename, ' - Help: '])
    disp('VeraSonics, Inc, (c) 2009.')
    disp('PortAudio driver interface to MatLab.')
    disp('Components:  audiostream.m and mexaudiostream.c (a mex file callable from MatLab) ')
    disp('>> audiostream.m will load a FIFO with audio chunks, and play on speakers as shown below.')
    disp('>> audiostream.m is a wrapper for mexaudiostream.c (also callable from command line).')
    disp('>> See mexaudiostream.c header for its methods documentation.')
    disp(' ')
    disp('Usage: no argument: display this Help Message.')
    disp('Usage: at least one argument:')
    disp('>> audiostream(''open''); % open/create the PortAudio stream object')
    disp('>> audiostream(''load'', audioChunk); % load a chunk, up to FIFO capacity, ')
    disp('    where audioChunk is interleaved stereo vector (class DOUBLE).')
    disp('    If audioChunk is Nx2 or 2xN, assumes stereo vectors (non-interleaved).')
    disp('>> audiostream(''play''); % start/resume play from FIFO data, and return control to matlab.')
    disp('>> audiostream(''pause'') % stop playing.')
    disp('>> audiostream(''close'') % delete the PortAudio stream object')
    disp('>> audiostream(''fifostate'') % retrieve state of audiostream FIFO in a structure. ')
    disp('>> % test functions: ')
    disp('>> audiostream(''test.fifo'') % test fifo functions ')
    disp('>> audiostream(''test.audio'') % test audio streaming with pitched siren.')
    disp(' ')
    disp('Build with: >> mex mexaudiostream ')
    disp(' ')
    disp(' Requires PortAudio library built in location: ')
    disp('     ~/Users/johnf0001/Documents/Programming/Verasonics/Matlab Simulator/portaudio/')
    disp(' To build PortAudio library, type in shell: ')
    disp(' > ./configure && make')
    disp(' - - - - - - - - - - - - - - - - - ')
    disp('Test history: ')
    disp('Built (Verasonics 1.5.1) and Verified on MacBook Pro, 2/13/2009' )
    disp('Migrated mexaudiostream.c to Verasonics environment 1.5.2+ , 3/26/2009')
    disp('*****************************************')
    return
end

persistent x %local audio chunk vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Kvai = length(varargin);kvai = 1; %use kvai and template below:

switch lower(method)
    case {'cleanup'}
        clear mexaudiostream

    case {'open'} %method %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ispc
          %DISABLED: may be bug (5/4/2013)  mexaudiostream(100); %directsound driver stream
            mexaudiostream(0);
            disp('Audio: open default stream ....') %directsound driver stream
        else
            mexaudiostream(0);
        end

    case {'close'} %method %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mexaudiostream(-1);

    case {'start','play'} %method %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if kvai<=Kvai, tempVar=varargin{kvai} ; %%%%%%%%%%%%%%
        else tempVar=[];
        end;kvai=kvai+1; PctOccThresh = tempVar; %in range [0.0 1.0]

        if isempty(PctOccThresh)
            mexaudiostream(1);
        else
            mexaudiostream(1,PctOccThresh);
        end

    case {'stop','pause'} %method %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mexaudiostream(2);

    case {'fifostate'} %method %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fifoState.occupancy = mexFifo('occupancy');
        fifoState.vacancy = mexFifo('vacancy');
        fifoState.capacity = mexFifo('capacity');
        fifoState.framecount = mexFifo('framecount');
        fifoState.isRunning = mexFifo('running');
        fifoState.FifoWasEmpty = mexFifo('wasempty');
        fifoState.FifoReadThreshFlt = mexFifo('readthresh');
        fifoState.minStartOccupancyFlt = mexFifo('minstartocc');

        varargout{1} = fifoState;

    case {'load'} %method %%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if ~isempty(varargin)
            %ensure row vector:
            [m,n]=size(varargin{1});
            if m==1
                %interleaved, is already row:
                x= varargin{1};
            elseif n==1
                %interleaved, is column:
                x = varargin{1}.';
            elseif m ==2
                %is stereo rows:
                x = varargin{1};
                x = x(:).'; %interleave
            elseif n ==2
                %is stereo columns:
                x = varargin{1}.';
                x = x(:).'; %interleave
            end

        else
            error('requires row vector input for audio chunk')
        end

        mexaudiostream(1000,x);

        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case { ...
            'test', 'test.audio' ...
            } %method %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %these commands in separate file:

        mexaudiostream(-1)
        clear mexaudiostream
        enableThrottling = 0;
        mexaudiostream(3)
        disp('opening stream: ')
        audiostream( 'open' )
        mexaudiostream(3)
        audioVolume = .25;
        started =0;
        floatsPerSample = 2;
        loopRecycleLimit = 80000;
        timingAlertThresh = .004;
        Nload = 256; % samples per DAC frame
        NFrameCyc = 500;
        statsbuff2(1:loopRecycleLimit) = struct('OccFlt' , -1, 'RelTime', -1 );
        statsbuff(1:NFrameCyc) = struct('timeSinceLastLoad',-1,'k',-1,'OccFlt',-1);
        y=getAudioTestDataChunk('siren.1',Nload,NFrameCyc,10);%init the signal generator

        NoccThreshFltsToLoad=8*Nload*floatsPerSample;
        NoccThreshFltsToStart = NoccThreshFltsToLoad;

        k=0;
        Twait = .002;
        fs=audiostream('fifostate');
        fs
        BuffSizeFlts =fs.capacity;
        NoccFracThresh=NoccThreshFltsToStart/BuffSizeFlts  %trigger to start
        kload = 0;
        pause(.5)
        reportedUnderflow = 0;
        stopCheck = 0;
        tic
        try
            while 1
                k=mod(k+1,loopRecycleLimit);

                NoccFlts = mexaudiostream(903);
                wasEmpty =  mexaudiostream(909);
                statsbuff2(k+1) = struct('OccFlt',NoccFlts,'RelTime',toc);
                if wasEmpty==0
                    %build one-shot trigger to report underflow
                    reportedUnderflow = 0;
                end
                if  started && (wasEmpty~=0) && (reportedUnderflow==0)
                    %flag and analyze an underflow:
                    timeSinceLastLoad = statsbuff(kload).timeSinceLastLoad;
                    kloadDisp = mod([-6:0]+kload,NFrameCyc)+1;
                    fprintf('Was empty, time since last load, load#: %f %i\n', timeSinceLastLoad , kload )
                    disp('[deltaT, kloop, Occ Flts] at load:')
                    diags = ...
                        [cat(2,statsbuff(kloadDisp).timeSinceLastLoad); ...
                        cat(2,statsbuff(kloadDisp).k); ...
                        cat(2,statsbuff(kloadDisp).OccFlt); ...
                        ];
                    diags(1,:)
                    diags(2:3,:)
                    disp('[Occ Flts, diff(Rel Time)] in loop check:')
                    kDisp = mod([-5:0]+k,loopRecycleLimit)+1;
                    diags2 =  [cat(2,statsbuff2(kDisp).OccFlt) ; cat(2,statsbuff2(kDisp).RelTime) ];
                    diags2(1,:)
                    diff(diags2(2,:))
                    %check for excessive time since last load:
                    if any(abs(diff(diags2(2,:)))>.02 & abs(diff(diags2(2,:)))<.9)
                        stopCheck=stopCheck+1;
                        disp('********************')
                        if stopCheck>0
                            switch 'flag'
                                case 'flag'
                                    disp('started,wasEmpty,reportedUnderflow,stopCheck')
                                    disp(num2str([started,wasEmpty,reportedUnderflow,stopCheck]))
                                case 'stop'
                                    audiostream('stop')
                                    disp([mfilename,'breaking ...'])
                                    keyboard
                                    audiostream( 'start' ,NoccFracThresh); %start the player
                                otherwise
                                    error('bad switch')
                            end
                        end
                    end
                    reportedUnderflow = 1; %notify that we've reported it
                end

                if started && enableThrottling
                    %if mod(k,200)==0, fprintf(' %i',k), end
                    %throttle the loading with a busywait
                    while NoccFlts>NoccThreshFltsToLoad
                        NoccFlts = mexaudiostream(903);
                        if ispc
                            zjunk=fft(rand(512,60));
                        else
                            pause(Twait);
                        end
                    end
                else
                end

                if NoccFlts< NoccThreshFltsToLoad
                    %load to the fifo:
                    kload = mod(kload+1,NFrameCyc);
                    if NoccFlts<(Nload*floatsPerSample)

                        fprintf('fifo is Low;\n')

                    end
                    if mod(kload,500)==0
                        fs=audiostream('fifostate');
                        fprintf('Loading Fifo, Occ now: %i floats\n',fs.occupancy),
                        fprintf('FifoWasEmpty: %i \n',fs.FifoWasEmpty),
                    end
                    y = audioVolume*getAudioTestDataChunk('siren.1');

                    %% Load - - - - - - - - - -
                    timeSinceLastLoad = toc;
                    statsbuff(kload+1) = ...
                        struct('timeSinceLastLoad',timeSinceLastLoad,'k',k,'OccFlt',NoccFlts);
                    audiostream( 'load' ,  y);
                    tic;
                end
                if k==0
                    disp(datestr(now))
                    %                     fs=audiostream('fifostate');
                    %                     fs
                end

                if ~started %& k>4,
                    audiostream( 'start' ,NoccFracThresh); %start the player
                    started = 1;
                end
                %                 if k<5,
                %                                         fs=audiostream('fifostate');
                % framecountvec(end+1)=fs.framecount;
                %                 end
            end

        catch ME

            disp('Error, closing audio stream: ')
            audiostream( 'stop' )
            audiostream( 'close' )

            throw(ME)

        end

    otherwise %method %%%%%%%%%%%%%%%%%%%%%%%%%%%%

        error('unknown method ')

end %method %%%%%%%%%%%%%%%%%%%%%%%%%%%%
end % main.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%private:
%functions % %%%%%%%%%%%
%functions % %%%%%%%%%%%
%functions % %%%%%%%%%%%
%functions % %%%%%%%%%%%
function nrmse = NRMSE(x,y)
nrmse = norm(x(:) - y(:))/norm(  x(:) );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = mexFifo ( comm   ,   pin  )
%mexFifo: gateway to   mexaudiostream.c  fifo  functions

switch lower(comm)
    case 'write'
        mexaudiostream( 902,    pin  );
    case 'read'
        N=pin;
        varargout{1}= mexaudiostream( 901  , N   );

    case 'clear'

        clear mexaudiostream
        status = [];

    case {'vacancy'} %%%%%%%
        status = mexaudiostream(904);

    case {'occupancy'} %%%%%%%
        status = mexaudiostream(903);

    case {'capacity'} %%%%%%%
        status = mexaudiostream(905);

    case {'framecount'} %%%%%%%
        status = mexaudiostream(906);
    case {'running'} %%%%%%%

        status = mexaudiostream(907);

    case {'samplerate'} %%%%%%%
        status = mexaudiostream(908);
    case {'wasempty'} %%%%%%%
        status = mexaudiostream(909);
    case {'readthresh'} %%%%%%%
        status = mexaudiostream(910);
    case {'minstartocc'} %%%%%%%
        status = mexaudiostream(911);


    otherwise
        error('bad switch ')
end
varargout{1} = status;
end
%%%%%%%%%%%%%%%%%%%%%%
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

%parse inputs:
Kvai = length(varargin);kvai = 1; %use kvai and template below:
if kvai<=Kvai, tempVar=varargin{kvai} ; %%%%%%%%%%%%%%
else tempVar=[];
end;kvai=kvai+1; NloadIn = tempVar;
if kvai<=Kvai, tempVar=varargin{kvai} ; %%%%%%%%%%%%%%
else tempVar=[];
end;kvai=kvai+1; NFrameCycIn = tempVar;
if kvai<=Kvai, tempVar=varargin{kvai} ; %%%%%%%%%%%%%%
else tempVar=[];
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
    PhzVecs = 2*pi*3*[0:NFrameCyc*Nload-1]/(Nload)+ ...
        (sin(2*pi*3*[0:NFrameCyc*Nload-1]/(NFrameCyc*Nload))*vibratoAmount)*2*pi;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


