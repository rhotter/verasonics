function [varargout] = rs_irrat(command,varargin)
%rs_irrat: irrational sample rate conversion
% Copyright 2001-2017 Verasonics, Inc.  All world-wide rights and remedies under all intellectual property laws and industrial property laws are reserved.  Verasonics Registered U.S. Patent and Trademark Office.
% usage:
% [yresampled,stateOut ] = rs_irrat(command, p1, p2 ... );
%
% specific signatures:
%
% [rsState ] = rs_irrat('initialize', FsIn, FsOut , nSamplesIn, anaResampFiltName );
%    default val of anaResampFiltName is 'rc.23'
% [y , stateOut ] = rs_irrat('', x, state );
% [y , stateOut ] = rs_irrat('resample', x, state );
%
% x = 1 x N input data, sampled at FsIn rate.
% y = resampled output data, at state.FsOut rate.
% state = state and parameter information.
%
% - - - - - - - - - -
%
% state = rs_irrat('init', FsIn, FsOut, paramSet);

%John Flynn 2/15/2009
% (c) 2009 Verasonics, Inc.

%
Kvai = length(varargin);kvai = 1; %use kvai and template below:

%debut/test setup: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%testCode = 'timevec';
testCode = '';
if ~isempty(testCode)
    disp([mfilename,':testCode:',testCode])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch command
    case {'init','initialize','init.state'} %command ----------------
        %initialize state - - - - - -  - - - - - - - - -

        %get inputs:
        if kvai<=Kvai, tempVar=varargin{kvai} ; %%%%%%%%%%%%%%
        else tempVar=[];
        end;kvai=kvai+1; FsIn = tempVar;
        %.
        if kvai<=Kvai, tempVar=varargin{kvai} ; %%%%%%%%%%%%%%
        else tempVar=[];
        end;kvai=kvai+1; FsOut = tempVar;
        %.
        if kvai<=Kvai, tempVar=varargin{kvai} ; %%%%%%%%%%%%%%
        else tempVar=[];
        end;kvai=kvai+1; N = tempVar;
        %.
        if kvai<=Kvai, tempVar=varargin{kvai} ; %%%%%%%%%%%%%%
        else tempVar=[];
        end;kvai=kvai+1; paramSet = tempVar;

        if isempty(paramSet)
            paramSet = 'rc.23';
        end

        switch paramSet
            case 'a.33'
                K = 33;
                B = K+1;
                kernelMethod = 'sinc.1';
                pctBW = 0.47680 ;
                kernelParams = pctBW;
            case 'a.13'
                K = 13;
                B = K+1;
                kernelMethod = 'sinc.1';
                pctBW = 0.47680 ;
                kernelParams = pctBW;
            case 'rc.23'
                K = 23;
                B = K+1;
                kernelMethod = 'rc.1';
                pctBW = .022;
                Tsym = 1.08 ;
                kernelParams = [pctBW , Tsym ] ;

            otherwise
                error('bad switch for param set ')
        end
        state.N = N;
        state.K = K; %kernel side samples not incl center
        state.Nkernel = 2*K+1;
        state.FsIn = FsIn;
        state.FsOut = FsOut;
        state.dt = 1/state.FsIn;
        state.dtau = 1/state.FsOut;
        state.Pad = B; %conversion window padding
        state.kernelMethod = kernelMethod;
        state.kernelParams = kernelParams;
        state.TconvWindow = N/FsIn;
        state.TrefFrame_Next = 0.0;
        state.tau_1Next = state.TrefFrame_Next ;
        state.xsave = zeros(1,2*B);
        state.outputSamplesTotal = 0; %interpolated samples done so far
        state.framecount = 0; %input frame count
        % - - - -
        varargout{1} = state;
        return ;

    case {'','resample'} %command ----------------
        %do nothing
        % now continue on
             %get inputs:

    case {'test'} %command ----------------
        % rs_irrat_test
        %

        FsOut = 2450;
        FsIn = 2400;
        anaResampFiltName = 'rc.23';

        Framesize = 98;Nframes = 1000;
        [rsState ] = rs_irrat('initialize', FsOut, FsIn, Framesize, anaResampFiltName );
        [b,a]=ellip(5,2,50,0.02*2);
        x = filtfilt(b,a,exp(i*randn(1,Nframes*Framesize)));
        X = reshape(x,Framesize,Nframes);
        clear yall

        Tk = 0;
        for k=1:Nframes
            tic;
            [y , rsState ] = rs_irrat('resample', X(:,k), rsState );
            Tk=Tk + toc;
            yall{k} = y;

        end
        yall = cat(2,yall{:});
        plot(yall);
        disp([mfilename,': exec Time = ',num2str(Tk/Nframes)])
        rsState
        return

    otherwise %command ----------------
        error('unknown command')
end %command ----------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%normal processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get inputs
if kvai<=Kvai, tempVar=varargin{kvai} ; %%%%%%%%%%%%%%
else tempVar=[];
end;kvai=kvai+1; xin = tempVar;
%.
if kvai<=Kvai, tempVar=varargin{kvai} ; %%%%%%%%%%%%%%
else tempVar=[];
end;kvai=kvai+1; state = tempVar;
%
if kvai<=Kvai, tempVar=varargin{kvai} ; %%%%%%%%%%%%%%
else tempVar=[];
end;kvai=kvai+1; tvMultMeth = tempVar;
if isempty(tvMultMeth)
    tvMultMeth = 1;
end

if size(xin,2)==1
    %orient to row:
    xin=xin.';
end
Ncheck = length(xin);
if state.N ~=Ncheck
    error('bad input frame size - does not match state.')
end
N=state.N;

% 1.  write input data to conversion buffer
xpad = [state.xsave , xin ];
state.xsave = [];
stateOut = state;

% 2.  Gen Input time vector (t)
B = state.Pad;
K = state.K;
tK = [-B+1:N+B]; %indices
t = tK*state.dt;

TconvWindow = N* state.dt ;

% 3.  Gen. Output time vector (tau)
tau_1 = state.tau_1Next ; % reference to current frame
tau_1Norm = tau_1*state.FsIn; %units: input samples

%   3.a Find M (number of output samples)
% they must fit in the window
M = floor((t(N+B) - tau_1 )*state.FsOut) + 1 ;

%   3.b output time vector (centers of kernel support)
tauNorm = tau_1Norm+ [0:M-1]*(state.dtau/state.dt); %ref to t(1) as zero time for each frame
[tauNorm_k] = round(tauNorm); %index
tauNorm_off = tauNorm- tauNorm_k; %normalized time offsets

% 4.  sample Kernel function
[hset,hsupportInd] = kernel(state,-tauNorm_off);

% 5.  T.V. convolution to gen. output samples
switch tvMultMeth
    case 1 %tvMultMeth %%%%%%%%%%%%%%%%%
        for m = 1:M
            y(m) = hset(m,:)* ...
                xpad(B+hsupportInd+tauNorm_k(m)).';
        end

    case 2 %tvMultMeth %%%%%%%%%%%%%%%%%

        %         for m=1:M,
        %             indset{m} = B+hsupportInd+tauNorm_k(m) ;
        %         end
        %        indset=cat(1,indset{:});

        Bh = repmat(B + hsupportInd,M,1);
        indset  =  Bh + repmat(tauNorm_k.',1,length(hsupportInd));
        XX = xpad(indset);
        YY  =    hset.*XX  ;
        y  = sum( YY , 2).';

    otherwise %tvMultMeth %%%%%%%%%%%%%%%%%
        error('bad switch')

end %tvMultMeth %%%%%%%%%%%%%%%%%

% 6.  Copy/shift input buffer.
stateOut.xsave = xpad(N+1:end);

% 7.  update state
stateOut.outputSamplesTotal = state.outputSamplesTotal + M;
stateOut.TrefFrame_Next = state.TrefFrame_Next + state.TconvWindow;
stateOut.framecount = state.framecount+1;
%updating/storing next frame's first output sample, assuming
%that any change in output sample rate begins at first sample of cursor frame:
stateOut.tau_1Next = tau_1  + M*state.dtau - state.TconvWindow;

%test/eval/debug
switch lower(testCode)
    case 'timevecs' %testCode
        commstr = 'length(stateOut.xsave)~= 2*B';
        disp(commstr),eval(commstr)


        subplot(2,2,1),  plot(tau,tau,'.r',t,t,'ob')
        subplot(2,2,2), plot(tK,tK,'.r', tauNorm,tauNorm,'o')

    otherwise
        %do nothing
end

% - - - - - - - - - - -

varargout{1} = y;
varargout{2} = stateOut;

end %main

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%functions %functions  %functions  %functions  %functions  %functions
%functions %functions  %functions  %functions  %functions  %functions
%functions %functions  %functions  %functions  %functions  %functions
%functions %functions  %functions  %functions  %functions  %functions
%functions %functions  %functions  %functions  %functions  %functions
%functions %functions  %functions  %functions  %functions  %functions
%functions %functions  %functions  %functions  %functions  %functions
%functions %functions  %functions  %functions  %functions  %functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h,ksupport] = kernel(state,tset, varargin)
%tset is center time
if size(tset,1)==1
    tset = tset(:);
end


switch lower(state.kernelMethod)
    case 'sinc.1' % - - - - - - - - - - - - - -
        %normalized time
        K =  state.K ; %one-sided window not incl centerpoint
        pctBW = state.kernelParams(1);

        ksupport = -K:K;
        tsupport =ksupport*(2*pctBW);
        [Tsupport,Tset] = meshgrid(tsupport,tset);
        h = sincvs(Tset+Tsupport);


    case 'rc.1'
        %%s = raisedCosine([-23:23]+.15,.022,1.08);freqz(s,1,1024,2400);

        K =  state.K ; %one-sided window not incl centerpoint
        pctBW = state.kernelParams(1);
        Tsym = state.kernelParams(2);
        ksupport = -K:K;
        J=2*K+1;
        Tgrid = [ones(length(tset),1),tset]*[ksupport;ones(1,J)];
        %same as this:
        %         [Tsupport,Tset] = meshgrid(ksupport,tset);
        %         Tgrid = Tset+Tsupport;

        h = raisedCosine(Tgrid, pctBW , Tsym );

    otherwise
        error('unknown kernel method ')
end

end % kernel gen.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x = raisedCosine(tvec,beta,T)
%raisedCosine: raisedcosine pulse shaper
%for definition see:
%     "Digital Communications", Proakis, third ed., page 546
%
%parameters:
% input:
%    tvec = vector of time samples
%    beta = excess bandwidth parameter
%    T = symbol interval
% output:
%    x = pulse samples at times specified in tvec

%john flynn 4/9/01

%note: this is generally not faster using interp1

tn = tvec.*(1.0 ./T);

x = sincvs(tn).*cos((pi*beta)*tvec)./ ...
    (1 - (2*beta*tn).^2);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
