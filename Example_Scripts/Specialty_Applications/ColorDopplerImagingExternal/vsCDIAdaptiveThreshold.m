function [img,pow , diags ] =vsCDIAdaptiveThreshold(ac1,inputParams)
% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use.
%
%function [img,pow,stateOut] =vsCDIAdaptiveThreshold(ac1,state)
%RunAcq.c Adaptive threshold/detection for color processing
% --------------------------------------------------------
%jaf 3/1/2012 TRANSCRIBED fron runacq.c "dopProc" routine,
% the preamble and section in runacq.c after wallfilter switch,
% branch "lowflow" in subversion matlab simulator.


%input params:
threshold = inputParams.threshold ; %
PowerEstMode = inputParams.PowerEstMode ;
%local params:
%[Doppler] = getRunAcqDopplerParams('hardCoded');
%[Doppler] = getRunAcqDopplerParams('special2');
[Doppler] = getRunAcqDopplerParams( 'computed' , inputParams.Npri );

powerCeilingFactor = Doppler.powerCeilingFactor;

%.
Npixels = prod(size(ac1));

%initial values and persistent variables definitions
maxPower = Doppler.maxPower;
minPwrInit = 1.0e+6;
minPwrCInit = 1.0e+7;
Ncolormap =  255;
Gain = sqrt(12/inputParams.Npri); %rel. to nominal ens length (12)

%state vars (will need to objectify for reentrant capability)
persistent minPwr  maxPwr maxPwrC minPwrC
if isempty(minPwr)
    minPwr = minPwrInit;
    maxPwr = 0.0;
    minPwrC = minPwrCInit;
end

DynRangeFraction = 0.9 ;

%initialize:
maxPwrM = maxPower; % multiply factor to convert minPwr to maxPwr
maxPwrC = minPwrC*maxPwrM;
minPwrC = min(minPwrC,minPwr);
[img,maxPwr_k,minPwr_k] = deal(zeros(size(ac1)));
Norm = Gain*Ncolormap/(maxPwrC - DynRangeFraction*minPwrC);   % compute constant to normalize power to 256.

pow = abs(ac1).^0.5 ; % apply some compression to power

%identify non-zero power
nzind = find(pow>0.0);

%adaptation speed
 adaptFactorFast =   1e-1; %runacq standard is
%adaptFactorFast =   2e-1;
adaptFactorSlow =  1e-4;

debugVarsString={};

%compute display dynamic range
%min power -----------------------------------
%tests
minPwr_kGT = pow(nzind)<minPwr; %test for minPwr too high
minPwr_kLT =  ~minPwr_kGT;      %minPwr too low
%adjusted minmap pixel locations
minPwr_kGTInd = find(minPwr_kGT);
minPwr_kLTInd = find(minPwr_kLT);
%directionally nonlinear AR smoothing:
delta = pow - minPwr;
rind = nzind(minPwr_kGTInd);NminAdaptFast = length(rind);
if NminAdaptFast>0,
    minPwr_k(rind) = minPwr + delta(rind)*adaptFactorFast;
    minPwrMFast = mean(minPwr_k(rind));
else
    minPwrMFast = 0;
end
rind = nzind(minPwr_kLTInd);NminAdaptSlow = length(rind);
if NminAdaptSlow>0,
    minPwr_k(rind) = minPwr + delta(rind)*adaptFactorSlow;
    minPwrMSlow= mean(minPwr_k(rind));
else
    minPwrMSlow = 0;
end
minPwr = ( NminAdaptFast*minPwrMFast + minPwrMSlow*NminAdaptSlow ) / ...
    (NminAdaptSlow + NminAdaptFast) ;

%compute max power -----------------------------------
%tests
maxPwr_kLT = pow(nzind)>maxPwr; %test for maxPwr too low
maxPwr_kGT =  ~maxPwr_kLT;      %maxPwr too high
%adjusted minmap pixel locations
maxPwr_kLTInd = find(maxPwr_kLT);
maxPwr_kGTInd = find(maxPwr_kGT);
%directionally nonlinear AR smoothing:
delta = pow - maxPwr;
rind = nzind(maxPwr_kLTInd);NmaxAdaptFast = length(rind);
if NmaxAdaptFast>0,
    maxPwr_k(rind) = maxPwr + delta(rind)*adaptFactorFast;
    maxPwrMFast = mean(maxPwr_k(rind));
else
    maxPwrMFast = 0;
end

rind = nzind(maxPwr_kGTInd);NmaxAdaptSlow = length(rind);
if NmaxAdaptSlow>0,
    maxPwr_k(rind) = maxPwr + delta(rind)*adaptFactorSlow;
    maxPwrMSlow = mean(maxPwr_k(rind));
else
    maxPwrMSlow = 0;
end
maxPwr = ( NmaxAdaptFast*maxPwrMFast + maxPwrMSlow*NmaxAdaptSlow ) / ...
    (NmaxAdaptSlow + NmaxAdaptFast) ;


if PowerEstMode,
    %COLOR POWER
    %test all pixels:
    flowCeiling = ( (powerCeilingFactor*maxPwrC)) ;
    flowFloor =  (  (minPwrC+(maxPwrC-minPwrC)*threshold) );
    %     [minPwr,maxPwr]
    %     [flowCeiling,flowFloor]
    ceilingTestPassed = pow<flowCeiling;
    PctCeilFail = sum(~ceilingTestPassed)/Npixels;
    floorTestPassed = pow>flowFloor;
    pltcTest =  ceilingTestPassed & floorTestPassed ;
    %find pixels passing:
    pltcInd = find(pltcTest);
    %find pixels failing:
    zInd = find(~pltcTest);

    %passed test:
    pow   = Norm*(pow - DynRangeFraction*minPwrC);
    %clip to colormap range:
    pow  = min( Ncolormap-1 , max(0,pow ) );

    %passed test:
    img(pltcInd) = pow(pltcInd);
    %failed test:
    img(zInd) = 0.0;

    diags.minPwr = minPwr;
    diags.minPwrC = minPwrC;
    diags.maxPwrC = maxPwrC;
    diags.Norm = Norm;
    diags.maxPwrMFast  = maxPwrMFast ;
    diags.maxPwrMSlow  = maxPwrMSlow ;
    diags.minPwrMFast  = minPwrMFast ;
    diags.minPwrMSlow  = minPwrMSlow ;
    diags.NmaxAdaptFast   = NmaxAdaptFast ;
    diags.NmaxAdaptSlow   =  NmaxAdaptSlow;
    diags.NminAdaptFast  = NminAdaptFast ;
    diags.NminAdaptSlow   =  NminAdaptSlow;
    diags.flowCeiling = flowCeiling;
    diags.pctFailedCeilingTest = PctCeilFail ;
    diags.flowFloor = flowFloor;

else

    %COLOR DOPPLER:
    %test all pixels:
    flowFloor  = (minPwrC + (maxPwrC-minPwrC)*threshold);
    pltcTest = pow >  flowFloor;

    %find pixels passing:
    pltcInd = find(pltcTest);
    %find pixels failing:
    zInd = find(~pltcTest);

    img(pltcInd) =  (Ncolormap-2)/(2*pi) * angle(ac1(pltcInd)); %(127/pi) * (+/- pi);

    %failed test:
    img(zInd) = 0.0;

    diags.minPwr = minPwr;
    diags.minPwrC = minPwrC;
    diags.maxPwrC = maxPwrC;
    diags.Norm = Norm;
    diags.maxPwrMFast  = maxPwrMFast ;
    diags.maxPwrMSlow  = maxPwrMSlow ;
    diags.minPwrMFast  = minPwrMFast ;
    diags.minPwrMSlow  = minPwrMSlow ;
    diags.NmaxAdaptFast   = NmaxAdaptFast ;
    diags.NmaxAdaptSlow   =  NmaxAdaptSlow;
    diags.NminAdaptFast  = NminAdaptFast ;
    diags.NminAdaptSlow   =  NminAdaptSlow;
    diags.flowCeiling = [];
    diags.pctFailedCeilingTest = 0.0 ;
    diags.flowFloor = flowFloor;

end


end %main  - adaptThresh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% functions %%%%%%%%% %%%%%%%%% functions %%%%%%%%% %%%%%%%%%
% functions %%%%%%%%% %%%%%%%%% functions %%%%%%%%% %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Doppler = getRunAcqDopplerParams(method,Npri);
switch method
    case 'hardCoded'
        Doppler.powerCeilingFactor = 3.5;%Doppler->powerCeilingFactor  HARDCODED in runAcq.c
        Doppler.maxPower =  50;      %Doppler->maxPower  HARDCODED in runAcq.c
    case 'special'
        Doppler.powerCeilingFactor = 40 ;
        % Doppler.powerCeilingFactor = 3.5;%Doppler->powerCeilingFactor  HARDCODED in runAcq.c
        Doppler.maxPower =  75;      %Doppler->maxPower  HARDCODED in runAcq.c
        %   Doppler.maxPower =  50;      %Doppler->maxPower  HARDCODED in runAcq.c
      case 'computed'
        dNpri = (Npri-12);
        Doppler.powerCeilingFactor = max(3.5,3.5 + dNpri/3) ;
        % Doppler.powerCeilingFactor = 3.5;%Doppler->powerCeilingFactor  HARDCODED in runAcq.c
        Doppler.maxPower =  max(20,50-dNpri);      %Doppler->maxPower  HARDCODED in runAcq.c
     otherwise
        error('param fetch method')
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
