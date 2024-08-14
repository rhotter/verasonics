function imgData = adaptiveThreshold(vel, pow, inputParams)
% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use.
%

persistent minPwr  maxPwr maxPwrC minPwrC

% Adaptation speed
ADAPTFACTORFAST = 1e-1; %runacq standard is
ADAPTFACTORSLOW = 1e-4;
MINPWRINIT = 1.0e+6;
MINPWRCINIT = 1.0e+7;
NCOLORMAP = 255;
DYNRANGEFRACTION = 0.9 ;

% Initial values and persistent variables definitions
maxPower = max(20, 50-(inputParams.Npri-12));
Gain = sqrt(12/inputParams.Npri); %rel. to nominal ens length (12)

if isempty(minPwr)
    minPwr = MINPWRINIT;
    maxPwr = 0.0;
    minPwrC = MINPWRCINIT;
end

% Initialize:
maxPwrM = maxPower; % multiply factor to convert minPwr to maxPwr
maxPwrC = minPwrC*maxPwrM;
minPwrC = min(minPwrC,minPwr);
imgData = zeros(size(vel));
maxPwr_k = zeros(size(vel));
minPwr_k = zeros(size(vel));
Norm = Gain*NCOLORMAP/(maxPwrC - DYNRANGEFRACTION*minPwrC);   % compute constant to normalize power to 256.

% Compute display dynamic range
nzind = find(pow>0.0);% identify non-zero power
minPwr_kGT = pow(nzind)<minPwr; %test for minPwr too high
minPwr_kLT =  ~minPwr_kGT;      %minPwr too low
%adjusted minmap pixel locations
minPwr_kGTInd = minPwr_kGT;
minPwr_kLTInd = minPwr_kLT;
%directionally nonlinear AR smoothing:
delta = pow - minPwr;
rind = nzind(minPwr_kGTInd);
NminAdaptFast = length(rind);
if NminAdaptFast>0
    minPwr_k(rind) = minPwr + delta(rind)*ADAPTFACTORFAST;
    minPwrMFast = mean(minPwr_k(rind));
else
    minPwrMFast = 0;
end
rind = nzind(minPwr_kLTInd);NminAdaptSlow = length(rind);
if NminAdaptSlow>0
    minPwr_k(rind) = minPwr + delta(rind)*ADAPTFACTORSLOW;
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
maxPwr_kLTInd = maxPwr_kLT;
maxPwr_kGTInd = maxPwr_kGT;
%directionally nonlinear AR smoothing:
delta = pow - maxPwr;
rind = nzind(maxPwr_kLTInd);NmaxAdaptFast = length(rind);
if NmaxAdaptFast>0
    maxPwr_k(rind) = maxPwr + delta(rind)*ADAPTFACTORFAST;
    maxPwrMFast = mean(maxPwr_k(rind));
else
    maxPwrMFast = 0;
end

rind = nzind(maxPwr_kGTInd);NmaxAdaptSlow = length(rind);
if NmaxAdaptSlow > 0
    maxPwr_k(rind) = maxPwr + delta(rind)*ADAPTFACTORSLOW;
    maxPwrMSlow = mean(maxPwr_k(rind));
else
    maxPwrMSlow = 0;
end
maxPwr = ( NmaxAdaptFast*maxPwrMFast + maxPwrMSlow*NmaxAdaptSlow ) / ...
    (NmaxAdaptSlow + NmaxAdaptFast) ;

if inputParams.PowerEstMode % COLOR POWER
    powerCeilingFactor = max(3.5, inputParams.Npri/3 - 0.5);
    flowCeiling = powerCeilingFactor*maxPwrC ;
    flowFloor = minPwrC + (maxPwrC-minPwrC)*inputParams.threshold;
    ceilingTestPassed = pow < flowCeiling;
    floorTestPassed = pow > flowFloor;
    showPixel =  find(ceilingTestPassed & floorTestPassed) ;
    
    % passed test:
    pow   = Norm*(pow - DYNRANGEFRACTION*minPwrC);
    % clip to colormap range:
    pow  = min( NCOLORMAP-1 , max(0,pow ) );
    
    imgData(showPixel) = pow(showPixel);
    
else % COLOR DOPPLER:
    flowFloor  = (minPwrC + (maxPwrC-minPwrC)*inputParams.threshold);
    showPixel = find(pow >  flowFloor);
    imgData(showPixel) =  vel(showPixel);
end

end