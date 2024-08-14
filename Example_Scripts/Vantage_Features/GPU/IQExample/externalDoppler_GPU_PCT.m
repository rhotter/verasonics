function ImgData = externalDoppler_GPU_PCT(varargin)

% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use.
%

persistent WF  %make persistent so it doesn't need to be recalculated with each call
IQData_g = vsGpuArray(single(complex(varargin{1}, varargin{2}))); % create complex array used by this function
dataSize = size(IQData_g);

Nr = dataSize(1);
Nc = dataSize(2);
if length(dataSize)==3
    Nf = 1;
    Np = dataSize(3);
elseif length(dataSize)==4
    Nf = dataSize(3);
    Np = dataSize(4);
end

%-- Initialize Wall Filter --%
if isempty(WF)||~existsOnGPU(WF) % if WallFilter has not yet been defined, do it here
    if (nargin == 3) % Either take it from the input
        WF = vsGpuArray(single(varargin{3})); 
    else % or generate it here
        ord = ceil(Np/8);
        WF = vsGpuArray(single(polyfilter(Np,ord)));
    end
end

%-- Processing --%
% reshape 4D data to 2D for faster processing
% assume only one frame to process with WF
IQD = reshape(IQData_g, Nr*Nc*Nf, Np);

% Wallfiltering
IQDwf = IQD*WF;

% lag-1 autocorrelation
ac1 = mean(IQDwf(:,2:end).*conj(IQDwf(:,1:end-1)),2)';

% Use adaptive thresholding
threshParams = evalin('base','threshParams');
pow = abs(ac1).^0.5 ; % apply some compression to power
vel = (253)/(2*pi) * angle(ac1);
ImgData_th = adaptiveThreshold( vel, pow, threshParams);

% Reshape linear data to 2D Image
ImgData = reshape(  double(gather(ImgData_th)), Nr, Nc, Nf );

return
