function RFout = externalFIR_GPU_PCT(RFin)
% RFout = externalFIR_GPU_PCT(RFin)
%
%   Simple 1-D FIR Filter on 2D data using Parallel Computing Toolbox
%   Test for comparison to Matlab CPU version and mex GPU version
%
%   Filter and RF data are transformed to frequency domain for a pointwise
%   multiplication and then transformed back to time domain.
%
% Example Test:  (signal with two frequency components, bandpass centered
% on higher frequency)
%
% t = 0.01:.01:2*pi;
% RFin1 = repmat(sin(t*94)*2^15,128,1)';
% RFin2 = repmat(sin(t*20)*2^15,128,1)';
% RFinC = RFin1 + RFin2;
% A=externalFIR_CPU_M(RFinC);
% plot(1:628, RFin1(1:628),1:628, RFin2(1:628),1:628, RFinC, 1:628, A(1:628)); %look at a single channel
% legend('Signal 1','Signal 2', 'Signal 1 + 2', 'Filtered')

persistent chs L N Nfft xzp H

if isempty(H)||~existsOnGPU(H) %initialize
    disp('initialize')
    N = size(RFin,1);%signal length
    chs = size(RFin,2);
    L = 12; %filter length (FIR)
    filterCoeffs = [-0.0282, 0.0539, 0.1229, -0.2134,  -0.2964, 0.3458, 0.3458, -0.2964, -0.2134, 0.1229, 0.0539, -0.0282];
    w = repmat(filterCoeffs, chs, 1);  %highpass filter
    Nfft = 2^(ceil(log2(L+N-1))); % or 2^nextpow2(L+M-1)
    hzp = [ w zeros(chs,Nfft-L) ]';
    H = vsGpuArray(fft(hzp)); % filter
end


%% -- Data In to GPU -- %%
% Zero pad the signal and impulse response:
if isempty(xzp)||~existsOnGPU(xzp)
    xzp = [ single(vsGpuArray(RFin)); zeros(Nfft-N,chs,'single')];
else
    xzp(1:N,:)=single(vsGpuArray(RFin));
end

%% -- Filter Operation -- %%
X = fft(xzp); % signal
Y = X .* H; % multiply the Frequency domain signal and filter
y = real(ifft(Y));

%% -- Data Out to CPU --%
RFout=int16(gather(y(L/2+(1:N),:)));  %get data from GPU
return
