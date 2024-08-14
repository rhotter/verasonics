% compileExampleDopplerCode.m
%
% compile external mex & mexcuda code
%   Requires an NVIDIA card, NVIDIA Drivers, and CUDA 11.0 SDK
%

examplePath = fileparts(mfilename('fullpath'));


%- 1. Compile mexcuda file -%

if (ispc)
    includePath = 'C:\ProgramData\NVIDIA Corporation\CUDA Samples\v11.0\common\inc';
elseif (isunix)
    includePath = '/usr/local/cuda-11.0/samples/common/inc';
end

fileNamePath = fullfile(examplePath, '/externalFIR_GPU_MEXCUDA.cu');
vsv.gpu.mexCuda(fileNamePath,'includePaths',includePath,'libs',{'cufft','cudart'});
