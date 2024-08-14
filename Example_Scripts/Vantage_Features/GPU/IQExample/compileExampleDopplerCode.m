% compileExampleDopplerCode.m
%
% compile external mex & mexcuda code
%   Requires an NVIDIA card, NVIDIA Drivers, and CUDA 11.0 SDK
%

examplePath = fileparts(mfilename('fullpath'));

%- 1. Compile mex file -%
fileNamePath = fullfile(examplePath, 'externalDoppler_CPU_MEX.c');
mexString = sprintf('mex %s -R2018a  -outdir %s', fileNamePath, examplePath);  % compile mex program in same directory as example script
disp('Compiling externalDoppler_CPU_MEX...');
eval(mexString)  %compile mex file

%- 2. Compile mex multithreading example
filePath = fileparts(mfilename('fullpath'));
srcPath = fullfile(filePath,'externalDoppler_CPU_MULTI_MEX.c'); % build a full file path for the .c code

if ispc()  %  Windows is more complicated to compile and requires non-native pthread files
    supportPath = fileparts(which('compileExtMexMultiThreads'));
    includesDir = fullfile(supportPath,'includes');
    libsDir = fullfile(supportPath,'libs');
    dllDir = fullfile(filePath,'dll');
    setenv('PATH',   [getenv('PATH')    ';' dllDir]);  % add .dll to path
    setenv('LIB',    [getenv('LIB')     ';' libsDir]);
    setenv('INCLUDE',[getenv('INCLUDE') ';' includesDir]);
    linkerString = '-lpthreadVC2';
    mexString = sprintf("mex %s -I%s -L%s %s -R2018a", srcPath, includesDir, libsDir, linkerString); % compilation string
else
    mexString = sprintf("mex %s -R2018a", srcPath); % compilation string
end
disp('Compiling externalDoppler_CPU_MULTI_MEX...');
eval(mexString)  %compile mex file

%- 2. Compile mexcuda file -%
fileNamePath = fullfile(examplePath, 'externalDoppler_GPU_MEXCUDA.cu');
if (ispc)
    includePath = 'C:\ProgramData\NVIDIA Corporation\CUDA Samples\v11.0\common\inc';
elseif (isunix)
    includePath = '/usr/local/cuda/samples/common/inc';
end
disp('Compiling externalDoppler_GPU_MEXCUDA...');
vsv.gpu.mexCuda(fileNamePath, 'includePaths', includePath);
