% Compile Multi-threading Example Code
%
%  A script to simplify the compilation of the external .mex multithreading
% example code.
%
% Copyright (C) 2001-2022, Verasonics, Inc.  All worldwide rights and
% remedies under all intellectual property laws and industrial property
% laws are reserved.
%

currentPath = fileparts(mfilename('fullpath'));  % get the directory for this file
codePath = fullfile(currentPath,'extmexMultiThreads.c'); % build a full file path for the .c code

if ispc()  %  Windows is more complicated to compile and requires non-native pthread files
    includesDir = fullfile(currentPath,'includes');
    libsDir = fullfile(currentPath,'libs');
    dllDir = fullfile(currentPath,'dll');
    addpath(dllDir);  % add .dll to path
    linkerString = '-lpthreadVC2';
    compileString = sprintf("mex %s -I%s -L%s %s -largeArrayDims", codePath, includesDir, libsDir, linkerString); % compilation string
else
    compileString = sprintf("mex %s -largeArrayDims", codePath); % compilation string
end
disp(compileString);
eval(compileString);