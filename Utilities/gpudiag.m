function gpudiag()
% Call the gpudiag executable for specific OS

% Copyright 2001-2022 Verasonics, Inc.  All world-wide rights and
% remedies under all intellectual property laws and industrial
% property laws are reserved.  Verasonics Registered U.S. Patent and
% Trademark Office.

if ispc()
    [status, result] = system(fullfile(getenv('VERASONICS_VPF_ROOT'), '\System\gpuDiag.exe'));
    
elseif isunix()
    [status, result] = system(fullfile(getenv('VERASONICS_VPF_ROOT'), './System/gpuDiag.linux64'));
    
elseif ismac()
    warning('GPUtoolkit is not supported for macOS');
    
else
    warning('Unable to detect system OS');
    
end

if((status==0)&&~isempty(result))
    disp(result)
else
    warning('gpuDiag command failed.  Try restarting Matlab to fix the problem.\nError reported is: \n\n%s', result);
end

end