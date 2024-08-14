function [exitStatus, resultString] = hwdiag(argumentsString)
% Call the HwDiag executable with the specified argments.

% Copyright 2001-2016 Verasonics, Inc.  All world-wide rights and
% remedies under all intellectual property laws and industrial
% property laws are reserved.  Verasonics Registered U.S. Patent and
% Trademark Office.

    switch computer
        case 'PCWIN64'
            exeExt = 'exe';
        case 'GLNXA64'
            exeExt = 'linux64';
        case 'MACI64'
            exeExt = 'mac64';
    end

    vpfRoot = getenv('VERASONICS_VPF_ROOT');
    exePath = fullfile(vpfRoot, 'System', ['HwDiag.' exeExt]);
    cmd = sprintf('"%s" %s', exePath, argumentsString);
    [exitStatus, resultString] = system(cmd);
end

