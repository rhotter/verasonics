function [ ] = showHardwareMemoryUsage()
% showHardwareMemoryUsage: Utility function for displaying the actual
% memory usage in the hardware system after a user script has been
% initialized and run by VSX - runAcq.  This utility queries the actual
% hardware event sequence that has been programmed into the hardware
% system, so it cannot be used if a hardware system is not present and the
% script was initialized for simulate-mode operation only.
%
% Copyright 2001-2019 Verasonics, Inc. Verasonics Registered U.S. Patent and Trademark Office.
%
% USAGE: Simply type "showHardwareMemoryUsage" at the Matlab command
% prompt.  This function has no input or output arguments.  The resulting
% memory usage report will be displayed in the Matlab Command Window.

% May 2020 VTS-343: updated to use new field Trans.HVMux.SHIAperture

disp(' ');

import com.verasonics.hal.hardware.*  % enables communication with Hal
% check to see if hw system is present
hwOpenResult = Hardware.openHardware(); % Does nothing if hardware is already open.
if(hwOpenResult ~= HwOpenResult.success)
    disp('Vantage hardware system not detected;');
    disp('  Hardware system memory usage not available.')
    return
elseif evalin('base', '~exist(''Trans'', ''var'') || ~exist(''VDAS'', ''var'') || VDAS == 0')
    disp('Hardware Event sequence structures not found in Matlab Workspace;');
    disp('  Hardware system memory usage not available.')
    return
end

HVMuxScript = evalin('base', 'isfield(Trans, ''HVMux'') && ~isempty(Trans.HVMux)');
if HVMuxScript
    SHIMemKB = Hardware.getInt(IntAttr.shiHvMuxRamInstalledInKb);
    HVMuxUsedKB = evalin('base', '(size(Trans.HVMux.SHIAperture, 1) + 1) * size(Trans.HVMux.SHIAperture, 2) / 1024');
end

CGDMemGB = Hardware.getInt(IntAttr.channelGroupRamInstalledInGb);
perCHCGDMemMB = CGDMemGB * 1024 / 32; % convert to per-channel, MB units

rx = Hardware.getLong(LongAttr.bytesUsedForReceiveData);
perChRcvMB = rx/(32 * 2^20); % convert to per-channel, MB units

desc = Hardware.getLong(LongAttr.bytesUsedByDescriptors);
perChDescMB = desc/(32 * 2^20); % convert to per-channel, MB units

sgl = Hardware.getLong(LongAttr.bytesUsedBySgl); % TODO: replace memory available with actual bytes used
perChSGLMB = sgl/(32 * 2^20); % convert to per-channel, MB units

seqMemMB = Hardware.getInt(IntAttr.sequencerRamInstalledInMb);

seq = Hardware.getLong(LongAttr.bytesUsedBySequenceInstructions);
SeqMB = seq/(2^20); % convert to MB units

tgc = Hardware.getLong(LongAttr.bytesUsedByTgcWaveforms);
tgcMB = tgc/(2^20); % convert to MB units

disp('Hardware System Memory Usage for the previously run acquisition sequence:');
disp('     Per-Channel Memory ("CGD memory") Usage in MB:');
disp(['         ', num2str(perCHCGDMemMB, '%2.3f'), ' MB total memory available']);
disp(['         ', num2str(perChRcvMB, '%2.3f'), ' MB used for Receive data buffers']);
disp(['         ', num2str(perChDescMB, '%2.3f'), ' MB used for TX-RX Descriptors']);
disp(['         ', num2str(perChSGLMB, '%2.3f'), ' MB space available for DMA Scatter-Gather Lists']);
disp('     Sequence Memory Usage for Event Sequence and TGC waveforms in MB:');
disp(['         ', num2str(seqMemMB, '%2.3f'), ' MB Event Sequence memory available']);
disp(['         ', num2str(SeqMB, '%2.3f'), ' MB used for Event Sequence Instructions']);
disp(['         ', num2str(tgcMB, '%2.3f'), ' MB used for TGC Waveforms']);
disp('     Note values in MB units defined as 1 MB = 2^20 or 1,048,576 bytes');
if HVMuxScript
    disp('     HVMux probe or UTA; UTA baseboard Memory Usage for HVMux programming in KB:');
    disp(['         ', num2str(SHIMemKB, '%2.1f'), ' KB HVMux programming memory available']);
    disp(['         ', num2str(HVMuxUsedKB, '%2.1f'), ' KB used for HVMux programming tables']);
    disp('     Note values in KB units defined as 1 KB = 2^10 or 1,024 bytes');
end
disp(' ');


