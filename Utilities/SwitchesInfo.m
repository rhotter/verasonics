function SwitchesInfo()
% SWITCHESINFO  Report PCI switches info
%
% Copyright 2001-2020 Verasonics, Inc.  Verasonics Registered U.S. Patent and Trademark Office.

% Aug 4 2020 VTS-1789 fix error reporting

import com.verasonics.hal.hardware.*

hwOpenResult = Hardware.openHardware(false);

if hwOpenResult ~= HwOpenResult.failedBecauseHostSwitchNotDetected && ...
        hwOpenResult ~= HwOpenResult.failedBecauseDeviceDriverRequiresReinstall && ...
        hwOpenResult ~= HwOpenResult.failedBecauseDeviceDriverNotFound
    Hardware.printPciSwitchesReport(false);
else
    disp('Failed to detect PCI switch on host adapter card;');
    disp('  cannot display any PCI Switches info.');
end