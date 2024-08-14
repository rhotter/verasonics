function [] = ExternalClockSelect()
%ExternalClockSelect Utility to determine whether system is using internal
%or external clock source, and change the source if desired
import com.verasonics.hal.hardware.*;
hwOpenResult = Hardware.openHardware(); % Does nothing if hardware is already open.
disp(' ');
    if(hwOpenResult == HwOpenResult.success)
        clockSource = Hardware.getBool(BoolAttr.isExternalClockBeingUsed);
        if clockSource
            disp('Vantage hardware system is using an external clock source.')
        else
            disp('Vantage hardware system is using its own internal clock source.')
        end
    else
        disp('Cannot access hardware system.');
        return
    end
disp(' ');
userReq = input('Enter any character to switch clock source, or just return to leave as-is: ', 's');
disp(' ');
if isempty(userReq)
    return
end
if clockSource
    disp('Switching to internal clock source.')
else
    disp('Switching to external clock source.')
end
swResult = Hardware.loadFpgaRuntime(true, clockSource);
if strcmp(swResult, 'Success')
    newClockSource = Hardware.getBool(BoolAttr.isExternalClockBeingUsed);
    if ~clockSource && ~newClockSource
        disp('External clock source not detected; still using internal clock.');
    elseif newClockSource
        disp('Success: system is now using an external clock source.')
    else
        disp('Success: system is now using its own internal clock source.')
    end
else
    disp('    An error occurred when trying to switch the clock source.');
    disp('    Try shutting down entire system to full power off state,');
    disp('      and then restarting to clear the error condition.');
    disp(' ');
    disp('    If you intend to select an external clock source,');
    disp('      make sure it is connected and active before you try to select it.')
    return
end
end

