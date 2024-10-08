function remoteMatlabTest()
activate
matlabVersion = version;
[ SysConfig, ~ ] = hwConfigCheck(1,1);
SWconfigFault = SysConfig.SWconfigFault;
HWconfigFault = SysConfig.HWconfigFault;
FPGAconfigFault = SysConfig.FPGAconfigFault;
VDAS = SysConfig.VDAS;
acqSlots = SysConfig.AcqSlots;
simLicense = SysConfig.simLicense;
reconLicense = SysConfig.reconLicense;
hwLicense = SysConfig.hwLicense;
multisystemRDMA = SysConfig.multisystemRDMA;

import com.verasonics.hal.hardware.*;
IOPanelType = Hardware.getIOPanelType;
isCopperOptical = strcmp(IOPanelType,'CopperOptical');
isBkpSwitchOptical = (PciSwitch.getPciSwitchCableType(PciSwitchId.bkp) == PciSwitchProgrammingType.optical);
isHaSwitchOptical = (PciSwitch.getPciSwitchCableType(PciSwitchId.host) == PciSwitchProgrammingType.optical);

resultStr = sprintf([...
    'MatlabVersion:     %s\n',...
    'SWconfigFault:     %i\n',...
    'HWconfigFault:     %i\n',...
    'FPGAconfigFault:   %i\n',...
    'VDAS:              %s\n',...
    'AcqSlots:          [%s]\n',...
    'IOPanelType:       %s\n',...
    'isCopperOptical:   %i\n',...
    'hcaOptical:        %i\n',...
    'ioPanelOptical:    %i\n',...
    'simLicense:        %i\n',...
    'reconLicense:      %i\n',...
    'hwLicense:         %i\n',...
    'multisystemRDMA:   %i\n'],...
    matlabVersion,...
    SWconfigFault,...
    HWconfigFault,...
    FPGAconfigFault,...
    num2str(VDAS),...
    num2str(acqSlots),...
    IOPanelType,isCopperOptical,isHaSwitchOptical,isBkpSwitchOptical,...
    simLicense,...
    reconLicense,...
    hwLicense,...
    multisystemRDMA);

disp(resultStr);