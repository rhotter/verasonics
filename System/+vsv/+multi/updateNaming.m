function [RDMAconfig, RDMAcontrol, RDMAsystem] = updateNaming(RDMAconfig, RDMAcontrol, RDMAsystem)

% Find and replace master/slave fields with primary/secondary fields

if(strcmp(RDMAconfig.rdmaRole,'master'))
     RDMAconfig.rdmaRole = 'primary';
elseif(strcmp(RDMAconfig.rdmaRole,'slave'))
     RDMAconfig.rdmaRole = 'secondary';
end

% -- Update RDMAconfig fields --%
if isfield(RDMAconfig,'numSlaves')
renameStructField(RDMAconfig, 'numSlaves', 'numSecondaries');       
end
if isfield(RDMAconfig,'numSlaveNodes')
renameStructField(RDMAconfig, 'numSlaveNodes', 'numSecondaryNodes');
end
if isfield(RDMAconfig,'masterVantageNode')
renameStructField(RDMAconfig, 'masterVantageNode', 'primaryVantageNode');
end

% -- Update RDMAcontrol fields --%
if isfield(RDMAcontrol,'masterDstBufNum')
renameStructField(RDMAcontrol,'masterDstBufNum','primaryDstBufNum');
end
if isfield(RDMAcontrol,'masterDstFrameNum')
renameStructField(RDMAcontrol,'masterDstFrameNum','primaryDstFrameNum');
end
if isfield(RDMAcontrol,'slaveSrcBufNum')
renameStructField(RDMAcontrol,'slaveSrcBufNum','secondarySrcBufNum');
end
if isfield(RDMAcontrol,'slaveSrcFrameNum')
renameStructField(RDMAcontrol,'slaveSrcFrameNum','secondarySrcFrameNum');
end
if isfield(RDMAcontrol,'slaveIndex')
renameStructField(RDMAcontrol,'slaveIndex', 'secondaryIndex');
end

%-- RDMAsystem --%
if isfield(RDMAcontrol,'masterDevName')
renameStructField(RDMAsystem, 'masterDevName', 'primaryDevName')
end
if isfield(RDMAcontrol,'masterIP')
renameStructField(RDMAsystem, 'masterIP', 'primaryIP')
end
if isfield(RDMAcontrol,'slaveDevName')
renameStructField(RDMAsystem, 'slaveDevName', 'secondaryDevName')
end
if isfield(RDMAcontrol,'slaveIP')
renameStructField(RDMAsystem, 'slaveIP', 'secondaryIP')
end

