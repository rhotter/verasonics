% rdmaSetup_5.m
%
% Settings and setup script for both standard and specialty system
% comment/uncomment depending on system configuration
%
warning('rdmaSetup_5:deprecated', 'rdmaSetup_5.m has been deprecated and should be replaced with vsv.multi.setup_5');
usingMultiSys = 1;
Resource.VDAS.dmaTimeout = 300*1000; % (ms) time software sequencer will wait for 'transferToHost'

%--- 5 system
%--- RDMAsystem 1---%
RDMAsystem(1).primaryDevName ='mlx5_3';
RDMAsystem(1).primaryIP = '10.10.0.1';
RDMAsystem(1).tcpPort = 5010;
RDMAsystem(1).secondaryDevName ='mlx5_0';
RDMAsystem(1).secondaryIP = '10.10.0.2';
%--- RDMAsystem 2---%
RDMAsystem(2).primaryDevName ='mlx5_0';
RDMAsystem(2).primaryIP = '10.10.1.1';
RDMAsystem(2).tcpPort = 5020;
RDMAsystem(2).secondaryDevName ='mlx5_0';
RDMAsystem(2).secondaryIP = '10.10.1.2';
%--- RDMAsystem 3---%
RDMAsystem(3).primaryDevName ='mlx5_1';
RDMAsystem(3).primaryIP = '10.10.2.1';
RDMAsystem(3).tcpPort = 5030;
RDMAsystem(3).secondaryDevName ='mlx5_0';
RDMAsystem(3).secondaryIP = '10.10.2.2';
%--- RDMAsystem 4---%
RDMAsystem(4).primaryDevName ='mlx5_2';
RDMAsystem(4).primaryIP = '10.10.3.1';
RDMAsystem(4).tcpPort = 5040;
RDMAsystem(4).secondaryDevName ='mlx5_0';
RDMAsystem(4).secondaryIP = '10.10.3.2';
%------------------------------%



%--- auto determine 'role' ---'
[status,result]=unix('ifconfig');
numSecondaries = length(RDMAsystem);
if (findstr(result,RDMAsystem(1).primaryIP))  %Primary
    RDMAconfig.rdmaRole = 'primary';
    RDMAconfig.nodeIndex=0;
    RDMAconfig.numSecondaries = length(RDMAsystem);
else %--- secondary ---%
    for a = 1:length(RDMAsystem)
        if (findstr(result,RDMAsystem(a).secondaryIP))
            RDMAconfig.rdmaRole = 'secondary';
            RDMAconfig.nodeIndex = a; %<--- for specialty
        end
    end
end

% if isempty(RDMAconfig.rdmaRole)
%     error('IP addresses are not configured for this sytem to be run in multisys')
% end
RDMAconfig.swSync = 1;
RDMAconfig.showTransferRate = 0;
RDMAconfig.initializeTimeoutms = 20000;
RDMAconfig.syncTimeoutms = 5000;
RDMAconfig.writeTimeoutms = 5000;
RDMAconfig.numVantageNodes = 4;
RDMAconfig.numSecondaryNodes = 4; %<--- 4 for specialty 3 for standard
RDMAconfig.primaryVantageNode = false; %<--- for specialty, the primary system isn't a vantage node
