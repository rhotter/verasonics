function [RDMAconfig, RDMAsystem]=generateConfig(numSystems)
%
% [RDMAconfig, RDMAsystem]=generateConfig(numSystems)
%
% Generate a configuration information for the multisystem
% 
% Inputs:   numSystems
%
% Outputs:  RDMAconfig
%           RDMAsystem
%

% Copyright 2001-2021 Verasonics, Inc.  All world-wide rights and remedies
% under all intellectual property laws and industrial property laws are
% reserved.  Verasonics Registered U.S. Patent and Trademark Office.

% Check arguments in and out
if nargin==0
    error('You must input the number of systems for which you wish to generate a config file.')
end
if nargout~=2
    error('You must have two output arguments')
end

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

switch numSystems
    case 5
        RDMAsystem=RDMAsystem;
    case 4
        RDMAsystem=RDMAsystem(2:4);
    case 3
        RDMAsystem=RDMAsystem(2:3);
    case 2
        RDMAsystem=RDMAsystem(2);
    case 1
       error('1 system configuration is not a multisystem.'); 
    otherwise
       error('%i system configuration is not yet supported.',numSystems);
end


%--- auto determine 'role' ---'
[~,result]=unix('ifconfig');
RDMAconfig.numSecondaries = length(RDMAsystem);
if (strfind(result,RDMAsystem(1).primaryIP))  %Primary
    RDMAconfig.rdmaRole = 'primary';
    RDMAconfig.nodeIndex=0;
else %--- secondary ---%
    for a = 1:length(RDMAsystem)
        if (strfind(result,RDMAsystem(a).secondaryIP))
            RDMAconfig.rdmaRole = 'secondary';
            RDMAconfig.nodeIndex = a;
        end
    end
end

% Default values 
RDMAconfig.swSync = 1;
RDMAconfig.initializeTimeoutms = 20000;
RDMAconfig.syncTimeoutms = 5000;
RDMAconfig.writeTimeoutms = 5000;

if numSystems<5
    RDMAconfig.numVantageNodes = numSystems;
    RDMAconfig.numSecondaryNodes = numSystems-1;
    RDMAconfig.primaryVantageNode = true;
else
    RDMAconfig.numVantageNodes = 4;
    RDMAconfig.numSecondaryNodes = 4;
    RDMAconfig.primaryVantageNode = false;
end