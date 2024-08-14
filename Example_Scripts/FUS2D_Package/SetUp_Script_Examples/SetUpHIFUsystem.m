% File name SetUpHIFUsystem.m:
%
% This file shows how to correctly program the HIFU system in the SetUp script.
% Compared with regular imaging script, comment lines with the prefix HIFU,
% identify the changes added to allow the script to exercise a TPC Profile 5
% HIFU transmit-only event. No image is created.
%
% Note: This script should be executed without connecting any transducer.
% It's mainly used for testing the communication with the external power
% supply for HIFU system.
%
% Last update 08/12/2020

clear all

SysConfig = hwConfigCheck(1);

% Define system parameters.
Resource.Parameters.numTransmit = 64;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 64;   % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
if contains(SysConfig.UTAname, '260-D')
    Resource.Parameters.Connector = [1,2];
else
    Resource.Parameters.Connector = 1;
end
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.fakeScanhead = 1;
Resource.Parameters.simulateMode = 0;
Resource.System.UTA = '260-M';
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% HIFU % The Resource.HIFU.externalHifuPwr parameter must be specified in a
% script using TPC Profile 5 with the HIFU option, to inform the system
% that the script intends to use the external power supply.  This is also
% to make sure that the script was explicitly written by the user for this
% purpose, and to prevent a script intended only for an Extended Transmit
% system from accidentally being used on the HIFU system.
Resource.HIFU.externalHifuPwr = 1;

% HIFU % The string value assigned to the variable below is used to set the
% port ID for the virtual serial port used to control the external HIFU
% power supply.  The port ID was assigned by the Windows OS when it
% installed the SW driver for the power supply; the value assigned here may
% have to be modified to match.  To find the value to use, open the Windows
% Device Manager and select the serial/ COM port heading.  If you have
% installed the driver for the external power supply, and it is connected
% to the host computer and turned on, you should see it listed along with
% the COM port ID number assigned to it.
Resource.HIFU.extPwrComPortID = 'COM5';

% HIFU % The system now supports two different commercial power supplies
% for HIFU: the AIM-TTI model QPX600DP and the Sorensen model XG40-38.
% These power supplies use different command formats for remote control
% through the USB-based virtual serial port, and so the power supply
% control funtion must be told which supply is present.  This is done in
% the setup script through the field Resource.HIFU.psType, which must be
% set to a string value of either 'QPX600DP' or 'XG40-38'.  If this field
% is not specified, a default value of 'QPX600DP' will be used.
Resource.HIFU.psType = 'QPX600DP'; % set to 'QPX600DP' to match supply being used

% Fake transducer definition because the transducer should not be connected
if isequal(1,getConnectorInfo) % the default return value will be 1 if the transducer is connected
    warning('Transducer is connected, please remove the transducer for testing the external power.');
end
Trans.name = 'custom';      % use custom to prevent the id check
Trans.units = 'wavelengths';
Trans.frequency = 2;                            % nominal frequency in MHz
Trans.Bandwidth = [0.7, 1.3] * Trans.frequency;  % 60% bandwidth default value
Trans.type = 0;                                 % assuming linear array (the array type doesn't matter if focal law and image Recon calculated by users).
Trans.numelements = Resource.Parameters.numTransmit;
Trans.connType = SysConfig.UTAtype(2);
Trans.id = hex2dec('0000'); % Dummy ID to be used with the 'fake scanhead' feature
Trans.ElementPos = zeros(Trans.numelements,4);
Trans.ElementSens = ones(1,101);
Trans.lensCorrection = 0;                           % this value in wavelength will be used in image reconstruction to correct the errors of pathlength caused by acoustic lens
Trans.impedance = 50;                               % this value will be used for TXEventCheck, a self-protection utility function preventing system from overload.
Trans.maxHighVoltage = 20;                          % set the max voltage from pulser, preventing damage of a probe and the system
Trans = computeUTAMux64(Trans);

% Set up Resources.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 4096;
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;  % this should match number of receive channels
Resource.RcvBuffer(1).numFrames = 10;

% Specify TW structure array.
TW.type = 'parametric';
TW.Parameters = [Trans.frequency,0.67,1e3,1];

% Specify TX structure array.
apod = ones(1,Trans.numelements);
apod(65:end) = 0;
TX(1).waveform = 1;            % use 1st TW structure.
TX(1).focus = 0;
TX(1).Apod = apod;
TX(1).Delay = zeros(1,Trans.numelements);
TX(1).aperture = 1;  % this line was added

% Specify SeqControl structure arrays.
%  - Jump back to start.
SeqControl(1).command = 'jump';
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';  % time between frames
SeqControl(2).argument = 20000;  % 20 ms
SeqControl(3).command = 'returnToMatlab';
SeqControl(4).command = 'triggerOut';

nsc = 5; % nsc is count of SeqControl objects

n = 1; % n is count of Events

% Initial TPC profile 5
Event(n).info = 'select TPC profile 5';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = nsc; % set TPC profile command.
n = n+1;
SeqControl(nsc).command = 'setTPCProfile';
SeqControl(nsc).argument = 5;
SeqControl(nsc).condition = 'immediate';
nsc = nsc + 1;

% Send HIFU pulse
for i = 1:Resource.RcvBuffer(1).numFrames
    Event(n).info = 'HIFU pulse';
    Event(n).tx = 1;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = [2,3,4];
    n = n+1;
end

Event(n).info = 'Jump back to first event';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 1;

% Save all the structures to a .mat file.
save('MatFiles/HIFUSystem');

return
