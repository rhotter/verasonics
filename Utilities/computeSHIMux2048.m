function [ HVMux ] = computeSHIMux2048()
% Define the HVMux structure for use with the Vantage 512 system with the
% SHI 2048-MUX captive SHI module installed.
%
% This initial developmental version is based on a direct copy of the
% existing function "computeUTAMux1024".  Once the actual hardware design
% and debug of the SHI 2048-MUX has been completed this function should be
% reviewed and updated to ensure parameter settings are correct.  Initial
% used of this function is just for software development and software
% integration testing, not for use with the actual V-512 HW system.
%
% Copyright 2001-2020 Verasonics, Inc. Verasonics Registered U.S. Patent
% and Trademark Office.
%
% May 2020 VTS-343 SHIAperture replaces old field "VDAS_Aperture"
% Initial version March 2020 for Vantage-512 integration using the 4.3.0
% software builds.
%
% Note the V-512 SHI 2048-MUX is intended only for use with the dynamic mux
% programming feature, and thus no predefined ApertureES or SHIAperture
% fields are provided. The computeUTA function for this module will
% define the element-to-channel mapping through the mux chips and SHI
% in the UTA.TransConnector array.

% The settings of 50 V. HV bias, 4.0 V. Logic rail, and 5 MHz clock rate
% are copied from the UTA 1024-MUX and may not be correct for the actual
% SHI 2048-MUX.
HVMux = struct('highVoltageRails', 50, ...
               'logicRail', 4, ...
               'clock', 5, ...
               'clockInvert', 0, ...
               'polarity', 0, ...
               'latchInvert', 0, ...
               'settlingTime', 4, ... % nominal value 4 usec
               'utaSF', 6, ...  % identifies this HVMux structure as the SHI 2048-MUX special feature
               'type', 'perEL', ...
               'ApertureES', [], ...
               'SHIAperture', []);

end
