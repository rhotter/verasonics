function [ HVMux ] = computeUTAMux1024()
% Define the HVMux structure for use with the Vantage 1024 mux-based UTA
% module "UTA 1024-MUX"
% Copyright 2001-2020 Verasonics, Inc. Verasonics Registered U.S. Patent and Trademark Office.
%
% May 2020 VTS-343 SHIAperture replaces old "VDAS_Aperture"
% Initial version July 2018 for 3.5-4.0 releases using dynamic Mux
% programming feature
%
% Note this Mux adapter is intended only for use with the dynamic mux
% programming feature, and thus no predefined ApertureES or SHIAperture
% fields are provided. The computeUTA function for this adapter module will
% define the element-to-channel mapping through the mux chips and adapter
% in the UTA.TransConnector array.

% The settings of 50 V. HV bias, 4.0 V. Logic rail, and 5 MHz clock rate
% are based on HW verification testing of the UTA 1024-MUX and use with the
% Vermon 2D array probes.
HVMux = struct('highVoltageRails', 50, ...
               'logicRail', 4, ...
               'clock', 5, ...
               'clockInvert', 0, ...
               'polarity', 0, ...
               'latchInvert', 0, ...
               'settlingTime', 4, ... % nominal value 4 usec
               'utaSF', 3, ...  % identifies this HVMux structure as the UTA 1024-MUX special feature
               'type', 'perEL', ...
               'ApertureES', [], ...
               'SHIAperture', []);

end
