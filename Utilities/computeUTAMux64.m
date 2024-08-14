function [ Trans ] = computeUTAMux64(TransIn)
% Define the HVMux structure for use with the Vantage 64 mux-based UTA
% module "UTA 260-MUX"
% Copyright 2001-2020 Verasonics, Inc. Verasonics Registered U.S. Patent and Trademark Office.
%
% May 2020 VTS-343 SHIAperture replaces old "VDAS_Aperture"
% July 2018 VTS-852: Add probe element mapping for Dynamic Mux Programming,
%    and error checking for probes with HVMux or probes with element to
%    connector mapping that is not one-to-one.
% July 2018 VTS-823: VDAS Aperture no longer specified here; created instead
%    through VsUpdate call to "computeHvMuxVdas Aperture"
if ~exist('TransIn', 'var')
    fprintf(2, ['\n\nERROR from computeUTAMux64: Trans Structure input argument is missing.\n' ...
    '    Also, return argument must be Trans, not Trans.HVMux.\n' ...
    '    You must modify your call to computeUTAMux64 to read as follows:\n' ...
    '           Trans = computeUTAMux64(Trans);  \n\n']);
    error(' ');
end
Trans = TransIn;
if isfield(Trans, 'HVMux')
    error('computeUTAMux64: UTA 260-Mux cannot be used with HVMux probes.');
end

% read in verbose level if it has been defined in script
try
    verbose = evalin('caller', 'Resource.Parameters.verbose');
catch
    verbose = 2;
end

Trans.HVMux = struct('highVoltageRails', 100, ...
               'logicRail', 3.3, ...
               'clock', 5, ...
               'clockInvert', 0, ...
               'polarity', 0, ...
               'latchInvert', 0, ...
               'settlingTime', 4, ... % nominal value 4 usec
               'utaSF', 2, ...  % identifies this HVMux structure as the UTA 260-MUX special feature
               'type', 'perEL', ...
               'ApertureES', [], ...
               'SHIAperture', []);

if Trans.numelements == 128
    % check if Transducer has one-to-one mapping from elements to connector
    if isequal(Trans.ConnectorES, (1:128)')
        % Create predefined set of "Tractor Tread" apertures. The Aperture
        % array consists of the 65 possible contiguous groups of 64
        % elements out of the 128.  The column number (aperture index) is
        % equal to the number of the first element signal in the aperture,
        % from 1 to 65.
        Trans.HVMux.ApertureES = zeros(128,65);
        for i = 1:65, Trans.HVMux.ApertureES((i:(i+63)),i) = (i:(i+63))'; end
    else
        if verbose > 1
            disp('STATUS NOTE UTA 260-MUX: Precomputed apertures are not available');
            disp('            for probes with non 1-1 element mapping.  Use');
            disp('            computeMuxAperture function to define your own apertures.');
        end
    end
elseif Trans.numelements == 98
%     % check if Transducer has one-to-one mapping from elements to connector
%     if isequal(Trans.ConnectorES, (1:96)')
%         % Create predefined set of "Tractor Tread" apertures. The Aperture
%         % array consists of the 33 possible contiguous groups of 64
%         % elements out of the 96.  The column number (aperture index) is
%         % equal to the number of the first element signal in the aperture,
%         % from 1 to 33.
%         Trans.HVMux.ApertureES = zeros(96,33);
%         for i = 1:33, Trans.HVMux.ApertureES((i:(i+63)),i) = (i:(i+63))'; end
%     else
%         if verbose > 1
%             disp('STATUS NOTE UTA 260-MUX: Precomputed apertures are not available');
%             disp('            for probes with non 1-1 element mapping.  Use');
%             disp('            computeMuxAperture function to define your own apertures.');
%         end
%     end
    Trans.HVMux.ApertureES = zeros(96,33);
    for i = 1:33, Trans.HVMux.ApertureES((i:(i+63)),i) = (i:(i+63))'; end
elseif Trans.numelements == 64
    % in this case a single precomputed aperture will be defined,
    % supporting all 64 elements.  First check to see if the
    % Trans.Connector mapping will allow all elements to be selected with a
    % single aperture:
    muxMap = [0:64, 1:64];
    [~,~,ActiveCh] = find(muxMap(Trans.ConnectorES + 1));
    if length(ActiveCh)~=length(unique(ActiveCh))
        if verbose > 1
            disp('STATUS NOTE UTA 260-MUX: A single precomputed aperture cannot be used');
            disp('            for the Trans.ConnectorES mapping of this 64 element probe.  Use');
            disp('            computeMuxAperture function to define your own apertures.');
        end
    else
        % all elements can be selected in a sigle aperture
        Trans.HVMux.ApertureES = Trans.ConnectorES;
        Trans.HVMux.type = 'preset'; % disable dynamic mux aperture generation
        if verbose > 1
            disp('STATUS NOTE UTA 260-MUX: A single precomputed aperture has been defined that');
            disp('            will access all elements of this 64 element probe.  Set the');
            disp('            aperture field to 1 in all of your TX and Receive structures.');
        end
    end
end

