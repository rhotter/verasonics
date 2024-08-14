function Trans = computeTrans(varargin)
% computeTrans - Assign array specifications for known transducer types
% Copyright 2001-2023 Verasonics Registered U.S. Patent and Trademark Office.
%
% Arguments:
%    One argument:
%        1st arg - Trans
%        if Trans is a structure, return the full Trans structure. If ID string, return transducer name only.
%    Two arguments:
%        1st arg should be transducer name as string.
%        2nd arg should be parameter desired (eg 'maxHighVoltage').
%        Returns value of parameter.
%
% Use built-in defaults, but allow prior specification of some attributes, such as Trans.frequency.
%
%
% Trans structure field definitions: (Refer to the Vantage Sequence
% Programming Manual for more details)
%
% Trans.name: string value representing probe name as used by the
%    Vantage software.  Note name is case sensitive and must be spelled
%    exactly as listed in the "KnownTransducers" array!
% Trans.id: Matlab double set to a 24 bit unsigned integer value
%    representing the probe eeprom ID code as defined and interpreted by
%    Verasonics.  Note that the id value is often listed or displayed as a
%    string of six hex digits, but the Trans.id value must always be set to
%    the numeric equivalent of that hex string!
% Trans.units: string variable that must be set to either 'mm' or
%    'wavelengths'; set by default to 'mm' if not specified by user before
%    calling computeTrans.  The units selection applies to the following
%    Trans structure fields: elementWidth, ElementPos, lensCorrection.
% Trans.type: integer-value Matlab double representing one of the probe
%    geometry types supported by the Vantage software.  Four types are
%    currently defined:
%        type = 0: Linear or phased array; all elements in a line along the
%            X-axis.  Y and Z coordinates of all elements are set to zero
%        type = 1: Curved linear array; all elements are in the X Z plane.
%            Y coordinate of all elements is zero
%        type = 2: Three-dimensional array; each element can have an
%            arbitrary X Y Z position, with an arbitrary orientation.  All
%            five entries in Trans.ElementPos are used ( X Y Z position
%            coordinates and azimuth and elevation orientation angles).  A
%            2D array with all elements in the Z=0 plane and oriented
%            straight ahead is a special case of type 2.
%        type = 3: Annular array; each element is a full-circle annular
%            ring with all elements concentric around the Z-axis through
%            the center of the array.  Each element can have an arbitrary
%            width and diameter, orientation, and offset in the Z
%            dimension.
%        type = 4: Row-Column 2D array 
% Trans.connType: integer-value Matlab double representing an index to one
%    of the probe connector types supported by the system, and identifying
%    the specific connector type and pinout used by the probe.
% Trans.frequency: Matlab double set to the nominal center frequency of the
%    probe in MHz.  The system software uses this value to define the
%    wavelength that will used by all system parameters that are specified
%    in "wavelength" units.  Receive demodulation frequency can be
%    specified independently of Trans.frequency, but since the default is
%    to set Receive.demodF = Trans.frequency, specifying a hardware
%    supported frequency here is helpful.
% Trans.Bandwidth: A [1 X 2] Matlab double array, set to the approximate
%    lower and upper bandwidth limits of the probe in MHz at a cutoff level
%    of -3 dB (one way) or -6 dB (round trip).
% Trans.numelements: Integer value Matlab double set to the number of
%    individual elements in the probe.
% Trans.ElementPos: An N X 5 Matlab double array, where N is equal to
%    Trans.numelements.  For Trans.type = 0, 1, or 2 the first three
%    entries in each row are the X Y Z coordinates of the associated
%    element and the last two entries are the orientation angles in azimuth
%    and elevation.:
%        ElementPos row:  [ X  Y  Z  azimuth  elevation ]
%    If elevation orientation is always zero, it can be
%    omitted and an N X 4 array can be used.  For Trans.type = 3 (annular
%    arrays) the first two entries in each row are the radius to the inner
%    and outer edges of the element's annular ring; the 3rd and 4th entries
%    are the radius and Z coordinates of the equal-area center of the
%    element; and the 5th entry is the element's angular orientation with
%    respect to the Z axis.
%        ElementPos row:  [ ri ro rc zc  azimuth ]
%    All position coordinates are in either mm or
%    wavelength units as selected by Trans.units.  All angular orientation
%    entries are in radians.
% Trans.elementWidth: Nominal width of individual elements, in either mm or
%    wavelength units as selected by Trans.units.
% Trans.elementLength: This optional field currently applies only to
%    row-column arrays, to specify the nominal length of individual
%    elements, in either mm or wavelength units as selected by Trans.units.
% Trans.spacing and Trans.spacingMm: Nominal center-to-center spacing of
%    individual elements; spacingMm is always in mm units and spacing is
%    always in wavelength units regardless of the setting of Trans.units.
% Trans.lensCorrection: Effective one-way path length through the lens that
%    separates the face of the transducer from the active element, in
%    either mm or wavelength units as selected by Trans.units.
% Trans.ElementSens: A 1 X 101 array of Matlab doubles representing the
%    relative off-axis sensitivity of an element as a function of angle
%    from -pi/2 to pi/2 in steps of pi/100.
% Trans.elevationApertureMm: Matlab double representing the elevation
%    aperture in mm. Applies only to "1D" arrays (Trans.type 0 or 1).
% Trans.elevationFocusMm: Matlab double representing the elevation
%    focus depth in mm. Applies only to "1D" arrays (Trans.type 0 or 1).
% Trans.maxHighVoltage: Matlab double representing the maximum allowed
%    transmit voltage setting for the probe in Volts (peak, not peak-to-peak).
% Trans.impedance: N X 2 complex double array, where the first entry in
%    each row specifies a frequency in MHz and the second entry is the
%    complex impedance in Ohms of the probe elements as seen by the system
%    at that frequency. N can be any positive integer representing the
%    number of impedance vs frequency pairs that have been specified.
% Trans.radius and Trans.radiusMm: Optional fields that apply only to
%    Trans.type values of 1, 2, or 3.  If specified for type 1 they give
%    the radius of a circular arc in the X-Z plane on which all elements of
%    the curved linear array are located.  If specified for type 2 or 3
%    they give the radius of a spherical surface upon which all elements
%    of the array are located.  radiusMm is always in mm units and radius
%    is always in wavelength units regardless of the setting of
%    Trans.units.

% revision history
% Mar 2023 VTS-5946 Correct specs for L11-5gH, C5-2gH, P4-2gH
% Apr 2021 VTS-2152 New L39-21gD Probe;VTS-2225 L10-4v support
% Oct 2020 VTS-1698 two versions of H-313 FUS probe
% Oct 2020 VTS-1970 Add P5-64vN and P5-128vN Imasonics NDE/NDT probes
% May 9, 2020 VTS-343 SHIAperture replaces VDAS_Aperture to allow 256-byte
%       HVMux programming tables.
% Aug 2019 4.2.0 VTS-1413 remove L10-4v
% Aug 2019 4.1.0 VTS-1351 Element bias change for L22-14vX and L22-14vX LF
% Jan 2019 VTS-1053 Freq. and BW corrections for 5 Verasonics probes
% Jan 2019 VTS-1049 Element bias change for L22-14vX and L22-14vX LF
% Aug 2018 VTS-738 UTA 256-Direct and Vermon Matrix array probe entries
% Jul 2018 VTS-843, 852 Dynamic HV Mux programming support
% Apr 2018 for 3.4.2 and 3.5, to fix bugs and update parameters for
% the FUS Elite probes, and for GE9LD.  See VTS-767, 769, 796.
% Jan 2018 to add FUS Elite probes, new Trans.type = 3 for annular arrays,
%   and updates to Trans structure field definitions in the comments
% Dec 13, 2017 to add L22-14vX and L22-14vX LF

% Known transducers and their corresponding ID, high voltage limit, and HVMux status. *** Do not modify these values. ***
KnownTransducers = {'L7-4',       '000250',  96,   0;...
                    'L10-5',      '00074C',  70,   1;...
                    'L10-4v',     '02ABA4',  75,   0;...
                    'L11-4v',     '01ABB4',  75,   0;...
                    'L11-5',      '000351',  96,   0;...
                    'L11-5v',     '05ABB5',  90,   0;...
                    'L11-5gH',    '15ABB5',  90,   0;... % VTS-2560 L11-5v but with UTA 408-GE connector
                    'L12-3v',     '01BBC3',  75,   1;...
                    'L12-5 38mm', '000755',  75,   1;...
                    'L12-5 50mm', '000B5B',  75,   1;...
                    'L22-8v',     '018A8A',  35,   1;...
                    'L22-14v',    '02AB18',  30,   0;...
                    'L22-14v LF',  '02AB17',  30,   0;...
                    'L22-14vX',   '02DB18',  45,   0;...
                    'L22-14vX LF','02DB17',  45,   0;...
                    'L35-16vX',   '02AB28',  40,   0;...
                    'L38-22v',    '03BB30',  35,   1;...
                    'L39-21gD',   '16AB30',  35,   0;... % VTS-2152 new ID and probe naming convention
                    'GE9LD'       '270212',  50,   0;...
                    '4DL7',       '02A4D7',  50,   0;... % no data
                    'GEC1-6D',    '27023E',  50,   0;...
                    'GE4CD',      '270200',  50,   0;...
                    'GEIC5-9D',   '270210',  50,   0;...
                    'GEL3-12D',   '270243',  50,   1;...
                    'GEM5ScD',    '27022E',  50,   0;...
                    'CL10-5',     '00034D',  96,   0;...
                    'CL15-7',     '00035C',  96,   0;...
                    'C4-2',       '0020D1',  96,   0;...
                    'C5-2',       '0020D9',  96,   0;...
                    'C5-2v',      '01AC52',  96,   0;...
                    'C5-2gH',     '15AC52',  96,   0;... % Humanscan curved array with UTA 408-GE connector VTS-2560
                    'C7-4',       '00224E',  96,   0;... % no data
                    'C8-4V',      '00228C',  96,   0;... % no data
                    'C8-5',       '0022DE',  96,   0;... % no data
                    'C9-5ICT',    '00228B',  96,   0;...
                    'P3-2',       '004428',  96,   0;... % no data
                    'P4-1',       '00483E',  96,   0;...
                    'P4-2',       '004439',  96,   0;...
                    'P4-2v',      '01AA42',  96,   0;...
                    'P4-2gH',     '15AA42',  96,   0;... % Humanscan phased array with UTA 408-GE connector VTS-2560
                    'P5-3',       '004529',  96,   0;... % no data
                    'P6-3',       '004D3B',  96,   0;...
                    'P7-4',       '00462A',  96,   0;...
                    'RC6gV',      '12A584',  96,   0;... % Vermon 128 X 128 6 MHz row-column array
                    'RC15gV',     '12A515',  30,   0;... % Vermon 80 X 80 15 MHz row-column array
                    'P5-64vN'     'FFFFFF',  96,   0;... % Imasonic linear I-PEX array 64-element 
                    'P5-128vN'    'FFFFFF',  96,   0;... % Imasonic linear I-PEX array 128-element
                    'IP-104',     '04AA52',  96,   0;... % FUS Elite imaging Probe, pairings 4, 5, 6
                    'IP-105',     '04AA73',  96,   0;... % FUS Elite imaging Probe, pairings 1, 2, 3
                    'H-101',      '04CD0B',  60,   0;... % FUS Elite 4 el HIFU Probe, pairing 2
                    'H-104',      '04CD05',  60,   0;... % FUS Elite 3 el HIFU Probe, pairing 1
                    'H-106',      '04CD22',  60,   0;... % FUS Elite 8 el HIFU Probe, pairing 3
                    'H-301',      '04C30B',  76.3, 0;... % FUS Elite 128 el HIFU Probe, pairing 5
                    'H-302',      '04C322',  96,   0;... % FUS Elite 128 el HIFU Probe, pairing 6
                    'H-313',      '04C304',  50,   0;... % FUS Elite 64 el HIFU Probe, pairing 4 serial numbers 02 and up
                    'H-313A',     '04C305',  50,   0;... % FUS Elite 64 el HIFU Probe, pairing 4 serial number 01 only
                    'Matrix1024-3',              '0D0400',  45,   0;... % Vermon 1024 element 3 MHz Matrix array
                    'Matrix1024-8',              '0D0400',  45,   0;... % Vermon 1024 element 8 MHz Matrix array
                    'UTA 1024-MUX',              '0D0400',  45,   1;... % no data
                    'UTA 256-Direct',            '0D0100',  45,   0;... % no data
                    'HIFU Adapter S.T.E',        '01FFF8',  96,   0;... % no data
                    'Adapter Embedded S.T.E',    '01FFFA',  96,   0;... % no data
                    '260 ZIF Backshell Kit',     '0BAD00',  96,   0;... % no data
                    '408 ZIF Backshell Kit',     '0BAD02',  96,   0;... % no data
                    '408-GE ZIF Backshell Kit',  '0BAD04',  96,   0;... % no data
                    '260 ZIF S.T.E. Fixture',    '01FFFC',  96,   0;... % no data
                    '408 ZIF S.T.E. Fixture',    '01FFFB',  96,   0;...   % no data
                    '260 ZIF Calibration STE',   '01FFF9',  96,   0};


% The probes listed above with the comment "no data" are only recognized in
% terms of the probe name and ID; there is no other Trans structure data
% provided by the computeTrans function for these probes.


switch nargin

    case 1
        Trans = varargin{1};
        if ~isstruct(Trans)  % if a structure is not provided as input, assume input is ID to translate into string.
            % input argument may be ID as a hex string, or the actual ID
            % numeric value
            if ischar(Trans)
                % convert hex string to number so it won't matter how many
                % leading zeros were provided
                probeID = hex2dec(Trans);
            else
                % not a string so presumably a numeric value
                probeID = round(Trans); % make sure we have an integer value
            end
            if probeID < 1 || probeID > (2^24 - 2)
                % not a valid ID value; must be in the range of a 24 bit
                % unsigned integer, and the values zero and all ones are
                % also excluded
                Trans = 'InvalidProbeID';
            else
                probeIDhex = num2str(probeID, '%06X');
                probenum = find(strcmpi(probeIDhex, KnownTransducers(:, 2)), 1, 'last');
                if isempty(probenum)
                    Trans = 'Unknown';
                else
                    Trans = KnownTransducers{probenum, 1}; % return the probe name (string value)
                end
            end
            return
        end
        if ~isfield(Trans,'name'), error('computeTrans: Trans.name must be provided in input structure.'); end
        probenum = find(strcmpi(Trans.name, KnownTransducers(:, 1)), 1);
        if isempty(probenum), error('computeTrans: Trans.name not recognized as known transducer.'); end
        speedOfSound = 1.540;  % default speed of sound in mm/usec
        verbose = 2;
        if evalin('base','exist(''Resource'',''var'')&&isfield(Resource,''Parameters'')')
            if evalin('base','isfield(Resource.Parameters,''speedOfSound'')')
                speedOfSound = evalin('base','Resource.Parameters.speedOfSound')/1000; % speed of sound in mm/usec
            end
            if evalin('base','isfield(Resource.Parameters,''verbose'')')
                verbose = evalin('base','Resource.Parameters.verbose');
            else
                verbose = 2; % default value- display warnings and status messages
            end
        end

        % check for user-specified units, and print warning message if not found
        if ~isfield(Trans,'units') || isempty(Trans.units)
            if verbose > 0
                fprintf(2, 'Warning: Trans.units not specified; selecting default units of mm.\n');
                fprintf(2, 'If script requires wavelength units, add an explicit definition of\n');
                fprintf(2, '"Trans.units = ''wavelengths'';" before calling computeTrans.\n');
            end
            Trans.units = 'mm';
        end
        if ~strcmp(Trans.units, 'mm') && ~strcmp(Trans.units, 'wavelengths')
            error('computeTrans: Unrecognized value for Trans.units.  Must be ''mm'' or ''wavelengths''.');
        end

        % if Trans.frequency value has already been specified, we will use
        % it as is.  VSX and update() will confirm the value matches the
        % A/D sample rate constraints, and will exit with an error message
        % to the user if not.  Therefore we do not need to validate the
        % Trans.frequency value here (and could not, since we don't know
        % intended use of 4/3 sampling or interleave, etc.).
        if isfield(Trans,'frequency')
            if isempty(Trans.frequency)
                % if empty, remove it so cases below will assign default frequency
                Trans = rmfield(Trans, 'frequency');
            end
        end
        % also allow user-specified Bandwidth to override the default:
        if isfield(Trans,'Bandwidth')
            if isempty(Trans.Bandwidth)
                % if empty, remove it so cases below will assign default
                % Bandwidth
                Trans = rmfield(Trans, 'Bandwidth');
            end
        end

        Trans.lensCorrection = 0; % specify default value, in case it is not set for a particular transducer;
        Trans.id = hex2dec(KnownTransducers(probenum, 2));

        switch KnownTransducers{probenum, 1}
            case 'L7-4'
                %% L7-4 HDI Probe
                if ~isfield(Trans,'frequency'), Trans.frequency = 5.208; end % nominal frequency in MHz
                % Vantage:  5.208 is closest supported frequency to 5 MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [4, 7]; end
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI connector=1
                Trans.numelements = 128;
                Trans.ConnectorES = (1:Trans.numelements)'; % 1:1 mapping from transducer elements to element signals at connector
                Trans.elementWidth = .250; % width in mm
                Trans.spacingMm = .298;   % Spacing between elements in mm.
                    Trans.elevationApertureMm = 7.5; % active elevation aperture in mm (estimate)
                    Trans.elevationFocusMm = 25; % nominal elevation focus depth from lens on face of transducer (estimate)
                Trans.ElementPos = zeros(Trans.numelements,5);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                Trans.lensCorrection = 0.887; % in mm units; was 3 wavelengths;
                Trans.impedance = [3 23.9-125i;  3.25 25.4-116i;  3.5 26-106i;  3.75 25.4-98.9i;  4 25.9-89.4i;  4.25 27-79.7i;...
                    4.5 32.8-72.6i;  4.75 39.2-66.2i;  5 46.1-69.6i;  5.25 46.5-72.4i;  5.5 41.9-71.6i;  5.75 43.2-69.8i;...
                    6 42.3-69.8i;  6.25 38.2-71i;  6.5 33.5-66.2i;  6.75 32-59.8i;  7 34.4-54.2i;  7.25 37.4-50.3i;...
                    7.5 42.3-48.2i;  7.75 47.8-47.9i;  8 53-51.3i];

            case 'L10-5'
                %% L10-5 HDI Probe
                if ~isfield(Trans,'frequency'), Trans.frequency = 7.813; end % nominal frequency in MHz
                % Vantage:  7.813 and 6.944 are closest supported frequencies to 7.5 MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = 7.5*[0.7, 1.3]; end % default assumed value of 60% of center frequency
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 192;
                Trans.elementWidth = .1729; % width in mm
                Trans.spacingMm = .1979;   % Spacing between elements in mm.
                    Trans.elevationApertureMm = 0; % active elevation aperture in mm (unknown)
                    Trans.elevationFocusMm = 0; % nominal elevation focus depth from lens on face of transducer (unknown)
                Trans.ElementPos = zeros(Trans.numelements,5);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                Trans.lensCorrection = 1.183; % in mm units; was 6 wavelengths;
                Trans.impedance = 51; % using default value for MUX probe
                Trans.HVMux = struct('highVoltageRails', 90, ...
                                     'logicRail', 10.5, ...
                                     'clock', 5, ...
                                     'clockInvert', 0, ...
                                     'polarity', 0, ...
                                     'latchInvert', 0, ...
                                     'type', 'perEL', ...
                                     'utaSF', 0, ...
                                     'settlingTime', 4, ... % nominal value 4 usec
                                     'ApertureES', zeros(Trans.numelements,65));
                % Add MuxMap table and length; these are used to create the
                % SHIAperture table for actually programming the HVMux
                % switches:
                [Trans.HVMux.MuxMap, Trans.HVMux.SHIAperLgth] = getMuxMap(Trans.name);
                for i = 0:64, Trans.HVMux.ApertureES(:,i+1) = [zeros(1,i),mod((i:i+127),128)+1,zeros(1,64-i)]'; end
                % create the HVMux.SHIAperture array to match
                % HVMux.ApertureES
                Trans = computeHvMuxSHIAperture(Trans);
                % define mapping through the HVMux switches from transducer
                % elements to element signals at connector:
                Trans.ConnectorES = [1:128 1:64]';

            case 'L10-4v'
                %% L10-4v Verasonics Probe
                % USER NOTE: This is a discontinued probe, and no example
                % scripts are provided for it in the current Vantage
                % software releases.  However, you can easily create
                % scripts for it by copying and modifying any L11-5v
                % example script since the L10-4v and L11-5v are virtually
                % identical in terms of center frequency and element pitch.
                % To create a script for the L10-4v, simply make a copy of
                % any L11-5v example script and then change the name of the
                % script to include "L10-4v" in the file name instead of
                % "L11-5v".  Then open the script with the Matlab editor,
                % and do a 'find and replace' to find all instances of
                % "L11-5v" and replace them with "L10-4v".  This will
                % update the call to computeTrans and the name of the .mat
                % file that will be created, and will also update any
                % comments, image display titles, etc. that may be present
                % in the script.
                if ~isfield(Trans,'frequency'), Trans.frequency = 7.5; end % nominal frequency in MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [Trans.frequency*(1-0.8/2), Trans.frequency*(1+0.8/2)]; end
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 128;
                Trans.ConnectorES = (1:Trans.numelements)'; % 1:1 mapping from transducer elements to element signals at connector
                Trans.elevationApertureMm = 4; % active elevation aperture in mm
                Trans.elevationFocusMm = 20; % nominal elevation focus depth from lens on face of transducer
                % Note elevation aperture and focus will not be converted
                % to wavelengths regardless of Trans.units
                Trans.elementWidth = 0.27; % element width in mm; assumes 50 micron kerf
                Trans.spacingMm = 0.30;   % Spacing between elements in mm.
                Trans.ElementPos = zeros(Trans.numelements,5);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                Trans.lensCorrection = 0.4; % velocities in m/sec; result in mm
                Trans.impedance = 50;
                if ~isfield(Trans,'maxHighVoltage'), Trans.maxHighVoltage = 50; end % much restricted than specified on data sheet

            case 'L11-4v'
                 %% L11-4v Verasonics Probe
                if ~isfield(Trans,'frequency'), Trans.frequency = 7.24; end % Center frequency in MHz from mfr's data sheet
                % Fractional bandwidth (-6 dB round trip) is 93.7% from mfr's data sheet
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [3.85, 10.63]; end
                % Bandwidth calculated from mfr's Center Frequency and fractional bandwidth
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 128;
                Trans.ConnectorES = (1:Trans.numelements)'; % 1:1 mapping from transducer elements to element signals at connector
                Trans.elementWidth = 0.270; % width in mm
                Trans.spacingMm = 0.300;   % Spacing between elements in mm.
                    Trans.elevationApertureMm = 5; % active elevation aperture in mm (estimae)
                    Trans.elevationFocusMm = 20; % nominal elevation focus depth from lens on face of transducer (spec)
                Trans.ElementPos = zeros(Trans.numelements,5);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                Trans.lensCorrection = 1.4785; % in mm units; was 5 wavelengths
                Trans.impedance = [ 3 17.4-123i;  3.25 18.6-110i;  3.5 20.5-100i;  3.75 20.4-91i;  4 23.3-83.3i;  4.25 21.1-78.3i;...
                    4.5 21.2-69.8i;  4.75 20.8-65.3i;  5 19.6-57.8i;  5.25 21.5-52.2i;  5.5 20.3-47.6i;  5.75 20.4-40.8i;...
                    6 21.3-37.2i;  6.25 20.4-31.5i;  6.5 22.7-25.8i;  6.75 23.1-22.8i;  7 23.5-17.3i;  7.25 26-13.8i;...
                    7.5 25.7-9.62i;  7.75 28.3-4.22i;  8 31.2-2.05i;  8.25 32.2+1.5i;  8.5 36.8+4.35i;  8.75 37.5+3.61i;...
                    9 38.7+7.27i;  9.25 39.7+7.01i;  9.5 38.5+11.5i;  9.75 42+15.4i;  10 42.9+18.2i;  10.25 47+24i;...
                    10.5 53.9+25.6i;  10.75 60.1+27.2i;  11 70.7+25.1i;  11.25 77.2+16.7i;  11.5 82.2+5.99i;  11.75 77-9.58i;...
                    12 65.4-15.8i;  12.25 54.3-17.7i;  12.5 44-15.6i;  12.75 35.2-10.6i;  13 29.1-3.17i;  13.25 26.2+4.47i;...
                    13.5 25.3+11i;  13.75 25.2+16.5i;  14 25.3+20.8i;  14.25 25+24.8i;  14.5 24.6+29.1i;  14.75 24.6+33.6i;...
                    15 25+38.1i;  15.25 25.6+42.4i;  15.5 26.5+46.5i;  15.75 27.5+50.7i;  16 28.9+54.9i];

            case 'L11-5'
                %% L11-5 HDI Probe
                if ~isfield(Trans,'frequency'), Trans.frequency = 7.813; end % nominal frequency in MHz
                % Vantage:  7.813 and 6.944 are closest supported frequencies to 7.5 MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = 7.5*[0.7, 1.3]; end % default assumed value of 60% of center frequency
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 128;
                Trans.ConnectorES = (1:Trans.numelements)'; % 1:1 mapping from transducer elements to element signals at connector
                Trans.elementWidth = 0.235; % width in mm
                Trans.spacingMm = 0.260;   % Spacing between elements in mm.
                    % Trans.elevationApertureMm = 'unknown; % active elevation aperture in mm (estimate)
                    % Trans.elevationFocusMm = 'unknown; % nominal elevation focus depth from lens on face of transducer (estimate)
                Trans.ElementPos = zeros(Trans.numelements,5);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                Trans.impedance = 50; % using default value

            case 'L11-5v'
                %% L11-5v Verasonics Probe
                if ~isfield(Trans,'frequency'), Trans.frequency = 7.6; end % as measured by Verasonics
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [4.68, 10.52]; end % 76.8% fractional bandwidth
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 128;
                Trans.ConnectorES = (1:Trans.numelements)'; % 1:1 mapping from transducer elements to element signals at connector
                Trans.elementWidth = 0.270; % width in mm
                    Trans.elevationApertureMm = 5; % (mm) active elevation aperture (spec)
                    Trans.elevationFocusMm = 18; % (mm) nominal elevation focus depth from lens on face of transducer (spec)
                Trans.spacingMm = 0.300;   % Spacing between elements in mm.
                Trans.ElementPos = zeros(Trans.numelements,5);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                Trans.lensCorrection = 0.6; % in mm units;
                Trans.impedance = 30;

            case 'L11-5gH'
                %% L11-5gH Verasonics Probe using the UTA 408-GE connector
                % This probe is identical to the L11-5v, but uses the UTA 408-GE
                % connector instead of the HDI connector.  Both versions
                % are manufactured by Humanscan.
                if ~isfield(Trans,'frequency'), Trans.frequency = 7.6; end % as measured by Verasonics
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [4.68, 10.52]; end % 76.8% fractional bandwidth
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 7; % UTA 408-GE connector
                Trans.numelements = 128;
                Trans.ConnectorES = (1:Trans.numelements)'; % 1:1 mapping from transducer elements to element signals at connector
                Trans.elementWidth = 0.270; % width in mm
                    Trans.elevationApertureMm = 5; % (mm) active elevation aperture (spec)
                    Trans.elevationFocusMm = 18; % (mm) nominal elevation focus depth from lens on face of transducer (spec)
                Trans.spacingMm = 0.300;   % Spacing between elements in mm.
                Trans.ElementPos = zeros(Trans.numelements,5);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                Trans.lensCorrection = 0.6; % in mm units;
                Trans.impedance = 30;

            case 'L12-3v'
                %% L12-3v Verasonics Probe
                if ~isfield(Trans,'frequency'), Trans.frequency = 7.54; end % Center frequency in MHz from mfr's data sheet
                % Fractional bandwidth (-6 dB round trip) is 93.1% from mfr's data sheet
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [4.03, 11.05]; end
                % Bandwidth calculated from mfr's Center Frequency and fractional bandwidth
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 192;
                Trans.elementWidth = .170; % width in mm
                Trans.spacingMm = .200;   % Spacing between elements in mm.
                    Trans.elevationApertureMm = 5; % active elevation aperture in mm (estimate)
                    Trans.elevationFocusMm = 20; % nominal elevation focus depth from lens on face of transducer (spec)
                Trans.ElementPos = zeros(Trans.numelements,5);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                Trans.lensCorrection = 1.183; % in mm units; was 6 wavelengths;
                Trans.impedance = 51; % using default value for MUX probe
                Trans.HVMux = struct('highVoltageRails', 90, ...
                                     'logicRail', 5.0, ...
                                     'clock', 5, ...
                                     'clockInvert', 0, ...
                                     'polarity', 0, ...
                                     'latchInvert', 0, ...
                                     'type', 'perEL', ...
                                     'utaSF', 0, ...
                                     'settlingTime', 4, ... % nominal value 4 usec
                                     'ApertureES', zeros(Trans.numelements,65));
                % Add MuxMap table and length; these are used to create the
                % SHIAperture table for actually programming the HVMux
                % switches:
                [Trans.HVMux.MuxMap, Trans.HVMux.SHIAperLgth] = getMuxMap(Trans.name);
                for i = 0:64, Trans.HVMux.ApertureES(:,i+1) = [zeros(1,i),mod((i:i+127),128)+1,zeros(1,64-i)]'; end
                % create the HVMux.SHIAperture array to match
                % HVMux.ApertureES
                Trans = computeHvMuxSHIAperture(Trans);
                % define mapping through the HVMux switches from transducer
                % elements to element signals at connector:
                Trans.ConnectorES = [1:128 1:64]';

            case 'L12-5 38mm'
                %% L12-5 38mm HDI Probe
                if ~isfield(Trans,'frequency'), Trans.frequency = 7.813; end % nominal frequency in MHz
                % Vantage:  7.813 and 6.944 are closest supported frequencies to 7.5 MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [5, 11]; end
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 192;
                Trans.elementWidth = .1729; % width in mm
                Trans.spacingMm = .1979;   % Spacing between elements in mm.
                    Trans.elevationApertureMm = 7.5; % active elevation aperture in mm (estimate)
                    Trans.elevationFocusMm = 20; % nominal elevation focus depth from lens on face of transducer (estimate)
                Trans.ElementPos = zeros(Trans.numelements,5);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                Trans.lensCorrection = 2.365; % in mm units; was 12 wavelengths;
                Trans.impedance = 51; % using default value for MUX probe
                Trans.HVMux = struct('highVoltageRails', 100, ...
                                     'logicRail', 10.5, ...
                                     'clock', 8, ...
                                     'clockInvert', 0, ...
                                     'polarity', 0, ...
                                     'latchInvert', 0, ...
                                     'type', 'perEL', ...
                                     'utaSF', 0, ...
                                     'settlingTime', 4, ... % nominal value 4 usec
                                     'ApertureES', zeros(Trans.numelements,65));
                % Add MuxMap table and length; these are used to create the
                % SHIAperture table for actually programming the HVMux
                % switches:
                [Trans.HVMux.MuxMap, Trans.HVMux.SHIAperLgth] = getMuxMap(Trans.name);
                for i = 0:64, Trans.HVMux.ApertureES(:,i+1) = [zeros(1,i),mod((i:i+127),128)+1,zeros(1,64-i)]'; end
                % create the HVMux.SHIAperture array to match
                % HVMux.ApertureES
                Trans = computeHvMuxSHIAperture(Trans);
                % define mapping through the HVMux switches from transducer
                % elements to element signals at connector:
                Trans.ConnectorES = [1:128 1:64]';

            case 'L12-5 50mm'
                %% L12-5 50mm HDI Probe
                if ~isfield(Trans,'frequency'), Trans.frequency = 7.813; end % nominal frequency in MHz
                % Vantage:  7.813 is closest supported frequency to 8.18 MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [5, 11]; end
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 256;
                Trans.elementWidth = .1703; % width in mm
                Trans.spacingMm = .1953;   % Spacing between elements in mm.
                    Trans.elevationApertureMm = 7.5; % active elevation aperture in mm (estimate)
                    Trans.elevationFocusMm = 20; % nominal elevation focus depth from lens on face of transducer (estimate)
                Trans.ElementPos = zeros(Trans.numelements,5);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                Trans.lensCorrection = 2.365; % in mm units; was 12 wavelengths;
                Trans.impedance = 51; % using default value for MUX probe
                Trans.HVMux = struct('highVoltageRails', 90, ...
                                     'logicRail', 10.5, ...
                                     'clock', 8, ...
                                     'clockInvert', 0, ...
                                     'polarity', 0, ...
                                     'latchInvert', 0, ...
                                     'type', 'perCH', ...
                                     'utaSF', 0, ...
                                     'settlingTime', 4, ... % nominal value 4 usec
                                     'ApertureES', zeros(Trans.numelements,129));
                % Add MuxMap table and length; these are used to create the
                % SHIAperture table for actually programming the HVMux
                % switches:
                [Trans.HVMux.MuxMap, Trans.HVMux.SHIAperLgth] = getMuxMap(Trans.name);
                for i = 0:128, Trans.HVMux.ApertureES(:,i+1) = [zeros(1,i),mod((i:i+127),128)+1,zeros(1,128-i)]'; end
                % create the HVMux.SHIAperture array to match
                % HVMux.ApertureES
                Trans = computeHvMuxSHIAperture(Trans);
                % define mapping through the HVMux switches from transducer
                % elements to element signals at connector:
                Trans.ConnectorES = [1:128 1:128]';

            case 'L22-8v'
                %% L22-8v Kolo CMUT probe
                if ~isfield(Trans,'frequency'), Trans.frequency = 15.625; end % nominal frequency in MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [8, 21.5]; end % estimate. 90% relative bandwidth
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 256;
                Trans.elementWidth = .0703; % width in mm
                Trans.spacingMm = 0.108;   % Spacing between elements in mm.
                    % Trans.elevationApertureMm = 'unknown; % active elevation aperture in mm (estimate)
                    % Trans.elevationFocusMm = 'unknown; % nominal elevation focus depth from lens on face of transducer (estimate)
                Trans.ElementPos = zeros(Trans.numelements,5);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                Trans.lensCorrection = 0; % in mm units; was 12 wavelengths;
                Trans.impedance = 5; % using default value for MUX probe
                if ~isfield(Trans,'maxHighVoltage'), Trans.maxHighVoltage = 35; end
                Trans.HVMux = struct('highVoltageRails', 90, ...
                                     'logicRail', 6.5, ... %10.5, ...
                                     'clock', 8, ...
                                     'clockInvert', 0, ...
                                     'polarity', 0, ...
                                     'latchInvert', 0, ...
                                     'type', 'perEL', ...
                                     'utaSF', 0, ...
                                     'settlingTime', 4, ... % nominal value 4 usec
                                     'ApertureES', zeros(Trans.numelements,129));
                % Add MuxMap table and length; these are used to create the
                % SHIAperture table for actually programming the HVMux
                % switches:
                [Trans.HVMux.MuxMap, Trans.HVMux.SHIAperLgth] = getMuxMap(Trans.name);
                for i = 0:128, Trans.HVMux.ApertureES(:,i+1) = [zeros(1,i),mod((i:i+127),128)+1,zeros(1,128-i)]'; end
                % create the HVMux.SHIAperture array to match
                % HVMux.ApertureES
                Trans = computeHvMuxSHIAperture(Trans);
                % define mapping through the HVMux switches from transducer
                % elements to element signals at connector:
                Trans.ConnectorES = [1:128 1:128]';

            case 'L22-14v'
                %% L22-14v Verasonics Probe
                if ~isfield(Trans,'frequency'), Trans.frequency = 15.625; end % nominal frequency in MHz
                % Note: use 15.625 MHz for 4X sampling (62.5 MHz sample rate),
                % or 18.75 MHz for 4/3 sampling (25.0 MHz sample rate, 50 MHz A/D rate)
                % Manufacturer specified center frequency is 18.0 MHz +/- 10%
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [14, 22]; end
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 128;
                Trans.ConnectorES = (1:Trans.numelements)'; % 1:1 mapping from transducer elements to element signals at connector
                    Trans.elevationApertureMm = 1.5; % active elevation aperture in mm
                    Trans.elevationFocusMm = 8; % nominal elevation focus depth from lens on face of transducer
                Trans.elementWidth = 0.08; % element width in mm; assumes 20 micron kerf
                Trans.spacingMm = 0.100;   % Spacing between elements in mm.
                Trans.ElementPos = zeros(Trans.numelements,5);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                % Lens Correction: From mfr data sheet, matching layers are 0.05 mm thick at 2145 m/sec
                % average velocity, and lens is 0.48 mm thick at 1147 m/sec.
                % Thus the net effective lens thickness in mm is given by the following
                % expression, which evaluates to 0.6804 mm for 1540 m/sec velocity
                Trans.lensCorrection = 1000 * speedOfSound * (0.05/2145 + 0.48/1147); % velocities in m/sec; result in mm
                Trans.impedance = [10.00, 11.16-54.28i; 11.00, 12.02-44.73i; 12.00, 14.38-39.52i; 13.00, 14.19-33.50i;...
                    14.00, 14.43-27.21i; 15.00, 16.01-21.87i; 16.00, 16.82-17.73i; 17.00, 17.81-13.12i;...
                    18.00, 18.65-8.77i; 19.00, 20.90-4.20i; 20.00, 23.60-1.48i; 21.00, 25.54+0.61i; 22.00, 27.02+1.89i;...
                    23.00, 26.97+2.95i; 24.00, 25.99+5.39i; 25.00, 25.24+9.19i];
                if ~isfield(Trans,'maxHighVoltage'), Trans.maxHighVoltage = 25; end % data sheet lists 30 Volt limit


            case 'L22-14v LF'
                %% L22-14v LF Verasonics Probe
                if ~isfield(Trans,'frequency'), Trans.frequency = 15.625; end % nominal frequency in MHz
                % Note: use 15.625 MHz for 4X sampling (62.5 MHz sample rate),
                % or 18.75 MHz for 4/3 sampling (25.0 MHz sample rate, 50 MHz A/D rate)
                % Manufacturer specified center frequency is 18.0 MHz +/- 10%
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [14, 22]; end
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 128;
                Trans.ConnectorES = (1:Trans.numelements)'; % 1:1 mapping from transducer elements to element signals at connector
                    Trans.elevationApertureMm = 3; % Long Focus active elevation aperture in mm
                    Trans.elevationFocusMm = 20; % Long Focus nominal elevation focus depth from lens on face of transducer
                Trans.elementWidth = 0.08; % element width in mm; assumes 20 micron kerf
                Trans.spacingMm = 0.100;   % Spacing between elements in mm.
                Trans.ElementPos = zeros(Trans.numelements,5);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                % Lens Correction: From mfr data sheet, matching layers are 0.05 mm thick at 2145 m/sec
                % average velocity, and lens is 0.48 mm thick at 1147 m/sec.
                % Thus the net effective lens thickness in mm is given by the following
                % expression, which evaluates to 0.6804 mm for 1540 m/sec velocity
                Trans.lensCorrection = 1000 * speedOfSound * (0.05/2145 + 0.48/1147); % velocities in m/sec; result in mm
                Trans.impedance = [10.00, 11.16-54.28i; 11.00, 12.02-44.73i; 12.00, 14.38-39.52i; 13.00, 14.19-33.50i;...
                    14.00, 14.43-27.21i; 15.00, 16.01-21.87i; 16.00, 16.82-17.73i; 17.00, 17.81-13.12i;...
                    18.00, 18.65-8.77i; 19.00, 20.90-4.20i; 20.00, 23.60-1.48i; 21.00, 25.54+0.61i; 22.00, 27.02+1.89i;...
                    23.00, 26.97+2.95i; 24.00, 25.99+5.39i; 25.00, 25.24+9.19i];
                if ~isfield(Trans,'maxHighVoltage'), Trans.maxHighVoltage = 25; end % data sheet lists 30 Volt limit


            case 'L22-14vX'
                %% L22-14vX Verasonics Probe
                if ~isfield(Trans,'frequency'), Trans.frequency = 15.625; end % nominal frequency in MHz
                % Note: use 15.625 MHz for 4X sampling (62.5 MHz sample rate),
                % or 18.75 MHz for 4/3 sampling (25.0 MHz sample rate, 50 MHz A/D rate)
                % Manufacturer specified center frequency is 18.0 MHz +/- 10%
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [14, 22]; end
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 128;
                Trans.ConnectorES = (1:Trans.numelements)'; % 1:1 mapping from transducer elements to element signals at connector
                    Trans.elevationApertureMm = 1.5; % active elevation aperture in mm
                    Trans.elevationFocusMm = 8; % nominal elevation focus depth from lens on face of transducer
                Trans.elementWidth = 0.08; % element width in mm; assumes 20 micron kerf
                Trans.spacingMm = 0.100;   % Spacing between elements in mm.
                Trans.ElementPos = zeros(Trans.numelements,5);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                % Lens Correction: From mfr data sheet, matching layers are 0.05 mm thick at 2145 m/sec
                % average velocity, and lens is 0.48 mm thick at 1147 m/sec.
                % Thus the net effective lens thickness in mm is given by the following
                % expression, which evaluates to 0.6804 mm for 1540 m/sec velocity
                Trans.lensCorrection = 1000 * speedOfSound * (0.05/2145 + 0.48/1147); % velocities in m/sec; result in mm
                Trans.impedance = [10.00, 11.16-54.28i; 11.00, 12.02-44.73i; 12.00, 14.38-39.52i; 13.00, 14.19-33.50i;...
                    14.00, 14.43-27.21i; 15.00, 16.01-21.87i; 16.00, 16.82-17.73i; 17.00, 17.81-13.12i;...
                    18.00, 18.65-8.77i; 19.00, 20.90-4.20i; 20.00, 23.60-1.48i; 21.00, 25.54+0.61i; 22.00, 27.02+1.89i;...
                    23.00, 26.97+2.95i; 24.00, 25.99+5.39i; 25.00, 25.24+9.19i];
                Trans.elBias = -20; % element DC bias voltage in Volts; see VTS-1351
                if ~isfield(Trans,'maxHighVoltage'), Trans.maxHighVoltage = 30; end % probe data sheet allows 30 V. max with no bias


            case 'L22-14vX LF'
                %% L22-14vX LF Verasonics Probe
                if ~isfield(Trans,'frequency'), Trans.frequency = 15.625; end % nominal frequency in MHz
                % Note: use 15.625 MHz for 4X sampling (62.5 MHz sample rate),
                % or 18.75 MHz for 4/3 sampling (25.0 MHz sample rate, 50 MHz A/D rate)
                % Manufacturer specified center frequency is 18.0 MHz +/- 10%
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [14, 22]; end
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 128;
                Trans.ConnectorES = (1:Trans.numelements)'; % 1:1 mapping from transducer elements to element signals at connector
                    Trans.elevationApertureMm = 3; % Long Focus active elevation aperture in mm
                    Trans.elevationFocusMm = 20; % Long Focus nominal elevation focus depth from lens on face of transducer
                Trans.elementWidth = 0.08; % element width in mm; assumes 20 micron kerf
                Trans.spacingMm = 0.100;   % Spacing between elements in mm.
                Trans.ElementPos = zeros(Trans.numelements,5);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                % Lens Correction: From mfr data sheet, matching layers are 0.05 mm thick at 2145 m/sec
                % average velocity, and lens is 0.48 mm thick at 1147 m/sec.
                % Thus the net effective lens thickness in mm is given by the following
                % expression, which evaluates to 0.6804 mm for 1540 m/sec velocity
                Trans.lensCorrection = 1000 * speedOfSound * (0.05/2145 + 0.48/1147); % velocities in m/sec; result in mm
                Trans.impedance = [10.00, 11.16-54.28i; 11.00, 12.02-44.73i; 12.00, 14.38-39.52i; 13.00, 14.19-33.50i;...
                    14.00, 14.43-27.21i; 15.00, 16.01-21.87i; 16.00, 16.82-17.73i; 17.00, 17.81-13.12i;...
                    18.00, 18.65-8.77i; 19.00, 20.90-4.20i; 20.00, 23.60-1.48i; 21.00, 25.54+0.61i; 22.00, 27.02+1.89i;...
                    23.00, 26.97+2.95i; 24.00, 25.99+5.39i; 25.00, 25.24+9.19i];
                Trans.elBias = -20; % element DC bias voltage in Volts; see VTS-1351
                if ~isfield(Trans,'maxHighVoltage'), Trans.maxHighVoltage = 30; end % probe data sheet allows 30 V. max with no bias


            case 'L35-16vX'
                %% L35-16vX Verasonics Probe
                if ~isfield(Trans,'frequency'), Trans.frequency = 28; end % nominal frequency in MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [20, 36]; end
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 128;
                Trans.ConnectorES = (1:Trans.numelements)'; % 1:1 mapping from transducer elements to element signals at connector
%                     Trans.elevationApertureMm = TBD; % active elevation aperture in mm
%                     Trans.elevationFocusMm = TBD; % nominal elevation focus depth from lens on face of transducer
                Trans.spacingMm = 0.069;   % Spacing between elements in mm.
                Trans.elementWidth = 0.8*Trans.spacingMm; % element width in mm
                Trans.ElementPos = zeros(Trans.numelements,5);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                % Lens Correction: TBD
                Trans.lensCorrection = 0.1; % velocities in m/sec; result in mm
                Trans.impedance = 50; % default value; actual impedance has not yet been measured
                Trans.elBias = -35; % element DC bias voltage in Volts
                if ~isfield(Trans,'maxHighVoltage'), Trans.maxHighVoltage = 40; end % data sheet lists 30 Volt limit


            case 'L38-22v'
                %% L38-22v Kolo CMUT probe
                if ~isfield(Trans,'frequency'), Trans.frequency = 30; end % nominal frequency in MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [22 38]; end
                Trans.type = 0;     % linear=0   Array geometry is linear (x values only).
                Trans.connType = 1; % HDI=1
                Trans.numelements = 256;
                Trans.elementWidth = 0.065; % width in mm
                Trans.spacingMm = 0.069;   % Spacing between elements in mm.
                    Trans.elevationApertureMm = 0; % active elevation aperture in mm (estimate)
                    Trans.elevationFocusMm = 8; % nominal elevation focus depth from lens on face of transducer (estimate)
                Trans.ElementPos = zeros(Trans.numelements,5);
                %   Set default element x positions (units in mm).
                Trans.ElementPos(1:256,1) = Trans.spacingMm*(-((256-1)/2):((256-1)/2));
                Trans.lensCorrection = 0.15; % in mm units
                Trans.impedance = 5;  % artificially low to prevent overdriving the array
                if ~isfield(Trans,'maxHighVoltage'), Trans.maxHighVoltage = 25; end % do not exceed 35 V!
                Trans.HVMux = struct('highVoltageRails', 90, ...
                                    'logicRail', 6.5, ... %10.5, ...
                                    'clock', 8, ...
                                    'clockInvert', 0, ...
                                    'polarity', 0, ...
                                    'latchInvert', 0, ...
                                    'type', 'perEL', ...
                                    'utaSF', 0, ...
                                    'settlingTime', 4, ... % nominal value 4 usec
                                    'ApertureES', zeros(Trans.numelements,129));
                % Add MuxMap table and length; these are used to create the
                % SHIAperture table for actually programming the HVMux
                % switches:
                [Trans.HVMux.MuxMap, Trans.HVMux.SHIAperLgth] = getMuxMap(Trans.name);
                for i = 0:128, Trans.HVMux.ApertureES(:,i+1) = [zeros(1,i),mod((i:i+127),128)+1,zeros(1,128-i)]'; end
                % create the HVMux.SHIAperture array to match
                % HVMux.ApertureES
                Trans = computeHvMuxSHIAperture(Trans);
                % define mapping through the HVMux switches from transducer
                % elements to element signals at connector:
                Trans.ConnectorES = [1:128 1:128]';

            case 'L39-21gD'
                %% L39-21gD Verasonics 30 MHz Linear array
                % VTS-2152 specifications provided by Daxsonics
                if ~isfield(Trans,'frequency'), Trans.frequency = 30; end % nominal frequency in MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [20 36]; end
                Trans.type = 0;     % linear=0   Array geometry is linear (x values only).
                Trans.connType = 7; % UTA 408-GE connector format
                Trans.numelements = 128;
                Trans.elementWidth = 0.0440; % width in mm, 80% of spacing
                Trans.spacingMm = 0.0550;   % Spacing between elements in mm.
                Trans.elevationApertureMm = 1.5; % active elevation aperture in mm
                Trans.elevationFocusMm = 6.5; % nominal elevation focus depth from lens on face of transducer
                Trans.ElementPos = zeros(Trans.numelements,5);
                %   Set element x positions (units in mm).
                Trans.ElementPos(1:128,1) = Trans.spacingMm*(-((128-1)/2):((128-1)/2));
                Trans.lensCorrection = 0.7480; % in mm units
                Trans.impedance = 100; % estimated value, to be Reviewed
                if ~isfield(Trans,'maxHighVoltage'), Trans.maxHighVoltage = 30; end  % max limit 40 from Daxsonics
                Trans.ConnectorES = [1  64   2  63   3  62   4  61   5  60   6  59   7  58   8  57 ...
                                     9  56  10  55  11  54  12  53  13  52  14  51  15  50  16  49 ...
                                    96  97  95  98  94  99  93 100  92 101  91 102  90 103  89 104 ...
                                    88 105  87 106  86 107  85 108  84 109  83 110  82 111  81 112 ...
                                    80 113  79 114  78 115  77 116  76 117  75 118  74 119  73 120 ...
                                    72 121  71 122  70 123  69 124  68 125  67 126  66 127  65 128 ...
                                    17  48  18  47  19  46  20  45  21  44  22  43  23  42  24  41 ...
                                    25  40  26  39  27  38  28  37  29  36  30  35  31  34  32  33]';

            case 'GE9LD'
                %% GE9LD GE probe
                % GE9LD Specification
                % 5.2083 MHz is closest supported 4X sampling frequency to
                % 5.3 MHz nominal value from GE
                if ~isfield(Trans,'frequency'), Trans.frequency = 5.2083; end % (results in 20.83 MHz sample rate)
                if ~isfield(Trans,'Bandwidth'), BW = 0.75; Trans.Bandwidth = 5.3*[1-BW/2, 1+BW/2]; end  % 75% relative bandwidth
                Trans.type = 0;     %linear array
                Trans.connType = 7; % =7 GE connector
                Trans.numelements = 192;
                Trans.spacingMm = .230; % spacing in mm.
                Trans.elementWidth = 0.9 * Trans.spacingMm; % width in mm (guess)
                   Trans.elevationApertureMm = 6.0; % active elevation aperture in mm
                   Trans.elevationFocusMm = 40; % nominal elevation focus depth from lens on face of transducer (estimate)
                % Set default element positions (units in mm).
                Trans.ElementPos = zeros(Trans.numelements,5);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                Trans.lensCorrection = 1; % in mm units; (guess)
                Trans.impedance = 50;    % value is set low - needs to be measured before using in profile 5
                Trans.maxHighVoltage = 50; % conservative value, max. voltage not known

                % Now add the Trans.ConnectorES mapping from the 192 elements
                % to the element signals at the connector.  For any probe
                % with more than 128 elements, the Vantage system uses the
                % GE pinout for a 256 element probe at the connector.  For
                % a 192 element probe GE centers the signals in the middle
                % of the 256, i.e. skips signals 1:32 and 225:256 so we
                % have to reflect that here:
                Trans.ConnectorES = (33:224)';

            case 'GEL3-12D'
                %% GEL3-12D GE probe
                if ~isfield(Trans,'frequency'), Trans.frequency = 6.25; end % 6.5 nominal frequency in MHz from GE data sheet (results in 50 MHz ADCrate)
                if ~isfield(Trans,'Bandwidth'), BW = 0.85; Trans.Bandwidth = 6.5*[1-BW/2, 1+BW/2]; end  % 85% relative bandwidth from GE data sheet
                Trans.type = 0;     % =0 Array geometry is linear (x values only).
                Trans.connType = 7; % GE 408 pin connector on UTA (GE "D" family of probes)
                Trans.numelements = 256;
                % define mapping through the HVMux switches from transducer
                % elements to element signals at connector:
                Trans.ConnectorES = [1:128 1:128]';
                Trans.spacingMm = 0.200;   % element pitch in mm from GE data sheet
                Trans.elementWidth = 0.9 * Trans.spacingMm; % width in mm (wild guess of 10% kerf; not from GE data sheet)
                    Trans.elevationApertureMm = 5; % active elevation aperture in mm (spec)
                    Trans.elevationFocusMm = 22; % nominal elevation focus depth from lens on face of transducer (spec)
                Trans.ElementPos = zeros(Trans.numelements,5);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                Trans.lensCorrection = 1; % in mm units  (empirically validated)
                Trans.impedance = 50; % totally artificial made-up value not from GE data sheet
                if ~isfield(Trans,'maxHighVoltage'), Trans.maxHighVoltage = 50; end

                % Add HVMux control structure; Voltage rail settings are based on notes
                % from Marc for the GE L3-12-D.  The clockInvert, polarity, and
                % latchInvert signals are not used by the GE L3-12-D so the values
                % assigned here are meaningless placeholders that will be ignored by
                % the UTA baseboard
                Trans.HVMux = struct('highVoltageRails', 100, ...
                                     'logicRail', 3.8, ...
                                     'clock', 10, ...
                                     'clockInvert', 0, ...
                                     'polarity', 0, ...
                                     'latchInvert', 0, ...
                                     'type', 'preset', ...
                                     'utaSF', 0, ...
                                     'settlingTime', 4, ... % nominal value 4 usec
                                     'ApertureES', zeros(Trans.numelements,129), ...
                                     'SHIAperture', zeros(2,129), ...
                                     'SHIAperType', 1);

                % Define Trans.HVMux.ApertureES for each of the 129 HVMux
                % apertures, to show the HVMux mapping of the selected
                % elements for that aperture to the 128 element signals at
                % the connector based on a 1:1 mapping to connector element
                % signals for the first 128 elements.  Note the
                % UTA.TransConnector mapping of connector element signals
                % to system channels will be added by VSX, to create the
                % overall element to channel mapping in
                % Trans.HVMux.Aperture
                for i = 0:128
                    Trans.HVMux.ApertureES(:,i+1) = [zeros(1,i), mod((i:i+127),128)+1, zeros(1,128-i)]';
                end

                % Defining the SHIAperture array:  The GE probe requires a 16 bit
                % control word to be sent to a CPLD in the probe, to tell it which
                % aperture to select through the HVMux chips.  Thus each
                % SHIAperture column will contain that 16 bit control word, with
                % the LS 8 bits in byte 2 and MS 8 bits in byte 1.  From the GE
                % documentation, the bit assignments are as follows in the control
                % word:
                % bits 15:12 Command, set to 0000 for legacy non-pipeline mode (9
                % usec load time)
                % bit ll Echo, set to 0 to disable echo back of the command received
                % bits 10:8 Row select (elevation aperture control) zero disables
                % all rows, 001 enables only the center row so we will assume that
                % is the value to use for L3-12 since it has only one elevation
                % row.
                % bits 7:0 Aperture select:  must be one less than the aperture
                % value we use in the script, since GE counts elements from zero
                % while we start at one.
                pldCmd = 0; % command value for bits 15:12
                pldEcho = 0; % value for bit 11
                pldRow = 1; % command valiue for bits 10:8
                % Set first byte to the 8 MSbits of the control word, which are
                % the same for all apertures since the only thing that changes is
                % the aperture select in bits 0:7
                MSbyte = 16*pldCmd + 8*pldEcho + pldRow;
                Trans.HVMux.SHIAperture(1, :) = MSbyte*ones(1, 129);
                % Set second byte to the aperture select value, one less than the
                % aperture index value we use which is in the range 1:129
                Trans.HVMux.SHIAperture(2, :) = 0:128;

            case 'CL10-5'
                %% CL10-5 HDI Probe
                if ~isfield(Trans,'frequency'), Trans.frequency = 6.25; end % nominal frequency in MHz (this is the closest allowed frequency)
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [4.5, 9]; end
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 128;
                Trans.ConnectorES = (1:Trans.numelements)'; % 1:1 mapping from transducer elements to element signals at connector
                Trans.spacingMm = 0.200;   % Spacing between elements in mm.
                Trans.elementWidth = Trans.spacingMm - 0.05; % width in mm
                    % Trans.elevationApertureMm = 'unknown; % active elevation aperture in mm (estimate)
                    % Trans.elevationFocusMm = 'unknown; % nominal elevation focus depth from lens on face of transducer (estimate)
                Trans.ElementPos = zeros(Trans.numelements,5);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                Trans.lensCorrection = 2.1; % in mm units; was 9 wavelengths;
                Trans.impedance = 50; % using default value


            case 'CL15-7'
                %% CL15-7 HDI Probe
                if ~isfield(Trans,'frequency'), Trans.frequency = 8.929; end % nominal frequency in MHz
                % Vantage:  8.929 is closest supported frequency to 9.0 MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = 9*[0.7, 1.3]; end % default assumed value of 60% of center frequency
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 128;
                Trans.ConnectorES = (1:Trans.numelements)'; % 1:1 mapping from transducer elements to element signals at connector
                Trans.elementWidth = 0.16; % width in mm
                Trans.spacingMm = 0.178;   % Spacing between elements in mm.
                    Trans.elevationApertureMm = 8; % active elevation aperture in mm (estimate)
                    Trans.elevationFocusMm = 15; % nominal elevation focus depth from lens on face of transducer (estimate)
                Trans.ElementPos = zeros(Trans.numelements,5);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                Trans.lensCorrection = 0.6899; % in mm units; was 4 wavelengths;
                Trans.impedance = 50; % using default value


            case 'C4-2'
                %% C4-2 HDI Probe
                if ~isfield(Trans,'frequency'), Trans.frequency = 2.976; end % nominal frequency in MHz
                % Vantage:  2.976 is closest supported frequency to 3.0 MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [2, 4]; end
                Trans.type = 1;     % Array geometry is curved linear (x and z values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 128;
                Trans.ConnectorES = (1:Trans.numelements)'; % 1:1 mapping from transducer elements to element signals at connector
                scanangle = 74.95 * (pi/180);    % degrees converted to radians
                radiusMm = 41.219;  % radius in mm.
                spacingMm = radiusMm * scanangle/(Trans.numelements-1); % spacing in mm.
                kerf = .050;   % guess (in mm)
                Trans.radiusMm = radiusMm; % radius in mm.
                Trans.spacingMm = spacingMm;  % Spacing in mm.
                Trans.elementWidth = (spacingMm - kerf);  % width in mm
                    % Trans.elevationApertureMm = 'unknown; % active elevation aperture in mm (estimate)
                    % Trans.elevationFocusMm = 'unknown; % nominal elevation focus depth from lens on face of transducer (estimate)
                firstangle = -(scanangle/2); %   first element angle = -0.65405 radians
                deltatheta = scanangle/(Trans.numelements-1);
                %   Set default element positions (units in mm).
                Trans.ElementPos = zeros(Trans.numelements,5);
                Angle = firstangle:deltatheta:-firstangle;
                Trans.ElementPos(:,1) = Trans.radiusMm*sin(Angle);
                Trans.ElementPos(:,2) = 0;
                Trans.ElementPos(:,3) = Trans.radiusMm*cos(Angle)-Trans.radiusMm;
                Trans.ElementPos(:,4) = Angle; % Orientation of element with respect to z axis.
                Trans.impedance = 50; % using default value

            case 'C5-2'
                %% C5-2 HDI Probe
                if ~isfield(Trans,'frequency'), Trans.frequency = 3.125; end % nominal frequency in MHz
                % Vantage:  3.125 is closest supported frequency to 3.2143 MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [2, 4.5]; end
                Trans.type = 1;     % Array geometry is curved linear (x and z values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 128;
                Trans.ConnectorES = (1:Trans.numelements)'; % 1:1 mapping from transducer elements to element signals at connector
                scanangle = 74.95 * (pi/180);    % degrees converted to radians
                radiusMm = 41.219;  % radius in mm.
                spacingMm = radiusMm * scanangle/(Trans.numelements-1); % spacing in mm.
                kerf = .050;   % guess (in mm)
                Trans.radiusMm = radiusMm; % radius in mm.
                Trans.spacingMm = spacingMm;  % Spacing in mm.
                Trans.elementWidth = (spacingMm - kerf);  % width in mm
                    Trans.elevationApertureMm = 0; % active elevation aperture in mm (unknown)
                    Trans.elevationFocusMm = 60; % nominal elevation focus depth from lens on face of transducer (estimate)
                firstangle = -(scanangle/2); %   first element angle = -0.65405 radians
                deltatheta = scanangle/(Trans.numelements-1);
                %   Set default element positions (units in mm).
                Trans.ElementPos = zeros(Trans.numelements,5);
                Angle = firstangle:deltatheta:-firstangle;
                Trans.ElementPos(:,1) = Trans.radiusMm*sin(Angle);
                Trans.ElementPos(:,2) = 0;
                Trans.ElementPos(:,3) = Trans.radiusMm*cos(Angle)-Trans.radiusMm;
                Trans.ElementPos(:,4) = Angle; % Orientation of element with respect to z axis.
                Trans.impedance = [ 1 19.4-255i;  1.25 20.1-186i;  1.5 22.4-136i;  1.75 29-94.8i;  2 43.6-62i;  2.25 64.8-48.9i;...
                    2.5 70.2-54.9i;  2.75 54.3-49.1i;  3 45.4-30.9i;  3.25 42.9-13.3i;  3.5 42.5+1.9i;  3.75 42.1+16.4i;...
                    4 41.5+31.1i;  4.25 41.2+46.3i;  4.5 42.3+62.7i;  4.75 44.8+78.7i;  5 45.2+93i;  5.25 40.9+114i;...
                    5.5 42.6+144i;  5.75 51.2+177i;  6 66.5+214i];

            case 'C5-2v'
                %% C5-2v Verasonics Probe
                if ~isfield(Trans,'frequency'), Trans.frequency = 3.57; end % Center frequency in MHz from mfr's data sheet
                % Fractional bandwidth (-6 dB round trip) is 79.6% from mfr's data sheet
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [2.15, 4.99]; end
                % Bandwidth calculated from mfr's Center Frequency and fractional bandwidth
                Trans.type = 1;     % Array geometry is curved linear (x and z values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 128;
                Trans.ConnectorES = (1:Trans.numelements)'; % 1:1 mapping from transducer elements to element signals at connector
                radiusMm = 49.57;  % radius in mm.
                spacingMm = .508; % spacing in mm.
                kerf = .048;   % in mm.
                scanangle = 128*spacingMm/radiusMm;    % arc length/radius
                Trans.radiusMm = radiusMm; % radius in mm.
                Trans.spacingMm = spacingMm;  % Spacing in mm.
                Trans.elementWidth = (spacingMm - kerf);  % width in mm
                    Trans.elevationApertureMm = 13.5; % active elevation aperture in mm (estimate)
                    Trans.elevationFocusMm = 60; % nominal elevation focus depth from lens on face of transducer (spec)
                deltatheta = spacingMm/radiusMm;
                firstangle = -(scanangle/2) + 0.5*deltatheta; % first element angle
                %   Set default element positions (units in mm).
                Trans.ElementPos = zeros(Trans.numelements,5);
                Angle = firstangle:deltatheta:-firstangle;
                Trans.ElementPos(:,1) = Trans.radiusMm*sin(Angle);
                Trans.ElementPos(:,2) = 0;
                Trans.ElementPos(:,3) = Trans.radiusMm*cos(Angle)-Trans.radiusMm;
                Trans.ElementPos(:,4) = Angle; % Orientation of element with respect to z axis.
                Trans.lensCorrection = 1.035; % in mm units; was 3 wavelengths;
                Trans.impedance = [ 1 39.2-216i;  1.25 41-162i;  1.5 41.1-118i;  1.75 47.5-87i;  2 50.8-69.7i;  2.25 44.7-49.5i;...
                    2.5 43.3-27i;  2.75 48.6-8.53i;  3 50.4+4.4i;  3.25 50.7+20.7i;  3.5 56.5+38.2i;  3.75 69.1+51.8i;...
                    4 83.2+53.8i;  4.25 91.9+45.7i;  4.5 81.7+34.6i;  4.75 66+40.1i;  5 56.6+47i;  5.25 42.4+57.2i;...
                    5.5 30.7+75i;  5.75 25.3+95.5i;  6 24.3+115i;  6.25 24.8+132i;  6.5 26+148i;  6.75 27.8+164i;  7 30+179i];

            case 'C5-2gH'
                %% C5-2gH Verasonics 3 MHz Curved Array from Humanscan
                if ~isfield(Trans,'frequency'), Trans.frequency = 3.3; end % Center frequency in MHz from mfr's data sheet
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [2.15, 4.46]; end % 70% fractional bandwidth (average at 6dB, taken from mfr's data sheet)
                % Bandwidth calculated from mfr's Center Frequency and fractional bandwidth
                Trans.type = 1;     % Array geometry is curved linear (x and z values only).
                Trans.connType = 7; % UTA 408-GE connector
                Trans.numelements = 128;
                Trans.radiusMm = 60;  % radius in mm.
                Trans.spacingMm = .475; % spacing in mm.
                kerf =  .1*Trans.spacingMm;   % in mm. (guess)
                scanangle = 128*Trans.spacingMm/Trans.radiusMm;    % arc length/radius
                Trans.elementWidth = (Trans.spacingMm - kerf);  % width in mm
                    Trans.elevationApertureMm = 13; % active elevation aperture in mm (spec)
                    Trans.elevationFocusMm = 60; % nominal elevation focus depth from lens on face of transducer (spec)
                deltatheta = Trans.spacingMm/Trans.radiusMm;
                firstangle = -(scanangle/2) + 0.5*deltatheta; % first element angle
                %   Set default element positions (units in mm).
                Trans.ElementPos = zeros(Trans.numelements,5);
                Angle = firstangle:deltatheta:-firstangle;
                Trans.ElementPos(:,1) = Trans.radiusMm*sin(Angle);
                Trans.ElementPos(:,2) = 0;
                Trans.ElementPos(:,3) = Trans.radiusMm*cos(Angle)-Trans.radiusMm;
                Trans.ElementPos(:,4) = Angle; % Orientation of element with respect to z axis.
                Trans.lensCorrection = computeLensCorrection(speedOfSound, 0.327, 2055, 0.74, 1050);
                Trans.ConnectorES = (1:Trans.numelements)'; % 1:1 mapping from transducer elements to element signals at connector
                Trans.impedance = 50;  % Needs to be measured before using in profile 5
                Trans.maxHighVoltage = 70;
                
           case 'GEC1-6D'
                %% GEC1-6D GE probe
                if ~isfield(Trans,'frequency'), Trans.frequency = 3.9; end % nominal frequency in MHz is 3.4
                if ~isfield(Trans,'Bandwidth'), BW = 0.95; Trans.Bandwidth = 3.4*[1-BW/2, 1+BW/2]; end
                Trans.type = 1;     % Array geometry is curved linear (x and z values only).
                Trans.connType = 7; % =7 GE connector
                Trans.numelements = 192;
                radiusMm = 56.8;  % radius in mm.
                spacingMm = .350; % spacing in mm.
                kerf = .9*spacingMm;   % in mm. (guess)
                scanangle = Trans.numelements*spacingMm/radiusMm;    % arc length/radius
                Trans.radiusMm = radiusMm; % radius in mm.
                Trans.spacingMm = spacingMm;  % Spacing in mm.
                Trans.elementWidth = (spacingMm - kerf);  % width in mm
                    Trans.elevationApertureMm = 11.5; % active elevation aperture in mm (estimate)
                    Trans.elevationFocusMm = 66; % nominal elevation focus depth from lens on face of transducer (estimate)
                deltatheta = spacingMm/radiusMm;
                firstangle = -(scanangle/2) + 0.5*deltatheta; % first element angle
                %   Set default element positions (units in mm).
                Trans.ElementPos = zeros(Trans.numelements,5);
                Angle = firstangle:deltatheta:-firstangle;
                Trans.ElementPos(:,1) = Trans.radiusMm*sin(Angle);
                Trans.ElementPos(:,2) = 0;
                Trans.ElementPos(:,3) = Trans.radiusMm*cos(Angle)-Trans.radiusMm;
                Trans.ElementPos(:,4) = Angle; % Orientation of element with respect to z axis.
                Trans.lensCorrection = 1.2; % in mm units;
                Trans.impedance = 50;    % Needs to be measured before using in profile 5
                if ~isfield(Trans,'maxHighVoltage'), Trans.maxHighVoltage = 50; end
                % Now add the Trans.ConnectorES mapping from the 192 elements
                % to the element signals at the connector.  For any probe
                % with more than 128 elements, the Vantage system uses the
                % GE pinout for a 256 element probe at the connector.  For
                % a 192 element probe GE centers the signals in the middle
                % of the 256, i.e. skips signals 1:32 and 225:256 so we
                % have to reflect that here:
                Trans.ConnectorES = (33:224)';

           case 'GE4CD'
                %% GE4CD GE probe
                if ~isfield(Trans,'frequency'), Trans.frequency = 3.125; end % nominal frequency = 3.2 MHz
                if ~isfield(Trans,'Bandwidth'), BW = 0.7; Trans.Bandwidth = 3.2*[1-BW/2, 1+BW/2]; end
                Trans.type = 1;     % 1= Array geometry is curved linear (x and z values only).
                Trans.connType = 7; % =7 GE connector
                Trans.numelements = 128;
                radiusMm = 60.;  % radius in mm.
                spacingMm = .478; % spacing in mm.
                kerf = 0.9*spacingMm;   % in mm. (guess)
                scanangle = 128*spacingMm/radiusMm;    % arc length/radius
                Trans.radiusMm = radiusMm; % radius in mm.
                Trans.spacingMm = spacingMm;  % Spacing in mm.
                Trans.elementWidth = (spacingMm - kerf);  % width in mm
                    Trans.elevationApertureMm = 13; % active elevation aperture in mm (spec)
                    Trans.elevationFocusMm = 67; % nominal elevation focus depth from lens on face of transducer (spec)
                deltatheta = spacingMm/radiusMm;
                firstangle = -(scanangle/2) + 0.5*deltatheta; % first element angle
                %   Set default element positions (units in mm).
                Trans.ElementPos = zeros(Trans.numelements,5);
                Angle = firstangle:deltatheta:-firstangle;
                Trans.ElementPos(:,1) = Trans.radiusMm*sin(Angle);
                Trans.ElementPos(:,2) = 0;
                Trans.ElementPos(:,3) = Trans.radiusMm*cos(Angle)-Trans.radiusMm;
                Trans.ElementPos(:,4) = Angle; % Orientation of element with respect to z axis.
                Trans.lensCorrection = 2.0; % in mm units;
                Trans.impedance = 50;  % Needs to be measured before using in profile 5
                if ~isfield(Trans,'maxHighVoltage'), Trans.maxHighVoltage = 50; end
                % Trans.ConnectorES mapping from probe elements to element
                % signals at the connector: Vantage system supports the GE
                % signal assignments for 128 element probes so for this
                % probe the mapping is one-to-one:
                Trans.ConnectorES = (1:128)';

            case 'C9-5ICT'
                %% C9-5ICT HDI Probe
                if ~isfield(Trans,'frequency'), Trans.frequency = 7.813; end % nominal frequency in MHz
                % Vantage:  7.813 and 6.944 are closest supported frequencies to 7.5 MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = 7.5*[0.7, 1.3]; end % default assumed value of 60% of center frequency
                Trans.type = 1;     % Array geometry is curved linear (x and z values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 128;
                scanangle = 146.677 * (pi/180);    % degrees converted to radians
                radiusMm = 8.511;  % radius in mm.
                spacingMm = radiusMm * scanangle/(Trans.numelements-1); % spacing in mm.
                kerf = .025;   % guess (in mm)
                    Trans.elevationApertureMm = 0; % active elevation aperture in mm (unknown)
                    Trans.elevationFocusMm = 0; % nominal elevation focus depth from lens on face of transducer (unknown)
                Trans.radiusMm = radiusMm; % radius in mm.
                Trans.spacingMm = spacingMm;  % Spacing in mm.
                Trans.elementWidth = (spacingMm - kerf);  % width in mm
                firstangle = -(scanangle/2); %   first element angle = -0.65405 radians
                deltatheta = scanangle/(Trans.numelements-1);
                %   Set default element positions (units in mm).
                Trans.ElementPos = zeros(Trans.numelements,5);
                Angle = firstangle:deltatheta:-firstangle;
                Trans.ElementPos(:,1) = Trans.radiusMm*sin(Angle);
                Trans.ElementPos(:,2) = 0;
                Trans.ElementPos(:,3) = Trans.radiusMm*cos(Angle)-Trans.radiusMm;
                Trans.ElementPos(:,4) = Angle; % Orientation of element with respect to z axis.
                Trans.impedance = 50; % using default value

            case 'GEIC5-9D'
                %% GEIC5-9D GE probe
                if ~isfield(Trans,'frequency'), Trans.frequency = 5.682; end % nominal frequency = 5.8 MHz
                if ~isfield(Trans,'Bandwidth'), BW = 0.75; Trans.Bandwidth = 5.8*[1-BW/2, 1+BW/2]; end
                Trans.type = 1;     % =1 Array geometry is curved (x,z values only).
                Trans.connType = 7; % =7 GE 408 pin connector on UTA (GE "D" family of probes)
                Trans.numelements = 192;
                radiusMm = 10.1;  % radius in mm.
                spacingMm = .138; % spacing in mm.
                kerf = spacingMm*0.2;   % in mm. (guess)
                    Trans.elevationApertureMm = 6; % active elevation aperture in mm (spec)
                    Trans.elevationFocusMm = 35; % nominal elevation focus depth from lens on face of transducer (spec)
                scanangle = Trans.numelements*spacingMm/radiusMm;    % arc length/radius
                Trans.radiusMm = radiusMm; % radius in mm.
                Trans.spacingMm = spacingMm;  % Spacing in mm.
                Trans.elementWidth = (spacingMm - kerf);  % width in mm
                deltatheta = spacingMm/radiusMm;
                firstangle = -(scanangle/2) + 0.5*deltatheta; % first element angle
                %   Set default element positions (units in mm).
                Trans.ElementPos = zeros(Trans.numelements,5);
                Angle = firstangle:deltatheta:-firstangle;
                Trans.ElementPos(:,1) = Trans.radiusMm*sin(Angle);
                Trans.ElementPos(:,2) = 0;
                Trans.ElementPos(:,3) = Trans.radiusMm*cos(Angle)-Trans.radiusMm;
                Trans.ElementPos(:,4) = Angle; % Orientation of element with respect to z axis.
                Trans.lensCorrection = 0.54; % in mm units
                Trans.impedance = 50; % totally artificial made-up value not from GE data sheet
                if ~isfield(Trans,'maxHighVoltage'), Trans.maxHighVoltage = 50; end
                % Now add the Trans.ConnectorES mapping from the 192 elements
                % to the element signals at the connector.  For any probe
                % with more than 128 elements, the Vantage system uses the
                % GE pinout for a 256 element probe at the connector.  For
                % a 192 element probe GE centers the signals in the middle
                % of the 256, i.e. skips signals 1:32 and 225:256 so we
                % have to reflect that here:
                Trans.ConnectorES = (33:224)';


            case 'GEM5ScD'
                %% GEM5ScD GE probe
                if ~isfield(Trans,'frequency'), Trans.frequency = 2.841; end % nominal frequency is 2.8 MHz from GE data sheet
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [1.7 4.2]; end % from frequency response plot on GE data sheet
                Trans.type = 0;  % specify 1.5D array as linear array with single row of elements.
                % Note this 1.5 D array cannot be steered in the Y direction, since the outer rows
                % of elements are connected together, allowing only the elevation focus to be adjusted.
                % For the Verasonics implementation, there is no dynamic focusing in the elevation
                % dimension - instead a dynamic aperture approach is used, with the near field using
                % acquisitions with the central row only, and the far field using all three rows.
                Trans.connType = 7; % GE 408 pin connector on UTA (GE "D" family of probes)
                Trans.numelements = 160;
                % Note the actual transducer element array is 80 X 3 for a total of 240 active
                % elements, but the pairs of side elements are connected together within the probe
                % so only 160 element signals are present at the probe connector.
                Trans.spacingMm = 0.270;   % element pitch in mm; total active aperture is 21.6 mm (80 * 0.270)
                Trans.elementWidth = 0.230; % width in mm from GE data sheet; kerf is 0.04 mm
                Trans.elevationApertureMm = 13; % active elevation aperture in mm from GE data sheet
                elevOffsetMm = 4.875; % distance in mm from center of the array to center of the side elements
                % total elevation aperture is 13 mm with a 6.5 mm center element and 3.25 mm side
                % element on each side, connected in parallel (and thus same total area; approx.
                % same impedance).  Center of a side element is therefore 4.875 mm from center of
                % the elevation aperture
                Trans.elevationFocusMm = 77; % nominal elevation focus depth from lens on face of
                % transducer, from GE spec sheet. Note that by adjusting the delay and weighting of
                % the side element rows relative to the center row, the elevation focus can be
                % adjusted over some range.
                Trans.ElementPos = zeros(Trans.numelements,5);
                % If treated as a full 2D array, Trans.type would be set to 2 and the ElementPos array
                % would require five entries for each element (X, Y, Z position coordinates plus angular
                % orientation in two dimensions).  When treated as a type 0 linear, all elements are
                % at z = 0 and facing straight ahead, so both angular directions are zero. In addition,
                % elements 81:160, which represent the outer row pairs, are specified to be located
                % at y = 0 (similar to the inner row), since in the far field the 3 element rows are
                % treated as if they were a single element.
                xpos = Trans.spacingMm*(-((Trans.numelements-2)/4):((Trans.numelements-2)/4));
                Trans.ElementPos(:,1) = [xpos, xpos]; % same xpos for both rows
                % If desiring to treat the GEMScD as a full 2D array (type 2), uncomment the lines below
                % to offset the outer row of elements in the y direction.  To provide zero delay focus
                % at elevationFocusMM,the outer rows can be offset in the z direction as well. This
                % allows computeTXDelays to place the elevation focus at the range focal point.
                % Trans.ElementPos(81:160,2) = elevOffsetMm; % add the y-offset for outer row
                % Trans.ElementPos(81:160,3) = Trans.elevationFocusMm - sqrt(elevOffsetMm^2 + Trans.elevationFocusMm^2);
                Trans.lensCorrection = 1; % in mm units  (wild guess; not from GE data sheet)
                Trans.impedance = 50; % totally artificial made-up value not from GE data sheet
                if ~isfield(Trans,'maxHighVoltage'), Trans.maxHighVoltage = 50; end

                % Notes on element numbering and element to channel mapping:  GE has defined element
                % signal numbers at the probe connector for 256 elements, with the convention that
                % probes using less than 256 elements will use the elements in the center of the
                % signal assignments for a 256 element array.  Thus for this 160 element probe the
                % GE element signal numbers at the connector will be 48 to 207; element signals 0:47
                % and 208:255 will not be used.

                % Note also that since Verasonics counts elements from 1 rather than zero, we will
                % be using elements signals 49 to 208 with the Verasonics numbering.

                % For the 1.5 D array of the M5ScD GE has mapped their
                % element numbers 1:16 to the first 16 center elements of
                % the 3-row array, and elements 17:32 to the first 16 outer
                % row element pairs.  This pattern repeats for each of the
                % five groups of 32 elements going across the array.  For
                % convenience in programming the system, we want to number
                % the center row as elements 1:80 and the outer rows as
                % elements 81:160.  Create Trans.ConnectorES to reflect that
                % mapping- the assignments below create the "center row,
                % outer row" mapping for the 160 elements, and then add an
                % offset of 48 to map them into GE's 256 element connector
                % signal assignments.
                Trans.ConnectorES = zeros(160, 1); % create the 160 element column vector
                Trans.ConnectorES(1:80, 1) = 48 + [1:16, 33:48, 65:72, 81:88, 97:112, 129:144]'; % center row of elements
                Trans.ConnectorES(81:160, 1) = 48 + [17:32, 49:64, 73:80, 89:96, 113:128, 145:160]'; % outer rows of elements

            case 'P4-1'
                %% P4-1 HDI Probe
                % The P4-1 is a 96 element array, so to use it with a 128 connector I/O system
                %    we have to define the connectivity in Trans.ConnectorES.
                if ~isfield(Trans,'frequency'), Trans.frequency = 2.5; end % nominal frequency in MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [1.5, 3.5]; end
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 96;
                Trans.ConnectorES = [023 022 021 041 024 042 046 043 045 044 047 018 017 048 013 020 ...
                                   019 014 015 016 049 050 054 051 053 052 009 055 056 011 012 005 ...
                                   006 007 008 010 004 003 002 001 040 039 038 037 033 034 035 036 ...
                                   093 094 095 096 092 091 090 089 128 127 126 125 119 121 122 123 ...
                                   124 117 118 073 074 120 077 076 078 075 079 080 113 114 115 110 ...
                                   109 116 081 112 111 082 085 084 086 083 087 105 088 108 107 106]';
                kerf = .050;   % guess (in mm)
                Trans.spacingMm = 0.2950;   % Spacing between elements in mm.
                Trans.elementWidth = (Trans.spacingMm - kerf); % width in mm
                    Trans.elevationApertureMm = 16; % active elevation aperture in mm (estimate)
                    Trans.elevationFocusMm = 80; % nominal elevation focus depth from lens on face of transducer (estimate)
                Trans.ElementPos = zeros(Trans.numelements,5);
                %   Set default element x positions (units in mm).
                Trans.ElementPos(1:96,1) = Trans.spacingMm*(-((96-1)/2):((96-1)/2));
                Trans.lensCorrection = 2.464; % in mm units; was 4 wavelengths;
                Trans.impedance = [ 1 33.9-390i;  1.25 43.1-286i;  1.5 63.9-217i;  1.75 80.1-185i;  2 74.1-159i;...
                    2.25 73.9-126i;  2.5 84-99.2i;  2.75 98-93.3i;  3 89.5-99.4i;  3.25 63.6-89.5i;  3.5 44.2-63.6i;...
                    3.75 32.7-34.5i;  4 25.8-3.09i;  4.25 24.1+29.5i;  4.5 26.3+62.6i;  4.75 33.2+96.9i;  5 45.2+130i];

            case 'P4-2'
                %% P4-2 HDI Probe
                % The P4-2 is a 64 element array, so to use it with a 128 connector I/O system
                %    we have to define the connectivity in Trans.ConnectorES.  Elements 32-63 are
                %    wired to connector inputs 97-128.
                if ~isfield(Trans,'frequency'), Trans.frequency = 2.5; end % nominal frequency in MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [2, 4]; end
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 64;
                Trans.ConnectorES = [1:32,97:128]';
                Trans.elementWidth = 0.2950; % width in mm
                Trans.spacingMm = 0.3200;   % Spacing between elements in mm.
                    Trans.elevationApertureMm = 0; % active elevation aperture in mm (unknown)
                    Trans.elevationFocusMm = 60; % nominal elevation focus depth from lens on face of transducer (estimate)
                Trans.ElementPos = zeros(Trans.numelements,5);
                %   Set default element x positions (units in mm).
                Trans.ElementPos(1:64,1) = Trans.spacingMm*(-((64-1)/2):((64-1)/2));
                Trans.lensCorrection = 3.08; % in mm units; was 5 wavelengths;
                Trans.impedance = [1.7, 98.0+91.0i; 2.0, 99.0; 2.5, 87.0+18.0i; 3.0, 77.0+16.0i; 3.5, 28.0+4.0i;...
                    4.0, 28.0+35.0i; 4.5, 42.0+83.0i; 5.0, 75.0+125.0i; 5.5, 128.0+159.0i; 6.0, 201.0+178.0i];

            case 'P4-2v'
                %% P4-2v Verasonics Probe
                % The P4-2 is a 64 element array, so to use it with the 128 channel UTA 260 connector
                %    we have to define the connectivity in Trans.ConnectorES.  Elements 1-64 are
                %    wired to connector Element Signals 33-96.
                if ~isfield(Trans,'frequency'), Trans.frequency = 2.72; end % Center frequency in MHz from mfr's data sheet
                % Fractional bandwidth (-6 dB round trip) is 74.2% from mfr's data sheet
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [1.71, 3.73]; end
                % Bandwidth calculated from mfr's Center Frequency and fractional bandwidth
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 64;
                Trans.ConnectorES = (33:96)';
                Trans.elementWidth = 0.250; % width in mm
                Trans.spacingMm = 0.300;   % Spacing between elements in mm.
                    Trans.elevationApertureMm = 14; % active elevation aperture in mm (estimate)
                    Trans.elevationFocusMm = 60; % nominal elevation focus depth from lens on face of transducer (spec)
                Trans.ElementPos = zeros(Trans.numelements,5);
                %   Set default element x positions (units in mm).
                Trans.ElementPos(1:64,1) = Trans.spacingMm*(-((64-1)/2):((64-1)/2));
                Trans.lensCorrection = 2.587; % in mm units; was 5 wavelengths;
                Trans.impedance = [ 1 33.8-366i;  1.25 49.4-252i;  1.5 64.1-196i;  1.75 75.5-140i;  2 99.5-118i;...
                    2.25 84.2-104i;  2.5 73.9-74.6i;  2.75 73.5-36i;  3 94.5-5.42i;  3.25 115-9.99i;  3.5 111-24.8i;...
                    3.75 74.9-14.8i;  4 66.9+17.5i;  4.25 70.2+22.4i;  4.5 37.7+37.8i;  4.75 26+75.3i;  5 24.8+107i;...
                    5.25 26.1+135i;  5.5 28.5+162i;  5.75 31.7+188i;  6 36.3+215i];


            case 'P4-2gH'
                %% P4-2gH Verasonics 3 MHz 96 El PA Probe from Humanscan
                % This probe uses UTA 408-GE connector interface.
                if ~isfield(Trans,'frequency'), Trans.frequency = 3.15; end % Center frequency in MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [1.81, 4.49]; end % 85% fractional bandwidth (average at 6dB, taken from mfr's data sheet)
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 7; % UTA 408-GE connector
                Trans.numelements = 96;
                Trans.ConnectorES = (1:96)';
                Trans.elementWidth = 0.180; % width in mm; estimated at 90% of spacing
                Trans.spacingMm = 0.200;   % Spacing between elements in mm.
                    Trans.elevationApertureMm = 14; % active elevation aperture in mm (estimate)
                    Trans.elevationFocusMm = 80; % nominal elevation focus depth from lens on face of transducer (spec)
                Trans.ElementPos = zeros(Trans.numelements,5);
                %   Set default element x positions (units in mm).
                Trans.ElementPos(1:96,1) = Trans.spacingMm*(-((96-1)/2):((96-1)/2));
                Trans.lensCorrection = computeLensCorrection(speedOfSound, 0.318, 2055, 0.670, 1050);
                Trans.impedance = 50;  % Needs to be measured before using in profile 5
                Trans.maxHighVoltage = 70;
                
                
            case 'P6-3'
                %% P6-3 HDI Probe
                if ~isfield(Trans,'frequency'), Trans.frequency = 4.464; end % nominal frequency in MHz
                % Vantage:  4.464 is closest supported frequency to 4.5 MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [3, 6]; end
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 128;
                Trans.ConnectorES = (1:Trans.numelements)'; % 1:1 mapping from transducer elements to element signals at connector
                Trans.elementWidth = (.218-.025); % width in mm
                Trans.spacingMm = 0.218;   % Spacing between elements in mm
                    % Trans.elevationApertureMm = 'unknown; % active elevation aperture in mm (estimate)
                    % Trans.elevationFocusMm = 'unknown; % nominal elevation focus depth from lens on face of transducer (estimate)
                Trans.ElementPos = zeros(Trans.numelements,5);
                %   Set default element x positions (units in mm).
                Trans.ElementPos(1:128,1) = Trans.spacingMm*(-((128-1)/2):((128-1)/2));
                Trans.lensCorrection = 1.380; % in mm units; was 4 wavelengths;
                Trans.impedance = 48; % Z @ 3.4 Mhz

            case 'P7-4'
                %% P7-4 HDI Probe
                % The P7-4 is a 64 element array, so to use it with a 128 connector I/O system
                %    we have to define the connectivity in Trans.ConnectorES.  Elements 32-64 are
                %    wired to connector inputs 97-128.
                if ~isfield(Trans,'frequency'), Trans.frequency = 5.208; end % nominal frequency in MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [4, 7]; end
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 64;
                Trans.ConnectorES = [1:32,97:128]';
                Trans.elementWidth = 0.1369; % width in mm
                Trans.spacingMm = 0.1711;   % Spacing between elements in mm.
                    % Trans.elevationApertureMm = 'unknown; % active elevation aperture in mm (estimate)
                    % Trans.elevationFocusMm = 'unknown; % nominal elevation focus depth from lens on face of transducer (estimate)
                Trans.ElementPos = zeros(Trans.numelements,5);
                %   Set default element x positions (units in mm).
                Trans.ElementPos(1:64,1) = Trans.spacingMm*(-((64-1)/2):((64-1)/2));
                Trans.lensCorrection = 0.7; % in mm units; was 5 wavelengths;
                Trans.impedance = 50;%[1.7, 98.0+91.0i; 2.0, 99.0; 2.5, 87.0+18.0i; 3.0, 77.0+16.0i; 3.5, 28.0+4.0i;...
%                     4.0, 28.0+35.0i; 4.5, 42.0+83.0i; 5.0, 75.0+125.0i; 5.5, 128.0+159.0i; 6.0, 201.0+178.0i];

            case 'RC6gV'
                %% RC6gV Vermon Row-Column array
                if ~isfield(Trans,'frequency'), Trans.frequency = 6.0; end % Center frequency in MHz from Vermon data sheet
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [3.6, 8.4]; end % 80% fractional bandwidth
                Trans.type = 4;     % Array geometry is 2D Row-Column array.
                Trans.connType = 7; % UTA 408-GE connector
                Trans.numelements = 256;
                Trans.spacingMm = 0.200; % element spacing in mm.
                Trans.ElementPos = zeros(Trans.numelements,5);
                Trans.ElementPos(1:128,1) = Trans.spacingMm*((-127/2):(127/2)); % in Wavelengths
                Trans.ElementPos(129:256,2) = Trans.spacingMm*((-127/2):(127/2));
                Trans.elementWidth = 0.175; % width in mm from Vermon data sheet
                Trans.elementLength = 25.6; % element Length in mm from Vermon data sheet
                Trans.lensCorrection = 0; % No lens
                Trans.impedance = 100;
                Trans.maxHighVoltage = 40;
                % The ConnectorES array below was derived from Brian's
                % original mapping for the first Vermon prototype using GE
                % connector; we have now verified it functions well with
                % both that probe and the second one from Vermon after
                % changing to GE connector
                Trans.ConnectorES = [ ...
               125   117   119   127    96   124   116   122    92   134   140   143   138   133   129   137 ...
               173   175   169   132   197   172   205   201   195   162   168   208   207   200   163   251 ...
               243    18   245   253   224   250   214   244   218   254    35    43    47    41    33     3 ...
                 9    10    12     6     5    38    73    77    42    67    66    40    76    78    72    34 ...
                94   123   115   121    90   126   118   120   128   130   135   141   142   136   131   167 ...
               171   144   139   165   193   199   203   174   170   164   161   206   176   202   166   256 ...
               249   241    17   247   255   252   246   242   248   222     2    37    45    14    39     1 ...
                 4     8    13    11     7    65    69    75    46    71    36    70    44    80    74    68 ...
                98   106   112   104    97   103   111   113   105   158   152   147   148   154   159   187 ...
               149   145   181   185   221   186   211   213   217   223   190   212   180   210   191   227 ...
               233    16   239   231   225   228   234   240   204   196    63    57    19    51    59    32 ...
                29    24    20    23    28    93    56    79    81    89    95    58    83    82    86    87 ...
               100   108   110   102   101   109   114   107    99   156   151   146   150   155   157   183 ...
               177   179   153   160   219   215   209   184   188   189   216   178   182   220   192   229 ...
               235    15   237   230   226   232   238   236   198   194    61    55    49    53    27    30 ...
                25    21    22    26    31    62    54    48    52    91    60    88    84    50    85    64 ]';  
            
            case 'RC15gV'
                %% RC15gV Vermon Row-Column array
                if ~isfield(Trans,'frequency'), Trans.frequency = 13.6; end % Center frequency in MHz from Vermon data sheet
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [8.16, 19.04]; end % 80% fractional bandwidth
                Trans.type = 4;     % Array geometry is 2D Row-Column array.
                Trans.connType = 7; % UTA 408-GE connector
                Trans.numelements = 160;
                Trans.ConnectorES = (1:160)'; % 1:1 mapping from transducer elements to element signals at connector
                Trans.spacingMm = 0.11; % element spacing in mm.
                Trans.ElementPos = zeros(Trans.numelements,5);
                Trans.ElementPos(1:80,1) = Trans.spacingMm*((-79/2):(79/2)); % in Wavelengths
                Trans.ElementPos(81:160,2) = Trans.spacingMm*((-79/2):(79/2));
                Trans.elementWidth = 0.10; % 0.9 * Trans.spacingMm
                Trans.elementLength = 8.8; % element Length in mm from Vermon data sheet
                Trans.lensCorrection = 0; % No lens
                Trans.impedance = 100; % Placeholder value; no data from Vermon
                Trans.maxHighVoltage = 30;

            case 'P5-64vN'
                %% P5-64vN Imasonics Linear Probe
                if ~isfield(Trans,'frequency'), Trans.frequency = 5.0; end % as specified by Imasonics
                if ~isfield(Trans,'Bandwidth'), BW = 0.6; Trans.Bandwidth = Trans.frequency*[1-BW/2, 1+BW/2]; end % 60% fractional bandwidth
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 10; % I-PEX connector
                Trans.id = -1; % Probe has no ID (overrides hex value used in 'KnownTransducers' array)
                Trans.numelements = 64;
                Trans.ConnectorES = (1:64)'; % 1:1 mapping from transducer elements to element signals at connector
                Trans.elementWidth = 0.400; % width in mm
                    Trans.elevationApertureMm = 10; % (mm) active elevation aperture (spec)
                    Trans.elevationFocusMm = 30; % (mm) nominal elevation focus depth from lens on face of transducer (spec)
                Trans.spacingMm = 0.500;   % Spacing between elements in mm.
                Trans.ElementPos = zeros(Trans.numelements,5);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                Trans.lensCorrection = 0.0; % in mm units;
                Trans.impedance = 50;

            case 'P5-128vN'
                %% P5-128vN Imasonics Linear Probe
                if ~isfield(Trans,'frequency'), Trans.frequency = 5.0; end % as specified by Imasonics
                if ~isfield(Trans,'Bandwidth'), BW = 0.6; Trans.Bandwidth = Trans.frequency*[1-BW/2, 1+BW/2]; end % 60% fractional bandwidth
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 10; % I-PEX connector
                Trans.id = -1; % Probe has no ID (overrides hex value used in 'KnownTransducers' array)
                Trans.numelements = 128;
                Trans.ConnectorES = (1:128)'; %  mapping from transducer elements to element signals at connector
                Trans.elementWidth = 0.400; % width in mm
                    Trans.elevationApertureMm = 10; % (mm) active elevation aperture (spec)
                    Trans.elevationFocusMm = 30; % (mm) nominal elevation focus depth from lens on face of transducer (spec)
                Trans.spacingMm = 0.500;   % Spacing between elements in mm.
                Trans.ElementPos = zeros(Trans.numelements,5);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                Trans.lensCorrection = 0.0; % in mm units;
                Trans.impedance = 50;

            case 'H-101'
                %% H-101 Sonic Concepts Probe
                % Sonic Concepts 4 ring annular HIFU probe, 1.1 MHz for FUS Elite 3000 3D pairing 2
                if ~isfield(Trans,'frequency'), Trans.frequency = 1.1; end % nominal frequency in MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [0.77, 1.43]; end
                Trans.type = 3;    % annular array
                Trans.connType = 1; % HDI connector
                Trans.numelements = 4;
                % ElementPos row:  [ ri ro rc zc  theta ]
                %    ri = distance from inner edge of ring to the Z axis (mm)
                %    ro = distance from outer edge of ring to the Z axis (mm)
                %    rc = distance from center of ring to the Z axis (mm)
                %    zc = distance from center of ring to the Z=0 plane (mm)
                %    theta = angle between the surface normal at center of the ring relative to the Z axis (rad)
                % A positive theta value means the element is pointed
                % toward the Z axis.  The "center" of the ring is defined
                % as the point that divides the ring into two subrings of equal area.
                %
                % ElementPos row:  [ ri     ro      rc      zc      theta ]
                Trans.ElementPos = [ 15.93	20.91	18.60	2.80	deg2rad(17.12);...
                                    20.91	24.86	22.98	4.33	deg2rad(21.32);...
                                    24.86	28.20	26.59	5.87	deg2rad(24.88);...
                                    28.20	31.04	29.67	7.40	deg2rad(28.00) ];
                Trans.spacingMm = 0; % dummy value
                Trans.radiusMm = 63.2; % radius in mm of the spherical probe surface
                Trans.maxHighVoltage = 50;  % arbitrary limit... user can change in setup script
                Trans.maxAvgPower = 700;    % Watts (Absolute Limit set by Sonic Concepts)
                Trans.impedance = 60; % average value near 1.1 MHz (at 1.1, Z=70 for quai-CW pulses)
                % Note each of the 4 elements is driven through a combiner
                % from 16 channels of the Vantage system.  Each row of the
                % Trans.ConnectorES array lists the channels driving one
                % element and effectively connected in parallel.
                Trans.ConnectorES = [ 65  69  70  71  73  74  75  76  77  78  79  80  84  85  90  91;...
                                    66  67  68  72  81  82  83  86  87  88  89  92  93  94  95  96;...
                                    97  98  99  101 102 103 107 110 114 116 117 118 119 121 122 126;...
                                    100 104 105 106 108 109 111  112  113  115  120  123  124  125  127  128];

            case 'H-104'
                %% H-104 Sonic Concepts Probe
                % Sonic Concepts 3 ring annular HIFU probe, 0.5 MHz for FUS Elite 3000 3D pairing 1
                if ~isfield(Trans,'frequency'), Trans.frequency = 0.5; end % nominal frequency in MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [0.35, 0.65]; end
                Trans.type = 3;    % annular array
                Trans.connType = 1; % HDI connector
                Trans.numelements = 3;
                Trans.radiusMm = 63.2; % radius in mm of the spherical probe surface
                % ElementPos row:  [ ri     ro      rc      zc      theta ]
                Trans.ElementPos = [ 16.07	21.92	19.24	3.00	deg2rad(17.72);...
                                    21.92	26.41	24.28	4.85	deg2rad(22.60);...
                                    26.41	30.05	28.30	6.69	deg2rad(26.60) ];
                Trans.spacingMm = 0; % dummy value
                Trans.maxHighVoltage = 50;  % arbitrary limit... user can change in setup script
                Trans.maxAvgPower = 550;    % Watts (Absolute Limit set by Sonic Concepts)
                Trans.impedance = 70; % nominal value @ 500 kHz
                % Note each of the 3 elements is driven through a combiner
                % from 16 channels of the Vantage system.  Each row of the
                % Trans.ConnectorES array lists the channels driving one
                % element and effectively connected in parallel.
                Trans.ConnectorES = [65  69  70  71  73  74  75  76  77  78  79  80  84  85  90  91;...
                                   66  67  68  72  81  82  83  86  87  88  89  92  93  94  95  96;...
                                   97  98  99 101 102 103 107 110 114 116 117 118 119 121 122 126];

            case 'H-106'
                %% H-106 Sonic Concepts Probe
                % Sonic Concepts 8 ring annular HIFU probe, 2.0 MHz for FUS Elite 3000 3D pairing 3
                if ~isfield(Trans,'frequency'), Trans.frequency = 2.0; end % nominal frequency in MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [1.4, 2.6]; end
                Trans.type = 3;    % annular array
                Trans.connType = 1; % HDI connector
                Trans.numelements = 8;
                % ElementPos row:  [ ri     ro      rc      zc      theta ]
                Trans.ElementPos = [ 15.81	18.56	17.25	2.40	deg2rad(15.84);...
                                    18.56	20.96	19.80	3.18	deg2rad(18.26);...
                                    20.96	23.10	22.06	3.97	deg2rad(20.43);...
                                    23.10	25.06	24.10	4.78	deg2rad(22.42);...
                                    25.06	26.86	25.98	5.59	deg2rad(24.27);...
                                    26.86	28.54	27.72	6.40	deg2rad(26.01);...
                                    28.54	30.12	29.35	7.23	deg2rad(27.67);...
                                    30.12	31.50	30.82	8.02	deg2rad(29.18) ];

                Trans.spacingMm = 0; % dummy value
                Trans.radiusMm = 63.2; % radius in mm of the spherical probe surface
                Trans.maxHighVoltage = 50;  % arbitrary limit... user can change in setup script
                Trans.maxAvgPower = 700;    % Watts (Absolute Limit set by Sonic Concepts)
                Trans.impedance = 50; % average value near 2 MHz.
                % Note each of the 8 elements is driven through a combiner
                % from 8 channels of the Vantage system.  Each row of the
                % Trans.ConnectorES array lists the channels driving one
                % element and effectively connected in parallel.
                Trans.ConnectorES = [ 65    69    70    71    75    76    78    79;...
                                    66    67    68    72    88    89    93    94;...
                                    97    98    99   101   102   103   107   110;...
                                   114   116   117   118   119   121   122   126;...
                                    73    74    77    80    84    85    90    91;...
                                    81    82    83    86    87    92    95    96;...
                                   100   104   105   106   108   109   111   112;...
                                   113   115   120   123   124   125   127   128];

            case 'H-301'
                %% H-301 Sonic Concepts Probe
                % Sonic Concepts 128 el, 1.1 MHz HIFU array for FUS Elite 3000 3D pairing 5
                if ~isfield(Trans,'frequency'), Trans.frequency = 1.1; end % nominal frequency in MHz
                % Vantage:  1.953 is closest supported frequency to 2.0 MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [0.77, 1.43]; end
                Trans.type = 2;    % Three dimensional array
                Trans.connType = 1; %  HDI 260-ZIF
                Trans.numelements = 128;
                Trans.elementWidth = 10.15; % width in mm
                Trans.spacingMm = 10.15/0.9; % in mm, wild guess
                % load H301 geometry
                arraygeom = computeHIFUgeometry(Trans);  % in mm
                Trans.ElementPos(1:128, 1:3) = arraygeom;
                Trans.radiusMm = 150; % Geometric Focus is 150 mm = 150 * (speedOfSound/Trans.frequence) in wavelength
                Trans.ElementPos(1:128, 4) = atan(arraygeom(:,1) ./ (Trans.radiusMm-arraygeom(:,3)) ); % AZ = atan(x/z)
                Trans.ElementPos(1:128, 5) = atan(arraygeom(:,2) ./ sqrt( arraygeom(:,1).^2 + (Trans.radiusMm-arraygeom(:,3)).^2 ) ); % EL = atan(y/sqrt(x^2+z^2)), where z has to be wrt geometric focus.
                Trans.lensCorrection = 1.0; % in mm units, value is a wild guess!!
                Trans.maxAvgPower = 14*128; % 1750;    % 14 Watts * 128 elements (Absolute Limit set in Test Report for Pair-05 SN-01 )
                Trans.impedance = 104;      % From SCI Test Report for Pair-05 SN-01
                Trans.maxHighVoltage = 75;  % sqrt(56*Trans.impedance) = 76.3 = sqrt(Pmax*Z)
                % Trans.ConnectorES same as for H-302
                Trans.ConnectorES = [ ...
                    49, 55, 52, 61, 17, 23, 29, 31,116,119,121,123,125,127, 84, 87,...
                    91, 93, 95, 62, 53, 56, 59, 50, 18, 21, 24, 27, 30, 32,114,117,...
                    120,122,124,126,128, 82, 85, 88, 90, 92, 94, 96, 51, 54, 58, 57,...
                    60, 63, 19, 20, 22, 25, 26, 28, 14,113,115,118,101,103,105,107,...
                    109,111, 81, 83, 86, 89, 71, 73, 75, 77, 79, 64, 35, 37, 39, 41,...
                    43, 45, 47, 1, 4, 6, 8, 10, 12, 16, 15, 97, 98, 99,100,102,...
                    104,106,108,110,112, 65, 66, 67, 68, 69, 70, 72, 74, 76, 78, 80,...
                    33, 34, 36, 38, 40, 42, 44, 46, 48, 2, 3, 5, 7, 9, 11, 13]';


            case 'H-302'
                %% H-302 Sonic Concepts Probe
                % Sonic Concepts 128 el, 2 MHz HIFU array for FUS Elite 3000 3D pairing 6
                if ~isfield(Trans,'frequency'), Trans.frequency = 1.953; end % nominal frequency in MHz
                % Vantage:  1.953 is closest supported frequency to 2.0 MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [1.4, 2.7]; end
                Trans.type = 2;    % Three dimensional array
                Trans.connType = 1; %  HDI 260-ZIF
                Trans.numelements = 128;
                Trans.elementWidth = 10.15; % width in mm
                Trans.spacingMm = 10.15/0.9; % in mm, wild guess
                Trans.elevationApertureMm = Trans.spacingMm;
                % load H302 geometry
                arraygeom = computeHIFUgeometry(Trans);  % in mm

                Trans.ElementPos = zeros(Trans.numelements,5);
                Trans.ElementPos(1:128, 1:3) = arraygeom;
                Trans.radiusMm = 150; % Geometric Focus is 150 mm = 150 * (speedOfSound/(1000*Trans.frequency)) in wavelength
                Trans.ElementPos(1:128, 4) = atan(arraygeom(:,1) ./ (Trans.radiusMm-arraygeom(:,3)) ); % AZ = atan(x/z)
                Trans.ElementPos(1:128, 5) = atan(arraygeom(:,2) ./ sqrt( arraygeom(:,1).^2 + (Trans.radiusMm-arraygeom(:,3)).^2 ) ); % EL = atan(y/sqrt(x^2+z^2)), where z has to be wrt geometric focus.
                Trans.lensCorrection = 1.0; % in mm units, value is a wild guess!!
                % --- from SCI Test Report for H-302 SN-05
                Trans.maxAvgPower = 1024;   % Watts
                Trans.impedance = 100;      % Ohms
                Trans.maxHighVoltage = 75;  % sqrt(100*Trans.impedance) = 100 = sqrt(Pmax*Z)
                % ---
                % If the HIFU transducer is connected to the connector 1
                % update on April 24 2018, follow the Test Report on H-302-B-02 2 MHz, SN: 005 12-11-2017
                Trans.ConnectorES = [ ...
                    49, 55, 52, 61, 17, 23, 29, 31,116,119,121,123,125,127, 84, 87,...
                    91, 93, 95, 62, 53, 56, 59, 50, 18, 21, 24, 27, 30, 32,114,117,...
                    120,122,124,126,128, 82, 85, 88, 90, 92, 94, 96, 51, 54, 58, 57,...
                    60, 63, 19, 20, 22, 25, 26, 28, 14,113,115,118,101,103,105,107,...
                    109,111, 81, 83, 86, 89, 71, 73, 75, 77, 79, 64, 35, 37, 39, 41,...
                    43, 45, 47, 1, 4, 6, 8, 10, 12, 16, 15, 97, 98, 99,100,102,...
                    104,106,108,110,112, 65, 66, 67, 68, 69, 70, 72, 74, 76, 78, 80,...
                    33, 34, 36, 38, 40, 42, 44, 46, 48, 2, 3, 5, 7, 9, 11, 13]';

            case {'H-313', 'H-313A'}
                %% H-313 Sonic Concepts Probe
                % Sonic Concepts 64 el, 0.5 MHz HIFU array for FUS Elite 3000 3D pairing 4
                % VTS-1698: There are two different versions of this probe,
                % one named H-313A and the other H-313.  The transducer
                % itself is identical in both versions, with differences in
                % the cable length and wiring:
                %    The H-313A (probe ID code 04C305) is a unique version
                %    for which only one was built, originally identified as
                %    H-313 serial number 01.  It has a shorter cable, does
                %    not have RoHs compliance, and has unique wiring from
                %    elements to the probe connector as listed below in
                %    Trans.ConnectorES
                %
                %    The H-313 (probe ID code 04C304) identifies ongoing
                %    production of the H-313 starting with serial number
                %    -02  It has a 4 meter cable, RoHs compliance, and
                %    revised wiring from elements to the probe connector as
                %    listed below in Trans.ConnectorES
                if ~isfield(Trans,'frequency'), Trans.frequency = 0.5; end % nominal frequency in MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [0.35, 0.65]; end
                Trans.type = 2;    % Three dimensional array
                Trans.connType = 1; % HDI 260-ZIF
                Trans.numelements = 64;
                Trans.elementWidth = 13.34; % Value from SCI
                Trans.spacingMm =  Trans.elementWidth/0.9; % Value estimated... not used
                % load H313 geometry
                arraygeom = computeHIFUgeometry(Trans);  % in mm
                arraygeom(:,2) = -arraygeom(:,2); % inverted geometry see (VTS-3223)
                Trans.ElementPos = zeros(Trans.numelements,5);
                Trans.ElementPos(1:64, 1:3) = arraygeom;
                Trans.radiusMm = 150; % Geometric Focus is 150 mm = 150 * (speedOfSound/Trans.frequence) in wavelength
                Trans.ElementPos(1:64, 4) = atan(arraygeom(:,1) ./ (Trans.radiusMm-arraygeom(:,3)) ); % AZ = atan(x/z)
                Trans.ElementPos(1:64, 5) = atan(arraygeom(:,2) ./ sqrt( arraygeom(:,1).^2 + (Trans.radiusMm-arraygeom(:,3)).^2 ) ); % EL = atan(y/sqrt(x^2+z^2)), where z has to be wrt geometric focus.
                Trans.lensCorrection = 1.0; %  in mm units, value is a wild guess!!
                Trans.maxAvgPower = 64*11;  % Watts (Absolute Limit set by Sonic Concepts)
                Trans.impedance = 60; % channel average Value @ 0.5 MHz
                Trans.maxHighVoltage = min(96, sqrt(56*64*abs(Trans.impedance)));  % actual limit based on Max pulsed Power 56W/ch

                if strcmp(KnownTransducers{probenum, 1}, 'H-313A')
                    % Element wiring for H-313A serial number 01
                    Trans.ConnectorES = [116, 122, 125, 128, 112, 106, 100,  98, ...
                                          93,  90,  88,  86,  84,  82,  77,  74, ...
                                          70,  68,  66, 115, 118, 121, 124, 127, ...
                                         111, 108, 105, 102,  99,  97,  95,  92, ...
                                          89,  87,  85,  83,  81,  79,  76,  73, ...
                                          71,  69,  67,  65, 114, 117, 119, 120, ...
                                         123, 126, 110, 109, 107, 104, 103, 101, ...
                                         113,  96,  94,  91,  80,  78,  75,  72]';
                else
                    % Element wiring for H-313 serial numbers 02 and up
                    % from Sonic Concepts "Test Report H-313 SNO2.pdf"
                    Trans.ConnectorES = [125, 119, 116, 113,  97, 103, 109, 111, ...
                                          84,  87,  89,  91,  93,  95,  68,  71, ...
                                          75,  77,  79, 126, 123, 120, 117, 114, ...
                                          98, 101, 104, 107, 110, 112,  82,  85, ...
                                          88,  90,  92,  94,  96,  66,  69,  72, ...
                                          74,  76,  78,  80, 127, 124, 122, 121, ...
                                         118, 115,  99, 100, 102, 105, 106, 108, ...
                                         128,  81,  83,  86,  65,  67,  70,  73]';
                end


            case 'IP-104'
                %% IP-104 Sonic Concepts Probe
                % Sonic Concepts 128 el, 3.5 MHz phased array for FUS Elite 3000 3D pairings 4, 5
                % note only the first 64 elements are used in pairing 4
                if ~isfield(Trans,'frequency'), Trans.frequency = 3.472; end % nominal frequency in MHz
                % Vantage:  3.472 is closest supported frequency to 3.5 MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [1.60, 4.70]; end % 89% of 3.5 MHz Fc, based on the test report
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI connector=1
                Trans.numelements = 128;
                Trans.ConnectorES = (1:128)';
                Trans.spacingMm = 0.22;
                Trans.elevationApertureMm = 13.5; % active elevation aperture in mm
                Trans.elevationFocusMm = 75; % nominal elevation focus depth from lens on face of transducer
                Trans.elementWidth = 0.9*Trans.spacingMm; % width in mm
                Trans.ElementPos = zeros(Trans.numelements,5);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                Trans.lensCorrection = 1.0; % in mm units, value is a wild guess!!
                Trans.impedance = 156.8;
                Trans.maxHighVoltage = 70;

            case 'IP-105'
                %% IP-105 Sonic Concepts Probe
                % Sonic Concepts 64 el, 5 MHz phased array for FUS Elite 3000 3D pairings 1, 2, 3
                if ~isfield(Trans,'frequency'), Trans.frequency = 5.208; end % nominal frequency in MHz
                % Vantage:  5.208 is closest supported frequency to 5 MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [2.89, 7.53]; end % 89% of 5 MHz Fc, based on the test report
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 64;
                Trans.ConnectorES = (1:64)';
                Trans.spacingMm = 0.150;   % Spacing between elements in mm.
                Trans.elevationApertureMm = 10; % active elevation aperture in mm (estimate)
                Trans.elevationFocusMm = 65; % nominal elevation focus depth from lens on face of transducer (estimate)
                Trans.elementWidth = 0.9*Trans.spacingMm; % width in mm
                Trans.ElementPos = zeros(Trans.numelements,5);
                %   Set default element x positions (units in mm).
                Trans.ElementPos(1:64,1) = Trans.spacingMm*(-((64-1)/2):((64-1)/2));
                Trans.lensCorrection = 1.5; % in mm units
                Trans.impedance = 129.5;
                Trans.maxHighVoltage = 70;

            case 'Matrix1024-3'
                %% Matrix1024-3   Vermon 1024 EL 3 MHz Matrix Array
                if ~isfield(Trans,'frequency'), Trans.frequency = 3.47; end % nominal frequency 3.5 MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = Trans.frequency*[0.7, 1.3]; end
                Trans.connType = 8;
                Trans.type = 2;             % Array geometry is 2D.
                Trans.numelements = 1024;
                Trans.impedance = 50;
                Trans.elementWidth = .275; % mm
                Trans.spacingMm = .3;        % Spacing between elements in mm.
                Trans.maxHighVoltage = 30;
                ConnMap = [128:-1:65, 1:1:64, 192:-1:129, 193:1:256];
                % ConnMap is the mapping for each of the four 256-channel
                % connector boards coming from the Vermon probe, from probe
                % elements Within each 256-element group to connector
                % signals for that group as defined by the connType 8 interface.
                Trans.ConnectorES = [ConnMap, ConnMap+256, ConnMap+512, ConnMap+768]';
                % Trans.ConnectorES is the overall mapping from probe elements to
                % the 1024 Element Signals as defined by the connType 8 interface.
                Trans.ElementPos = getVermonMatrixElpos(Trans.spacingMm); % ElementPos in mm units

            case 'Matrix1024-8'
                %% Matrix1024-8   Vermon 1024 EL 8 MHz Matrix Array
                if ~isfield(Trans,'frequency'), Trans.frequency = 7.81; end % nominal frequency 8 MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = Trans.frequency*[0.7, 1.3]; end % 69% rel. bw
                Trans.connType = 8;
                Trans.type = 2;             % Array geometry is 2D.
                Trans.numelements = 1024;
                Trans.impedance = 50;
                Trans.elementWidth = .275; % mm
                Trans.spacingMm = .3;        % Spacing between elements in mm.
                Trans.maxHighVoltage = 30;
                ConnMap = [128:-1:65, 1:1:64, 192:-1:129, 193:1:256];
                % ConnMap is the mapping for each of the four 256-channel
                % connector boards coming from the Vermon probe, from probe
                % elements Within each 256-element group to connector
                % signals for that group as defined by the connType 8 interface.
                Trans.ConnectorES = [ConnMap, ConnMap+256, ConnMap+512, ConnMap+768]';
                % Trans.ConnectorES is the overall mapping from probe elements to
                % the 1024 connector signals as defined by the connType 8 interface.
                Trans.ElementPos = getVermonMatrixElpos(Trans.spacingMm); % ElementPos in mm units


            otherwise
                Trans.type = [];
                Trans.frequency = [];
                Trans.spacingMm = [];
                Trans.elementWidth = [];
                Trans.ElementSens = [];
                Trans.lensCorrection = [];
                Trans.ElementPos = [ 0 0 0 0 ];
                if verbose > 2
                    disp(' ');
                    disp(['computeTrans Status: Data not available for ', Trans.name]);
                    disp('Trans structure must be provided in user script.');
                    disp(' ');
                end
        end

        % Set a conservative value for maxHighVoltage, if not already defined
        if ~isfield(Trans,'maxHighVoltage'), Trans.maxHighVoltage = 50; end

        % Now convert all units as required, based on Trans.units
        scaleToWvl = Trans.frequency/speedOfSound; % conversion factor from mm to wavelengths

        % regardless of units, always provide spacing and radius in
        % wavelengths if they have been defined
        if isfield(Trans, 'spacingMm') && ~isempty(Trans.spacingMm)
            Trans.spacing = Trans.spacingMm * scaleToWvl;   % Spacing between elements in wavelengths.
        end
        if  isfield(Trans, 'radiusMm') && ~isempty(Trans.radiusMm)
            Trans.radius = Trans.radiusMm * scaleToWvl; % convert radiusMm to wavelengths
        end

        % define Trans.ElementSens based on Trans.elementWidth, but only if
        % user has not already defined it; assign default elementWidth if
        % it doesn't exist or is empty
        if ~isfield(Trans,'ElementSens')
            % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
            if ~isfield(Trans,'elementWidth') || isempty(Trans.elementWidth)
                % create default value of zero if not assigned in case
                % statements above (zero implies the element is a point
                % source)
                Trans.elementWidth = 0;
            end
            Theta = (-pi/2:pi/100:pi/2);
            Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
            % note at this point elementWidth is in mm, so we have to
            % convert to wavelengths for the ElementSens calculation
            eleWidthWl = Trans.elementWidth * scaleToWvl;
            if eleWidthWl < 0.01
                % avoid the divide by zero for very small values (in this
                % case the sinc function will be extremely close to 1.0 for
                % all Theta, so we only need the cos term)
                Trans.ElementSens = abs(cos(Theta));
            else
                Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
            end
        end


        if strcmp(Trans.units, 'wavelengths')
            % convert all mm unit variables to wavelengths.  Note columns 4
            % and 5 of the ElementPos array are angles in radians, and do
            % not require units conversion
            % VTS-2610 added support for optional Trans.elementLength
            if  isfield(Trans, 'elementLength') && ~isempty(Trans.elementLength)
                Trans.elementLength = Trans.elementLength * scaleToWvl;
            end
            Trans.elementWidth = Trans.elementWidth * scaleToWvl;
            Trans.ElementPos(:,1) = Trans.ElementPos(:,1) * scaleToWvl;
            Trans.ElementPos(:,2) = Trans.ElementPos(:,2) * scaleToWvl;
            Trans.ElementPos(:,3) = Trans.ElementPos(:,3) * scaleToWvl;
            if Trans.type == 3
                % for type 3 annular arrays, the fourth column is also a distance, not an angle
                Trans.ElementPos(:,4) = Trans.ElementPos(:,4) * scaleToWvl;
            end
            Trans.lensCorrection = Trans.lensCorrection * scaleToWvl;
        end

    case 2
        % Two inputs provided - 1st input is Trans.name, 2nd is parameter desired (currently only
        % 'maxHighVoltage' allowed).
        nameStr = varargin{1};
        if ischar(nameStr)
            probenum = find(strcmpi(nameStr, KnownTransducers(:, 1)), 1);
            % if 1st input is not a recognized transducer name, set probenum to
            % zero to flag as unrecognized transducer
            if isempty(probenum)
                probenum = 0; % special case of zero will trigger assignment of default value
            end
            Param = varargin{2};
            switch Param
                case 'maxHighVoltage'
                    if probenum == 0
                        % unrecognized transducer name, so assign default
                        % maxHighVoltage limit at hw max value
                        Trans = 96;
                    else
                        % return maxHighVoltage from KnownTransducers array
                        Trans = KnownTransducers{probenum, 3};
                    end
                case 'HVMux'
                    if probenum == 0
                        % unrecognized transducer name, so assign default
                        % of no HVMux
                        Trans = 0;
                    else
                        % return HVMux status from KnownTransducers array
                        Trans = KnownTransducers{probenum, 4};
                    end
                otherwise
                    error('computeTrans: Unrecognized parameter in 2nd argument.');
            end
        else
            error('computeTrans: When called with 2 inputs, 1st input must be transducer name.');
        end

    otherwise
        error('computeTrans: computeTrans accepts one or two input arguments.');

end
end  % End of the computeTrans function

%% Subfunctions defined within computeTrans
function ELpos = getVermonMatrixElpos(spacing)
    % Create element position array for Vermon 1024 Matrix probes
    % Define a square array EM (Element Map) that shows the location of
    % each element (numbered 1 to 128) on the physical grid of transducer
    % elements.  Locations that represent either an element that does not exist
    % or an element that is not connected are identified by setting the element
    % number to zero.
    GSRow = 35;
    GSCol = 32;
    EM = zeros(GSRow,GSCol);  % define the array and fill it with zeros
    % Physical element array order defined in a for loop with three row gaps
    for ii = 1:GSRow
        if ii==1; jj=0; end
        if ii==9 || ii==18 || ii==27
            continue;
        end
        jj=jj+1;
        EM(ii,:) = (jj-1)*GSCol+(1:32);
    end
    % Next use EM array to assign xy position for each element in Trans.ElementPos array
    ELpos = zeros(1024,5,'double');  % define array and fill with zeros
    ELpos(:,1:2) = -1;  % fill x,y positions with -1, to allow test for
    for ii=1:GSRow
        for jj=1:GSCol
            if(EM(ii,jj)>0)   % note positions are relative to center of array
                ELpos(EM(ii,jj),1) = (jj-1)*spacing - ((GSCol-1)*spacing)/2;
                ELpos(EM(ii,jj),2) = ((GSRow-1)*spacing)/2 - (ii-1)*spacing;
            end
        end
    end
end

function lensCorrectionMm = computeLensCorrection(speedOfSound, matchingLayerThicknessMm, matchingLayerVelocity, lensThicknessMm, lensVelocity)
    % Compute the lens correction
    % matchingLayerVelocity and lensVelocity in m/sec, 
    % matchingLayerThicknessMm and lensThicknessMm in mm
    % speedOfSound in km/s
    % lensCorrectionMm in mm
    
    lensCorrectionMm = 1000 * speedOfSound * (matchingLayerThicknessMm/matchingLayerVelocity + lensThicknessMm/lensVelocity);
    
end

