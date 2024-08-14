function [ RcvProfile ] = computeRcvProfile( RcvProfile )
% computeRcvProfile: Generate complete RcvProfile structure from partial
% user input.  Check for validity of user-defined values, assign defaults
% to all fields undefined by user, and create the VDASRcvProfile output
% array to be sent to the HW system descriptor tables by the HAL.
% Copyright 2001-2019 Verasonics, Inc.  Verasonics Registered U.S. Patent and Trademark Office.
%

% the input argument is a predefined RcvProfile structure or a row vector
% array of multiple structures of any length.  The output argument is the
% same structure or array of structures, with all values filled in.

% For each of the defined fields in the RcvProfile structure that can be
% set by the user (identified here with a section heading such as "%%
% AntiAliasCutoff"), computeRcvProfile first checks to see if it has been
% set by the user and if not the default value is assigned.  Then the field
% value is translated into the register bit settings used within the AFE
% chip.  These settings are placed in the "VDASRcvProfile" array, which is
% accessed by the HAL to define the descriptors in the HW system used to
% actually program the AFE devices at run time.




% First step: get the receive index circuit configuration from the base
% workspace.  This "RXindex" value identifies which type of AFE chip is
% present in the system; allowed values of some of the RcvProfile fields
% are different for the different AFE devices.  Three different types of
% AFE devices are supported by the Vantage system SW:
    % AFE5808A: used in early standard-frequency systems.  Does not support
    % the high freuency option.
    % PAFE5812: used in early high-frequency systems, based on an
    % unreleased "pre-production" version of the AFE5812 chip from Texas
    % Instruments.
    % AFE5812: Released version of the AFE5812 from TI, will be included in
    % all new Vantage systems of all configuratons once the inventory of
    % AFE5808A and PAFE5812 chips has been exhausted.

% The programmable registers in the AFE chip all are set to a default value
% of all zeros at power-up.  TI has defined the register encoding for most
% bit fields such that a value of zero represents the typical default value
% most applications would want. Therefore we only need to write to those
% few registers where we use a value different than the TI default, or
% where there is a field that we wish to make user-programmable thorugh the
% setRcvProfile command.  The number of registers that need to be written
% to varies with the different versions of the AFE device so the
% VDASProfile output array is initialized uniquely for each AFEtype in the
% case statements below.
%
% NOTE:  Some AFE registers are not allowed to be set.  Attempting to
% do so will give an error from the HAL.  An example is the registers
% used for ADC DC Offset Calibration, which is performed and set by the HAL
% so the HAL does not allow those registers to be specified in an AFE
% register set to prevent the HAL register values from being overwritten.
% Those registers are 13, 15, 17, 19, 22, 25, 27, 29, and 31.

% Revised Sep 1 2019 VTS-1416 add TX, RX index values for acquisition
%     module with active clamp


RXindex = evalin('base', 'Resource.SysConfig.RXindex');
switch RXindex
    case {1, 3}
        % HW configuration has AFE 5808A chips installed, AAF cutoff
        % can go from 10 to 30 MHz
        AFEtype = 0; % AFE5808A
        numRegisters = 8; % number of registers to be written to for this AFEtype
        % create VDASProfile array for output through Hal to CGD descriptor
        % memory
        VDASProfile = zeros(numRegisters, 2);
        % Fill in the register addresses in the first column
        VDASProfile(:, 1) = [1; 3; 21; 33; 51; 52; 53; 59];
    case {2, 4}
        %  HW configuration has pre-production PAFE5812 chips installed,
        %  AAF cutoff can go from 10 to 60 MHz
        AFEtype = 1; % pre-released version of AFE5812
        numRegisters = 8; % number of registers to be written to for this AFEtype
        % create VDASProfile array for output through Hal to CGD descriptor
        % memory
        VDASProfile = zeros(numRegisters, 2);
        % Fill in the register addresses in the first column
        VDASProfile(:, 1) = [1; 3; 21; 33; 51; 52; 53; 59];
    case {5, 6}
        %  HW configuration has production AFE5812 chips installed, AAF
        %  cutoff can go from 5 to 50 MHz
        AFEtype = 2; % production version of AFE5812
        numRegisters = 9; % number of registers to be written to for this AFEtype
        % create VDASProfile array for output through Hal to CGD descriptor
        % memory
        VDASProfile = zeros(numRegisters, 2);
        % Fill in the register addresses in the first column
        VDASProfile(:, 1) = [1; 3; 21; 33; 51; 52; 53; 59; 61];
    otherwise
        error('computeRcvProfile: Unrecognized RXindex value of %d.\n', RXindex);
end

% determine whether the system is High Frequency, since some parameters
% have different defaults for that configuration
if evalin('base', 'Resource.SysConfig.TXindex == 2 || Resource.SysConfig.TXindex == 4 || Resource.SysConfig.TXindex == 6')
    % this is a high frequency system
    highFreqSys = 1;
else
    % standard or low frequency configurations
    highFreqSys = 0;
end


numProfiles = size(RcvProfile, 2);
% the number of unique RcvProfile structures the user has defined in their script.  If
% the script does not define the RcvProfile structure at all, VSX will
% create one before calling computeRcvProfile (and thus it will be
% populated with the default values Verasoncis has defined for all fields)



% Populate the VDASProfile array with the 'hard coded' register settings
% that are not user-programmable through fields in the RcvProfile
% structure.  This same pre-populated VDASProfile array will be used as the
% starting point for each of the numProfiles structures defined by the
% user.
for i = 1:numRegisters % row index for the VDASProfile array
    registerAddress = VDASProfile(i, 1); % AFE register address for this row
    switch registerAddress
        case 1
            % Register contents at address 1
            % Set bit 13 to select external A/D reference
            VDASProfile(i, 2) = bitset(VDASProfile(i, 2), 14);
            % note matlab counts from one not zero for bitset argument
        case 3
            % Register contents at address 3
            % Set bit 15 to select external A/D reference
            VDASProfile(i, 2) = bitset(VDASProfile(i, 2), 16);
            % note matlab counts from one not zero for bitset argument
        case 52
            % Register contents at address 52
            % Set bits 5 and 8 as shown to enable the direct control of the
            % resistor switches used to set LNA input impedance.  Then the
            % RcvProfile.LnaZinSel field will be used to set bits 4:0 as
            % required for one of the 32 choices available to the user.
            VDASProfile(i, 2) = bin2dec( '1 0 0 1 0 0 0 0 0');

            % Set bits 10:9 to [0 1] to select 1.5 V. LNA input clamp
            VDASProfile(i, 2) = bitset(VDASProfile(i, 2), 10);
            % note matlab counts from one not zero for bitset argument
        otherwise
            % for all registers other than those identified in the case
            % statements above, nothing needs to be written since the
            % hard-coded settings desired by Verasonics are the same as the
            % defaults defined by TI, i.e. what we get if we don't write to
            % the register at all- all bits set to zero, since the
            % VDASProfile array was created by filling it with all zeros.
    end
end

% paste hard-coded values into each profile structure, and add the
% user-defined values
statusMsg = 0; % flag to avoid repeating same message for every profile
for i = 1:numProfiles
    newVDASProfile = VDASProfile;

    %% PowerDownA2D
    if ~isfield(RcvProfile(i),'PowerDownA2D') || isempty(RcvProfile(i).PowerDownA2D)
        RcvProfile(i).PowerDownA2D = 0; % default value is powered up
    end
    switch RcvProfile(i).PowerDownA2D
        case 0 % Power Up A/D
            % set D0 and D10 of reg adr 1 to zero, so don't need to do anything
        case 1 % Partial power-down of A/D
            % set D10 of reg adr 1 to one
            newVDASProfile(1, 2) = bitset(newVDASProfile(1, 2), 11);
        case 2 % Complete power-down of A/D
            % set D0 of reg adr 1 to one
            newVDASProfile(1, 2) = bitset(newVDASProfile(1, 2), 1);
        otherwise
            error('RcvProfile: PowerDownA2D must be 0, 1, or 2.');
    end

    %% PowerDownPreamp
    if ~isfield(RcvProfile(i),'PowerDownPreamp') || isempty(RcvProfile(i).PowerDownPreamp)
        RcvProfile(i).PowerDownPreamp = 0; % default value is powered up
    end
    switch RcvProfile(i).PowerDownPreamp
        case 0 % Power Up A/D
            % set D14 and D15 of reg adr 53 to zero, so don't need to do anything
        case 1 % Partial power-down of A/D
            % set D14 of reg adr 53 to one
            newVDASProfile(7, 2) = bitset(newVDASProfile(7, 2), 15);
        case 2 % Complete power-down of A/D
            % set D15 of reg adr 53 to one
            newVDASProfile(7, 2) = bitset(newVDASProfile(7, 2), 16);
        otherwise
            error('RcvProfile: PowerDownPreamp must be 0, 1, or 2.');
    end

    %% DCsubtract
    if ~isfield(RcvProfile(i),'DCsubtract') || isempty(RcvProfile(i).DCsubtract)
        RcvProfile(i).DCsubtract = 'on'; % default value
    end
    switch RcvProfile(i).DCsubtract
        case 'off'
            % set D8 of reg adr 3 to zero, so don't need to do anything
        case 'on'
            % set D8 of reg adr 3 to one
            newVDASProfile(2, 2) = bitset(newVDASProfile(2, 2), 9);
        otherwise
            error('RcvProfile: DCsubtract must be ''on'' or ''off''.');
    end

    %% AntiAliasCutoff
    % available cutoff frequencies differ for each AFEtype.  The lowest
    % available cutoff that is equal to or above the frequency set by the
    % user will be applied.  The RcvProfile.AntiAliasCutoff value will be
    % changed to the actual setting being used, and user will be informed
    % if the value has changed and verbose >1
    if ~isfield(RcvProfile(i),'AntiAliasCutoff') || isempty(RcvProfile(i).AntiAliasCutoff)
        % if not spoecified by user, set default based to the upper limit of
        % Trans.Bandwidth.
        cutoffRequested = evalin('base', 'Trans.Bandwidth(2)');
        setbyuser = 0;
    else
        % user set a value, so copy that into cutoffRequested
        cutoffRequested = RcvProfile(i).AntiAliasCutoff;
        setbyuser = 1;
    end

    switch AFEtype
        case 0 % AFE5808A
            % AFE5808A or equivalent, with only four settings
            if cutoffRequested <= 10
                RcvProfile(i).AntiAliasCutoff = 10;
                code = 8; % 10 MHz cutoff setting
            elseif cutoffRequested <= 15
                RcvProfile(i).AntiAliasCutoff = 15;
                code = 0; % 15 MHz cutoff setting
            elseif cutoffRequested <= 20
                RcvProfile(i).AntiAliasCutoff = 20;
                code = 4; % 20 MHz cutoff setting
            else
                % for any requested cutoff greater than 20
                RcvProfile(i).AntiAliasCutoff = 30;
                code = 6; % 30 MHz cutoff setting
            end
        case 1 % PAFE5812 (early pre-production version)
            % PAFE5812 adds 40 and 60 MHz cutoff settings
            if cutoffRequested <= 10
                RcvProfile(i).AntiAliasCutoff = 10;
                code = 8; % 10 MHz cutoff setting
            elseif cutoffRequested <= 15
                RcvProfile(i).AntiAliasCutoff = 15;
                code = 0; % 15 MHz cutoff setting
            elseif cutoffRequested <= 20
                RcvProfile(i).AntiAliasCutoff = 20;
                code = 4; % 20 MHz cutoff setting
            elseif cutoffRequested <= 30
                RcvProfile(i).AntiAliasCutoff = 30;
                code = 6; % 30 MHz cutoff setting
            elseif cutoffRequested <= 40
                RcvProfile(i).AntiAliasCutoff = 40;
                code = 5; % 40 MHz cutoff setting
            else
                % for any requested cutoff greater than 40
                RcvProfile(i).AntiAliasCutoff = 60;
                code = 7; % 60 MHz cutoff setting
            end
        case 2 % AFE5812 (released version)
            % AFE5812 adds 5, 35, and 50 MHz cutoff settings
            if cutoffRequested <= 5
                RcvProfile(i).AntiAliasCutoff = 5;
                % 5 MHz AAF setting requires the 10 MHz setting in register
                % 51, and also setting bit 14 of register address 61 (for
                % all other AAF settings this bit is left in its default
                % zero state).
                code = 8; % 10 MHz cutoff setting for register 51
                % set D14 of reg adr 61 to one to enable 5 MHz AAF
                newVDASProfile(9, 2) = bitset(newVDASProfile(9, 2), 15);
            elseif cutoffRequested <= 10
                RcvProfile(i).AntiAliasCutoff = 10;
                code = 8; % 10 MHz cutoff setting
            elseif cutoffRequested <= 15
                RcvProfile(i).AntiAliasCutoff = 15;
                code = 0; % 15 MHz cutoff setting
            elseif cutoffRequested <= 20
                RcvProfile(i).AntiAliasCutoff = 20;
                code = 4; % 20 MHz cutoff setting
            elseif cutoffRequested <= 30
                RcvProfile(i).AntiAliasCutoff = 30;
                code = 6; % 30 MHz cutoff setting
            elseif cutoffRequested <= 35
                RcvProfile(i).AntiAliasCutoff = 35;
                code = 5; % 40 MHz cutoff setting
            else
                % for any requested cutoff greater than 35
                RcvProfile(i).AntiAliasCutoff = 50;
                code = 7; % 50 MHz cutoff setting
            end
    end
    % map the 4-bit code to bits 3:0 of register adr 51
    newVDASProfile(5, 2) = code + newVDASProfile(5, 2);
    if RcvProfile(i).AntiAliasCutoff ~= cutoffRequested
        % notify the user of change in value, if verbose>1
        if evalin('base', 'Resource.Parameters.verbose > 1') && setbyuser
            fprintf('\ncomputeRcvProfile status: RcvProfile(%i).AntiAliasCutoff has been changed from %g to %g MHz,\n',...
                i, cutoffRequested, RcvProfile(i).AntiAliasCutoff);
            fprintf('        the closest setting supported by the hardware system with RXindex value of %i.\n', RXindex);
        end
    end

    %% PgaInputClamp
    % note that this field is only supported with production AFE5812.
    % Attempting to enable this clamp on any other configuration will cause
    % an error.
    if AFEtype == 2
        if ~isfield(RcvProfile(i),'PgaInputClamp') || isempty(RcvProfile(i).PgaInputClamp)
            RcvProfile(i).PgaInputClamp = 'off'; % default value
        end
        switch RcvProfile(i).PgaInputClamp
            case 'off'
                % set D13 of reg adr 61 to zero, so don't need to do anything
            case 'on'
                % set D13 of reg adr 61 to one
                newVDASProfile(9, 2) = bitset(newVDASProfile(9, 2), 14);
            otherwise
                error('RcvProfile: PgaInputClamp must be ''on'' or ''off''.');
        end
    elseif isfield(RcvProfile(i),'PgaInputClamp') && strcmpi(RcvProfile(i).PgaInputClamp, 'on')
        error('RcvProfile: PgaInputClamp cannot be used with Resource.SysConfig.RXindex value of %i.\n', RXindex);
    end

    %% PgaOutputClamp
    if ~isfield(RcvProfile(i),'PgaOutputClamp') || isempty(RcvProfile(i).PgaOutputClamp)
        RcvProfile(i).PgaOutputClamp = 2; % default value
    end
    switch RcvProfile(i).PgaOutputClamp
        case 0
            code = 4; % off state
        case 1
            code = 2; % 0 db compression threshold
        case 2
            code = 0; % -2 db compression threshold
        case 3
            switch AFEtype
                case 0 % AFE5808A
                    code = 1; % -4 db compression threshold
                case 1 % PAFE5812 (early pre-production version)
                    error('RcvProfile: PgaOutputClamp must be 0, 1, or 2 for the Resource.SysConfig.RXindex value of %i.\n', RXindex);
                case 2 % AFE5812 (released version)
                    % For released 5812 this is a -6 dB compression
                    % threshold, requiring register 61 bit 15 to be set
                    % while bits 7:5 of register 51 are 000.  This
                    % compression level is only available when PGA gain is
                    % set to 24 dB.
                    code = 0; % clear bits 7:5 of adr 51
                    % set D15 of reg adr 61 to one to enable -6 dB level
                    newVDASProfile(9, 2) = bitset(newVDASProfile(9, 2), 16);
            end
        otherwise
            error('RcvProfile: PgaOutputClamp must be 0, 1, 2, or 3.');
    end
    % map the 3-bit code to bits 7:5 of register adr 51
    newVDASProfile(5, 2) = 32*code + newVDASProfile(5, 2);


    %% PgaGain
    if ~isfield(RcvProfile(i),'PgaGain') || isempty(RcvProfile(i).PgaGain)
        RcvProfile(i).PgaGain = 24; % default value
    end
    switch RcvProfile(i).PgaGain
        case 24
            % set D13 of reg adr 51 to zero, so don't need to do anything
        case 30
            % set D13 of reg adr 51 to one, but not allowed when using -6dB
            % compression level
            if AFEtype == 2 && RcvProfile(i).PgaOutputClamp == 3
                error('RcvProfile: PgaGain must be 24 dB when using -6 dB PgaOutputClamp.');
            end
            newVDASProfile(5, 2) = bitset(newVDASProfile(5, 2), 14);
        otherwise
            error('RcvProfile: PgaGain must be 24 or 30.');
    end


    %% PgaHPF
    if ~isfield(RcvProfile(i),'PgaHPF') || isempty(RcvProfile(i).PgaHPF)
        RcvProfile(i).PgaHPF = 80; % default value
    end
    switch RcvProfile(i).PgaHPF
        case 80 % 80 KHz
            % This is the default state with reg 51[4] set to zero, so
            % don't need to do anything
        case 0   % integrator disabled; will pass DC in this state
            % Need to set bit 4 of register 51 to disable the PGA
            % integrator:
            newVDASProfile(5, 2) = 16 + newVDASProfile(5, 2);
        otherwise
            error('RcvProfile: PgaHPF must be 0(disabled) or 80 KHz.');
    end


    %% digHighPassFilter
    if ~isfield(RcvProfile(i),'digHighPassFilter') || isempty(RcvProfile(i).digHighPassFilter)
        RcvProfile(i).digHighPassFilter = 0; % default state is disabled
    end
    switch RcvProfile(i).digHighPassFilter
        case 0
            % disabled; set bits 4:0 to zero, so don't need to do anything
        case 2:10
            % set D0 to a one and D4:D1 to the specified value in registers
            % at address 21 and 33
            newVDASProfile(3, 2) = newVDASProfile(3, 2) + ...
                2*RcvProfile(i).digHighPassFilter + 1;
            newVDASProfile(4, 2) = newVDASProfile(4, 2) + ...
                2*RcvProfile(i).digHighPassFilter + 1;
        otherwise
            error('RcvProfile: digHighPassFilter value must be 0 or in the range 2:10.');
    end


    %% LnaGain
    if ~isfield(RcvProfile(i),'LnaGain') || isempty(RcvProfile(i).LnaGain)
        RcvProfile(i).LnaGain = 18; % default value
    end
    if AFEtype == 2
        % production version of AFE5812 lowest gain is 15 dB not 12
        switch RcvProfile(i).LnaGain
            case 12
                code = 2;
                RcvProfile(i).LnaGain = 15;
                if evalin('base', 'Resource.Parameters.verbose > 1')
                    disp('RcvProfile status: Changing RcvProfile LnaGain from 12 to 15,');
                    disp(' the closest supported value for this HW system configuration.');
                end
            case 15
                code = 2;
            case 18
                code = 0;
            case 24
                code = 1;
            otherwise
                error('RcvProfile: LnaGain must be 15, 18, or 24 for the Resource.SysConfig.RXindex value of %i.\n', RXindex);
        end
    else
        % AFE5808A or PAFE5812
        switch RcvProfile(i).LnaGain
            case 12
                code = 2;
            case 15
                code = 2;
                RcvProfile(i).LnaGain = 12;
                if evalin('base', 'Resource.Parameters.verbose > 1') && statusMsg == 0
                    disp('RcvProfile status: Changing RcvProfile LnaGain from 15 to 12,');
                    disp(' the closest supported value for this HW system configuration.');
                    statusMsg = 1; % don't repeat this message again
                end
            case 18
                code = 0;
            case 24
                code = 1;
            otherwise
                error('RcvProfile: LnaGain must be 12, 18, or 24 for the Resource.SysConfig.RXindex value of %i.\n', RXindex);
        end
    end
    % map the 2-bit code to bits 14:13 of register adr 52
    newVDASProfile(6, 2) = 8192*code + newVDASProfile(6, 2);


    %% LnaHPF
    if ~isfield(RcvProfile(i),'LnaHPF') || isempty(RcvProfile(i).LnaHPF)
        RcvProfile(i).LnaHPF = 150; % default value
    end
    switch RcvProfile(i).LnaHPF
        case 200 % 200 KHz: set reg 59[3,2] to [1 0]
            code = 2;
        case 150 % 150 KHz: set reg 59[3,2] to [1 1]
            code = 3;
        case 100 % 100 KHz: set reg 59[3,2] to [0 0]
            code = 0;
        case 50  %  50 KHz: set reg 59[3,2] to [0 1]
            code = 1;
        case 0   % integrator disabled; will pass DC in this state
            code = 0;
            % also need to set bit 12 of register 52 to disable the LNA
            % integrator:
            newVDASProfile(6, 2) = 4096 + newVDASProfile(6, 2);
        otherwise
            error('RcvProfile: LnaHPF must be 0(disabled), 50, 100, 150, or 200 KHz.');
    end
    % map the 2-bit code to bits 3:2 of register adr 59
    newVDASProfile(8, 2) = 4*code + newVDASProfile(8, 2);


    %% LnaZinSel
    % Selects feedback resistors to set LNA input impedance; actual
    % impedance is also a function of the LNA gain setting.  The
    % "LnaZinSel" variable must be an integer in the range 0 to 31, and
    % sets a relative impedance level where 0 gives the lowest available
    % impedance and 31 gives the highest.
    if ~isfield(RcvProfile(i),'LnaZinSel') || isempty(RcvProfile(i).LnaZinSel)
        % unspecified, so set a default value based on system configuration
        if highFreqSys
            % this is a high frequency system so set default of the high-Z
            % state to disable highpass filter.  To use the highpass
            % filter, user must explicitly set LnaZinSel in the setup
            % script.
            RcvProfile(i).LnaZinSel = 31; % default value is max
        else
            % not a high frequency system so set default of minimum Zin.
            RcvProfile(i).LnaZinSel = 0; % default value is min
        end
    end
    if RcvProfile(i).LnaZinSel >= 0 && RcvProfile(i).LnaZinSel <= 31
        % Map the 5-bit impedance index value to bits 4:0 of register adr 52.
        % The bits in register 52(4:0) control switches to select feedback
        % resistors, but in reverse order from their significance to the
        % composite resistance. So we remap the bit positions using the
        % array Zindexmap such that input impedance increases monotonically
        % with increasing values of LnaZinSel as provided by the user in
        % the RcvProfile structure.

        % Check for the VTS-622 issue of potential instability in a high
        % frequency system with LnaGain = 24 dB and LnaZinSel less than 4
        if highFreqSys && RcvProfile(i).LnaGain == 24 && RcvProfile(i).LnaZinSel < 4
            % Avoid the instability problem by changing LnaZinSel to at
            % least 4
            RcvProfile(i).LnaZinSel = 4;
            if evalin('base', 'Resource.Parameters.verbose > 0')
                % Send warning message to user
                disp('WARNING: High Frequency system with LnaGain=24 and LnaZinSel<4 may cause instability artifacts in receive data.');
                fprintf('    RcvProfile(%d).LnaZinSel is being increased to 4 to avoid this issue.\n', i);
            end
        end
        Zindexmap = [31 15 23  7 27 11 19  3 29 13 21  5 25  9 17  1 ...
                     30 14 22  6 26 10 18  2 28 12 20  4 24  8 16  0];
        newVDASProfile(6, 2) = Zindexmap(round(RcvProfile(i).LnaZinSel)+1) + newVDASProfile(6, 2);
    else
        error('RcvProfile: LnaZinSel must be in the range 0:31.');
    end

    RcvProfile(i).VDASRcvProfile = newVDASProfile;

end


end

