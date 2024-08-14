function [PSerror, Imeas] = extPwrCtrl(command, valueIn)

% script for remote control of AIM-TTI QPX600DP Power Supply.
% Copyright 2001-2017 Verasonics, Inc.  All world-wide rights and remedies under all intellectual property laws and industrial property laws are reserved.  Verasonics Registered U.S. Patent and Trademark Office.
%
%
% AIM-TTI provides a USB driver for Windows that will be installed when the
% Windows OS sees the unrecognized USB device connected.  This driver will
% create a virtual .com serial port that we can access from matlab using
% the matlab serial features.  The serial COM port ID number assigned by
% the OS will then have to be entered in the "serial" command listed below.
% The user can specify this ID value from their SetUp script, using the
% parameter "Resource.HIFU.extPwrComPortID".
%
% Since AIM-TTI does not provide a USB driver for MAC OS, we may be able to
% do the same thing by buying a commercial USB-to-RS232 interface that
% would connect directly to the RS232 serial port on the power supply, and
% uses a generic driver that would work on the MAC and thus would result in
% a serial port device that matlab could communicate with.  This approach
% has not yet been tested at Verasonics.
%
% The commands used to control the QPX600DP (strings of ASCII characters)
% are identical whether sent directly to the RS232 port or through the
% 'virtual serial port' using the USB physical interface.
%
% See "Serial Port I/O" in matlab help for definitions, examples, etc.:
    % matlab getting started example:
        % s = serial('COM1');
        % set(s,'BaudRate',4800);
        % fopen(s);
        % fprintf(s,'*IDN?')
        % out = fscanf(s);
        % fclose(s)
        % delete(s)
        % clear s
    % (but note that with the USB-based virtual serial port, the Baudrate
    % parameter becomes a meaningless don't care.)

% Matlab function "instrfind" is useful to find and clean up serial port
% objects that were not closed, deleted, and cleared successfully.

PSerror = 2;  % return value.  Non zero indicates error- will be set to zero if function completes successfully
% Possible return values for PSerror:
    %  0: function completed without errors
    % -2: commanded voltage exceeds voltage limits
    % -1: unrecognized command string
    %  1: initialization of external power supply failed
    %  2: abnormal return from the function (Failure to open serial port, for
    %     example)

Imeas = 0;  % measured current return value; default to zero if not actually checked

% define and open the serial port with each function call, and then close
% and delete it before returning (leaving it open doesn't seem to work
% well, or maybe I just don't know how to do that).  But closing and
% deleting it with every call has another advantage:  if VSX exits
% abnormally due to an error or other fault condition, we won't have a
% dangling serial port that could prevent the system from communicating
% with the power supply the next time VSX is run.

% TBD issue:  Seems like we should be able to use the "idvalid" command to
% determine if port was not closed successfully from a previous call,
% rather than generating an error by trying to open something that is
% already open.  But as it stands now, this function will just return with
% an error, and the user will then have to manually clear and delete the serial
% port object before they can run the system again.

% The input variable "valueIn" is used to convey the desired Voltage
% setting in units of Volts.  Before sending this value to the power
% supply, logic in this function will check it against voltage limits set
% by the system, and will also create an appropriate current limit value.
% If the actual load current exceeds the current limit setting, the power
% supply automatically switches to constant-current mode and will reduce
% the output voltage to maintain that current.

% Another important consideration is that the QPX600DP specification claims
% a 2 millisecond response time for a load change from 5% to 95% of maximum
% load current.  Since the VDAS system typically will draw a load current
% that varies dramatically over time scales much less than 2 msec, the
% power supply cannot be expected to follow those changes.  Instead, the
% power supply will produce a more-or-less constant output current over
% those time scales, and the push capacitor installed within the VDAS
% system will absorb the short term veriations in actual transmit load
% current.

% First, check for Resource.HIFU structure and determine power supply type,
% number of them, and whether supply outputs are connected in series or
% parallel
if evalin('base', 'exist(''Resource'', ''var'')')
    Resource = evalin('base', 'Resource');
    if ~isfield(Resource, 'HIFU')
        fprintf(2, 'extPwrCtrl: Could not find Resource.HIFU structure. \n');
        return
    end
else
    fprintf(2, 'extPwrCtrl: Could not find Resource structure. \n');
    return
end

% identify power supply type, or default to QPX600DP
if ~isfield(Resource.HIFU, 'psType') || isempty(Resource.HIFU.psType) || strcmp(Resource.HIFU.psType, 'QPX600DP')
    % select QPX600DP by default, if psType is undefined
    psID = 1;
    numOutputs = 2;
    PSVmax = 60; % Maximum voltage per output
elseif strcmp(Resource.HIFU.psType, 'XG40-38')
    psID = 2;
    numOutputs = 1;
    PSVmax = 40; % Maximum output voltage
else
    fprintf(2, 'extPwrCtrl: Unrecognized power supply type in Resource.HIFU.psType. \n');
    disp('Supported types are "QPX600DP" and "XG40-38".');
    return
end

% determine number of power supplies being used.
if ~isfield(Resource.HIFU, 'externalHifuPwr')
    fprintf(2, 'extPwrCtrl: Could not find Resource.HIFU.externalHifuPwr. \n');
    return
else
    numPS = Resource.HIFU.externalHifuPwr;
end

if psID == 1 && numPS~=1 && numPS~=2
    % If two, the series/ parallel setting applies to all four individual
    % outputs across both supplies for QPX600DP
    fprintf(2, 'extPwrCtrl: Unsupported value for Resource.HIFU.externalHifuPwr; must be 1 or 2. \n');
    return
elseif psID == 2 && numPS~=1
    % More than one supply is not currently supported for XG40-38
    fprintf(2, 'extPwrCtrl: Multiple XG40-38 power supplies not currently supported. \n');
    return
end


if (isfield(Resource.HIFU,'extPwrConnection') && strcmp(Resource.HIFU.extPwrConnection, 'series'))
    % we have a series connection specified by the setup script
    PSVmax = min(100, PSVmax * numOutputs * numPS); % Maximum operating output voltage (Volts) from multiple outputs in series
    vScale = 1/(numOutputs * numPS); % scaling of voltage setpoint per supply for series connection
    iScale = 1.0; % scaling of current setpoint per supply for series connection
    psMode = '3'; % series connection mode for power supply
else
    % default to parallel if series not explicitly listed in setup script
    % Previously set PSVmax is maximum operating output voltage (Volts) from each output
    vScale = 1.0; % scaling of voltage setpoint per supply for parallel connection
    iScale = 1/(numOutputs * numPS); % scaling of current setpoint per supply for parallel connection
    psMode = '4'; % parallel connection mode for power supply
end

% second step is to precondition the voltage and current limit values:

% Determine the currently active maximum voltage limit (min of all the
% following):
try % catch any errors if Trans or TPC undefined
    VLim = evalin('base', 'Trans.maxHighVoltage');
    TPC = evalin('base', 'TPC'); % get the TPC structure to check its limits
    VLim = min(VLim, min(TPC(5).maxHighVoltage, TPC(5).highVoltageLimit));
    VLim = min(VLim, PSVmax); % PSVmax is maximum allowed by power supply configuration
catch
    fprintf(2, 'extPwrCtrl: Could not find required TPC or Trans structure values. \n');
    return
end

if valueIn > VLim % if over the limit, we will return with an error code
    fprintf(2, 'External power supply voltage command exceeds allowed limits. \n');
    PSerror = -2;
    return
end

% Now scale up the commanded voltage level to allow for IR drop in power
% supply cable and system wiring:

Rdist = 0.02;
% Rdist is estimated net resistance of the wiring between power supply and
% acquisition boards, based on empirical measurement on a typical system

if isfield(TPC, 'currentLimit') && ~isempty(TPC(5).currentLimit)
    maxI = TPC(5).currentLimit;
else
    maxI = 1; % dummy value for setting Vdrop, if TPC current limit not set yet
end

Vdrop = Rdist * maxI * valueIn / VLim;
% scale TPC(5).currentLimit to reflect actual voltage setting for calculating Vdrop.

PSVset = min(PSVmax, valueIn + Vdrop);
% Add Vdrop to commanded voltage to get desired power supply output
% setting, but don't exceed PSVmax

Vset = vScale * max(PSVset, 1.6); % Vset is the desired output Voltage setting per output
% To avoid TPC problems in the VDAS system,
% never set the voltage below 1.6 Volts. (instead, disable the transmit
% pulser outputs when you want to go all the way down to zero.)

% Define Max output current limit for power supply, per output
switch psID
    case 1 % QPX600DP
        Imax = min(50, 900/Vset);
        % QPX600D max output is 600 Watts or 50 Amps per output.  We add 50% to
        % power rating for the maximum limit we will set.  This ensures power
        % supply will be at maximum output 'unreg' state before it hits our limits.
        % Thus at very high power levels we are not really imposing a current limit
        % at all; but we set the current limit parameter to protect from excessive
        % fault currents when operating well below the power supply's maximum
        % output.
    case 2 % XG40-38
        Imax = 38; % fixed maximum output current at any voltage setting
end



% Now set up the serial port to talk to the power supply:
Resource = evalin('base', 'Resource'); % Check for port ID in the HIFU structure
if (isfield(Resource,'HIFU')&&isfield(Resource.HIFU,'extPwrComPortID')&&~isempty(Resource.HIFU.extPwrComPortID))
    portID = Resource.HIFU.extPwrComPortID;
else
    fprintf(2, 'extPwrCtrl: Required field Resource.HIFU.extPwrComPortID not found or empty. \n');
    return
end
% This value depends on the particular computer and configuration in which
% the VDAS system and external supply have been installed.

if numPS == 2
    portID2 = Resource.HIFU.extPwrComPortID2;
    if strcmp(portID, portID2)
        fprintf(2, 'Second supply must have a different com port ID. \n');
        PSerror = -2;
        return
    end
end

for PSnum = 1:numPS
    if PSnum == 1
        portIDset = portID;
    else
        portIDset = portID2;
    end

    extPwrComPort = serial(portIDset);
    switch psID
        case 1 % QPX600DP
            fopen(extPwrComPort);
        case 2 % XG40-38
            waitTime = 0.08;
            set(extPwrComPort,'BaudRate',9600);
            set(extPwrComPort,'StopBits',1);
            set(extPwrComPort,'DataBits',8);
            set(extPwrComPort,'Terminator','CR');
            set(extPwrComPort,'Timeout',2);
            fopen(extPwrComPort);
            fprintf(extPwrComPort,'*ADR 1'); pause (waitTime)
%             fprintf(extPwrComPort,'*IDN?'); pause (waitTime)
%             out = fscanf(extPwrComPort);
%             disp(out);
    end



    switch upper(command)
        case('INIT')
            % Set the voltage, set the current limit, and enable the
            % outputs
            Iset = 2*iScale; % 2 Amps maximum total current for initial state
            switch psID
                case 1 % QPX600DP
                    % Note the QPX600DP actually has two independent outputs, which can
                    % be configured for either serial or parallel connection to
                    % effectively provide one output with higher power.  This function
                    % configures the supply for parallel operation, as selected from
                    % the user's setup script
                    %
                    % The INIT command first disables the outputs, then sets mode to  4
                    % or 3 (outputs in parallel or series), and then queries for the
                    % actual mode to make sure we are communicating
                    fprintf(extPwrComPort,['OPALL 0 \n CONFIG ', psMode, ' \n CONFIG? \n']);
                    PSresponse = fscanf(extPwrComPort);
                    % Note: the reply from this query evidently contains some ascii
                    % characthers folloing the mode ID number.  So we just check the
                    % first character in the response string.
                    if ~strcmp(PSresponse(1), psMode) % check if actually in mode requested
                        fprintf(2,'extPwrCtrl: Failed to initialize external power supply.\n');
                        PSerror = 1;
                    else
                        % No problem setting mode, so go ahead and set voltage and
                        % enable the outputs.
                        fprintf(extPwrComPort,['V1 ', num2str(Vset), ' \n',...
                            'I1 ', num2str(Iset), '  \n OPALL 1 \n']);
                        PSerror = 0; % clear the error flag
                    end
                case 2 % XG40-38
                    fprintf(extPwrComPort,[':CURRent:LEVel ', num2str(Iset)]); pause (waitTime)
                    fprintf(extPwrComPort,[':VOLTage:LEVel ', num2str(Vset)]);pause (waitTime)
                    fprintf(extPwrComPort,':OUTP ON');pause (waitTime)
                    % see if power supply is responding
                    fprintf(extPwrComPort,'*IDN?'); pause (waitTime)
                    out = fscanf(extPwrComPort);
                    if ~strcmp(out(1:15), 'AMETEK,XG 40-38') % check if ID is what we expect
                        fprintf(2,'extPwrCtrl: Failed to initialize external power supply.\n');
                        PSerror = 1;
                    else
                        PSerror = 0; % clear the error flag
                    end
            end

        case('SETV')
            % set the voltage and set the current limit, and then read the
            % measured output current
%             ILim = TPC(5).currentLimit * iScale;  THis field may not
%             exist? (see lines 177-181)
            ILim = maxI * iScale;
            HVLim = VLim * vScale; % scale the voltage limit to per-output values
            Iset = max(2*iScale, (0.3 + 0.7*Vset/HVLim)*ILim);  % set current limit in proportion to hv/maxhv, but no less than 2 Amps
            Iset = min(1.5*Iset, Imax); % add 50% margin for current limit, but not more than power supply max set previously

            switch psID
                case 1 % QPX600DP
                    fprintf(extPwrComPort,['V1 ', num2str(Vset), ' \n',...
                        'I1 ', num2str(Iset), ' \n']);
                    % the commented out lines below would query for actual power supply
                    % output current
            %             'I1O? \n']);
            %         Istring = fscanf(extPwrComPort);
            %         Imeas = str2num(Istring(1:4));
            %         fprintf(extPwrComPort,'I2O? \n'); % read current from both outputs
            %         Istring = fscanf(extPwrComPort);
            %         Imeas = Imeas + str2num(Istring(1:4)); % sum and return the total
                    PSerror = 0; % clear the error flag
                case 2 % XG40-38
                    fprintf(extPwrComPort,[':CURRent:LEVel ', num2str(Iset)]);pause (waitTime)
                    fprintf(extPwrComPort,[':VOLTage:LEVel ', num2str(Vset)]);pause (waitTime)
                    fprintf(extPwrComPort,':OUTP ON');pause (waitTime)
                    PSerror = 0; % clear the error flag
            end

        case('CLOSE')
            % the close command sets voltage back to minimum and disables the
            % outputs before exiting VSX
            Iset = iScale * 0.2; % 0.2 Amp total for the idle state
            switch psID
                case 1 % QPX600DP
                    fprintf(extPwrComPort,['V1 ', num2str(Vset), ' \n',...
                        'I1 ', num2str(Iset), ' \n OPALL 0 \n']);
                    PSerror = 0; % clear the error flag
                case 2 % XG40-38
                    fprintf(extPwrComPort,[':CURRent:LEVel ', num2str(Iset)]);pause (waitTime)
                    fprintf(extPwrComPort,[':VOLTage:LEVel ', num2str(Vset)]);pause (waitTime)
                    fprintf(extPwrComPort,':OUTP OFF');pause (waitTime)
                    PSerror = 0; % clear the error flag
            end
        otherwise
            fprintf(2,['Command "',command, '" not recognized by extPwrCtrl.\n']);
            PSerror = -1;
    end
    % not matter what happened, close delete and clear the port before
    % returning
    fclose(extPwrComPort);
    delete(extPwrComPort);
    clear extPwrComPort
end
