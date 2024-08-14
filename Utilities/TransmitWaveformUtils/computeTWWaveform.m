function [Wvfm2Wy, peak, numsamples, rc, TWout] = computeTWWaveform(TWin)
%
% Copyright 2001-2020 Verasonics, Inc. Verasonics Registered U.S. Patent and Trademark Office.
%
% Compute transmit waveform from TW.type, and associated variables for a
% specific type.
%
% This function returns the five variables identified in the function
% definition, based on the TW.type waveform definition provided by the user
% in the TW structure provided as an input argument.  Note that this
% function does not modify the TW structure in the base workspace; it is up
% to the caller to write the returned variables back to the TW structure as
% appropriate. This same function supports use by the computeTXPD and
% showTXPD functions, which do not write the results back to the TW
% structure at all.
%
% The output variables are in three categories, listed below with
% definitions for each variable:
%
% 1. Transmit Waveform definition for use by the computeTXPD function:
    % Wvfm2Wy:  a sequence of time-domain samples representing a
        % simulation of the waveform that would be produced by the transmit
        % HW after passing through the transducer's frequency response twice.
        % This is an N X 1 array of doubles, with a sample rate fixed at
        % 250MHz and N set as large as is needed to represent the actual
        % transmit waveform at that sample rate.
    % numsamples: the length N of the Wvfm2Wy vector as defined above
    % peak: the delay from the start of the Waveform to the peak
        % amplitude of the waveform envelope, in units of wavelengths of
        % Trans.frequency.
% 2. Error status to notify the calling function if this function had an
% error condition (such as missing or out-of-range input variable value).
    % rc:  return code- a value of zero if there were no errors, or a non-
        % zero value if an error occured (missing or out-of-range input
        % variable, etc.).  The returned value will be an index into the
        % TWin structure array, identifying the individual TW structure
        % within which the error was found.
% 3. Optional output structure TWout, which is a copy of the TWin input structure
% with the simulation variables (Wvfm2Wy, numsamples, peak) added, plus
% the waveform definition variable "States" for
% use by the HW system, in the States format used to program the
% transmit waveform generator.  TWout is used by the update function to
% replace the TW structure in the base workspace.
    % States: N x 2 array of doubles defining the actual transmit
        % waveform that will be produced by the HW.  Refer to system
        % documentation for detailed definition.

%% Revision History
% April 2020 VTS-1365 remove indirect "StatesIndex" mapping; add "trilevel"
%    TW.type
% Sep 1 2019 VTS-1416 TX, RX index values for Active Clamp Acq module
% June 10, 2019 changes & bug fixes for 4.1.0 release
% January 23, 2018 fix problems with 'function' and 'sampled' waveform types.
% May 6, 2016 to correct error in copy of user-specified IR2wy into IR1wy (VTS-3050
% Nov. 23, 2015 no longer exit with an error if the peak output
%    current limit or transformer flux limit is exceeded.  Instead, a warning
%    message will be sent to the user and Trans.maxHighVoltage (and
%    TPC(n).maxHighVoltage as well, if present) will be reduced to the level
%    required to allow the script to execute.
% Feb. 23, 2015 enforce minimum pulse width of 3 sysClk cycles
%    and maximum of 175 sysClk cycles (0.7 usec).  Parametric and envelope
%    waveform algorithms also modified to apply this restriction, so they
%    won't create an invalid waveform that will then result in an error.


%% Initialization
% Check to see if the full TW structure is requested by the function call
% and if so, check for VDASupdates to see if the arbwave tables are
% needed.
if nargout == 5
 returnTW = 1;
else
 returnTW = 0;
end

rc = 0;  % default return code to success; change it later if an error is encountered.
TWout = TWin; % start with a copy of the input structure

% Set dummy output values, in case function returns with an error before
% the actual values are created.
Wvfm2Wy = 0;
peak = 0;
numsamples = 1;


Resource = evalin('base', 'Resource');
% check for verbose or assign default value
if ~isfield(Resource, 'Parameters') || ~isfield(Resource.Parameters,'verbose') || isempty(Resource.Parameters.verbose)
    Resource.Parameters.verbose = 2; % default setting if field doesn't exist
end
if isfield(Resource,'VDAS') && isfield(Resource.VDAS, 'sysClk') && ~isempty(Resource.VDAS.sysClk)
    sysClk = Resource.VDAS.sysClk;  % system clock rate being used within Vantage HW system
else
    sysClk = 250;  % default to 250 MHz if Resource.VDAS value not present
end

% Get Trans structure parameters, to use in creating simulation waveforms
Trans = evalin('base','Trans');

% Verify frequency is specified and quit with an error if not
if ~isfield(Trans, 'frequency') || isempty(Trans.frequency)
    fprintf(2,'Error (computeTWWaveform): Trans.frequency must be specified in user script.\n')
    rc=1;
    return
end
rxFc=Trans.frequency; % Center frequency in MHz for receive processing
% Note that at this point we do not know whether the value of
% Trans.frequency meets the system clock rate constraints, and just use it
% as is.  If it is not a supported value and the HW system is actually
% going to be used, update(Receive) will declare an error later and will
% prevent the script from running.  For simulation only, the clock rate
% constraints are not applied

% check for a value of Trans.impedance and set default if not present
if ~isfield(Trans,'impedance') || isempty(Trans.impedance)
    % default value of 5 Ohms will allow system to run, but may excessively
    % restrict transmit voltage limit.  If VSX is being used, it will have
    % prompted the user to supply a more accurate value.
    Trans.impedance = 5;
end

% Transducer Impulse Response:  If the IR1wy and IR2wy arrays exist in the
% Trans structure, use them as is. If only one exists, use it to synthesize
% the other.  If neither exists, synthesize both from Trans.Bandwidth.  If
% Bandwidth is also missing, create it from Trans.bandwidth or a default 60
% percent bandwidth if bandwidth is also not present.  If the IR array(s)
% had to be generated, save the Trans structure back to the base workspace
% so they won't have to be re-created again on subsequent calls.
if ~isfield(Trans,'IR1wy') || isempty(Trans.IR1wy) || ~isfield(Trans,'IR2wy') || isempty(Trans.IR2wy)
    if isfield(Trans,'IR2wy') && ~isempty(Trans.IR2wy)
        % Two way response is present so just copy it into one way and
        % we're done.  The complexity of doing a more accurate synthesis of
        % one-way response from two-way is not worth the only slightly
        % improved simulation accuracy
        Trans.IR1wy = Trans.IR2wy;
    elseif isfield(Trans,'IR1wy') && ~isempty(Trans.IR1wy)
        % we have the one way response, so convolve it with itself to
        % create the two way.
        Trans.IR2wy=filter(Trans.IR1wy,1,Trans.IR1wy); % convolve with itself for 2-way response
    else
        % neither impulse response exists, so synthesize both from
        % Trans.Bandwidth
        if ~isfield(Trans,'Bandwidth') || isempty(Trans.Bandwidth)
            % create Bandwidth array from bandwidth scalar, or create
            % default at 60 percent bandwidth
            if isfield(Trans,'bandwidth') && ~isempty(Trans.bandwidth)
                Trans.Bandwidth = Trans.bandwidth * [-0.5, +0.5] + rxFc;
            else
                % default transducer bandwidth 60% of receive processing center freq.
                Trans.Bandwidth = rxFc * [0.7, 1.3];
            end
        end
        %
        % Create the one-way impulse response with a 2nd order butterworth
        % IIR bandpass filter, with bandwidth set by Trans.Bandwidth; note
        % the "butter" function uses relative frequency units where 1 is
        % the Nyquist limit.
        [bflt,aflt]=buttervs(2, Trans.Bandwidth*2/sysClk);
        % set number of samples to save in impulse response, by truncating it at a
        % period of 3.5 times the period of the net bandwidth
        numIRsamp = ceil(3.5*sysClk/(Trans.Bandwidth(2)-Trans.Bandwidth(1)));

        Trans.IR1wy = impz(bflt,aflt,numIRsamp,sysClk); % get 1-way impulse response

        Trans.IR2wy = filter(bflt,aflt,Trans.IR1wy); % convolve with itself for 2-way response
    end
    % save the new arrays we created back to the base workspace, so they
    % won't have to be re-created on future calls to computeTWWaveform
    assignin('base', 'Trans', Trans);
end

%% Define Limits

% Get frequency range option status, based on TXindex value from
% acquisition modules
if isfield(Resource,'SysConfig') && isfield(Resource.SysConfig, 'TXindex') && ~isempty(Resource.SysConfig.TXindex)
    % If TXindex value is available, use it to set frequency limits
    TXindex = Resource.SysConfig.TXindex; % will use this index later to find source impedance
    switch TXindex
        case 1
            % standard frequency configuration- no frequency options
            % installed, DA2319 xfmt
            minTXfreq = 0.4; % minimum allowed nominal frequency
            maxTXfreq = 23; % maximum allowed nominal frequency
            fluxLimit = 25; % volt-usec flux limit for transmit transformer
            % flux decay time constant is L/R (primary inductance 300 uH/
            % (FET on resistance 2.35 Ohms plus primary resistance 0.78
            % Ohms)), scaled by 4/3 for margin due to tolerance variations
            % etc.
            fluxDecayTimeConst = (4/3) * 300/(2.35 + 0.78); % exponential decay time constant in usec.
            maxPulseDur = 175; % maximum allowed active pulse duration in periods of sysClk (0.7 usec)
            invert = 0; % polarity of waveform definition is the same at probe connector
        case 2
            % high frequency configuration, orig PWB1010-1L xfmr
            minTXfreq = 1.9; % minimum allowed nominal frequency
            maxTXfreq = 42; % maximum allowed nominal frequency
            fluxLimit = 4.2; % volt-usec flux limit for transmit transformer
            % flux decay time constant is L/R (primary inductance 95 uH/
            % (FET on resistance 2.35 Ohms plus primary resistance 0.20
            % Ohms)), scaled by 4/3 for margin due to tolerance variations
            % etc.
            fluxDecayTimeConst = (4/3) * 95/(2.35 + 0.20); % exponential decay time constant in usec.
            maxPulseDur = 175; % maximum allowed active pulse duration in periods of sysClk (0.7 usec)
            invert = 1; % correct for inversion of output polarity at transformer
        case 3
            % low frequency configuration, MSD7342-105ML_1000 transformer
            % total DC resistance 15.6 Ohms, 7.6 uH leakage inductance
            % (est. typical)
            minTXfreq = 0.04; % minimum allowed nominal frequency
            maxTXfreq = 1.7; % maximum allowed nominal frequency
            fluxLimit = 317; % volt-usec flux limit for transmit transformer
            % flux decay time constant is L/R (primary inductance 1000 uH/
            % (FET on resistance 2.35 Ohms plus primary resistance 7.8
            % Ohms)), scaled by 4/3 for margin due to tolerance variations
            % etc.
            fluxDecayTimeConst = (4/3) * 1000/(2.35 + 7.8); % exponential decay time constant in usec.
            maxPulseDur = 2500; % maximum allowed active pulse duration in periods of sysClk (10.0 usec)
            invert = 0; % polarity of waveform definition is the same at probe connector
        case 4
            % high frequency configuration, new PWB1010L xfmr
            minTXfreq = 1.9; % minimum allowed nominal frequency
            maxTXfreq = 42; % maximum allowed nominal frequency
            fluxLimit = 6.25; % volt-usec flux limit for transmit transformer
            % flux decay time constant is L/R (primary inductance 780 uH/
            % (FET on resistance 2.35 Ohms plus primary resistance 0.32
            % Ohms)), scaled by 4/3 for margin due to tolerance variations
            % etc.
            fluxDecayTimeConst = (4/3) * 780/(2.35 + 0.32); % exponential decay time constant in usec.
            maxPulseDur = 175; % maximum allowed active pulse duration in periods of sysClk (0.7 usec)
            invert = 1; % correct for inversion of output polarity at transformer
        case 5
            % standard frequency configuration- Active Clamp acq. module,
            % no frequency options installed, DA2319 xfmt
            minTXfreq = 0.4; % minimum allowed nominal frequency
            maxTXfreq = 23; % maximum allowed nominal frequency
            fluxLimit = 25; % volt-usec flux limit for transmit transformer
            % flux decay time constant is L/R (primary inductance 300 uH/
            % (FET on resistance 2.35 Ohms plus primary resistance 0.78
            % Ohms)), scaled by 4/3 for margin due to tolerance variations
            % etc.
            fluxDecayTimeConst = (4/3) * 300/(2.35 + 0.78); % exponential decay time constant in usec.
            maxPulseDur = 175; % maximum allowed active pulse duration in periods of sysClk (0.7 usec)
            invert = 0; % polarity of waveform definition is the same at probe connector
        case 6
            % high frequency configuration, Active Clamp, new PWB1010L xfmr
            minTXfreq = 1.9; % minimum allowed nominal frequency
            maxTXfreq = 42; % maximum allowed nominal frequency
            fluxLimit = 6.25; % volt-usec flux limit for transmit transformer
            % flux decay time constant is L/R (primary inductance 780 uH/
            % (FET on resistance 2.35 Ohms plus primary resistance 0.32
            % Ohms)), scaled by 4/3 for margin due to tolerance variations
            % etc.
            fluxDecayTimeConst = (4/3) * 780/(2.35 + 0.32); % exponential decay time constant in usec.
            maxPulseDur = 175; % maximum allowed active pulse duration in periods of sysClk (0.7 usec)
            invert = 1; % correct for inversion of output polarity at transformer
        case 7
            % low frequency configuration, Active Clamp, MSD7342-105ML_1000 transformer
            % total DC resistance 15.6 Ohms, 7.6 uH leakage inductance
            % (est. typical)
            minTXfreq = 0.04; % minimum allowed nominal frequency
            maxTXfreq = 1.7; % maximum allowed nominal frequency
            fluxLimit = 317; % volt-usec flux limit for transmit transformer
            % flux decay time constant is L/R (primary inductance 1000 uH/
            % (FET on resistance 2.35 Ohms plus primary resistance 7.8
            % Ohms)), scaled by 4/3 for margin due to tolerance variations
            % etc.
            fluxDecayTimeConst = (4/3) * 1000/(2.35 + 7.8); % exponential decay time constant in usec.
            maxPulseDur = 2500; % maximum allowed active pulse duration in periods of sysClk (10.0 usec)
            invert = 0; % polarity of waveform definition is the same at probe connector
        otherwise
            fprintf(2, ['Error (computeTWWaveform): Unrecognized TXindex value "', Resource.SysConfig.TXindex, '".\n']);
            rc=1;
            return
    end
else
    % TXindex not available in base workspace, so assign default of full
    % range (may be in simulation mode, or doing a TXPD evaluation in setup
    % script)
    minTXfreq = 0.04; % minimum allowed nominal frequency for low frequency configuration
    maxTXfreq = 42; % maximum allowed nominal frequency for high frequency option
    fluxLimit = 100; % dummy flux limit value to prempt any flux limit restrictions
    fluxDecayTimeConst = 50; % 50 usec dummy value
    TXindex = 1; % default to standard frequency transformer for source impedance estimate
    maxPulseDur = 2500; % maximum allowed active pulse duration in periods of sysClk (10 usec)
    invert = 0; % polarity of waveform definition is the same at probe connector
end

% Minimum pulse duration limits required by the transmitter (currently
% common for all TXindex values)
minPulseDur = 3; % minimum allowed active pulse duration in periods of sysClk (12 nsec)
pulseRecover = 3; % minimum zero-state interval needed to recharge gate driver AC coupling


%% Process each TW struccture
for i=1:size(TWout, 2)
    % Check for required user-specified type field.
    if ~isfield(TWout(i),'type')  || isempty(TWout(i).type)
        fprintf(2, 'Error (computeTWWaveform): Required TW.type attribute is missing.\n');
        rc=i;
        return
    end

    % Set a flag to indicate that a States array has been generated by the
    % TW.type being used.  TW types that do not provide a States array for
    % HW operation will set this flag to zero
    statesArrayExists = 1;

    switch TWout(i).type
        case 'parametric'
            %%        TW.type = 'parametric'
            if ~isfield(TWout(i),'Parameters')
                fprintf(2, 'Error (computeTWWaveform): TW.Parameters required with TW.type=''parametric''.\n');
                rc=i;
                return
            end

            % assign default values if needed
            if ~isfield(TWout(i),'equalize') || isempty(TWout(i).equalize)
                TWout(i).equalize = 1;
            else
                % make sure its an integer
                TWout(i).equalize = round(TWout(i).equalize);
            end

            % check for valid equalize value
            if TWout(i).equalize > 2 || TWout(i).equalize < 0
                fprintf(2, 'Error (computeTWWaveform): unrecognized value for TW.equalize.\n');
                rc=i;
                return
            end

            % now make a copy of the Parameters array, for range checking
            % and translation into Gen3 system clock periodunits
            G3Param = TWout(i).Parameters;

            % check for valid values of TW "D" parameter (which is the same
            % for either V1 or Gen3 versions of the parameters array).
            if any((G3Param(:,4)~=1) & (G3Param(:,4)~=-1))
                fprintf(2, 'Error (computeTWWaveform): Invalid value for TW Parametric Parameter D.  Must be one or minus one.\n');
                rc=i;
                return
            end

            % now check the "TWB" parameter value to see if we are doing V1
            % backward compatible translation or not.  For a V1 Parameters
            % array TWB will always be greater than 1, while for Gen3 it
            % will be within the range 0:1.
            if G3Param(1,2) > 1.1
                % This is a V1 Parameters array, so convert it to the
                % G3 clock rate after applying V1 range checks

                % check for twClock field and assign default value if not
                % found or empty
                if ~isfield(TWout(i),'twClock') || isempty(TWout(i).twClock)
                    TWout(i).twClock = 180;
                end

                % check for extendBL fiels and assign default if needed
                if ~isfield(TWout(i),'extendBL') || isempty(TWout(i).extendBL)
                    TWout(i).extendBL = 0;
                end

                % Check for supported value of twClock
                if TWout(i).twClock ~= 180 && TWout(i).twClock ~= 90 && TWout(i).twClock ~= 45
                    % we have an error condition due to an illegal value for
                    % twClock.
                    fprintf(2, 'Error (computeTWWaveform): unrecognized TW.twClock value of %d MHz.\n', TWout(i).twClock);
                    rc=i;
                    return
                end

                % First check for valid values of A, B, C, D and modify B and C
                % to adjust B scaling and add equalization pulses per VDAS HW spec
                if any((G3Param(:,1)<5) | (G3Param(:,1)>63))
                    fprintf(2, 'Error (computeTWWaveform): Invalid value for V1-format TW Parametric Parameter A.  Must be in the range [5 : 63].\n');
                    rc=i;
                    return
                end

                G3Param(:,2) = G3Param(:,2) + G3Param(:,2) .* (G3Param(:,1)>31);  % B scaled up by two if MSB of A is one.
                if any((G3Param(:,2)<2) | (G3Param(:,2)>G3Param(:,1)))
                    fprintf(2, 'Error (computeTWWaveform): Invalid value for V1-format TW Parametric Parameter B.  Must be in the range [2 : A] after rescaling.\n');
                    rc=i;
                    return
                end

                if any((G3Param(:,3)<0) | (G3Param(:,3)>31))
                    fprintf(2, 'Error (computeTWWaveform): Invalid value for V1-format TW Parametric Parameter C.  Must be in the range [0 : 31].\n');
                    rc=i;
                    return
                end

                % scale up the "C" parameter by 64 if extendBL is true
                G3Param(:,3) = G3Param(:,3) * (1 + 63*TWout(i).extendBL);


                % now scale up the A and B parameters for the 250 MHz clock
                % rate used for Gen 3 transmit
                % first find relative pulse width using V1 integer values
                G3Param(:,2) = G3Param(:,2) ./ G3Param(:,1); % real number in range 0:1

                % next scale A to 250 MHz and round
                G3Param(:,1) = round(G3Param(:,1) * sysClk/TWout(i).twClock);

                % finally, find a B value using relative pulse width, and
                % then round to integer
                G3Param(:,2) = round(G3Param(:,2) .* G3Param(:,1));

            else
                % We have a gen3 version of the Parameters specification,
                % so check for illegal values using the Gen3 limits:

                % TW "A" can be in the range from minTXfreq to maxTXfreq
                if any((G3Param(:,1)<minTXfreq) | (G3Param(:,1)>maxTXfreq)) && Resource.Parameters.verbose > 1
                    % note operation outside of this frequency range is not
                    % considered an error; just results in a status message
                    % at verbose 2 or greater
                    disp('computeTWWaveform Status: TW Parametric Parameter A is not within the system''s');
                    disp(['operating frequency range of [' num2str(minTXfreq, '%3.2f') ' : ' num2str(maxTXfreq, '%3.2f') '] MHz.']);
                end

                % TW "B" can be in the range from 0 to 1
                if any((G3Param(:,2)<0) | (G3Param(:,2)>1))
                    fprintf(2, 'Error (computeTWWaveform): Invalid value for TW Parametric Parameter B.  Must be in the range [0 : 1].\n');
                    rc=i;
                    return
                end

                % TW "C" can be in the range from 0 to 2e7 (10 second burst
                % duration at 1 MHz, or 2 seconds at 5 MHz)
                if any((G3Param(:,3)<0) | (G3Param(:,3)>2e7))
                    fprintf(2, 'Error (computeTWWaveform): Invalid value for TW Parametric Parameter C.  Must be in the range [0 : 2e7].\n');
                    rc=i;
                    return
                end

                % scale the G3Param array to integer 250 MHz clock periods
                G3Param(:,1) = round(sysClk ./ (2*G3Param(:,1))); % find the period of one half-cycle
                G3Param(:,2) = round(G3Param(:,2) .* G3Param(:,1)); % absolute pulse width

                % now check for illegal pulse widths of 1 or 2 and force to
                % 0 or 3 respectively
                for ii = 1:size(G3Param, 1)
                    if G3Param(ii, 2) == 1
                        G3Param(ii, 2) = 0;
                    elseif G3Param(ii, 2) == 2
                        G3Param(ii, 2) = 3;
                    end
                end


                % update the Parameters array with the quantized values
                TWout(i).Parameters(:, 1) = sysClk./(2*G3Param(:,1)); % actual burst freq. MHz
                TWout(i).Parameters(:, 2) = G3Param(:,2)./G3Param(:,1); % actual relative pulse width


            end

            % Now we have a qualified G3-format parameters array in
            % G3Param, and are ready to build the simulation waveform and
            % States array output for the HW.

            % Generate the States output array
            numCh = size(G3Param, 1);
            TWout(i).States = zeros(13, 2, numCh); % array long enough for any parametric waveform
            TWout(i).States(:, 1, :) = 30; % fill it with end waveform commands
            eq = TWout(i).equalize;
            for chnum = 1:numCh % do it for every row in TW.Parameters
                A = G3Param(chnum, 1); % get the ABCD values
                B = G3Param(chnum, 2);
                C = G3Param(chnum, 3);
                D = G3Param(chnum, 4);

                %% Parametric to States encoding

                % translate TW.Parameters into States array
                %   This function translates the V1-style ABCD parameters into a Gen3
                %   States array.

                % The A, B, C, D values are assumed to be in units of 250 MHz clock
                % periods, and already checked as being within range.

                % The eq input value is TW.equalize, specifying whether equalization pulses
                % are enabled or not.  Regardless of the eq value, waveform always starts
                % 1/4 of the way into the first half-cycle, preceeding the first
                % non-equalization half-cycle (leading edge of longest equalization pulse).

                    % eq = 0 disables equalization pulses. First 3/4 half-cycle outputs
                    % nothing.

                    % eq = 1 creates centered equalization pulses, centered on midpoint of
                    % each full half-cycle.

                    % eq = 2 creates truncated equalization pulses, formed by taking a full
                    % half-cycle and truncating the leading/ trailing half of the
                    % equalization pulses.  First equalization pulse always starts at
                    % center of first half-cycle and thus 1/4 halfcycle after start of
                    % waveform.

                States = zeros(1,2); % start out with a do-nothing waveform
                stIndx = 1; % row index of the States array

                if eq && B>4
                    % disable equalization pulses if full pulse is 4 or less, since minimum
                    % equalization pulse is 3
                    Peq1 = ceil(B/2); % will be 3 for minimum B of 5
                    Seq1 = -D; % sign of first eq pulse
                    if mod(C,2) == 0
                        % even number of half-cycles so equalization pulses have opposite
                        % sign and same duration so they will sum to zero
                        Peq2 = Peq1;
                        Seq2 = D;
                    else
                        % odd number of half-cycles so equalization pulses have same sign
                        % and duration must sum to B (but we violate this rule if B=5 and
                        % force both equalization pulses to be 3)
                        Peq2 = max(floor(B/2), 3);
                        Seq2 = -D;
                    end
                end

                if C == 0
                    % for zero half-cycles, output null array of end
                    % command and quit
                    States = [30, 0];
                else
                    % There is an actual waveform to generate
                    if (eq == 0) || (B<5)
                        % no equalization pulses; zero state in first row provides the added delay
                        % needed from start of waveform
                        States(stIndx, :) = [0, round(.75*A - (A-B)/2)]; % initial zero state replacing eq pulse
                        stIndx = stIndx + 1;
                    elseif eq == 2
                        % equalization pulses are same as a full half-cycle with the leading edge of
                        % first one and trailing edge of last one chopped off.  Always allow
                        % for one half of a full half-cycle duration of the first equalization pulse
                        Z1 = round(.75*A - Peq1 - (A-B)/2);
                        Z2 = A-B;
                        Z3 = A-B;
                        States(stIndx, :) = [0, Z1];
                        stIndx = stIndx + 1;
                        States(stIndx, :) = [Seq1, Peq1];
                        stIndx = stIndx + 1;
                        States(stIndx, :) = [0, (Z2-(A-B))];
                        stIndx = stIndx + 1;
                    elseif eq == 1
                        % equalization pulses centered on each half-cycle interval
                        Z1 = round(A/4 - Peq1/2);
                        Z2 = round(A - Peq1/2 -B/2);
                        Z3 = round(A - Peq2/2 - B/2);
                        States(stIndx, :) = [0, Z1];
                        stIndx = stIndx + 1;
                        States(stIndx, :) = [Seq1, Peq1];
                        stIndx = stIndx + 1;
                        States(stIndx, :) = [0, (Z2-(A-B))];
                        stIndx = stIndx + 1;
                    end

                    pSgn = D; % sign for first full pulse
                    if mod(C,2) == 1 || C==2
                        % C is odd so add the first half cycle prior to the repeating loop of
                        % full cycles
                        States(stIndx, :) = [0, (A-B)];
                        stIndx = stIndx + 1;
                        States(stIndx, :) = [pSgn, B];
                        stIndx = stIndx + 1;
                        pSgn = -pSgn;
                    end

                    if C<4
                        % there will be no repeating loop so just fill in the remaining pulses
                        for sti = 2:C
                            States(stIndx, :) = [0, (A-B)];
                            stIndx = stIndx + 1;
                            States(stIndx, :) = [pSgn, B];
                            stIndx = stIndx + 1;
                            pSgn = -pSgn;
                        end
                    else
                        % create a repeating loop for remaining pulses
                        % start with the loop command
                        States(stIndx, :) = [10, floor(C/2)]; % loop command with loop count
                        stIndx = stIndx + 1;
                        States(stIndx, :) = [0, (A-B)]; % first zero state in full cycle
                        stIndx = stIndx + 1;
                        States(stIndx, :) = [pSgn, B]; % first pulse
                        stIndx = stIndx + 1;
                        States(stIndx, :) = [0, (A-B)]; % second zero state in full cycle
                        stIndx = stIndx + 1;
                        States(stIndx, :) = [-pSgn, B]; % second pulse
                        stIndx = stIndx + 1;
                        States(stIndx, :) = [20, 1]; % loop end command with loop level
                        stIndx = stIndx + 1;
                    end

                    if eq && B>4
                        % fill in last eq pulse
                        States(stIndx, :) = [0, Z3];
                        stIndx = stIndx + 1;
                        States(stIndx, :) = [Seq2, Peq2];
                    end
                    TWout(i).States(1:size(States,1), 1:size(States,2), chnum) = States;
                end
            end

        case 'envelope'
            %%        TW.type = 'envelope'

            % The envelope waveform definition requires three variables to
            % specify the waveform:
                % TW.envNumCycles is the duration of the overall waveform
                % in whole cycles.  Each individual cycle of the waveform
                % will be a symmetric pair of positive and negative
                % half-cycle pulses of the exact same duration, ensuring no
                % DC or even harmonic content in the resulting waveform.
                % TW.envNumCycles is a scalar double which must have an
                % integer value of 1 or greater. For a per-channel envelope
                % waveform definition, all channels must have the same
                % number of cycles and thus TW.envNumCycles will still be a
                % scalar

                % TW.envFrequency is a row-vector which must be of length
                % equal to TW.envNumCycles.  Each entry is a double
                % representing the frequency for the associated cycle in
                % MHz, with an allowed range of 0.5 to 32 MHz.  For a
                % per-channel waveform a separate row should be specified
                % for each channel, so TW.envFrequency will become a 2D
                % array of size (number of channels X envNumCycles).

                % TW.envPulseWidth is a row-vector which must be of length
                % equal to TW.envNumCycles.  Each entry is a double
                % representing the relative pulse width for each half-cycle
                % of the associated transmit cycle.  The allowed range of
                % each entry is from 0.0 to 1.0.  A negative value can be
                % specified to invert the polarity of the associated cycle
                % (the first half-cycle will be positive if the
                % TW.envPulseWidth value is positive, or negative if the
                % TW.envPulseWidth value is negative).  For a
                % per-channel waveform a separate row should be specified
                % for each channel, so TW.envFrequency will become a 2D
                % array of size (number of channels X envNumCycles).
                
            if ~isfield(TWout(i),'envNumCycles')
                fprintf(2, 'Error (computeTWWaveform): TW.envNumCycles required with TW.type=''envelope''.\n');
                rc=i;
                return
            end
            if ~isfield(TWout(i),'envFrequency')
                fprintf(2, 'Error (computeTWWaveform): TW.envFrequency required with TW.type=''envelope''.\n');
                rc=i;
                return
            end
            Frequency = TWout(i).envFrequency;
            
            if ~isfield(TWout(i),'envPulseWidth')
                fprintf(2, 'Error (computeTWWaveform): TW.envPulseWidth required with TW.type=''envelope''.\n');
                rc=i;
                return
            end
            PulseWidth = TWout(i).envPulseWidth;

            % check range limits
            if TWout(i).envNumCycles < 1 || TWout(i).envNumCycles > 1e4
                fprintf(2, 'Error (computeTWWaveform): TW.envNumcycles must be in range 1:1e4 with TW.type=''envelope''.\n');
                rc=i;
                return
            end
            if Resource.Parameters.verbose > 1 && (min(Frequency(:)) < minTXfreq || max(Frequency(:)) > maxTXfreq)
                % note operation outside of this frequency range is not
                % considered an error; just results in a status message
                % at verbose 2 or greater
                disp('computeTWWaveform Status: TW.envFrequency is not within the system''s');
                disp(['operating frequency range of [' num2str(minTXfreq, '%3.2f') ' : ' num2str(maxTXfreq, '%3.2f') '] MHz.']);
            end
            if min(PulseWidth(:)) < -1.0 || max(PulseWidth(:)) > 1.0
                fprintf(2, 'Error (computeTWWaveform): envelope pulse width must be in range -1.0 to 1.0 with TW.type=''envelope''.\n');
                rc=i;
                return
            end
            if size(Frequency,2) ~= TWout(i).envNumCycles || size(PulseWidth, 2) ~= TWout(i).envNumCycles
                fprintf(2, 'Error (computeTWWaveform): size of envelope PulseWidth and Frequency arrays must equal TW.envNumcycles.\n');
                rc=i;
                return
            end
            numEnvCh = max( size(Frequency, 1), size(PulseWidth, 1)) ; % number of channels to be defined
            if numEnvCh > 1 && (size(Frequency,1) ~= Resource.Parameters.sizeApod || size(PulseWidth, 1) ~= Resource.Parameters.sizeApod)
                fprintf(2, 'Error (computeTWWaveform): envelope PulseWidth and Frequency channel count must equal Resource.Parameters.sizeApod.\n');
                rc=i;
                return
            end
            


            % synthesize the States array
            TWout(i).States = zeros(4*TWout(i).envNumCycles + 1, 2, numEnvCh); % four states rows per cycle, plus end command
            TWout(i).States(:, 1, :) = 30; % fill entire array with end commands
            
            for chnum = 1:numEnvCh
                % fill in States array for this channel
                prevZ3 = 0;
                for cyclenum = 1:TWout(i).envNumCycles
                    period = round(sysClk/(2*Frequency(chnum, cyclenum))); % period of one half-cycle in G3 system clock periods
                    pw = round(PulseWidth(chnum, cyclenum) * period); % pulse width for this cycle
                    abspw = abs(pw); % find magnitude to use in duration calculations
                    if abspw == 1
                        % force a pulsewidth of 1 to zero since minimum pulse
                        % width is 3
                        abspw = 0;
                        pw = 0;
                    elseif abspw == 2
                        % force pulsewidth of 2 up to 3
                        abspw = 3;
                        pw = sign(pw)*3; %preserve the sign but set magnitude to 3
                    end
                    TWout(i).States(4*(cyclenum-1)+2, :, chnum) = [sign(pw), abspw];
                    TWout(i).States(4*(cyclenum-1)+4, :, chnum) = [-sign(pw), abspw];
                    TWout(i).States(4*(cyclenum-1)+1, :, chnum) = [0, round((period-abspw)/2 + prevZ3)];
                    TWout(i).States(4*(cyclenum-1)+3, :, chnum) = [0, period-abspw];
                    prevZ3 = (period-abspw)/2;
                end
            end


        case 'trilevel'
            %%        TW.type = 'trilevel'
            if ~isfield(TWout(i),'TriLevelWvfm')
                fprintf(2, 'Error (computeTWWaveform): TW.TriLevelWvfm required with TW.type=''trilevel''.\n');
                rc=i;
                return
            end
            % for this type the user specifies an arbitrary waveform as a
            % TriLevelWaveform vector array of samples (with values +1, 0,
            % or -1) at the 250 MHz system clock sample rate.

            % Call the "convertTriLvlToStates" utility function for
            % converting the TriLevelWvfm array to the States equivalent.
            % This function will check for and report any errors found in
            % the TriLevelWvfm input array
            [ TWout(i).States, result ] = convertTriLvlToStates( TWout(i).TriLevelWvfm );
            if ~strcmp(result, 'Success')
                % error from convertTriLvlToStates so report it and quit
                fprintf(2, 'Error (computeTWWaveform): convertTriLvlToStates reported the following error:\n');
                fprintf(2, ['    ', result, '.\n']);
                rc=i;
                return
            end


        case 'pulseCode'
            %%        TW.type = 'pulseCode'
            if ~isfield(TWout(i),'PulseCode')
                fprintf(2, 'Error (computeTWWaveform): TW.PulseCode required with TW.type=''pulseCode''.\n');
                rc=i;
                return
            end
            % for this type the user specifies an arbitrary waveform in the
            % pulseCode format and all we do here is call the
            % "convertPulseCodeToStates" utility function for converting a
            % PulseCode array to the States equivalent.  This function will
            % check for and report any errors found in the PulseCode input
            % array
            [ TWout(i).States, result ] = convertPulseCodeToStates(TWout(i).PulseCode);
            if ~strcmp(result, 'Success')
                % error from convertPulseCodeToStates so report it and quit
                fprintf(2, 'Error (computeTWWaveform): convertPulseCodeToStates reported the following error:\n');
                fprintf(2, ['    ', result, '.\n']);
                rc=i;
                return
            end

        case 'states'
            %%        TW.type = 'states'
            if ~isfield(TWout(i),'States')
                fprintf(2, 'Error (computeTWWaveform): TW.States required with TW.type=''states''.\n');
                rc=i;
                return
            end
            % for this type the user specifies an arbitrary waveform in the
            % States format and all we do here is pass it on to the
            % subsequent States processing functions, which will verify
            % it is in a valid format.





        case 'function'
            %%        TW.type = 'function'
            % TW.type = 'function' is for use in simulation only.  The
            % Wvfm2Wy, peak, and numsamples variables will be returned
            % but no States array.
            statesArrayExists = 0; % this will bypass all States processing below
            if ~isfield(TWout(i),'frequency') || isempty(TWout(i).frequency)
                fprintf(2, 'Error (computeTWWaveform):TW.frequency not provided for TW.type =''function''.\n');
                rc=i;
                return
            end
            if ~isfield(TWout(i),'Descriptors') || isempty(TWout(i).Descriptors)
                fprintf(2, 'Error (computeTWWaveform): TW.Descriptors required with TW.type=''function''.\n');
                rc=i;
                return
            end
            if ~isfield(TWout(i),'numsamples') || isempty(TWout(i).numsamples)
                fprintf(2, 'Error (computeTWWaveform):TW.numsamples at 250MHz not provided for TW.type =''function''.\n');
                rc=i;
                return
            end
            numsamples = TWout(i).numsamples; % individual variable return value
            T = 0:(TWout(i).numsamples-1);
            Window = TWout(i).Descriptors(1) - TWout(i).Descriptors(2)*cos(T*2*pi/TWout(i).numsamples) + ...
                TWout(i).Descriptors(3)*cos(2*T*2*pi/TWout(i).numsamples);
            Wvfm2Wy = (TWout(i).Descriptors(4)*Window.*sin(2*pi*(T/250)*TWout(i).frequency))';
            Wvfm1Wy = Wvfm2Wy; % return a duplicate for users of 1-way waveform

        case 'sampled'
            %%        TW.type = 'sampled'
            % TW.type = 'sampled' is for use in simulation only.  The
            % Wvfm2Wy, peak, and numsamples variables will be returned
            % but no States array.
            statesArrayExists = 0; % this will bypass all States processing below
            if ~isfield(TWout(i),'numsamples') || isempty(TWout(i).numsamples)
                fprintf(2, 'Error (computeTWWaveform): TW.numsamples at 250MHz not provided for TW.type =''sampled''.\n');
                rc=i;
                return
            end
            % For the 'sampled' type, TW.Waveform at 250 samples per wavelength of Trans.frequency must be provided.
            if ~isfield(TWout(i),'Waveform') || isempty(TWout(i).Waveform)
                fprintf(2, 'Error (computeTWWaveform): TW.Waveform required with TW.type =''sampled''.\n');
                rc=i;
                return
            end
            if size(TWout(i).Waveform,2)==1
                Wvfm2Wy = TWout(i).Waveform;
            else
                Wvfm2Wy = TWout(i).Waveform';
            end
            Wvfm1Wy = Wvfm2Wy; % return a duplicate for users of 1-way waveform
            numsamples = TWout(i).numsamples;
            if numsamples ~= length(TWout(i).Waveform)
                fprintf(2, 'computeTWWaveform: TW.numsamples not equal to length of TW.Waveform for TW.type =''sampled''.\n');
                rc=i;
                return
            end
        otherwise
            fprintf(2, ['Error (computeTWWaveform): Unrecognized TW.type of "', TWout(i).type, '".\n'])
            rc=i;
            return
    end % end of switch cases over all TW.type values



    if statesArrayExists
        %% Check States array for errors and consolidate rows
        % At this point a States array has been generated.  Now before
        % doing anything with it, check for illegal pulse durations,
        % remove null rows, and combine rows where possible

        % first, verify the required array size of 2 columns
        if size(TWout(i).States, 2) ~= 2
            fprintf(2, 'Error (computeTWWaveform): TW.States array must have two columns.\n');
            rc=i;
            return
        end

        StatesIn = TWout(i).States; % copy the input array to be processed

        % Test for illegal pulse durations less than minPulseDur
        AbsStates = abs(StatesIn); % eliminate all negative values for duration checks
        activePulses = AbsStates(:, 1, :) == 1;
        AbsPulses = activePulses .* AbsStates(:, 2, :);
        shortPulses = (AbsPulses > 0 & AbsPulses < minPulseDur);
        if any(shortPulses(:))
            fprintf(2, '\nError: (computeTWWaveform)  TW(%d) States contains active pulses\n', i);
            fprintf(2, 'of duration less than allowed minimum of %d system clock periods.\n', minPulseDur);
            rc=i;
            return
        end

        StatesOut = zeros(size(StatesIn));
        StatesOut(:, 1, :) = 30; % put a "waveform end" command in every row
        % create the output States array of same size as input array and
        % filled with waveform end commands.  In addition to combining any
        % redundant rows found in StatesIn, we will eliminate null rows
        % (state with zero duration) and make sure any inactive rows
        % (beyond end of waveform) will be filled with end commands as part
        % of the process.
        numch = size(StatesIn, 3); % 3rd index selects channels
        % we have to apply the row-combining algorithm independently to each
        % individual channel in the input array
        if numch == 1
            TWout(i).perChWvfm = 0;
        elseif numch == Resource.Parameters.sizeApod
            TWout(i).perChWvfm = 1;
        else
            error('computeTWWaveform: number of waveforms must be either 1 or Resource.Parameters.sizeApod.');
        end
        numRows = size(StatesIn, 1); % 1st index selects rows

        maxOutRows = 1; % this variable will track the largest number of rows needed over all channels
        for chnum = 1:numch
            StatesCh = zeros(numRows, 2);
            outRow = 0; % number of the output row being generated
            for rownum = 1:numRows
                if AbsStates(rownum, 1, chnum) > 1
                    % this is a command row so copy it into output array as
                    % is
                    outRow = outRow + 1;
                    StatesCh(outRow,:) = StatesIn(rownum,:,chnum);
                    if StatesIn(rownum, 1, chnum) == 30
                        % this is the end of the waveform command, so
                        % break out of the rownum loop and leave the rest of this
                        % channel's StatesCh array filled with zeros
                        break
                    else
                        % for any other command, we have now moved it into
                        % StatesCh so continue to next row in the for loop
                        continue
                    end
                end
                % now we know this row is a state in a waveform segment; if
                % it has duration zero remove it
                if StatesIn(rownum,2,chnum) == 0
                    % just skip it and go on to the next row
                    continue
                end
                % now we have an active state to process
                if outRow == 0
                    % have not yet found a row that does anything, but this row
                    % does so copy it over as the first one.
                    outRow = outRow + 1;
                    StatesCh(outRow,:) = StatesIn(rownum,:,chnum);
                else
                    % at least one output row has been generated, so we can
                    % compare the previous row to the new one being
                    % processed from StatesIn
                    if StatesIn(rownum, 1,chnum) ~= StatesCh(outRow,1)
                        % this is a new State with different level than
                        % previous state, so increment outRow and copy it
                        % over
                        outRow = outRow + 1;
                        StatesCh(outRow,:) = StatesIn(rownum,:,chnum);
                    else
                        % this is a continuation of the same level as in
                        % previous row, so sum the durations into it.
                        StatesCh(outRow,2) = StatesCh(outRow,2) + StatesIn(rownum,2,chnum);
                    end
                end
            end
            % done generating waveform; add an end command if one is not
            % already present
            if ~isequal([30 0], StatesCh(outRow, :))
                outRow = outRow + 1;
                StatesCh(outRow,:) = [30, 0];
            end
            maxOutRows = max(maxOutRows, outRow); % update max row count if needed
            % StatesCh has been created; now copy it into the corresponding
            % row of StatesOut
            StatesOut(1:outRow, :, chnum) = StatesCh(1:outRow, :);
        end

        TWout(i).States = zeros(maxOutRows,2,numch); % redefine output array with correct size
        TWout(i).States = StatesOut((1:maxOutRows),:,:); % copy in the active part of StatesOut



        %% Select Channel for limits, and check other limits
        % A single simulation waveform is used for all channels.  Therefore
        % if States contains unique waveforms for each channel, we first
        % have to find the one to use.  For this purpose select the
        % waveform with longest cumulative active transmit output time;
        % this same reference waveform will be used to check transmit
        % transformer core saturation limit and DC content, and to define
        % the 'estimated average frequency' for this TW structure.

        % note the code below assumes States waveform has already been
        % pruned of zero-duration states and contiguous rows at same level


        refStchnum = 1;
        refCumOnTime = 0;
        refBdur = 0;
        TWout(i).CumOnTime = zeros(numch, 1); % total of all active pulse durations, in usec
        TWout(i).Numpulses = zeros(numch, 1);

        for ch = 1:size(TWout(i).States, 3)
            StatesCh = TWout(i).States(:, :, ch);
            % eliminate all negative values to simplify finding active
            % pulses
            AbsStatesCh = abs(StatesCh);
            cumOnTimeCh = zeros(1, 5);
            numPulsesCh = zeros(1, 5);
            BdurCh = zeros(1, 5);
            rptCount = zeros(1, 5);
            rptLvl = 1;
            % step through States array to find and account for repeat
            % commands etc.
            for rownum = 1:size(StatesCh, 1)
                if AbsStatesCh(rownum, 1) > 1
                    % this is a command row; decode the command
                    if AbsStatesCh(rownum, 1) == 30
                        % end command, so drop out of rownum loop
                        break
                    elseif AbsStatesCh(rownum, 1) == 10
                        % start of a repeating loop.  increment loop level
                        % and reset accumulators for this level to zero and
                        % record the loop count
                        rptLvl = rptLvl + 1;
                        cumOnTimeCh(rptLvl) = 0;
                        numPulsesCh(rptLvl) = 0;
                        BdurCh(rptLvl) = 0;
                        rptCount(rptLvl) = AbsStatesCh(rownum, 2);
                    elseif AbsStatesCh(rownum, 1) == 20
                        % end of a repeating loop.  Decrement loop level,
                        % and sum into the lower level accumulators the
                        % contribution from this loop and repeat count
                        rptLvl = rptLvl - 1;
                        cumOnTimeCh(rptLvl) = cumOnTimeCh(rptLvl) + rptCount(rptLvl+1)*cumOnTimeCh(rptLvl+1);
                        numPulsesCh(rptLvl) = numPulsesCh(rptLvl) + rptCount(rptLvl+1)*numPulsesCh(rptLvl+1);
                        BdurCh(rptLvl) = BdurCh(rptLvl) + rptCount(rptLvl+1)*BdurCh(rptLvl+1);
                    else
                        % unrecognized command so throw error and quit
                        error('computeTWWaveform: unrecognized command in States array for TW(%d).\n')
                    end
                else
                    % this row is a waveform state so process it
                    BdurCh(rptLvl) = BdurCh(rptLvl) + AbsStatesCh(rownum, 2);
                    if AbsStatesCh(rownum, 1) == 1
                        % this is an active pulse not a zero state
                        cumOnTimeCh(rptLvl) = cumOnTimeCh(rptLvl) + AbsStatesCh(rownum, 2);
                        numPulsesCh(rptLvl) = numPulsesCh(rptLvl) + 1;
                    end
                end
            end
            TWout(i).CumOnTime(ch) = cumOnTimeCh(1)/sysClk;
            TWout(i).Numpulses(ch) = numPulsesCh(1);
            if TWout(i).CumOnTime(ch) > refCumOnTime
                refCumOnTime = TWout(i).CumOnTime(ch);
                refStchnum = ch;
                refBdur = BdurCh(1)/sysClk;
            end
        end
        % refStchnum is now the index into the States array identifying the
        % reference channel
        TWout(i).refchnum = refStchnum;
        TWout(i).Bdur = refBdur;
        TWout(i).refRelPulseWidth = refCumOnTime/refBdur;

        %% Maximum effective pulse duration and estimated average frequency

        % Find max pulse duration and estimate a center frequency for the
        % reference waveform.  Count polarity transitions in that waveform,
        % add one and divide by two to estimate number of cycles, and then
        % divide by refBdur to estimate center frequency. Note freuqency is
        % a scalar for that particular channel; frequency is not estimated
        % on a per-channel basis.

        RefChStates = TWout(i).States(:, :, refStchnum);
        % eliminate all negative values to simplify finding active
        % pulses and commands
        AbsRefChStates = abs(RefChStates);
        cumPulse = zeros(1, 5); % used to accumulate pulse duration over multiple pulses of same polarity
        maxPulse = 0; % longest effective pulse duration we've found so far
        % accumulate transistion count through entire waveform, to determine
        % estimatedAverageFrequency
        cumXsn = zeros(1, 5); % initial value for cumulative transition count
        prevSign = zeros(1, 5); % will be the sign from last active pulse
        rptflg = zeros(1, 5);
        rptCount = zeros(1, 5);
        rptLvl = 1;
        % step through States array to find and account for repeat
        % commands etc.
        for rownum = 1:size(RefChStates, 1)
            if AbsRefChStates(rownum, 1) > 1
                % this is a command row; decode the command
                if AbsRefChStates(rownum, 1) == 30
                    % end command, so drop out of rownum loop
                    break
                elseif AbsRefChStates(rownum, 1) == 10
                    % start of a repeating loop.  increment loop level
                    % and reset accumulators for this level to zero and
                    % record the loop count
                    rptLvl = rptLvl + 1;
                    cumXsn(rptLvl) = 0;
                    cumPulse(rptLvl) = 0;
                    prevSign(rptLvl) = 0; % identifies first pulse within each segment
                    rptflg(rptLvl) = 1; % indicates repeats will apply to cumPulse
                    rptCount(rptLvl) = AbsRefChStates(rownum, 2);
                elseif AbsRefChStates(rownum, 1) == 20
                    % end of a repeating loop.  Decrement loop level,
                    % and sum into the lower level accumulators the
                    % contribution from this loop and repeat count
                    if rptflg(rptLvl) == 1
                        % this loop never reset polarity and never reset AC
                        % coupling so we have to increase duration by the
                        % repeat count
                        cumPulse(rptLvl) = cumPulse(rptLvl) * rptCount(rptLvl);
                        maxPulse = max(maxPulse, cumPulse(rptLvl));
                    end
                    % for transition count, if this loop has an odd number
                    % of transitions there will be one more each time it
                    % loops back
                    oddXsns = mod(cumXsn(rptLvl), 2); % will be one if cumXsn is odd
                    cumXsn(rptLvl) = (cumXsn(rptLvl) + oddXsns) * rptCount(rptLvl) - oddXsns; % total transitions in this loop
                    rptLvl = rptLvl - 1;
                    % pass cumPulse and prevSign down to the lower repeat
                    % level we are returning to, and add cumXsn
                    cumPulse(rptLvl) = cumPulse(rptLvl+1);
                    prevSign(rptLvl) = prevSign(rptLvl+1);
                    cumXsn(rptLvl) = cumXsn(rptLvl) + cumXsn(rptLvl+1);
                else
                    % unrecognized command so throw error and quit
                    disp(['Unrecognized command ID of ', num2str(AbsRefChStates(rownum, 1))]);
                    error('computeTWWaveform: unrecognized command in States array for TW(%d).\n', i);
                end
            else
                % this row is a waveform state so process it
                if AbsRefChStates(rownum, 1) == 0
                    % this is a zero level state; will it reset cumPulse?
                    if AbsRefChStates(rownum, 2) >= pulseRecover
                        % zero state duration is "pulseRecover" sysClk
                        % periods or more and thus will reset AC coupling so we
                        % can clear the cumPulse history and rptflg for
                        % this repeat level
                        cumPulse(1:rptLvl) = 0;
                        rptflg(rptLvl) = 0;
                    end
                else
                    % this is an active pulse not a zero state
                    if prevSign(rptLvl) == 0
                        % this is first active pulse at this repeat level
                        prevSign(rptLvl) = RefChStates(rownum, 1);
                        cumPulse(rptLvl) = RefChStates(rownum, 2);
                        if rptLvl > 1 && prevSign(rptLvl) == prevSign(rptLvl-1)
                            % cumPulse will accumulate from prior to the
                            % repeat
                            cumPulse(rptLvl) = cumPulse(rptLvl) + cumPulse(rptLvl-1);
                            maxPulse = max(maxPulse, cumPulse(rptLvl));
                        elseif rptLvl > 1 && prevSign(rptLvl) ~= prevSign(rptLvl-1)
                            % polarity change at start of this loop
                            cumXsn(rptLvl-1) = cumXsn(rptLvl-1) + 1;
                        end
                    elseif RefChStates(rownum, 1) ~= prevSign(rptLvl)
                        % new pulse of different polarity so restart
                        % cumPulse
                        cumPulse(1:rptLvl) = 0;
                        cumPulse(rptLvl) = RefChStates(rownum, 2);
                        prevSign(rptLvl) = RefChStates(rownum, 1);
                        rptflg(rptLvl) = 0;
                        cumXsn(rptLvl) = cumXsn(rptLvl) + 1; % increment transition count
                    else
                        % another pulse of same polarity so accumulate
                        cumPulse(rptLvl) = cumPulse(rptLvl) + RefChStates(rownum, 2);
                    end
                    maxPulse = max(maxPulse, cumPulse(rptLvl));
                end
            end
        end


        % all done so save the max pulse duration that was found in output
        % TW structure after converting to usec
        TWout(i).maxPulseusec = maxPulse/sysClk;

        if maxPulse > maxPulseDur
            fprintf(2, '\nError: (computeTWWaveform)  TW(%d) maximum pulse duration of %d usec exceeds %d usec limit.\n', i, TWout(i).maxPulseusec, maxPulseDur/sysClk);
            rc=i;
            return
        end

        % add one to cumXsn to get effective number of pulses; divide by two to
        % get effective number of cycles in the waveform; divide by refBdur to
        % get "average" frequency in MHz for the waveform.
        if refBdur>0
            TWout(i).estimatedAvgFreq = (cumXsn(1) + 1) / (2*refBdur);
        else
            TWout(i).estimatedAvgFreq = 0;
        end

        % check for frequency within supported range, and report a warning if
        % not:
        if Resource.Parameters.verbose && (TWout(i).estimatedAvgFreq > 0 && (TWout(i).estimatedAvgFreq < minTXfreq || TWout(i).estimatedAvgFreq > maxTXfreq))
            disp(['WARNING(computeTWWaveform): estimated average TW waveform frequency of ', num2str(TWout(i).estimatedAvgFreq, '%2.2f'), ' MHz']);
            disp(['    exceeds operating frequency range of ', num2str(minTXfreq, '%3.2f'), ' to ', num2str(maxTXfreq, '%3.2f'), ' MHz.']);
            disp('    Transmit performance may be significantly degraded outside this range.');
        end


        %% Estimate current using load impedance and frequency

        % find load impedance at estimated burst frequency
        if size(Trans.impedance, 1) == 1
            % only one row, so use impedance value from that row
            Zload = Trans.impedance(1, end);
        else
            % we have an impedance array with multiple rows
            if TWout(i).estimatedAvgFreq <= Trans.impedance(1, 1)
                % burst frequency lower than first impedance frequency so use
                % first value as is
                Zload = Trans.impedance(1, 2);
            elseif TWout(i).estimatedAvgFreq >= Trans.impedance(end, 1)
                % burst frequency greater than last impedance frequency so use
                % last value as is
                Zload = Trans.impedance(end, 2);
            else
                % find first row with higher frequency
                indxF = find(Trans.impedance(:, 1)>=TWout(i).estimatedAvgFreq, 1, 'first');
                % interpolate between the two values around TW.estimatedAvgFreq
                scale = (TWout(i).estimatedAvgFreq - Trans.impedance(indxF-1, 1))/(Trans.impedance(indxF,1) - Trans.impedance(indxF-1, 1));
                Zload = (1-scale)*Trans.impedance(indxF-1, 2) +...
                    scale*Trans.impedance(indxF, 2);
            end
        end
        TWout(i).Zload = Zload; % add the Zload value to TW structure



        % set transmit output source impedance, depending on which
        % transformer is being used (based on TXindex value), and the
        % estimated center frequency
        txFreq = TWout(i).estimatedAvgFreq; % make expressions easier to read
        switch TXindex % from earlier definition we know TXindex is in range 1:7
            case 1
                % standard frequency configuration- no frequency options
                % installed, DA2319 xfmr
                % Zsource is the equivalent of 8 Ohms in series with 1.4 uH
                TWout(i).Zsource = complex(8, 2*pi*txFreq*1.4);
            case 2
                % high frequency configuration, orig PWB1010-1L xfmr
                % Zsource is a function of the ferrite bead impedance since
                % leakage inductance of the transformer is negligible.
                % Bead impedance is approximated using following model:
                if txFreq < 10
                    TWout(i).Zsource = complex(7, 1.91*txFreq^0.98);
                elseif txFreq < 20
                    TWout(i).Zsource = complex(7, 18 + 1.2*(txFreq-10));
                else % frequencies above 20 MHz
                    TWout(i).Zsource = complex(7 + 1.37*(txFreq-20), 30 + 0.37*(txFreq-20));
                end
            case 3
                % low frequency configuration,  MSD7342-105ML_1000 transformer
                % total XFMR DC resistance 15.6 Ohms, 7.6 uH leakage inductance
                % (est. typical)
                TWout(i).Zsource = complex(22, 2*pi*txFreq*7.7);
            case 4
                % high frequency configuration, new PWB1010L xfmr
                % Zsource is a function of the ferrite bead impedance since
                % leakage inductance of the transformer is negligible.
                % Bead impedance is approximated using following model:
                txFreq = TWout(i).estimatedAvgFreq; % make expressions easier to read
                if txFreq < 10
                    TWout(i).Zsource = complex(7, 1.91*txFreq^0.98);
                elseif txFreq < 20
                    TWout(i).Zsource = complex(7, 18 + 1.2*(txFreq-10));
                else % frequencies above 20 MHz
                    TWout(i).Zsource = complex(7 + 1.37*(txFreq-20), 30 + 0.37*(txFreq-20));
                end
            case 5
                % standard frequency configuration- Active Clamp, no
                % frequency options installed, DA2319 xfmr. Zsource is the
                % equivalent of 8 Ohms in series with 1.4 uH
                TWout(i).Zsource = complex(8, 2*pi*txFreq*1.4);
            case 6
                % high frequency configuration, Active Clamp, new PWB1010L xfmr
                % Zsource is a function of the ferrite bead impedance since
                % leakage inductance of the transformer is negligible.
                % Bead impedance is approximated using following model:
                txFreq = TWout(i).estimatedAvgFreq; % make expressions easier to read
                if txFreq < 10
                    TWout(i).Zsource = complex(7, 1.91*txFreq^0.98);
                elseif txFreq < 20
                    TWout(i).Zsource = complex(7, 18 + 1.2*(txFreq-10));
                else % frequencies above 20 MHz
                    TWout(i).Zsource = complex(7 + 1.37*(txFreq-20), 30 + 0.37*(txFreq-20));
                end
            case 7
                % low frequency configuration, Active Clamp, MSD7342-105ML_1000 transformer
                % total XFMR DC resistance 15.6 Ohms, 7.6 uH leakage inductance
                % (est. typical)
                TWout(i).Zsource = complex(22, 2*pi*txFreq*7.7);
        end

        % now find total source plus load impedance, power factor, etc.
        Ztotal = TWout(i).Zsource + Zload;
        TWout(i).chIpk1V = (4/pi*sin(TWout(i).refRelPulseWidth*pi/2))/abs(Ztotal);
        % chIpk1V is amplitude of sine wave at estimatedAverageFrequency
        % based on relative pulse width of the reference waveform if it is
        % a tri-level square wave at that frequency.  Scale by HV in Volts
        % to get peak amplitude of current wavform; multiply by 1/sqrt(2)
        % to get RMS current.
        %% Check limits

        % Check transmit current amplitude against peak current limit
        chIpkLimit = 2; % peak current limit in Amps, for both high frequency and standard configurations
        chIpk = TWout(i).chIpk1V * Trans.maxHighVoltage;
        if chIpk > chIpkLimit
            % reduce maxHighVoltage as needed to stay within the limit, and
            % display warning message to the user
            newMaxHV = floor(chIpkLimit/TWout(i).chIpk1V);
            Trans.maxHighVoltage = newMaxHV; % update Trans with the new limit
            assignin('base', 'Trans', Trans); % and update Trans in base workspace as well
            % also update TPC limits if needed
            if evalin('base', 'exist(''TPC'', ''var'')')
                TPC = evalin('base', 'TPC');
                for tpcnum = 1:length(TPC)
                    if ~isempty(TPC(tpcnum).maxHighVoltage)
                        TPC(tpcnum).maxHighVoltage = min(TPC(tpcnum).maxHighVoltage, newMaxHV);
                    end
                end
                assignin('base', 'TPC', TPC);
            end
            if Resource.Parameters.verbose
                fprintf(2, '\nWarning: (computeTWWaveform)  TW(%d) transmit waveform peak output current estimate\n', i);
                fprintf(2, 'of %2.2f Amps exceeds maximum limit of %2.1f Amps.\n', chIpk, chIpkLimit);
                fprintf(2, 'Trans.maxHighVoltage is being reduced to %2.0f Volts to stay within the limit.\n\n', newMaxHV);
            end
        end

        % Evaluate and check flux limit
        % Using the reference channel States array, find the peak
        % magnitude of the integral over all individual active pulses.  This
        % will be used to limit transmit voltage to the transmitter's
        % volt-seconds maximum flux limit.  This same measure is a good
        % indication of excessive DC content at any point in the waveform, even
        % if DC content of overall waveform is zero.

        StatesFluxLim = TWout(i).States(:, :, refStchnum); % this is States array to be used for the limit check
        AbsStFlxLim = abs(StatesFluxLim);
        wvfmIntegral = 0;
        % this variable will represent the integrated pulse durations as we
        % step through the waveform.  It is initialized to zero at the start.

        wvfmIntegralPk = [0, 0];
        % will represent the maximum and minimum level of integrated pulse
        % durations over the waveform.
        % fluxDecayTimeConst is exponential decay time constant in usec.
        decayPerClk = exp(-1/(fluxDecayTimeConst * sysClk)); % decay per system clock period
        cumZeroDur = zeros(1, 5);
        loopDur = zeros(1, 5);
        loopIntegral = zeros(1, 5);
        for rownum = 1:size(StatesFluxLim, 1)
            if AbsStFlxLim(rownum,1) > 1
                % this is a command row
                if AbsStFlxLim(rownum, 1) == 30
                    % this is the end of the waveform command, so
                    % break out of the rownum loop
                    break
                elseif AbsStFlxLim(rownum, 1) == 10
                    % start of a repeating loop.  increment loop level
                    % and reset accumulators for this level to zero and
                    % record the loop count
                    rptLvl = rptLvl + 1;
                    cumZeroDur(rptLvl) = 0;
                    loopDur(rptLvl) = 0;
                    loopIntegral(rptLvl) = 0;
                    rptCount(rptLvl) = AbsStFlxLim(rownum, 2);
                elseif AbsStFlxLim(rownum, 1) == 20
                    % end of a repeating loop.
                    if rptCount(rptLvl) > 1
                        % if this loop repeats more than once, add the integral over the
                        % additional repeats (will be zero if all pulses in
                        % the loop sum to zero).
                        % First apply decay to previous value of
                        % wvfmIntegral, for the net effect of the total
                        % duration of the repeating loop (remaining repeats
                        % times total duration of the states in the row):
                        wvfmIntegral = wvfmIntegral * decayPerClk^(loopDur(rptLvl) * (rptCount(rptLvl)-1));
                        wvfmIntegral = wvfmIntegral + (loopIntegral(rptLvl) * (rptCount(rptLvl)-1));
                        if wvfmIntegral > wvfmIntegralPk(1)
                            wvfmIntegralPk(1) = wvfmIntegral; % update max if this is the biggest yet
                        elseif wvfmIntegral < wvfmIntegralPk(2)
                            wvfmIntegralPk(2) = wvfmIntegral; % update min if this is the smallest yet
                        end
                    end
                    % Decrement loop level
                    rptLvl = rptLvl - 1;
                else
                    % unrecognized command so throw error and quit
                    error('computeTWWaveform: unrecognized command in States array for TW(%d).\n')
                end
            else
                % this row is a waveform state so process it
                if AbsStFlxLim(rownum,1) == 0
                    % this is a zero state; add duration to cumZeroDur and
                    % go on to next State
                    cumZeroDur(rptLvl) = cumZeroDur(rptLvl) + AbsStFlxLim(rownum,2);
                    loopDur(rptLvl) = loopDur(rptLvl) + AbsStFlxLim(rownum,2);
                else
                    % this is an active pulse.
                    % first apply decay to previous value of wvfmIntegral, for
                    % the current pulse and any previous zero states:
                    wvfmIntegral = wvfmIntegral * decayPerClk^(cumZeroDur(rptLvl) + AbsStFlxLim(rownum,2));
                    wvfmIntegral = wvfmIntegral + StatesFluxLim(rownum, 1) * AbsStFlxLim(rownum,2); % add the signed pulse to the integral
                    if wvfmIntegral > wvfmIntegralPk(1)
                        wvfmIntegralPk(1) = wvfmIntegral; % update max if this is the biggest yet
                    elseif wvfmIntegral < wvfmIntegralPk(2)
                        wvfmIntegralPk(2) = wvfmIntegral; % update min if this is the smallest yet
                    end
                    loopDur(rptLvl) = loopDur(rptLvl) + AbsStFlxLim(rownum,2);
                    loopIntegral(rptLvl) = loopIntegral(rptLvl) + StatesFluxLim(rownum, 1) * AbsStFlxLim(rownum,2); % integral within this loop
                    cumZeroDur(rptLvl) = 0; % reset duration of preceding zero states
                end
            end
        end
        TWout(i).integralPkUsec = wvfmIntegralPk / sysClk; % convert integral peak value to units of microseconds.
        TWout(i).wvfmIntegral = wvfmIntegral; % final value at end of waveform represents DC content

        maxpk = max(abs(TWout(i).integralPkUsec));
        if maxpk > 0
            TWout(i).fluxHVlimit = ceil(fluxLimit/maxpk) + 1; % Maximum transmit voltage allowed by the flux limit
        else
            TWout(i).fluxHVlimit = 101;
        end

        %% Synthesize simulation waveform

        % We will use the States array for the channel with largest
        % cumulative pulse width to synthesize a time waveform for
        % simulation, or use the user-specified channel selection if it was
        % provided:

        if isfield(TWout(i), 'userSimChNumEnable') &&  ~isempty(TWout(i).userSimChNumEnable) && TWout(i).userSimChNumEnable == 1
            % user says they have set simulation channel number to use
            if ~isfield(TWout(i), 'simChNum') || isempty(TWout(i).simChNum)
                fprintf(2, 'Error (computeTWWaveform): User did not specify TW.simChNum.\n')
                rc=i;
                return
            end
            % copy in the user-specified value and use it
            simStChNum = TWout(i).simChNum;
        else
            % use the refStchnum waveform we already found with max cum pulse width
            simStChNum = refStchnum;
        end

        % if all transmitters disabled we generate a do-nothing simulation
        % waveform: one sample set to zero.
        if refCumOnTime == 0
            LVL(1).TriLvlWvfm_Sim = 0; % Null output waveform if no active pulses
        else
            WvfmStates = TWout(i).States(:, :, simStChNum);
            % transducer impulse response and transmit waveform are
            % generated at the system clock rate.
            absWvfmStates = abs(WvfmStates);
            clear LVL;
            LVL.TriLvlWvfm_Sim = []; % TriLvlWvfm_Sim will be the square wave waveform used for simulation
            LVL = repmat(LVL, 1, 5);
            maxRepeat = 20; % Don't allow waveform segment repeats longer than 20 cycles
            for rownum = 1:size(WvfmStates, 1) % step through each row of WvfmStates array
                if absWvfmStates(rownum,1) > 1
                    % this is a command row
                    if absWvfmStates(rownum, 1) == 30
                        % this is the end of the waveform command, so
                        % break out of the rownum loop
                        break
                    elseif absWvfmStates(rownum, 1) == 10
                        % start of a repeating loop.
                        rptLvl = rptLvl + 1;
                        LVL(rptLvl).TriLvlWvfm_Sim = [];
                        rptCount(rptLvl) = min(maxRepeat, absWvfmStates(rownum, 2));
                    elseif absWvfmStates(rownum, 1) == 20
                        % end of a repeating loop.
                        % Decrement loop level and add repeating segment
                        % from loop just ended
                        rptLvl = rptLvl - 1;
                        LVL(rptLvl).TriLvlWvfm_Sim = [LVL(rptLvl).TriLvlWvfm_Sim; repmat(LVL(rptLvl+1).TriLvlWvfm_Sim, rptCount(rptLvl+1), 1)];
                    else
                        % unrecognized command so throw error and quit
                        error('computeTWWaveform: unrecognized command in States array for TW(%d).\n')
                    end
                else
                    % this row is a waveform state so process it
                    thisState = WvfmStates(rownum, 1) * ones(WvfmStates(rownum, 2), 1);
                    LVL(rptLvl).TriLvlWvfm_Sim = [LVL(rptLvl).TriLvlWvfm_Sim; thisState];
                end
            end
        end % done with generating the simulation waveforms from the States array

        TriLvlWvfm_Sim = LVL(1).TriLvlWvfm_Sim;

        % Pad the waveform with zeros, to allow for convolution with
        % impulse response
        TriLvlWvfm_Sim = [TriLvlWvfm_Sim; zeros(length(Trans.IR2wy), 1)];

        % Generate the transmit waveform for simulation:  convolve the
        % square wave "TriLvlWvfm_Sim" as defined by the States array with
        % the transducer 2-way impulse response:
        Wvfm2Wy=filter(Trans.IR2wy,1,TriLvlWvfm_Sim);

        % find the length of the resampled waveform
        numsamples=size(Wvfm2Wy,1);

        % repeat the above steps to define the one-way simulated transmit
        % waveform:
        Wvfm1Wy=filter(Trans.IR1wy,1,TriLvlWvfm_Sim);

        % If TW output structure is required, add the separate output variables to the TWout(i) structure
        if returnTW
            TWout(i).TriLvlWvfm_Sim = TriLvlWvfm_Sim;
            TWout(i).Wvfm2Wy = Wvfm2Wy;
            TWout(i).Wvfm1Wy = Wvfm1Wy;
            TWout(i).simChNum = simStChNum;
        end
    end % this ends the processing done only if a States array has been generated




    %% Generate peak value and apply nonlinearity if specified
    if isfield(TWout(i), 'userPeakEnable') && ~isempty(TWout(i).userPeakEnable) && TWout(i).userPeakEnable==1
        % user says they have set the peak value to be used
        if ~isfield(TWout(i), 'peak') || isempty(TWout(i).peak)
            fprintf(2, 'Error (computeTWWaveform): User did not specify TW.peak.\n')
            rc=i;
            return
        end
        peak = TWout(i).peak; % copy in the user-specified value
    else
        % peak not set by user, so find the location of the peak point
        % on the waveform, by going to the middle of segment that is
        % within 3 dB of the peak.  This will take us to the middle of
        % a very long burst, even if the absolute peak is near one end
        % of this segment.
        W2=abs(hilbertvs(Wvfm2Wy(:)));
        Wvfpk=max(W2);
        k=1;
        Wvfmaxstart=k;
        while (W2(k)<.7*Wvfpk)
            Wvfmaxstart=k;
            k=k+1;
        end
        k=numsamples;
        Wvfmaxend=k;
        while (W2(k)<.7*Wvfpk)
            Wvfmaxend=k;
            k=k-1;
        end
        peak=Trans.frequency*((Wvfmaxstart+Wvfmaxend)/2-1)/250;
    end

    % --- Add nonlinear transformation here ---
    if isfield(TWout(i),'nonLinear') && ~isempty(TWout(i).nonLinear) && TWout(i).nonLinear~=0
        % the "nonlin" function is invoked only if the optional
        % TW.nonLinear variable exists and is set to a non-zero value
        Wvfm2Wy = nonlin( Wvfm2Wy, TWout(i).nonLinear);
    end

    % If TW output structure is required, add the separate output variables to the TWout(i) structure
    if returnTW
        TWout(i).Wvfm2Wy = Wvfm2Wy;
        TWout(i).Wvfm1Wy = Wvfm1Wy;
        TWout(i).numsamples = numsamples;
        TWout(i).peak = peak;
    end


end % end of the TW(i) for loop over each TW structure i

end  % end of the computeTWWaveform function


%% Heuristic Nonlinearity Transformation (for simulation only)
function [out] = nonlin (input, nl)
    % NONLIN is simple nonlinear transformation that mimics acoustic propagation through water heuristically, by effectively "accelerating"
    % compressional phases and decelerating low pressure phases. The algorithm expands (contracts) time
    % during positive (negative) waves, and resamples the waveform at the new times to achieve the desired distortion.
    % NONLIN resamples the 'input' waveform, originally sampled at 'times', advancing or delaying the samples depending on the value of
    % 'input'. For positive values, samples are fewer in density (advancing the wave) and for negative values the
    % samples are more dense (retarding the waveform). The nonlinearity parameter 'nl' defines the scale of the effect, with
    % nl=0 producing no distortion.  Empirically, values of 0 < nl < 1 produce reasonable results.

    Resource = evalin('base', 'Resource');
    if ~isfield(Resource, 'Parameters') || ~isfield(Resource.Parameters,'verbose') || isempty(Resource.Parameters.verbose)
        Resource.Parameters.verbose = 2; % default setting if field doesn't exist
    end
    if Resource.Parameters.verbose
        if nl<0 || nl>1
            disp( ['WARNING --- computeTWWaveform/nonlin: The nonlinearity coefficient nl= ' num2str(nl) ' produces poor results for (nl<0 || nl>1)'])
        end
    end

    times = (1:length(input))'; % simply an integer count of time points
    newtimes = times;
    out = input;%  .* exp(input/5);    % add an asymmetrical amplitude scaling
    for ii=2:length(input)
        deltat = ( 1 - 10*nl*out(ii) );  % the maximum value of TW.Wvfm2Wy is approximately 1
        newtimes(ii) = times(ii-1) + deltat; % modify the original sample time for each value
    end

    out = interp1 (times, out, newtimes); % now resample the original waveform at the new time points
end
