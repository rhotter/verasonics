function VsUpdate(SeqObj)
% Copyright 2001-2020 Verasonics, Inc.  Verasonics Registered U.S. Patent and Trademark Office.
%
% Updates dependent or uninitialized parameters of the referenced sequence
% object.  If VDASupdates is true, also creates/ updates all associated
% VDAS variables for the referenced object.
%   Allowed inputs: 'TW', 'TX','Receive','SeqControl'.
%   If 'TW' is specified, TX must be updated as well since some TX
%   variables are dependent on the referenced TW.  This function does not
%   do that automatically, but there is logic in VSX to ensure
%   VsUpdate(TX) is always called after a VsUpdate(TW) call.  The VSX logic
%   also avoids redundant TX calls, such as if the user placed both TW and
%   TX in their 'update&Run' command structure.
%   If 'Receive' is specified, the ReconInfo.Aperture values are
%   also updated.

% Revision history:
% May 14 2020 VTS-1365 eliminate TW.StatesIndex per-channel waveform mapping
% May 2020 VTS-1617 auto enable of new synchronous DMA scheme
% Feb 2020 VTS-1607 new TX waveform Converter
% Sep 1 2019 VTS-1416 support for "ActiveClamp" acq. module
% Apr. 2019 VTS-1147 Partial Recon support
% June 2018 VTS-827, 843 Dynamic HVMux Programming; use of
%                  Trans.Connector and HVMux.Aperture
% May 2018 VTS-804 allowing both all-element and active element Apod arrays
%          VTS-808 eliminate Resource.Parameters.numTransmit and
%                  numRcvChannels as required fields
% 8-19-2016 - Changes for 3.2.0 release including VTS-392.

% Get structures needed from base.
VDASupdates = evalin('base','VDASupdates');
Resource = evalin('base','Resource');
if VDASupdates
    Cch2Vch = Resource.VDAS.Cch2Vch;
end

if isfield(Resource,'VDAS') && isfield(Resource.VDAS, 'sysClk') && ~isempty(Resource.VDAS.sysClk)
    sysClk = Resource.VDAS.sysClk;  % system clock rate being used within Vantage HW system
else
    sysClk = 250;  % default to 250 MHz if Resource.VDAS value not present
end

% get "usingMultiSys" as initialized by VSX to indicate whether we are
% running a script for the Volume Imaging Package
usingMultiSys = evalin('base', 'usingMultiSys');

% get other variables that will be used for either TX/TW or Receive updates
Trans = evalin('base','Trans');
updateTrans = 0; % indicates Trans has been updated so return it to base workspace

numBoards = evalin('base', 'numBoards');

switch SeqObj
    %% TX, TW structure updates
    case 'TW'
        TW = evalin('base','TW');

        % update TW structure using computeTWWaveform.
        % Check each TW structure for correct attributes by calling
        % computeTWWaveform, which will also generate waveform vector for
        % simulation and TW.peak value for reconstruction
        [~, ~, ~, rc, TW] = computeTWWaveform(TW);
        if rc~=0, error('VsUpdate: error in computeTWWaveform for TW(%d).\n',rc); end

        assignin('base', 'TW', TW);

    case 'TX'
        % start of the TX update function; will be called by runAcq for
        % updates to either TW or TX
        TX = evalin('base','TX');
        TW = evalin('base','TW');

        % find lowest maxHighVoltage from Trans and TPC if it exists
        lowestMaxHV = Trans.maxHighVoltage;
        if evalin('base', 'exist(''TPC'', ''var'')')
            TPC = evalin('base','TPC');
            if isfield(TPC, 'maxHighVoltage')
                % TPC structure exists and has a maxHighVoltage field so
                % find the smallest one
                for tpcnum = 1:length(TPC)
                    if ~isempty(TPC(tpcnum).maxHighVoltage) && TPC(tpcnum).inUse
                       lowestMaxHV = min(lowestMaxHV, TPC(tpcnum).maxHighVoltage);
                    end
                end
            end
        end

        % Check for existence of required variables and array sizes.  This
        % is repeated on every VsUpdate(TX) call since the user may have
        % modified TX variables or created new entries in the TX structure
        % array.

        %% TX per-structure initialization loop
        for i = 1:size(TX,2)
            % check for presence and size of Apod and Delay
            if ~isfield(TX(i),'Apod')||isempty(TX(i).Apod)
                error('VsUpdate: ''Apod'' field is missing or empty in TX(%d).\n',i);
            end
            % Verify correct size of TX.Apod array
            if size(TX(i).Apod,2) ~= Resource.Parameters.sizeApod
                error('VsUpdate: Size of TX(%d).Apod does not match Resource.Parameters.sizeApod.\n',i);
            end
            if ~isfield(TX(i),'Delay')||isempty(TX(i).Delay)
                error('VsUpdate: ''Delay'' field is missing or empty in TX(%d).\n',i);
            end
            if size(TX(i).Delay,2) ~= Resource.Parameters.sizeApod
                error('VsUpdate: Size of TX(%d).Delay does not match Resource.Parameters.sizeApod.\n',i);
            end
            % For HVMux probe, verify that a valid aperture has been
            % selected and that it is compatible with the probe/UTA
            % capabilities and Apod definition
            if isfield(Trans,'HVMux')
                numApertures = size(Trans.HVMux.Aperture, 2); % number of Aperture tables that exist
                if ~isfield(TX(i), 'aperture') || isempty(TX(i).aperture)
                    error('VsUpdate(TX): TX(%d).aperture field has not been specified.\n', i);
                elseif TX(i).aperture > numApertures || TX(i).aperture < 1
                    error('VsUpdate(TX): TX(%d).aperture field points to nonexistent Trans.HVMux.Aperture.\n', i);
                end
                % A valid aperture index exists, so copy the identified
                % HVMux.Aperture into TXAperture
                TXAperture = Trans.HVMux.Aperture(:, TX(i).aperture); % Aperture being used by this TX
                % We have a valid Aperture table; check for elements in
                % parallel
                [~,~,ActiveCh] = find(TXAperture);
                if length(ActiveCh)~=length(unique(ActiveCh))
                    % This Aperture will select elements in parallel; check
                    % to see if that is allowed and issue an error if not
                    if ~strcmp(Trans.HVMux.type, 'perEL')
                        error('VsUpdate: TX(%d).aperture enables elements in parallel, which is not supported by HVMux being used.\n',i);
                    elseif Resource.Parameters.sizeApod < Trans.numelements
                        % active-element Apod cannot select paralleled
                        % elements
                        error('VsUpdate: TX(%d).aperture enables elements in parallel, which cannot be used with active-element Apod.\n',i);
                    else
                        % force the paralleled elements to use the Apod and
                        % Delay entries from the lowest element number of
                        % those in parallel
                        [~,FirstIndices] = ismember(TXAperture, TXAperture);
                        % For channels that appear more than once in
                        % newAperture, FirstIndices will contain the
                        % element number of the first element using that
                        % channel. Apply the Apod and Delay entries for
                        % that element to all other elements in parallel
                        % with it.
                        TX(i).Apod = TX(i).Apod(FirstIndices);
                        TX(i).Delay = TX(i).Delay(FirstIndices);
                    end
                end
            else
                % not an HVMux probe or UTA; get TXAperture from
                % Trans.Connector
                TXAperture = Trans.Connector; % Aperture mapping for this probe
            end
            if Resource.Parameters.sizeApod == Trans.numelements
                % All-element Apod is being used; now check to see if
                % TXAperture includes all nonzero entries in Apod
                LApod = (TX(i).Apod ~= 0)'; % logical column array of active elements in Apod
                LAperture = (TXAperture ~= 0); % logical column array of active elements in TXAperture
                if min(LAperture(:,1) - LApod) < 0 && Resource.Parameters.verbose > 1
                    % Apod has active entries that are not in Aperture so
                    % warn the user if verbose level 2 or greater
                    fprintf('\nVsUpdate(TX) Status: TX(%d).Apod contains active entries that are\n', i);
                    fprintf('     not supported by Trans.Connector or HVMux.Aperture.\n');
                    fprintf('     These entries will be ignored but will be left intact\n');
                    fprintf('     for potential use by Recon Processing.\n\n');
                end
            end
            if VDASupdates % don't add VDAS parameters unless requested
                if ~isfield(TX(i),'VDASBistDriveEnable')||isempty(TX(i).VDASBistDriveEnable)
                    TX(i).VDASBistDriveEnable = 0;
                end
                TWnum = TX(i).waveform; % index to TW structure used with this TX
                %% TX.Apod and TX.Delay checks and VDAS remapping

                % Sort TX.Apod into TX.VDASApod.
                numParallel = 1; % default value
                if isfield(Trans,'HVMux')
                    % for HVMux transducer, mapping is set by selected
                    % aperture
                    [ActiveEL, ~, Aper2Chnl] = find(Trans.HVMux.Aperture(:,TX(i).aperture));
                else
                    % not an HVMux transducer, so mapping is set by
                    % Trans.Connector (note VSX will have created a default
                    % Connector array with 1:1 mapping if none was provided by
                    % setup script).
                    numParallel = size(Trans.Connector, 2);
                    % if numParallel is greater than one, channels are
                    % being paralleled to drive one element
                    [ActiveEL, ~, Aper2Chnl] = find(Trans.Connector(:, 1));
                    if numParallel > 1
                        for idx = 2:numParallel
                            [~, ~, newCol] = find(Trans.Connector(:, idx));
                            Aper2Chnl = [Aper2Chnl, newCol];
                        end
                    end
                end

                % Note array size checks in VSX will already have thrown an
                % error and quit if script calls for an active aperture that is
                % incompatible with the number of boards in the system.

                % Use Cch2Vch to map connector channels into VDAS channels
                TX(i).VDASApod = zeros(1,Resource.VDAS.numTransmit);
                if TX(i).VDASBistDriveEnable && isfield(TX(i),'VDASBistDriveSelect') && ~isempty(TX(i).VDASBistDriveSelect)
                    % BIST driver is being used and we have a
                    % BistDriveSelect field, so use it to enable the
                    % specified BIST drivers and ignore TX.Apod
                    BDS = TX(i).VDASBistDriveSelect;
                    if size(BDS, 2) == 1 && numBoards > 1
                        BDS = repmat(BDS, 1, numBoards);
                    elseif size(BDS, 2) ~= numBoards
                        error('VsUpdate(TX): Length of TX.VDASBistDriveSelect does not match number of boards.');
                    end
                    for n = 1:numBoards
                        bdOffset = (n-1)*64;
                        if BDS(n) == 1 || BDS(n) == 3
                            TX(i).VDASApod(1 + bdOffset) = 1;
                        end
                        if BDS(n) == 2 || BDS(n) == 3
                            TX(i).VDASApod(33 + bdOffset) = 1;
                        end
                    end
                else
                    % Map the Apod array into VDASApod
                    for idx = 1:numParallel
                        % VTS-804 active element vs all element Apod?
                        if Resource.Parameters.sizeApod == Trans.numelements
                            % need to use only the active entries in TX Apod
                            TX(i).VDASApod(Cch2Vch(Aper2Chnl(:, idx))) = TX(i).Apod(ActiveEL);
                        else
                            % active aperture Apod so use all entries
                            TX(i).VDASApod(Cch2Vch(Aper2Chnl(:, idx))) = TX(i).Apod;
                        end
                    end
                end

                %% Generate TX.VDASStates from TW.States
                % and read in other TW variables needed
                if ~isfield(TW(TWnum),'States') || isempty(TW(TWnum).States)
                    error(['VSUpdate: ''States'' field is missing in TW(%d).  ' ...
                           'This may be a simulation-only waveform.\n'], TWnum);
                end
                StatesIn = TW(TWnum).States;
                perChWvfm = TW(TWnum).perChWvfm;
                numCHin = size(StatesIn,3);
                if perChWvfm
                    % VTS-6351 Select only the TW per channel waveforms
                    % that are active in this TX as defined by TX.Apod
                    TWCumOnTime = zeros(Resource.VDAS.numTransmit, 1);
                    TWNumpulses = zeros(Resource.VDAS.numTransmit, 1);
                    for idx = 1:numParallel
                        % VTS-804 active element vs all element Apod?
                        if Resource.Parameters.sizeApod == Trans.numelements
                            % need to use only the active entries in TX Apod
                            TWCumOnTime(Cch2Vch(Aper2Chnl(:, idx))) = TW(TWnum).CumOnTime(ActiveEL);
                            TWNumpulses(Cch2Vch(Aper2Chnl(:, idx))) = TW(TWnum).Numpulses(ActiveEL);
                        else
                            % active aperture Apod so use all entries
                            TWCumOnTime(Cch2Vch(Aper2Chnl(:, idx))) = TW(TWnum).CumOnTime;
                            TWNumpulses(Cch2Vch(Aper2Chnl(:, idx))) = TW(TWnum).Numpulses;
                        end
                    end
                else
                    % single waveform for all elements; CumOnTime and
                    % Numpulses will be scalar values.
                    TWCumOnTime = TW(TWnum).CumOnTime;
                    TWNumpulses = TW(TWnum).Numpulses;
                end
                % Create TX.VDASStates as null waveform for all VDAS
                % channels; actual States array will later be copied in for
                % all active channels
                TX(i).VDASStates = zeros(size(StatesIn, 1), 2, Resource.VDAS.numTransmit);
                TX(i).VDASStates(:, 1, :) = 30; % fill with end commands for unused channels


                % There are three approaches to generating the per-channel
                % arbwave tables used to program the HW:
                    % 1. If the user provides a per-channel array of States waveforms in
                    % the TW structure, they will be used as-is to generate the
                    % transmit waveforms with no PWM or polarity inversion
                    % based on TX.VDASApod (but TX.VDASApod will still be
                    % interpreted as a logical array, to enable/ disable
                    % transmit on a per-channel basis).  The arbwave tables
                    % will have already been provided in TW; TX just needs
                    % to reorder them into the VDAS mapping, and disable
                    % (fill with zeros) channels with zero in VDASApod.

                    % 2. If there is a single waveform in TW.States and all
                    % TX.VDASApod values are either 1 or zero, then the same
                    % transmit waveform will be used on all channels with
                    % TX.VDASApod interpreted as a logical array, to enable/
                    % disable transmit on a per-channel basis.  Arbwave
                    % table will have already been built in TW; TX just
                    % replicates it for all active channels.

                    % 3. If there is a single waveform in TW.States and
                    % TX.VDASApod contains values other than 1 or zero, then a
                    % unique waveform will be generated for each channel based
                    % on TW.States but with the polarity set by the sign of
                    % the corresponding entry in TX.VDASApod and the relative
                    % pulse width scaled by the magnitude of that entry.
                    % In this case, per-channel arbwave tables must be
                    % built here for TX using the scaled States arrays.

                if perChWvfm
                    %% per-channel TW.States array
                    TX(i).perChWvfm = 1; % unique waveform for each channel
                    % This is case 1 above: a unique waveform has been
                    % predefined by the user for each active channel.
                    % Verify correct size of TW.States array
                    if numCHin ~= Resource.Parameters.sizeApod
                        error('VsUpdate: Number of channels defined in TW structure must equal Resource.Parameters.sizeApod.');
                    end

                    % map TW.States (which has been copied into StatesIn) into TX.VDASStates
                    TX(i).CumOnTime = zeros(Resource.VDAS.numTransmit, 1);
                    TX(i).Numpulses = zeros(Resource.VDAS.numTransmit, 1);
                    for idx = 1:numParallel
                        TX(i).CumOnTime = TWCumOnTime;
                        TX(i).Numpulses = TWNumpulses;
                    end
                    TX(i).CumOnTime = TX(i).CumOnTime .* logical(TX(i).VDASApod)'; % zero out the inactive channels
                    TX(i).Numpulses = TX(i).Numpulses .* logical(TX(i).VDASApod)';
                    for elin = 1:length(ActiveEL) %count through the # of active elements
                        for idx = 1:numParallel
                            if logical(TX(i).VDASApod(Cch2Vch(Aper2Chnl(elin, idx))))  %if the channel is apodized...
                                TX(i).VDASStates(:, :, Cch2Vch(Aper2Chnl(elin, idx))) = StatesIn(:, :, ActiveEL(elin)); %... copy the Cstates
                            end
                        end
                    end
                else
                    % single States array for all active channels.
                    LApod = logical(TX(i).VDASApod); % convert to logical array
                    % Apply TWNumpulses to all active channels per TXApod
                    TX(i).Numpulses = TWNumpulses * LApod;
                    % and scale TWCumOnTime by TXApod
                    TX(i).CumOnTime = TWCumOnTime * abs(TX(i).VDASApod);

                    % Test VDASApod for values other than zero or one:
                    if any(TX(i).VDASApod~=0 & TX(i).VDASApod~=1)
                        %%%%%%%%%%%% PWM pulse width scaling needs to be applied %%%%%%%%%%%%%%%%%
                        TX(i).perChWvfm = 1; % unique waveform for each channel
                        invertWvfm = sign(TX(i).VDASApod); % get the sign (or zero state) of each Apod entry
                        scaleWvfm = abs(TX(i).VDASApod); % get the magnitude for PWM scaling
                        if max(scaleWvfm)>1
                            error('VsUpdate(TX): TX(%d) Apod array contains entries greater than one.', i);
                        end
                        % create StatesOut as the same size as StatesIn,
                        % plus enough additional rows to allow insertion of
                        % a zero state in between every active pulse and at
                        % every command
                        numCmds = sum(abs(StatesIn(:, 1, 1)) > 1); % total number of commands in the array
                        numPulses = sum(abs(StatesIn(:, 1, 1)) == 1); % total number of active pulses
                        numRowsOut = 2*numCmds + 2*numPulses + 2; % number of rows to create in StatesOut
                        numRowsIn = size(StatesIn, 1); % number of rows in StatesIn
                        % create initial VDASStates array full of null
                        % waveforms
                        TX(i).VDASStates = zeros(size(StatesIn, 1), 2, Resource.VDAS.numTransmit);
                        TX(i).VDASStates(:, 1, :) = 30; % fill with end commands for unused channels


                        % now create the StatesOut array for each channel
                        % to incorporate the PWM scaling, add it to
                        % TW.States if it is unique, and fill in
                        % TX.VDAStwStatesIndex pointers to the scaled
                        % waveforms
                        maxOutRow = 0; % track the longest array actually generated
                        for jj = 1:Resource.VDAS.numTransmit
                            StatesOut = zeros(numRowsOut, 2);
                            StatesOut(:, 1) = 30; % puts an end command in every row
                            % loop to create States array for each active
                            % channel
                            % skip this channel if it is disabled (zero entry
                            % in TX.VDASApod) and also disable if pulse width
                            % scaling is less than 0.2
                            if TX(i).VDASApod(jj) == 0 || invertWvfm(jj) == 0 || scaleWvfm(jj) < 0.2
                                continue
                            end
                            % create unscaled output array with polarity
                            % inversion applied
                            StatesInP(1:numRowsIn, :) = StatesIn(:, :, 1);
                            if invertWvfm(jj) == -1
                                StatesInP(:, 1) = -StatesInP(:, 1); % apply the inversion
                            end
                            outRowNum = 0;
                            for rownum = 1:numRowsIn
                                % loop to apply scaling to each pulse in
                                % the States array
                                if abs(StatesInP(rownum, 1)) > 1
                                    % this is a command row so just copy it
                                    % over
                                    outRowNum = outRowNum + 1;
                                    StatesOut(outRowNum, :) = StatesInP(rownum, :);
                                    if abs(StatesInP(rownum, 1)) == 30
                                        % this is an end command, so break
                                        % out of the rownum loop at this
                                        % point
                                        break
                                    end
                                elseif StatesInP(rownum, 1) == 0
                                    % this is a zero state; copy it over
                                    % unless it follows a previous zero
                                    % state then just add duration
                                    if outRowNum > 0 && StatesOut(outRowNum, 1) == 0
                                        % previous row is a zero state so
                                        % just add duration
                                        StatesOut(outRowNum, 2) = StatesOut(outRowNum, 2) + StatesInP(rownum, 2);
                                    else
                                        % copy over the new zero state
                                        outRowNum = outRowNum + 1;
                                        StatesOut(outRowNum, :) = StatesInP(rownum, :);
                                    end
                                else
                                    % this is an active pulse to be scaled
                                    % by "scaleWvfm" value; first find
                                    % number of zero clocks to add on each
                                    % side of the pulse
                                    newPw = scaleWvfm(jj)*StatesInP(rownum, 2);
                                    if newPw < 2
                                        newPw = 0;
                                    elseif newPw < 3
                                        newPw = 3;
                                    else
                                        newPw = round(newPw);
                                    end
                                    if newPw == StatesInP(rownum, 2)
                                        % pulse width will not be reduced
                                        % so just copy the row over
                                        % unchanged
                                        outRowNum = outRowNum + 1;
                                        StatesOut(outRowNum, :) = StatesInP(rownum, :);
                                    elseif newPw == 0
                                        % pulse has been reduced to zero so
                                        % change it to a zero state, or add
                                        % to duration of previous zero
                                        % state
                                        if outRowNum > 0 && StatesOut(outRowNum, 1) == 0
                                            % previous row is a zero state so
                                            % just add duration
                                            StatesOut(outRowNum, 2) = StatesOut(outRowNum, 2) + StatesInP(rownum, 2);
                                        else
                                            % copy over the new zero state
                                            outRowNum = outRowNum + 1;
                                            StatesOut(outRowNum, :) = [0, StatesInP(rownum, 2)];
                                        end
                                    else
                                        % There is a valid pulse width
                                        % scaling to be applied so create
                                        % the dZ values, add dZ1
                                        % to previous zero state (or create
                                        % one if needed), then add the
                                        % scaled active pulse row, then add
                                        % dZ2 zero state row.
                                        dZ1 = ceil((StatesInP(rownum, 2) - newPw)/2);
                                        dZ2 = floor((StatesInP(rownum, 2) - newPw)/2);
                                        if outRowNum > 0 && StatesOut(outRowNum, 1) == 0
                                            % previous row is a zero state so
                                            % just add duration
                                            StatesOut(outRowNum, 2) = StatesOut(outRowNum, 2) + dZ1;
                                        elseif dZ1 > 0
                                            % copy over the new zero state
                                            outRowNum = outRowNum + 1;
                                            StatesOut(outRowNum, :) = [0, dZ1];
                                        end
                                        % scale the pulse duration and copy
                                        % over the pulse row
                                        outRowNum = outRowNum + 1;
                                        StatesOut(outRowNum, :) = [StatesInP(rownum, 1), newPw];
                                        % add the following dZ2 zero state
                                        % if not zero
                                        if dZ2 > 0
                                            outRowNum = outRowNum + 1;
                                            StatesOut(outRowNum, :) = [0, dZ2];
                                        end
                                    end
                                end
                            end
                            % done processing the scaled states array for
                            % this channel.  Now add an end command if one
                            % is not already present
                            if StatesOut(outRowNum, 1) ~= 30
                                outRowNum = outRowNum + 1;
                                StatesOut(outRowNum, :) = [30, 0];
                            end
                            maxOutRow = max(maxOutRow, outRowNum); % find the longest StatesOut that has been created

                            if maxOutRow > numRowsIn
                                % need to add some end commands to lengthen
                                % the StatesIn array and also VDASStates
                                EndFillSt = zeros((maxOutRow - numRowsIn), 2, 1);
                                EndFillSt(:, 1, :) = 30;
                                StatesIn = [StatesIn; EndFillSt];
                                EndFill = zeros((maxOutRow - numRowsIn), 2, Resource.VDAS.numTransmit);
                                EndFill(:, 1, :) = 30;
                                TX(i).VDASStates = [TX(i).VDASStates; EndFill];
                                numRowsIn = maxOutRow;
                            end

                            % Now write StatesOut just created to its
                            % position in VDASStates
                            TX(i).VDASStates(1:outRowNum, :, jj) = StatesOut(1:outRowNum, :);
                        end
                    else
                        %%%%%%%%%%%%%%% No per-channel PWM scaling %%%%%%%%%%%%%%%%%%%%%%
                        TX(i).perChWvfm = 0; % Single waveform will be applied to all active channels
                        TX(i).VDASApod = LApod; % force Apod to be a logical
                        % copy the shared waveform into all active channels
                        for Vchnum = 1:Resource.VDAS.numTransmit
                            if LApod(Vchnum)
                                TX(i).VDASStates(:, :, Vchnum) = StatesIn;
                            end
                        end
                    end
                end % finished processing the TW waveform(s) into TX per-channel with VDAS mapping

                % Determine total TXbus power consumption for this burst,
                % to set profile 5 required threshold level
                Bdur = TW(TWnum).Bdur;
                Ztot = TW(TWnum).Zsource + TW(TWnum).Zload;

                CperBd = 18; % TXbus capacitance in uF on each Acq board
                VdroopMax = 7; % Max droop in Volts of TXbus for replenish stage protection error
                ImaxTPC = 3; % max current the TPC will supply during the burst
                % find total chage available from TPC for an imaging burst,
                % for this waveform burst duration:
                QtpcOut = CperBd * numBoards * VdroopMax + Bdur * ImaxTPC;

                % find per-channel peak current for the fundamental
                % term of Fourier series, based pn per-channel cumOnTIme
                PerChIpk1V = 4*sin(TX(i).CumOnTime*pi/(2*Bdur))/(pi*abs(Ztot));

                % find total TXBus load current over all channels during
                % this burst, including effect of power factor
                TXbusIburst1V = cos(angle(Ztot)) * sum(PerChIpk1V) / sqrt(2);
                TX(i).imgProfileMaxHv = QtpcOut / (TXbusIburst1V * Bdur);

                % Set the extended burst flag if this TX would exceed
                % imaging supply capacity
                if TX(i).imgProfileMaxHv < lowestMaxHV
                    TX(i).sysExtendBL = 1;
                elseif ~isfield(TX(i), 'sysExtendBL') || isempty(TX(i).sysExtendBL)
                    % if user had already set it to 1, leave it alone
                    TX(i).sysExtendBL = 0;
                end

            end % finished with VDASupdates- based calculations for this TX

        end % end of the loop over all TX structures


        %% TX: Create TX Descriptor tables using convertStatesToDescriptor Function
        assignin('base','TX',TX); % save what we have so far, in case of an error
        if VDASupdates
            TX = computeDescriptors(TX);
        end

        assignin('base','TX',TX);
        assignin('base','TW',TW);
        if updateTrans
            % changes have been made to Trans, so return it as well
            assignin('base','Trans',Trans);
        end

    %% Receive Structure updates
    case 'Receive'
        % first, check if a Receive structure exists and just return if not
        if evalin('base','~exist(''Receive'', ''var'')')
            return
        end
        Receive = evalin('base','Receive');

        if isempty(Receive)
            error('VsUpdate: Receive structure is empty.');
        elseif length(Receive) > 1 && length(Receive(:)) == size(Receive, 1)
            % was defined as a column vector; convert to row and warn user if verbose
            if Resource.Parameters.verbose
                fprintf(2, 'VSX warning: Receive structure array was defined as a column vector;\n');
                fprintf(2, '     converting to required row vector format.\n');
            end
            Receive = Receive';
        elseif length(Receive(:)) > 1 && length(Receive(:)) ~= size(Receive, 2)
            error('VSX: Multidimensional Receive structure array not supported.  It must be a row vector.');
        end

        % Get values needed from base workspace
        LEsys = evalin('base','LEsys'); % Flag to indicate whether we enforce the Vantage 64 rcv mux constraint

        % Check for missing required fields.
        if ~isfield(Receive,'Apod'), error('VsUpdate: Receive.Apod field missing.'); end
        if ~isfield(Receive,'startDepth'), error('VsUpdate: Receive.startDepth field missing.'); end
        if ~isfield(Receive,'endDepth'), error('VsUpdate: Receive.endDepth field missing.'); end
        if ~isfield(Receive,'TGC'), error('VsUpdate: Receive.TGC field missing.'); end
        if ~isfield(Receive,'bufnum'), error('VsUpdate: Receive.bufnum field missing.'); end
        if ~isfield(Receive,'framenum'), error('VsUpdate: Receive.framenum field missing.'); end
        if ~isfield(Receive,'acqNum'), error('VsUpdate: Receive.acqNum field missing.'); end
        if ~isfield(Receive,'mode'), error('VsUpdate: Receive.mode field missing.'); end

        % Add fields that are optional to specify. If not provided, set to empty.
        Flds = {'sampleMode','decimSampleRate','demodFrequency','samplesPerWave','ADCRate'};
        for j=1:length(Flds), if ~isfield(Receive,Flds{j}), [Receive.(Flds{j})] = deal([]); end, end

        % Add empty valued dependent fields (overwrite any user defined values).
        Flds = {'decimFactor','quadDecim','startSample','endSample'};
        for j=1:length(Flds), [Receive.(Flds{j})] = deal([]); end

        % Remove obsolete fields
        if isfield(Receive,'interleave')
            fprintf(2,'VsUpdate: Receive.interleave attribute no longer supported; removing from Receive structure.\n');
            Receive = rmfield(Receive,'interleave');
        end
        warnPl = 0; % flag to disable multiple warning messages for parallel elements in Mux aperture

        % Check for empty or incorrect required attributes in each receive structure.
        Flds = {'Apod','startDepth','endDepth','TGC','bufnum','framenum','acqNum','mode'};
        for i = 1:size(Receive,2)
            for j=1:length(Flds)
                if isempty(Receive(i).(Flds{j})), error(['VsUpdate: Receive(%d).',Flds{j},' empty.\n'],i); end
            end
            % Verify proper size of Receive.Apod array
            if Resource.Parameters.sizeApod ~= size(Receive(i).Apod,2)
                error(['VsUpdate: Receive(%d).Apod array is not sized correctly. It should be a row\n',...
                       'vector with the same length as Resource.Parameters.sizeApod, %d.\n'],i,Resource.Parameters.sizeApod);
            end
            % For HVMux probe using Dynamic Mux Programming, verify that an
            % aperture has been selected and that it supports the Apod
            % array
            if isfield(Trans,'HVMux') && Resource.Parameters.sizeApod == Trans.numelements
                LApod = (Receive(i).Apod ~= 0)'; % logical column array of active elements in Apod
                numApertures = size(Trans.HVMux.Aperture, 2); % number of Aperture tables that already exist
                if ~isfield(Receive(i), 'aperture') || isempty(Receive(i).aperture)
                    error('VsUpdate(Receive): Receive(%d).aperture field has not been specified.\n', i);
                elseif Receive(i).aperture > numApertures || Receive(i).aperture < 1
                    error('VsUpdate(Receive): Receive(%d).aperture field points to nonexistent Trans.HVMux.Aperture.');
                else
                    % a valid aperture index exists, confirm that the
                    % identified column of HVMux.Aperture includes all
                    % nonzero entries in Apod
                    RcvAperture = Trans.HVMux.Aperture(:, Receive(i).aperture);
                    LAperture = (RcvAperture ~= 0);
                    if min(LAperture - LApod) < 0
                        % Apod has active entries that are not in Aperture
                        error('VsUpdate(Receive): Receive(%d).Apod is not supported by Aperture index %d.\n', i, Receive(i).aperture);
                    end
                end
                % We have a valid Aperture table; check for elements in
                % parallel
                [~,~,ActiveCh] = find(RcvAperture);
                if length(ActiveCh)~=length(unique(ActiveCh))
                    % This Aperture will select elements in parallel; check
                    % to see if that is allowed and issue an error if not
                    if ~strcmp(Trans.HVMux.type, 'perEL')
                        error('VsUpdate: Receive(%d).Apod enables elements in parallel, which is not supported by HVMux being used.\n',i);
                    else
                        % force the paralleled elements to use the Apod
                        % entry from the lowest element number of those in
                        % parallel
                        [~,FirstIndices] = ismember(RcvAperture, RcvAperture);
                        % For channels that appear more than once in
                        % newAperture, FirstIndices will contain the
                        % element number of the first element using that
                        % channel. Apply the Apod entry for that element to
                        % all other elements in parallel with it.
                        Receive(i).Apod = Receive(i).Apod(FirstIndices);
                        if Resource.Parameters.verbose>0 && warnPl == 0
                            % issue a warning message since parallel
                            % elements on receive may be a mistake
                            fprintf(2, ['VsUpdate Warning: Receive(%d).Apod enables elements in parallel through the HVMux switches.\n', ...
                                'This warning message will not be repeated if additional parallel-Apod arrays are found.\n'], i);
                            warnPl = 1;
                        end
                    end
                end
            end
            % Check for sampleMode or samplesPerWave attributes provided; if not found, set defaults
            if isempty(Receive(i).sampleMode)
                if ~isempty(Receive(i).samplesPerWave) % for backwards compatibility
                    switch Receive(i).samplesPerWave
                        case 4
                            Receive(i).sampleMode = 'NS200BW';
                        case 2
                            Receive(i).sampleMode = 'BS100BW';
                        case 4/3
                            Receive(i).sampleMode = 'BS67BW';
                        case 1
                            Receive(i).sampleMode = 'BS50BW';
                        otherwise
                            Receive(i).sampleMode = 'custom';
                            Receive(i).decimSampleRate = Receive(i).samplesPerWave * Trans.frequency;
                    end
                else
                    Receive(i).sampleMode = 'NS200BW'; % default sampleMode if not provided.
                end
            end
        end
        clear Flds

        %% Compute table of allowed decimSampleRates and their corresponding decimation values.
        SR = zeros(100,2); % Compute 100 entries from 125 MHz to 1.25 MHz.
        % These are the allowed sample rates that meet the constraints of
        % A/D rate in range 10 to 62.5 MHz (divisor in range 4 to 25)
        % followed by a decimation factor in range 1 to 8. decimSampleRates
        % of 71.43, 83.33, 100, and 125 MHz are special cases that imply 2X
        % interleaved sampling, and thus the actual sample rate will be one
        % half of this value, with interleaving of the RF data from two
        % separate acquisitions.  Note transmit must be shifted by half of
        % the A/D sample period between the two acquisitons to achieve the
        % correct interleave.  For the 71.43 and 100 MHz sample rates this
        % would require transmit delay shifts (in periods of the 250 MHz
        % system clock) of 3.5 and 2.5 respectively- so some small
        % performance compromise will result when using the closest
        % realizable integer value delays.
        k = 1;
        for i = [2, 2.5, 3, 3.5, 4:200]
            for j=8:-1:1 % find highest decimation factor that can be used
                d = i/j;
                if (d - floor(d)) < .0001 && d > 3
                    break % use this j if d is an integer > 3
                    % otherwise continue until j=1 and d=i
                end
            end
            if (d <= 25) % This is a valid A/D rate, so add it to the list
                SR(k,1) = sysClk/i;
                SR(k,2) = j;
                k = k+1;
            end
        end

        % Loop over all Receives to set various attributes.
        warn = 1; % used to limit warnings to one time
        for i = 1:size(Receive,2)
            %% Compute ADRate, decimFactor, demodFrequency and quadDecim
            % check for special case of 'custom' sampleMode
            if strcmp(Receive(i).sampleMode,'custom')
                if isempty(Receive(i).ADCRate)
                    % ADC rate not specified so check for decimSampleRate
                    if ~isempty(Receive(i).decimSampleRate)
                        % Find closest realizable decimSampleRate, and default ADCRate with its associated decimFactor
                        [~,j] = min(abs(SR(:,1)-Receive(i).decimSampleRate));
                        j = max(5, j); % custom sample rates cannot exceed 62.5 MHz
                        Receive(i).ADCRate = SR(j,1)*SR(j,2);  % set the derived ADCRate
                        Receive(i).decimSampleRate = SR(j,1);  % set new decimSampleRate
                        Receive(i).decimFactor = SR(j,2);      % set associated decimFactor
                    else
                        error('VsUpdate(Receive): No sample rate specified for Receive sampleMode of ''custom''.');
                    end
                else
                    % ADC rate has been set; force it to a legitimate value
                    Receive(i).ADCRate = sysClk/round(sysClk/Receive(i).ADCRate); % set to nearest submultiple of sysClk
                    Receive(i).ADCRate = min(62.5, max(10, Receive(i).ADCRate)); % force to the range 10:62.5 MHz
                    if ~isempty(Receive(i).decimSampleRate)
                        % decimSampleRate has been set, so adjust it to
                        % nearest integer submultiple of ADCRate and set
                        % decimFactor to that integer
                        Receive(i).decimFactor = round(Receive(i).ADCRate/Receive(i).decimSampleRate);
                        Receive(i).decimSampleRate = Receive(i).ADCRate/Receive(i).decimFactor;
                    else
                        % decimSampleRate has not been set, so assign a
                        % default decimFactor of 1 so decimSampleRate will
                        % be set equal to ADCRate.
                        Receive(i).decimFactor = 1;
                        Receive(i).decimSampleRate = Receive(i).ADCRate;
                    end
                end
            else
                % Set target decimSampleRate - if provided, use it, otherwise use (4 or 4/3)*Trans.frequency.
                % Note that when sampleMode is not 'custom', any
                % user-specified value of ADCRate will be ignored and
                % overwritten with the derived value
                if isempty(Receive(i).decimSampleRate) && isempty(Receive(i).demodFrequency)
                    % No sample rate has been set, so assign a default
                    % value based on sampleMode and Trans.frequency
                    if strcmp(Receive(i).sampleMode,'BS67BW')
                        Receive(i).decimSampleRate = (4/3)*Trans.frequency;
                    else
                        Receive(i).decimSampleRate = 4*Trans.frequency;
                    end
                elseif isempty(Receive(i).decimSampleRate)
                    % demodFrequency has been set but not decimSampleRate
                    % so derive it from the value provided for
                    % demodFrequency
                    if strcmp(Receive(i).sampleMode,'BS67BW')
                        Receive(i).decimSampleRate = (4/3)*Receive(i).demodFrequency;
                    else
                        Receive(i).decimSampleRate = 4*Receive(i).demodFrequency;
                    end
                end
                % Find closest realizable ADCRate with its associated decimFactor
                [~,j] = min(abs(SR(:,1)-Receive(i).decimSampleRate));
                Receive(i).ADCRate = SR(j,1)*SR(j,2);  % set the derived ADCRate
                Receive(i).decimSampleRate = SR(j,1);  % set new decimSampleRate
                Receive(i).decimFactor = SR(j,2);
                if Receive(i).ADCRate > 62.5
                    Receive(i).ADCRate = Receive(i).ADCRate/2;
                    if ~strcmp(Receive(i).sampleMode,'NS200BWI')
                        error(['VsUpdate: The Receive(%d).ADCRate specified or derived from 4*Trans.frequency\n',...
                               'of %2.3f exceeds the sample rate of the system. Receive.sampleMode must be set\n',...
                               'to ''NS200BWI'' for interleaved sampling.  ADCRate will be set to %2.3f MHz and\n',...
                               'Two acquisitions with successive acqNums are required with an offset in TX.Delay\n',...
                               'for the first. Interleaving will be performed by reconstruction with a ReconInfo\n',...
                               'that references the 1st acqNum.\n'],i,2*Receive(i).ADCRate,Receive(i).ADCRate);
                    end
                end
            end
            if strcmp(Receive(i).sampleMode,'BS67BW')
                demodFrequency = Receive(i).decimSampleRate * 3/4;
            else
                demodFrequency = Receive(i).decimSampleRate/4;
            end
            if ~isempty(Receive(i).demodFrequency)
                if abs(Receive(i).demodFrequency-demodFrequency)/Receive(i).demodFrequency > .05 % is error > 5%?
                    if (warn==1)&&Resource.Parameters.verbose>0
                        fprintf('VsUpdate: The specified Receive(%d).demodFrequency value could not be\n',i);
                        fprintf('    matched to within 5 percent.  The closest match of %2.3f MHz\n',demodFrequency);
                        fprintf('    is being used. This warning message is only printed once.\n');
                        warn = 0;
                    end
                end
            end
            Receive(i).demodFrequency = demodFrequency;
            % Set quadDecim and samplesPerWave attributes. The quadDecim factor is the quadrature
            %   decimation performed after the Inputfilter.
            switch Receive(i).sampleMode
                case 'NS200BW'
                    Receive(i).quadDecim = 1;
                    Receive(i).samplesPerWave = Receive(i).decimSampleRate/Trans.frequency;
                case 'NS200BWI'
                    Receive(i).quadDecim = 1;
                    Receive(i).samplesPerWave = Receive(i).decimSampleRate/(2*Trans.frequency);
                case 'BS100BW'
                    Receive(i).quadDecim = 2;
                    Receive(i).samplesPerWave = Receive(i).decimSampleRate/(2*Trans.frequency);
                case 'BS67BW'
                    Receive(i).quadDecim = 1;
                    Receive(i).samplesPerWave = Receive(i).decimSampleRate/Trans.frequency;
                case 'BS50BW'
                    Receive(i).quadDecim = 4;
                    Receive(i).samplesPerWave = Receive(i).decimSampleRate/(4*Trans.frequency);
                case 'custom'
                    Receive(i).quadDecim = 1;
                    Receive(i).samplesPerWave = Receive(i).decimSampleRate/Trans.frequency;
               otherwise
                   error('VsUpdate: Unrecognized Receive(%d).sampleMode.\n',i);
            end

            % Provide default InputFilter and LowPassCoef filter coefficient arrays
            % if the user has not specified them.
            %   If running with the hardware, the default filters will get set from
            %   the decimFactor and quadDecim values; otherwise they will default
            %   to pass through. If the user defines filter coefficients, they will
            %   only get set at initialization (a runtime update of Receive will
            %   not change the filter coefficients if they already exist).  If user
            %   changes Receive structure while system is running and wants a new
            %   set of default coefficients to be selected, they can set the
            %   LowPassCoef or InputFilter array to empty.

            % Listed below are default Gen3 LPF filter coef, for
            % 23 coefficients and 16 bit coefficient resolution.
            DefaultLPFcoef =[...
                [ +0.0000  +0.0000  +0.0000  +0.0000  +0.0000  +0.0000  +0.0000  +0.0000  +0.0000  +0.0000  +0.0000  +1.0000];...
                [ -0.0004  +0.0063  +0.0038  -0.0139  -0.0137  +0.0234  +0.0359  -0.0328  -0.0868  +0.0398  +0.3105  +0.4559];...
                [ -0.0057  -0.0005  +0.0115  +0.0197  +0.0095  -0.0208  -0.0497  -0.0411  +0.0285  +0.1444  +0.2545  +0.2997];...
                [ +0.0110  +0.0072  -0.0042  -0.0202  -0.0334  -0.0337  -0.0128  +0.0311  +0.0914  +0.1538  +0.2008  +0.2183];...
                [ -0.0006  -0.0058  -0.0130  -0.0191  -0.0192  -0.0084  +0.0161  +0.0531  +0.0972  +0.1394  +0.1699  +0.1810];...
                [ -0.0137  -0.0183  -0.0194  -0.0148  -0.0030  +0.0163  +0.0419  +0.0710  +0.1001  +0.1248  +0.1414  +0.1472];...
                [ -0.0212  -0.0198  -0.0139  -0.0031  +0.0124  +0.0315  +0.0528  +0.0744  +0.0942  +0.1101  +0.1205  +0.1240];...
                [ -0.0192  -0.0144  -0.0060  +0.0060  +0.0208  +0.0376  +0.0552  +0.0722  +0.0873  +0.0992  +0.1067  +0.1093];...
                [ +0.0030  +0.0034  -0.0093  -0.0068  +0.0214  +0.0106  -0.0441  -0.0141  +0.0932  +0.0166  -0.3135  +0.4820]; ...
            ];

            if ~isfield(Receive(i),'LowPassCoef')||isempty(Receive(i).LowPassCoef)
                % There is no LowPassCoef array, so generate it from the
                % precomputed defaults
                if isfield(Receive(i),'decimFactor') && ~isempty(Receive(i).decimFactor)
                    if ~strcmp(Receive(i).sampleMode,'BS67BW') % if sampleMode not 67% BW
                        Receive(i).LowPassCoef = DefaultLPFcoef(Receive(i).decimFactor,:);
                    elseif Receive(i).decimFactor == 2
                        % special case of 8/3 A/D rate followed by
                        % subsampling by 2 for 4/3 sampling; implement
                        % highpass filter from row 9
                        Receive(i).LowPassCoef = DefaultLPFcoef(9,:);
                    else
                        % 4/3 sampling with no subsampling after A/D so
                        % we have to use the pass-through filter
                        Receive(i).LowPassCoef = DefaultLPFcoef(1,:);
                    end
                else
                    Receive(i).LowPassCoef = DefaultLPFcoef(1,:);
                end
            else
                % user-defined filter exists; if of V1 size expand it to be
                % compatible with Gen3
                if size(Receive(i).LowPassCoef, 2) == 6
                    Receive(i).LowPassCoef = [zeros(1,6), Receive(i).LowPassCoef];
                end
            end

            % Set up the bandpass filter coefficient array.
            % Gen3 default BPF coefficients: 41 taps, 21 symmetric coefficients with center
            % tap last.  V1 had 21 taps with even taps set to zero, and thus used 6
            % coefficients.

            DefaultBPFcoef = zeros(8, 21); % default coef array of 8 filters, 21 taps each
            DefaultBPFcoef(:, 21) = ones(8, 1); % set center tap to 1.0 for a 'pass through' response

            % BPF #1 (4 samples per wavelength, default coef row 1):
            DefaultBPFcoef(1, :) = [ -0.00113 +0.00000 -0.00116 +0.00000 +0.00549 +0.00000 +0.00720 ...
                                     +0.00000 -0.01419 +0.00000 -0.02640 +0.00000 +0.02606 +0.00000 ...
                                     +0.07816 +0.00000 -0.03671 +0.00000 -0.30786 +0.00000 +0.54108 ];

            % BPF #2 (2 samples per wavelength, default coef row 2):
            DefaultBPFcoef(2, :) = [ +0.00034 +0.00000 +0.00244 +0.00000 -0.00629 +0.00000 -0.00333 ...
                                     +0.00000 +0.02188 +0.00000 -0.00897 +0.00000 -0.04745 +0.00000 ...
                                     +0.06076 +0.00000 +0.07294 +0.00000 -0.30048 +0.00000 +0.41632 ];

            % BPF #3 (1 sample per wavelength, default coef row 4):
            DefaultBPFcoef(4, :) = [ -0.00162 +0.00000 +0.00568 +0.00000 -0.01065 +0.00000 +0.01349 ...
                                     +0.00000 -0.00858 +0.00000 -0.00955 +0.00000 +0.04312 +0.00000 ...
                                     -0.08841 +0.00000 +0.13550 +0.00000 -0.17130 +0.00000 +0.18463 ];

            % BPF #4 (4/3 samples per wavelength aliased, default coef row 5):
            DefaultBPFcoef(5, :) = [ -0.00159 +0.00000 -0.00549 +0.00000 -0.01157 +0.00000 -0.02066 ...
                                     +0.00000 -0.03275 +0.00000 -0.04721 +0.00000 -0.06281 +0.00000 ...
                                     -0.07785 -0.00003 -0.09039 -0.00003 -0.09875 -0.00003 +0.89832 ];
            % BPF #5 (interleave, default coef row 6):
            DefaultBPFcoef(6, :) = [ +0.00000 +0.00214 +0.00000 -0.00409 +0.00000 +0.00693 +0.00000 ...
                                     -0.01093 +0.00000 +0.01654 +0.00000 -0.02457 +0.00000 +0.03665 ...
                                     +0.00000 -0.05713 +0.00000 +0.10217 +0.00000 -0.31735 +0.50067];

            if isfield(Receive(i),'InputFilter') && ~isempty(Receive(i).InputFilter)
                % There is a Receive.InputFilter array so use as is unless it has multiple rows.
                if size(Receive(i).InputFilter,1) ~= 1 % if multiple rows, check to see if all are the same
                    if max(Receive(i).InputFilter,[],1) == min(Receive(i).InputFilter,[],1)
                        % all rows identical, so copy row 1 into a new single-row
                        % vector
                        bpcoef = Receive(i).InputFilter(1,:);
                        clear Receive(i).InputFilter
                        Receive(i).InputFilter = bpcoef;
                        clear bpcoef
                    else
                        error('VSX: Receive(%d).InputFilter has multiple rows that are not identical.\n',i);
                    end
                end
                if size(Receive(i).InputFilter,2) == 6
                    % looks like a V1 filter provided by user so expand it to Gen3 size
                    BPF = zeros(1, 21);
                    for ntap = 1:6
                        BPF(1,9 + 2*ntap) = Receive(i).InputFilter(:, ntap);
                    end
                    Receive(i).InputFilter = BPF;
                    clear BPF
                end
                % Now we have a single-row array; ready to convert to format
                % used by HAL
            else
                % No Input Filter array exists, so create default from
                % the DefaultBPFcoef array using sampleMode to
                % select the appropriate response
                switch Receive(i).sampleMode
                    case 'BS50BW'
                        % 4X sampling with 50% BW QBS subsampling
                        Receive(i).InputFilter = DefaultBPFcoef(4,:);
                    case 'BS100BW'
                        % 4X sampling with 100% BW QBS subsampling
                        Receive(i).InputFilter = DefaultBPFcoef(2,:);
                    case 'NS200BW'
                        % 4X uniform sampling 200% BW
                        Receive(i).InputFilter = DefaultBPFcoef(1,:);
                    case 'BS67BW'
                        % 4/3 aliased sampling
                        Receive(i).InputFilter = DefaultBPFcoef(5,:);
                    case 'NS200BWI'
                        % interleave sampling
                        Receive(i).InputFilter = DefaultBPFcoef(6,:);
                    otherwise
                        % any other value for sampleMode means default
                        % to the pass-through filter
                        Receive(i).InputFilter = DefaultBPFcoef(8,:);
                end
            end
        end
        clear d SR demodFrequency

        %% Set startSample and endSample attributes to correspond to startDepth and endDepth attributes.
        % - Force numRcvSamples to be a multiple of 128.  Round up to find nearest multiple,
        %   if the user specified endDepth doesn't come within a few samples.
        i = 1;
        n = size(Receive,2);
        while i <= n
            % Set values for 1st Receive of frame
            if Receive(i).acqNum ~= 1, error('VsUpdate: Receive(%d).acqNum for start of frame not equal to 1.\n',i); end
            maxAcqNum = 1;  % maxAcqNum keeps track of the largest acqNum for the frame.
            nSmpls = 2*(Receive(i).endDepth - Receive(i).startDepth) * Receive(i).samplesPerWave;
            if abs(round(nSmpls/128) - nSmpls/128) < .01
                numRcvSamples = 128*round(nSmpls/128);
            else
                numRcvSamples = 128*ceil(nSmpls/128);
            end
            if numRcvSamples == 0
                numRcvSamples = 128;
                if Resource.Parameters.verbose>1
                    fprintf('VsUpdate: Setting Receive(%d).endDepth to acquire a minimum of 128 samples.\n',i);
                end
            end
            Receive(i).endDepth = Receive(i).startDepth + (numRcvSamples/2)/Receive(i).samplesPerWave;
            Receive(i).startSample = 1;
            Receive(i).endSample = numRcvSamples;
            j = numRcvSamples + 1; % j is the next startSample to use
            i = i+1;
            while (i <= n) && (Receive(i).bufnum == Receive(i-1).bufnum) && (Receive(i).framenum == Receive(i-1).framenum)
                nSmpls = 2*(Receive(i).endDepth - Receive(i).startDepth) * Receive(i).samplesPerWave;
                if abs(round(nSmpls/128) - nSmpls/128) < .01
                    numRcvSamples = 128*round(nSmpls/128);
                else
                    numRcvSamples = 128*ceil(nSmpls/128);
                end
                if numRcvSamples == 0
                    numRcvSamples = 128; % minimum no. of samples is 128
                    if Resource.Parameters.verbose>1
                        fprintf('VsUpdate: Setting Receive(%d).endDepth to acquire a minimum of 128 samples.\n',i);
                    end
                end
                Receive(i).endDepth = Receive(i).startDepth + (numRcvSamples/2)/Receive(i).samplesPerWave;
                if Receive(i).acqNum == maxAcqNum + 1
                    maxAcqNum = maxAcqNum + 1;
                    Receive(i).startSample = j;
                    Receive(i).endSample = j + numRcvSamples - 1;
                    j = Receive(i).endSample + 1;
                elseif Receive(i).acqNum <= maxAcqNum  % If acqNum has been previously used, overwrite previous acquisition.
                    k = i-1;  % search backwards for a previous Receive with the same acqNum
                    while (k > 0)
                        if Receive(k).acqNum == Receive(i).acqNum && ...
                                Receive(k).bufnum == Receive(i).bufnum && ...
                                Receive(k).framenum == Receive(i).framenum
                            break
                        end
                        k = k - 1;
                    end
                    if k == 0, error('VsUpdate: Out of order Receive(%d).acqNum which doesn''t match previous acqNum.\n',i); end
                    Receive(i).startSample = Receive(k).startSample;
                    Receive(i).endSample = Receive(i).startSample + numRcvSamples - 1;
                    if Receive(i).endSample > Receive(k).endSample
                        error('VsUpdate: Receive(%d).endSample for overwrite acqNum larger than in previous Receive.\n',i);
                    end
                else
                    error('VsUpdate: Receive(%d).acqNum is not in sequential order.\n',i);
                end
                i = i+1;
            end
        end

        if VDASupdates
            % Set attributes in each receive structure.
            for i = 1:size(Receive,2)
                % Build an aperture to connector channel map. For an HVMux
                %   transducer, the aperture is defined in Trans.HVMux.Aperture.
                %   For a non-HVMux, it is defined by the Trans.Connector array.
                if isfield(Trans,'HVMux')
                    % for HVMux transducer verify the required TX.aperture field
                    % is present
                    if (~isfield(Receive(i),'aperture'))||(isempty(Receive(i).aperture))||(Receive(i).aperture<=0)
                        error('VsUpdate: ''aperture'' field required for HVMux transducer is missing or empty in Receive(%d).\n',i);
                    end
                    [ActiveEL, ~, Aper2Chnl] = find(Trans.HVMux.Aperture(:,Receive(i).aperture));
                else
                    % use first column of Trans.Connector
                    [ActiveEL, ~, Aper2Chnl] = find(Trans.Connector(:,1));
                end

                if Resource.Parameters.sizeApod == Trans.numelements
                    % For all-element Apod, force all entries
                    % that are not part of the active aperture to zero
                    ApodMask = zeros(1, Trans.numelements);
                    ApodMask(ActiveEL) = 1;
                    Receive(i).Apod = ApodMask .* Receive(i).Apod;
                end

                %%  Receive.Apod
                if any(Receive(i).Apod < -4.0) || any(Receive(i).Apod > 3.997)
                    error('VsUpdate: Receive(%d).Apod value exceeds allowed range of [-4.0:3.997].\n',i);
                end

                % Apply receive mux constraint to 32LE and 64LE system configurations.
                if LEsys && numBoards == 1
                    if nnz(Receive(i).Apod) > 32
                        error('VsUpdate: Receive(%d).Apod is attempting to enable more than 32 channels on Vantage 32LE system.\n', i);
                    end
                elseif LEsys && numBoards == 2
                    if nnz(Receive(i).Apod) > 64
                        error('VsUpdate: Receive(%d).Apod is attempting to enable more than 64 channels on Vantage 64LE system.\n', i);
                    end
                end
            end
        end % end of the VDAS updates (if enabled) for Receive

        % Return updated Receive object to base.
        assignin('base','Receive',Receive);
        if updateTrans
            % changes have been made to Trans, so return it as well
            assignin('base','Trans',Trans);
        end



        %% If 'ReconInfo' structures exist, update ReconInfo.Aperture
        if evalin('base','exist(''ReconInfo'', ''var'')')
            ReconInfo = evalin('base','ReconInfo');
            % The ReconInfo.Aperture field is treated as a "dependent
            % variable" that will always be overwritten here based on the
            % mapping specified by Trans.Connector or Trans.HVMux.Aperture.
            % However, we will not do this if this is a multisystem
            % configuration and the user script has already defined
            % ReconInfo.Aperture
            if ~(usingMultiSys && isfield(ReconInfo(1),'Aperture'))
                for i = 1:size(ReconInfo,2)
                    % always update ReconInfo.Aperture, to make sure if matches
                    % any changes to Receive.Apod
                    % If we have an HVmux aperture definition, build a connector I/O channel to element map.
                    if isfield(Receive(ReconInfo(i).rcvnum),'aperture') && ...
                            ~isempty(Receive(ReconInfo(i).rcvnum).aperture)
                        n = Receive(ReconInfo(i).rcvnum).aperture;
                        [Indices,~,Aper2Chnl] = find(Trans.HVMux.Aperture(:,n));
                        Chnl2Ele(Aper2Chnl) = Indices-1;
                    else
                        % Not using HVMux, so set aperture mapping from
                        % first column of Trans.Connector
                        [Indices,~,Aper2Chnl] = find(Trans.Connector(:, 1));
                        Chnl2Ele(Aper2Chnl) = Indices-1;
                    end

                    % Now generate the actual ReconInfo.Aperture array
                    % With new HAL CG mapping, Recon numchannels should always match
                    % the size of the receive buffer:
                    % initialize Aperture array to the disabled channel state
                    ReconInfo(i).Aperture(1:Resource.Parameters.numTransmit) = -1;
                    % VTS-804 only map the active elements since sizeApod may
                    % be all elements
                    if Resource.Parameters.sizeApod == Trans.numelements
                        for j = 1:length(Indices)
                            if Receive(ReconInfo(i).rcvnum).Apod(Indices(j)) ~= 0
                                ReconInfo(i).Aperture(Aper2Chnl(j)) = Chnl2Ele(Aper2Chnl(j));
                            end
                        end
                    else
                        for j = 1:length(Indices)
                            if Receive(ReconInfo(i).rcvnum).Apod(j) ~= 0
                                ReconInfo(i).Aperture(Aper2Chnl(j)) = Chnl2Ele(Aper2Chnl(j));
                            end
                        end
                    end
                end  % end of ReconInfo(i) for loop
                % Return updated ReconInfo object to base.
                assignin('base','ReconInfo', ReconInfo);
            end
        end  % done updating ReconInfo


    case 'SeqControl'
        SeqControl = evalin('base', 'SeqControl');
        % VTS-1617 add check for old synchronous DMA method, and
        % automatically enable new method
        enableNewDMA = 0; % flag to indicate new DMA scheme should be enabled
        for i = 1:length(SeqControl)
            % ensure condition field exists.
            if ~isfield(SeqControl(i), 'condition')
                SeqControl(i).condition = [];
            end

            % condition field of timeToNext* commands defaults to 'report'
            if strcmp(SeqControl(i).command, 'timeToNextAcq') || ...
                     strcmp(SeqControl(i).command, 'timeToNextEB')
                if ~isfield(SeqControl(i), 'condition') || ...
                        isempty(SeqControl(i).condition)
                    SeqControl(i).condition = 'report';
                end
            end

            % condition field of loopCnt and loopTst commands default to
            % 'counter1'
            if strcmp(SeqControl(i).command, 'loopCnt') || ...
                    strcmp(SeqControl(i).command, 'loopTst')
                if ~isfield(SeqControl(i), 'condition') || ...
                        isempty(SeqControl(i).condition)
                    SeqControl(i).condition = 'counter1';
                end
            end

            % condition field of triggerOut command defaults to 'syncNone'
            if strcmp(SeqControl(i).command, 'triggerOut')
                if ~isfield(SeqControl(i), 'condition') || ...
                        isempty(SeqControl(i).condition)
                    SeqControl(i).condition = 'syncNone';
                end
            end

            % condition field of multiSysSync command defaults to 'normal'
            if strcmp(SeqControl(i).command, 'multiSysSync')
                if ~isfield(SeqControl(i), 'condition') || ...
                        isempty(SeqControl(i).condition)
                    SeqControl(i).condition = 'normal';
                end
            end

            % VTS-1617 check for old TTH sync DMA wait for processing
            % condition and replace it with new scheme
            if strcmp(SeqControl(i).command, 'transferToHost')
                if isfield(SeqControl(i), 'condition') && strcmp(SeqControl(i).condition, 'waitForProcessing')
                    enableNewDMA = 1; % need to enable new DMA scheme and disable old one
                    SeqControl(i).condition = [];
                    SeqControl(i).argument = [];
                end
            end
        end

        if enableNewDMA
            % enable the new DMA scheme with Resource.Parameters
            evalin('base', 'Resource.Parameters.waitForProcessing = 1;');
            if Resource.Parameters.verbose > 1
                % notify the user
                disp('VsUpdate Status: transferToHost SeqControl command with waitForProcessing condition has been found;');
                disp('    Setting the field Resource.Parameters.waitForProcessing to enable new synchronous DMA scheme,');
                disp('    and removing condition and argument from the transferToHost Command.');
            end
        end

        assignin('base', 'SeqControl', SeqControl);
    otherwise  %% check for unsupported call to VsUpdate
        error('Unrecognized sequence object in ''VsUpdate''.');
end
end

function TX = computeDescriptors(TX)
    % This is the functionality to be transferred into Sequence Load
    % per VTS-388
    % The TX.VDASStates array contains a States array for every channel,
    % with a null array (first row set to [30 0] for inactive channels.
    % This routine will create a TX.Descriptor table for every active
    % channel, and a null Descriptor (all zeros) for inactive channels.

    % To avoid re-creating the identical descriptor over and over when only
    % a single waveform is being used on all active channels, the field
    % TX.perChWvfm can be used.  If it is true, search TX.VDASStates for an
    % active waveform, convert it to the single active Descriptor, and then
    % just copy that descriptor into all active channels.

    % We have to determine whether the transmit HW requires inversion of
    % the States waveform definition, to produce the correct polarity at
    % the transmit output.  This decision is based on the TXindex value
    % from the acquisition boards
    TXindex = evalin('base', 'Resource.SysConfig.TXindex');
    switch TXindex
        case {1, 5}
            % standard frequency transmit DA2319 xfmr
            invert = 0;
        case 2
            % High Frequency original PWB1010-1L xfmr
            invert = 1;
        case {3, 7}
            % Low frequency transmit
            invert = 1;
        case {4, 6}
            % High frequency current PWB1010L transformer
            invert = 1;
        otherwise
            % assume no inversion
            invert = 0;
    end

    % VTS-1365 updates:  Go through the TX structures, and create only one
    % descriptor if TX(i).perChWvfm is false and then apply it to all
    % channels. Otherwise, create a unique Descriptor for each active
    % channel.
    for i = 1:length(TX)
        numch = size(TX(i).VDASApod, 2);
        if TX(i).perChWvfm
            % Initialize the per-channel arbwave table with all zeros.
            % This automatically programs all disabled channels to the
            % inactive state.  Since we don't know the length of the
            % descriptors that will be generated, define the table
            % initially at max length.

            maxLength = 1024; % Maximum arbwave length supported by the system
            FullDescriptor = uint16(zeros(numch, maxLength));
            actualLength = 4; % to keep track of the longest length we actually get
            % Initial length of 4 ensures the null waveform for inactive
            % channels will be a descriptor array of four zeros, which is a
            % legitimate null waveform descriptor for both active clamp and
            % pre-active clamp system configurations
            for chnum = 1:numch
                if TX(i).VDASStates(1, 1, chnum) == 30
                    % this is a null waveform so just continue to next
                    % channel (Descriptor has already been filled with
                    % zeros, which results in a null waveform for the HW system)
                    continue
                else
                    States = TX(i).VDASStates(:, :, chnum);
                    Descriptor = States2Descriptor(States, invert, i);
                    chLength = size(Descriptor, 2); % how long is this channel's arbwave?
                    actualLength = max(actualLength, chLength);
                    FullDescriptor(chnum,1:chLength) = Descriptor; % copy in the arbwave table for this channel
                end
            end
            % now populate the output arbwave table using the
            % actual length
            TX(i).Descriptor = FullDescriptor(:, 1:actualLength);
        else
            % PerChWvfm is false so create the single descriptor that will
            % be used for all active channels; first have to find an active
            % channel
            activeCHnum = find(TX(i).VDASStates(1, 1, :) ~= 30, 1);
            if isempty(activeCHnum)
                % no active waverorms were found so provide a null waveform
                singleDescriptor = States2Descriptor([30 0], invert, i);
            else
                singleDescriptor = States2Descriptor(TX(i).VDASStates(:, :, activeCHnum), invert, i);
            end
            % Create per-channel descriptor table initialized to a null
            % waveform for every channel
            TX(i).Descriptor = uint16(zeros(numch, size(singleDescriptor, 2)));
            for chnum = 1:numch
                if TX(i).VDASStates(1, 1, chnum) == 30
                    % this is a null waveform so just continue to next
                    % channel (Descriptor has already been filled with
                    % zeros, which results in a null waveform for the HW system)
                    continue
                else
                    TX(i).Descriptor(chnum, :) = singleDescriptor;
                end
            end
        end
    end
end

function Descriptor = States2Descriptor(States, invert, txnum)
    if invert
        % invert the entire waveform
        States(:, 1) = -States(:, 1);
    end
    % VTS-1607 use the new waveform Converter
    [status, Descriptor, ~] = convertStatesToDescriptor(0, States);
    if ~isequal(status, 'Success')
        % convertStatesToDescriptor has reported an error condition;
        % report it and quit
        disp(['convertStatesToDescriptor error: ', status]);
        error('VsUpdate(TX): see above error from convertStatesToDescriptor for TX(%d).\n', txnum);
    end
end
