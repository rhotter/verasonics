function varargout = ProbeImpedanceCheck(testFreq, Trans, connSel)
% Wrapper function for making approximate measurements of probe impedance.
% This is the function a user invokes to make the probe impedance
% measurement.  This function assembles the parameters needed to define the
% desired measurement, and then calls "ProbeImpedanceCheckAcq" to make the
% actual measurements through the hardware system.
%
% Copyright 2001-2021, Verasonics, Inc. All worldwide rights and remedies under
% all intellectual property laws and industrial property laws are reserved.
%
% The optional output argument "ProbeZdata" is a structure containing
% fields defining the system configuration and the measurement results, as
% listed below.  If called without an output argument, the ProbeZdata
% structure will still be created, and written to the Matlab base workspace.
%    ProbeZdata.sys: Name of the system configuration used for the
%         measurement
%    ProbeZdata.UTA: Name of the UTA module being used
%    ProbeZdata.probe: Name of the probe being tested, or "unknown" if not
%         identified in the Trans structure being referenced
%    ProbeZdata.testFreq: Frequency in MHz at which measurements were made
%    ProbeZdata.MeasZperCh: Estimated impedance magnitude in Ohms as seen
%         by the system at the probe connector, for each element as
%         identified by the Trans structure.
%    ProbeZdata.medianMeasZ: Median of the per-element values in MeasZperCh
%    ProbeZdata.apparentOpens: Element numbers of any elements that appear
%         to be open
%    ProbeZdata.apparentShorts: Element numbers of any elements that appear
%         to be shorted to ground (measured impedance less than 10 Ohms)
%
% If the optional input arguments are specified, they will be used to
% define the test.  If they are not specified, this function will look for
% a known probe connected to the system and if found, use that probe's
% Trans structure to define the test.  If a known probe cannot be
% identified, this function will create a default Trans structure using
% connector 1 of the UTA that is present on the system.  For more
% information refer to the document "Probe Impedance Check User Guide".

    %% Get HW system configuration

    SysConfig = hwConfigCheck(0);

    VDAS = SysConfig.VDAS; % will be one if HW system is present, zero if not
    if VDAS == 0 || SysConfig.SWconfigFault || SysConfig.HWconfigFault || SysConfig.FPGAconfigFault
        if nargout == 0
            fprintf(2, 'HW system not detected or in a fault condition. This Utility cannot be used.\n')
            return
        else
            ProbeZdata.sys = 'HW system not detected or in a fault condition';
            varargout{1} = ProbeZdata;
            return
        end
    end
    
    if strcmp(SysConfig.Frequency,'LF')
        if nargout == 0
            fprintf(2, 'Low Frequency configuration is not supported; will be added in a future release.\n')
            return
        else
            ProbeZdata.sys = 'Low Frequency configuration is not yet supported by this utility';
            varargout{1} = ProbeZdata;
            return
        end
    end

    ProbeZdata.sys = SysConfig.System; % HW system configuration
    
    if SysConfig.UTA == 1
        % identify the UTA module that is present
        if SysConfig.UTAtype(4) == 2 || SysConfig.UTAtype(4) == 3
            % UTA is 260-MUX or 1024-MUX
            if nargout == 0
                fprintf(2, 'HVMux UTA modules are not supported; will be added in a future release.\n')
                return
            else
                ProbeZdata.UTA = 'HVMux UTA modules are not yet supported by this utility';
                varargout{1} = ProbeZdata;
                return
            end
        end
        ProbeZdata.UTA = SysConfig.UTAname;
    end
    
    numConn = SysConfig.UTAtype(3);
    if numConn > 1
        % specify probe connector of set of connectors to be used but only if
        % more than 1 are available
        if nargin ~= 3
            % not defined by user so default to connector 1
            connSel = 1;
        end
        ProbeZdata.connSel = connSel; % user-provided value
    else
        connSel = 1;
    end
    
    % If Trans structure not provided, look for recognized probe at
    % available connectors and use it to define Trans; if no identifiable
    % probe is found create a default Trans structure for the test.
    if nargin < 2
        Trans = []; % this indicates we didn't find an identifiable probe
        % Check connector status to see if a probe can be identified
        [ connected, ~, probeId ] = getConnectorInfo();
        for connNum = 1:size(connected, 2)
            if connected(connNum) == 1 && probeId(connNum) > 0
                % something with a probe ID is connected so see if
                % computeTrans recognizes it
                probeName = computeTrans(probeId(connNum));
                if ~isempty(probeName) && ~strcmpi(probeName, 'Unknown')
                    Trans.name = probeName;
                    Trans.units = 'mm';
                    Trans = computeTrans(Trans);
                    connSel = connNum;
                    break % leave the for loop since a recognized probe has been found
                end
            end
        end
        if isempty(Trans)
            % Cannot find a recognizable probe, so look in base workspace
            % to see if a Trans structure exists.  If so, use it otherwise
            % make up a default Trans structure
            if evalin('base', 'exist(''Trans'', ''var'')') && evalin('base', 'isstruct(Trans)')
                Trans = evalin('base', 'Trans');
            else
                % create default Trans structure
                Trans.name = 'unknown';
            end
        end
    end
    % get probe name, numelements, ConnectorES or set defaults
    if ~isfield(Trans, 'name') || isempty(Trans.name)
        Trans.name = 'unknown';
    end
    ProbeZdata.probe = Trans.name;
    if ~isfield(Trans, 'numelements') || isempty(Trans.numelements)
        % get default from computeUTA
        UTA  = computeUTA( SysConfig.UTAtype, connSel);
        Trans.numelements = UTA.numCh;
        Trans.ElementPos = zeros(Trans.numelements, 4);
    end
    if ~isfield(Trans, 'lensCorrection') || isempty(Trans.lensCorrection)
        Trans.lensCorrection = 0; % dummy non-empty value
    end
    if nargin == 0 || isempty(testFreq) || testFreq == 0
        % set testFreq to default value
        if isfield(Trans,'frequency') && ~isempty(Trans.frequency)
            testFreq = Trans.frequency;
            Trans.frequency = 15.625;
        else
            testFreq = 5; % default value if no clues given
            Trans.frequency = 15.625;
        end
    elseif nargin > 0
        Trans.frequency = 15.625;
    end
    if ~isfield(Trans, 'spacing') || isempty(Trans.spacing)
        Trans.spacing = 0.5; % dummy half-wavelength spacing
        Trans.spacingMm = Trans.spacing * 1.54 / Trans.frequency;
        Trans.elementWidth = 0.45;
    end

    if isfield(Trans, 'HVMux')
        if nargout == 0
            fprintf(2, 'HVMux probes are not supported; will be added in a future release.\n')
            return
        else
            ProbeZdata.probe = 'HVMux probes are not yet supported by this utility';
            varargout{1} = ProbeZdata;
            return
        end
    end
    
    %% find FreqList index for closest match to testFreq
    % create the list of supported frequencies
    numFreq = 30; % number of frequencies supported
    Freqin = 2.^((-5:24)/5); % 30 steps from 0.5 to 31.25 MHz, approx 15% per step
    FreqList = zeros(1, numFreq);
    % Convert to the actual TX frequencies we wil use
    Fs = 125;
    for freqnum = numFreq:-1:1
        div = round(Fs/Freqin(freqnum));
        if div > 9
            div = div/2;
            Fs = Fs/2;
        end
        FreqList(freqnum) = Fs/div;
    end
    
    fIndxH = find(FreqList>testFreq, 1); % index to first frequency greater than testFreq
    if isempty(fIndxH) || fIndxH > 23
        % testFreq higher than 10 MHz (index 23) not supported for now
        txFreqIndx = 23;
    elseif fIndxH == 1
        % testFreq lower than first entry so use first one
        txFreqIndx = 1;
    else
        % in between two entries so find the closest
        if testFreq/FreqList(fIndxH-1)<FreqList(fIndxH)/testFreq
            % lower one is closer
            txFreqIndx = fIndxH-1;
        else
            % higher one is closer
            txFreqIndx = fIndxH;
        end
    end
    testFreq = FreqList(txFreqIndx);
    ProbeZdata.testFreq = testFreq;
%     ProbeZdata.TransIn = Trans;  % for debug purposes only
    assignin('base', 'ProbeZdat', ProbeZdata);
    
    %% Call the Acquisition function, that will run VSX to get the RMS data
    [ProbeRMSdata] = ProbeImpedanceCheckAcq(txFreqIndx, Trans);

    % Following commented-out lines are for debug only, providing intermediate results
    % not needed in normal use of the tool
%     ProbeZdata.TransOut = ProbeRMSdata.Trans;
%     ProbeZdata.medianRMSlvl = ProbeRMSdata.medianRMSlvl;
%     ProbeZdata.RMSlevel = ProbeRMSdata.ESrmsTot;

    % estimated impedance will be derived from these measured values for
    % each transducer element
    RMSlevel = ProbeRMSdata.ESrmsTot;
    
    
    %% calculate impedance from measured RMS level

    % ProbeZcalRMSlevel will contain four test values at each frequency, with open
    % circuit value scaled by 0.8 to represent highest impedance that will not
    % be reported as "apparent open", and a fifth value that is the resulting
    % highest non-open estimated impedance value: [0, 50, 100, 0.8*open, maxR]
    
    if strcmp(SysConfig.Frequency,'HF')
        % Data for high frequency system with 732-109 Acq boards, CGD file ID 2    
        ProbeZcalRMSlevel = 1.0e+04 * [ ...
            0.0193    0.0756    0.1231    1.0438    0.1071; ...
            0.0165    0.0822    0.1336    1.0284    0.0970; ...
            0.0140    0.0927    0.1550    1.0509    0.0819; ...
            0.0113    0.1012    0.1751    1.0447    0.0689; ...
            0.0108    0.1056    0.1837    1.0386    0.0647; ...
            0.0110    0.1097    0.1930    1.0224    0.0598; ...
            0.0117    0.1133    0.1996    1.0009    0.0564; ...
            0.0110    0.1035    0.1854    0.8597    0.0512; ...
            0.0121    0.0964    0.1735    0.7341    0.0464; ...
            0.0114    0.0889    0.1604    0.6353    0.0432; ...
            0.0117    0.0836    0.1511    0.5454    0.0392; ...
            0.0119    0.0790    0.1441    0.4579    0.0341; ...
            0.0119    0.0736    0.1318    0.3707    0.0305; ...
            0.0108    0.0645    0.1141    0.2750    0.0262; ...
            0.0109    0.0609    0.1080    0.2343    0.0234; ...
            0.0106    0.0563    0.0967    0.1923    0.0218; ...
            0.0109    0.0538    0.0919    0.1613    0.0191; ...
            0.0129    0.0599    0.1013    0.1471    0.0155; ...
            0.0136    0.0573    0.0882    0.1139    0.0142; ...
            0.0139    0.0573    0.0851    0.0992    0.0126; ...
            0.0150    0.0574    0.0853    0.0880    0.0105; ...
            0.0153    0.0523    0.0698    0.0678    0.0094; ...
            0.0186    0.0560    0.0616    0.0592    0.0079; ...
                 0         0         0         0         0; ...
                 0         0         0         0         0; ...
                 0         0         0         0         0; ...
                 0         0         0         0         0; ...
                 0         0         0         0         0; ...
                 0         0         0         0         0; ...
                 0         0         0         0         0];

    elseif strcmp(SysConfig.Frequency,'SF') || strcmp(SysConfig.Frequency,'HIFU')
        % Data for standard frequency system with 732-04 Acq boards, CGD file ID 1  
        ProbeZcalRMSlevel = 1.0e+04 * [ ...
            0.0222    0.1957    0.3499    1.1405    0.0356; ...
            0.0229    0.1990    0.3554    1.1521    0.0355; ...
            0.0243    0.2098    0.3747    1.1939    0.0348; ...
            0.0247    0.2095    0.3643    1.1926    0.0368; ...
            0.0246    0.2017    0.3577    1.1850    0.0365; ...
            0.0250    0.2014    0.3478    1.1861    0.0386; ...
            0.0250    0.1933    0.3388    1.1704    0.0386; ...
            0.0230    0.1696    0.2953    1.0341    0.0394; ...
            0.0218    0.1533    0.2643    0.9118    0.0392; ...
            0.0206    0.1400    0.2391    0.8192    0.0393; ...
            0.0197    0.1308    0.2237    0.7438    0.0380; ...
            0.0191    0.1231    0.2087    0.6659    0.0367; ...
            0.0187    0.1155    0.1945    0.5782    0.0343; ...
            0.0167    0.1027    0.1724    0.4525    0.0301; ...
            0.0165    0.0988    0.1674    0.3968    0.0267; ...
            0.0157    0.0933    0.1511    0.3387    0.0262; ...
            0.0159    0.0934    0.1458    0.2982    0.0245; ...
            0.0160    0.0883    0.1343    0.2493    0.0225; ...
            0.0172    0.0883    0.1229    0.2098    0.0226; ...
            0.0185    0.0883    0.1167    0.1929    0.0234; ...
            0.0208    0.0886    0.1142    0.1805    0.0230; ...
            0.0212    0.0807    0.0953    0.1602    0.0323; ...
            0.0289    0.0864    0.0900    0.1542    0.0988; ...
                 0         0         0         0         0; ...
                 0         0         0         0         0; ...
                 0         0         0         0         0; ...
                 0         0         0         0         0; ...
                 0         0         0         0         0; ...
                 0         0         0         0         0; ...
                 0         0         0         0         0];

    else
        if nargout == 0
            fprintf(2, 'Unrecognized/unsupported Acquisition Module Type.\n')
            return
        else
            ProbeZdata.sys = 'Unrecognized Acquisition Module Type';
            varargout{1} = ProbeZdata;
            return
        end
    end         
         
    
    testRin = [0 50 100]; % Load resistance values used for testing
    
    Zinterp = griddedInterpolant(ProbeZcalRMSlevel(txFreqIndx, 1:3), testRin(:));
    ProbeZdata.MeasZperCh = Zinterp(RMSlevel);
    ProbeZdata.medianMeasZ = Zinterp(ProbeRMSdata.medianRMSlvl);
    ProbeZdata.apparentOpens = find(ProbeZdata.MeasZperCh > ProbeZcalRMSlevel(txFreqIndx, 5));
    ProbeZdata.apparentShorts = find(ProbeZdata.MeasZperCh < 10);
    

    %% finish and return
    if nargout == 0
        % write ProbeZdata to base workspace since no return argument specified
        assignin('base', 'ProbeZdata', ProbeZdata);
        
        % plot the results
        numEL = length(ProbeZdata.MeasZperCh);
        figure('Position',[20 200 800 400],'Name','Probe impedance for each Element Signal');
        probeZplotRange = 200;
        probeZplothndl = axes('XLim',[1 numEL], 'YLim',[0 probeZplotRange], 'NextPlot','replacechildren');
        
        % identify apparent shorts
        numShorts = length(ProbeZdata.apparentShorts);
        if numShorts == 0
            shortsString = 'No apparent shorted elements ( < 10 Ohms)';
        elseif numShorts == 1
            shortsString = ['Apparent Short ( < 10 Ohms) at Element ', num2str(ProbeZdata.apparentShorts)];
        elseif numShorts < 9
            shortsString = ['Apparent Shorts ( < 10 Ohms) at Elements ', num2str(ProbeZdata.apparentShorts)];
        else
            shortsString = [num2str(numShorts), ' Elements appear to be shorted ( < 10 Ohms)'];
        end
        
        % identify apparent opens
        numOpens = length(ProbeZdata.apparentOpens);
        opstr = [' ( > ', num2str(round(ProbeZcalRMSlevel(txFreqIndx, 5))), ' Ohms)'];
        if numOpens == 0
            opensString = ['No apparent Open elements ', opstr];
        elseif numOpens == 1
            opensString = ['Apparent Open ', opstr, ' at Element ', num2str(ProbeZdata.apparentOpens)];
        elseif numOpens < 9
            opensString = ['Apparent Opens ', opstr, ' at Elements ', num2str(ProbeZdata.apparentOpens)];
        else
            opensString = [num2str(numOpens), ' Elements appear to be open', opstr];
        end
            
        
        plot(probeZplothndl, 1:numEL, ProbeZdata.MeasZperCh);
        title(probeZplothndl, {['Estimated Element Impedance for Probe:   ', ProbeZdata.probe, '   measured at ', num2str(ProbeZdata.testFreq,'%.2f'), ' MHz']; ...
            ['Using ', ProbeZdata.sys, ' system, with ', ProbeZdata.UTA]; ...
            ['Median impedance ', num2str(ProbeZdata.medianMeasZ,'%.0f'), ' Ohms.']; ...
            shortsString; ......
            opensString});
        xlabel(probeZplothndl, 'Element Number');
        ylabel(probeZplothndl, 'Impedance, Ohms');
        
        
        
    else
        % ProbeImpedanceCheck was called with a return argument, so just
        % return the results structure and do not create the plot
        varargout{1} = ProbeZdata;
    end
end

