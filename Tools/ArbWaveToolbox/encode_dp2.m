% encode_dp.m
%
% Encoder scripting example using MLSE algorithm with and
% impulse response shortening method.
% Examples:
%   LFM waveform encoding
%   DPSS waveforms encoding
% Gain search is employed by repeatedly invoking the encoder with
% different gains.  The gain giving the best encoding error is chosen as
% optimal.
% Results are saved in a mat file.

% John Flynn 4/18/2017 (c) 2017 Verasonics, Inc
% John Flynn updated 9/5/2018 (c) 2017 Verasonics, Inc
% John Flynn updated 8/13/2018 (c) 2017 Verasonics, Inc

if ~exist( 'testsig' ,'var')
    testsig = [];
end
if ~exist( 'nGainSearch' ,'var')
    nGainSearch = [];
end
if ~exist( 'gainSetup' ,'var')
    gainSetup = '';
end
if ~exist( 'functionalityTestMode' ,'var')
    functionalityTestMode = 0;
end

if isempty(functionalityTestMode)
    functionalityTestMode = 0 ;
end
if functionalityTestMode>0,
    disp([mfilename,': using functionalityTestMode :',num2str(functionalityTestMode)])
end
switch functionalityTestMode
    case {1,'fast'}
        %not effective result - only to verify functionality
        nGainSearch = 4;
        gainSetup = 15.5;
        testsig='dpss.L11-5.wb.1';
    case {2,'fast.2'}
        %sub optimal - only searching near a previously found optimum.
        nGainSearch = 2;
        gainSetup = 16.241; %
        testsig='dpss.L11-5.wb.1';
    case {3,'search.3'}
        %somewhat optimal - searching near a previously found optimum.
        nGainSearch = 30;
        gainSetup = 16.241; %
        testsig='dpss.L11-5.wb.1';
    case {4,'search.4'}
        %fast - not optimal.
        nGainSearch = 2;
        gainSetup = 16.241; %
        testsig = 'dpss.L22-14.wb.1' ;

    case 0
        %nothing - normal operation
    otherwise
        error('bad switch: functionalityTestMode')
end
if isempty(testsig)
    %Default test / example design signal:
    %testsig = 'lfm.L7-4.1'  ;
    %  testsig = 'dpss.L7-4.1' ;
    testsig = 'dpss.L22-14.wb.1' ;
    disp([mfilename,': using default "testsig" :',testsig])
end

if ~exist( 'figureCleanup' ,'var')
    figureCleanup = [];
end
if isempty(figureCleanup)
    figureCleanup = 0;  %default
end
if figureCleanup
    disp([mfilename,': figure cleanup enabled (figureCleanup).'])
end

disp([mfilename,':testsig:',testsig])

switch testsig,

    case {'dpss.L7-4.1','dpss.L7-4.wb','dpss.L7-4.nb' ,...
            'dpss.L11-5.wb.1', ...
            'cex-awe.L11-5v.1', ...
            'cex-awe.L11-5v.1a', 'cex-awe.L11-5v.1b', ...
            'dpss.L22-14.wb.1'  ,  'dpss.L22-14.wb.2' , 'dpss.L22-14.wb.3' , ...
            'dpss.L22-14.wb.4' ,   'dpss.L22-14.wb.4b' ,  'dpss.L22-14.wb.4c' , 'dpss.L22-14.wb.4d' , ...
            'dpss.L22-14.wb.5' ...
            'lfm.L7-4.50cyc.1' ...
            'file-ir.L22-14.c#18' ,  'file-ir.L22-14.d#18' ...
            'lfm.L7-4.1' ,'lfm.L7-4.2' ,'lfm.L7-4.3' ...  %LFM waveforms
            } %testsig %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %test signal:  discrete prolate spheroidal sequence function
        disp([mfilename,': encode for specific transducer to produce Slepian sequence waveform...'])
        disp([mfilename,': testsig: ',testsig])

        showFigures = 0;

        if isempty(nGainSearch)
            nGainSearch = 200;  %default, may override below
        end
        %setup gain/scaling search for specific reference signal:
        vitEvalIRSource = 'VAmodel'; %default/legacy - use same IR as VA encoding model
        % vitEvalIRSource = 'offset:3'; %use different IR as VA encoding model (three elements over)
        % vitEvalIRSource = 'offset:-10'; %use different IR as VA encoding model (-10 elements over)
        %vitEvalIRSource = 'offset:102'; %use different IR as VA encoding model (102 elements over)
        ...
            %vitEvalIRSource = 'VAmodel'; %default/legacy - use same IR as VA encoding model
        ...
            % vitEvalIRSource = 'offset:16'; %use different IR as VA encoding model (three elements over)
        % vitEvalIRSource = 'offset:102'; %use different IR as VA encoding model (three elements over)
        %vitEvalIRSource = 'SOME OTHER ELEMENT OR PROBE ...';

        if isempty(gainSetup)
            gainSetup = 15.1 ;  % fri 9/7
            gainSetup = 15.2 ;  % sat 9/8
            gainSetup = 15.3 ;  % sun 9/9
            gainSetup = 15.4 ;  % mon 9/10
            gainSetup = 15.5 ;  % mon 9/10   FAST TEST
            gainSetup = 15.51 ;  % tue 9/11  FAST TEST
            gainSetup = 15.52 ;  % tue 9/11
            gainSetup = 15.5 ;  %      FAST TEST
            gainSetup = 15.511 ;  % tue 9/11  FAST TEST
            gainSetup = 15.21 ;  % tue 9/15
            gainSetup = 15.22 ;  % tue 9/15 and ~9/24
            gainSetup = 15.23 ;  % mon 9/24 zoom to non-optimal point
            %  gainSetup = 15.6 ;  % mon 9/10
            %     reading 2019 format IR file:
            % gainSetup = 16.23 ;  % sun 6/7/19 quick
            % gainSetup = 16.24 ;  % sun 6/09/19 full search
            gainSetup = 16.25 ;  % tue 6/11/19 half search
        end

        voltRat = 0.01; %[may override] use empty matrix to force automatic computation
        switch gainSetup
            case 1 %gainSetup .....................
                meanVoltTrial=0.6;
                %  [gainStart ,gainStop,nGainSearch] = deal( meanVoltTrial-.4, meanVoltTrial+.6 ,nGainSearch); %  L7-4 voltrat=0.01
                [gainStart ,gainStop,nGainSearch] = deal( meanVoltTrial-.3, meanVoltTrial+.3 ,nGainSearch); %  L7-4 voltrat=0.01

            case {2,'hi-low.1'} %gainSetup .....................
                [gainStart ,gainStop,nGainSearch]=deal( 1 ,   2 ,2);
                meanVoltTrial = mean([gainStart ,gainStop]);
            case {3,'hi-low.3'} %gainSetup .....................
                [gainStart ,gainStop,nGainSearch]=deal( 2 ,   8 ,2);
                meanVoltTrial = mean([gainStart ,gainStop]);
            case {3.1,'hi-low.4','amplitude modulation example'} %gainSetup .....................
                [gainStart ,gainStop,nGainSearch]=deal(0.5, 1 , 3);
                meanVoltTrial = mean([gainStart ,gainStop]);
            case {3.2,'hi-low.4b','amplitude modulation example search'} %gainSetup .....................
                [gainStart ,gainStop,nGainSearch]=deal(0.4 ,   1.1 , 5);
                meanVoltTrial = mean([gainStart ,gainStop]);
            case {4,'single.1'} %gainSetup .....................
                [gainStart ,gainStop,nGainSearch]=deal( 0.885 ,   0.885 ,1);
                meanVoltTrial = mean([gainStart ,gainStop]);
            case {5,5.1,'search.1'} %gainSetup .....................
                [gainStart ,gainStop,nGainSearch]=deal( 0.4 ,   1.5 ,nGainSearch);
                meanVoltTrial = mean([gainStart ,gainStop]);
            case {15.1,'search.15a'} %gainSetup .....................
                %  nGainSearch = 3; %temp debug
                %     nGainSearch = 400;

                [gainStart ,gainStop,nGainSearch]=deal( 0.001 ,   1.5 ,nGainSearch);
                meanVoltTrial = mean([gainStart ,gainStop]);
            case {15.2,'search.15b0'} %gainSetup .....................
                %  nGainSearch = 3; %temp debug
                nGainSearch = 400;

                [gainStart ,gainStop,nGainSearch]=deal( 0.001 ,  0.9 ,nGainSearch);
                meanVoltTrial = mean([gainStart ,gainStop]);
            case {15.21,'search.15b1'} %gainSetup .....................
                %  nGainSearch = 3; %temp debug
                nGainSearch = 100;

                [gainStart ,gainStop,nGainSearch]=deal( 0.01 ,  2.0 ,nGainSearch);
                meanVoltTrial = mean([gainStart ,gainStop]);
            case {15.22,'search.15b2'} %gainSetup .....................
                %  nGainSearch = 3; %temp debug
                nGainSearch = 200;

                [gainStart ,gainStop,nGainSearch]=deal( 0.001 ,  0.25 ,nGainSearch);
                meanVoltTrial = mean([gainStart ,gainStop]);

            case {15.23,'search.15b3'} %gainSetup .....................
                %  nGainSearch = 3; %temp debug
                nGainSearch = 3;
                %zoom in secondary local optimum (not the global):
                [gainStart ,gainStop,nGainSearch]=deal( 0.008 ,  0.011 ,nGainSearch);
                meanVoltTrial = mean([gainStart ,gainStop]);
            case {15.3,'search.15c'} %gainSetup .....................
                %  nGainSearch = 3; %temp debug
                nGainSearch = 200;

                [gainStart ,gainStop,nGainSearch]=deal( 0.001 ,  0.3 ,nGainSearch);
                meanVoltTrial = mean([gainStart ,gainStop]);
            case {15.4,'search.15d'} %gainSetup .....................
                %  nGainSearch = 3; %temp debug
                nGainSearch = 60;

                [gainStart ,gainStop,nGainSearch]=deal( 0.02 ,  0.06 ,nGainSearch);
                meanVoltTrial = mean([gainStart ,gainStop]);
            case {15.5,'search.15e','search.15e0'} %gainSetup .....................
                if isempty(nGainSearch)
                    nGainSearch = 3; %temp debug
                    % nGainSearch = 60;
                end
                [gainStart ,gainStop,nGainSearch]=deal( 0.03 ,  0.05 ,nGainSearch);
                meanVoltTrial = mean([gainStart ,gainStop]);
            case {15.51,'search.15e1'} %gainSetup .....................
                % nGainSearch = 3; %temp debug
                nGainSearch = 60;

                [gainStart ,gainStop,nGainSearch]=deal(0.032 ,   0.052 , nGainSearch );

                meanVoltTrial = mean([gainStart ,gainStop]);
            case {15.511,'search.15e1.fast'} %gainSetup .....................
                nGainSearch = 3;

                [gainStart ,gainStop,nGainSearch]=deal(0.032 ,   0.052 , nGainSearch );

                meanVoltTrial = mean([gainStart ,gainStop]);

            case {15.52,'search.15e2'} %gainSetup .....................
                nGainSearch = 3; %temp debug
                nGainSearch = 100;

                [gainStart ,gainStop,nGainSearch]=deal(0.030 ,   0.07 , nGainSearch );

                meanVoltTrial = mean([gainStart ,gainStop]);

            case {15.6,'search.15f'} %gainSetup .....................
                nGainSearch = 400;

                [gainStart ,gainStop,nGainSearch]=deal( 0.01 ,  0.1 ,nGainSearch);
                meanVoltTrial = mean([gainStart ,gainStop]);

            case {16.23,'search.16a'} %gainSetup .....................
                %  nGainSearch = 3; %temp debug
                voltRat = 1; %[may override] use empty matrix to force automatic computation
                nGainSearch = 3;
                %zoom in secondary local optimum (not the global):
                [gainStart ,gainStop,nGainSearch]=deal(1 ,  10 ,nGainSearch);
                meanVoltTrial = mean([gainStart ,gainStop]);
            case {16.24,'search.16b'} %gainSetup .....................
                %  nGainSearch = 3; %temp debug
                voltRat = 1; %[may override] use empty matrix to force automatic computation
                if isempty(nGainSearch)
                    nGainSearch = 20;
                end
                %zoom in secondary local optimum (not the global):
                [gainStart ,gainStop,nGainSearch]=deal(6 ,  20 ,nGainSearch);
                meanVoltTrial = mean([gainStart ,gainStop]);
            case {16.241,'search.16b1'} %gainSetup .....................
                %  nGainSearch = 3; %temp debug
                voltRat = 1; %[may override] use empty matrix to force automatic computation
                if isempty(nGainSearch)
                    nGainSearch = 30;
                end
                %zoom in secondary local optimum (not the global):
                [gainStart ,gainStop,nGainSearch]=deal(6 ,  20 ,nGainSearch);
                meanVoltTrial = mean([gainStart ,gainStop]);
            case {16.25,'search.16c'} %gainSetup .....................
                %  nGainSearch = 3; %temp debug
                voltRat = 1; %[may override] use empty matrix to force automatic computation
                nGainSearch = 14;
                %zoom in secondary local optimum (not the global):
                [gainStart ,gainStop,nGainSearch]=deal(13 ,  26 ,nGainSearch);
                meanVoltTrial = mean([gainStart ,gainStop]);


            otherwise  %gainSetup .....................
                error('bad switch')
        end  %gainSetup .....................
        gainSearchVec = linspace(gainStart ,gainStop,nGainSearch); %
        FtxClk = 250;

        %note:  some test signal methods below will normalize
        % the test signal before applying the following scaling:
        signalLevel=2; %test signal scaling

        %user settings: .....................

        %make the reference/design signal: - - - - - - - - -
        % discrete prolate spheroid sequence imaging pulse
        refSigPad = 'short.pulse.3'; %default - may be overwritten
        prePostAmbMeth = 'zPadTail';  %DEFAULT
        padMeth = prePostAmbMeth;  %DEFAULT
        %     padMeth = 'none';%DEFAULT
        %     padMeth = 'zPadBoth';%DEFAULT
        RedLenFact = 2; %(length adj.) for linear inv. post-proc to viterbi  %DEFAULT
        switch testsig
            case {'dpss.L7-4.nb','dpss.L7-4.1'} %testsig-----------------------
                xdcr='L7-4'; %
                time_halfbandwidth = 28.0; %narrower bandwidth
                FcProbe=5.8;
                seq_length = 512;
                num_seq = 2*(2.5)-1;
                regFact = 0.001; %larger value allows more stable solution but reduces accuracy.
                symbolType = 7; %defines Viterbi model: see "viterbi_symbol2.m"
            case {'dpss.L7-4.wb'}%testsig-----------------------
                xdcr='L7-4'; %
                time_halfbandwidth = 40.0; %wider bandwidth
                FcProbe=5.8;
                seq_length = 512;
                num_seq = 2*(2.5)-1;
                regFact = 0.001; %larger value allows more stable solution but reduces accuracy.
                symbolType = 7; %defines Viterbi model: see "viterbi_symbol2.m"
                %.............................
            case { 'cex-awe.L11-5v.1' , ...
                    'cex-awe.L11-5v.1a', 'cex-awe.L11-5v.1b' ... %half tx ap
                    } %testsig---------------
                if ~exist('refSigFilename','var')
                    refSigFilename=   'analogWavformDesign';
                    disp([mfilename,':using default: refSigFilename',refSigFilename])
                end
                refSigPad = 'short.pulse.1'; %temp %%%%%%%%%%%%%%%%%%%%%%%

                xdcr='L11-5v'; %

                %define Fc of probe so we can make the modulation band for our reference waveform.
                timeSupport = 'long';
                regFact = 0.001; %larger value allows more stable solution but reduces accuracy.

                symbolType = 41; %defines Viterbi model: see "viterbi_symbol2.m"

                %marshall the design waveform, PRI and TX channel
                if ~exist('testSigFileDefinition','var')
                    disp([mfilename,':using default: testSigFileDefinition',testSigFileDefinition])
                    testSigFileDefinition = struct;
                    testSigFileDefinition.filename =refSigFilename;
                    AWDFile = load(refSigFilename);
                    %   NEED TO PARSE FROM "testsig"
                    testSigFileDefinition.TXchannel = [];
                    testSigFileDefinition.PRI = [];
                    NTXCH = length(AWDFile.sigCEXAWE);
                    NPRI = length(AWDFile(1).NumPRI);

                end

            case {'dpss.L11-5.wb.1'} %testsig---------------
                %%% Trans.frequency = 7.8125 ; %default
                %      refSigPad = 'imaging.pulse.1';
                refSigPad = 'short.pulse.1' %temp %%%%%%%%%%%%%%%%%%%%%%%

                xdcr='L11-5v'; %

                %define Fc of probe so we can make the modulation band for our reference waveform.
                NtxClocksPerCycle = 32;
                FcProbe=1/(NtxClocksPerCycle*(1/FtxClk)); %MHz
                timeSupport = 'long';
                %     timeSupport = 'short';  %temp %%%%%%%%%%%%%%%%%%%%%%%

                switch timeSupport
                    case {'short'}
                        %   seq_length = 160; time_halfbandwidth = 10.0; %wider bandwidth
                        %    seq_length = 120; time_halfbandwidth = 6.0; %magenta wider bandwidth
                        %   seq_length = 120; time_halfbandwidth =20.0; %green wider bandwidth
                        seq_length = 120; time_halfbandwidth = 3.0; %yellow wider bandwidth
                    case {'long','standard'}
                        seq_length = 512; time_halfbandwidth = 50.0; %wider bandwidth
                    case 'medium'
                        seq_length = 256; time_halfbandwidth = 25.0; %wider bandwidth
                    otherwise
                        error('bad switch')
                end
                num_seq = 2*(2.5)-1;
                regFact = 0.001; %larger value allows more stable solution but reduces accuracy.

                symbolType = 41; %defines Viterbi model: see "viterbi_symbol2.m"
                %    symbolType = 42; %%temp %%%%%%%%%%%%%%%%%%%%%%%
                %..............................
            case {'dpss.L22-14.wb.1'}%testsig-----------------------
                xdcr='L22-14vX'; %
                time_halfbandwidth = 120.0; %wider bandwidth
                %define Fc of probe so we can make the modulation band for our reference waveform.
                NtxClocksPerCycle = 14;
                FcProbe=1/(NtxClocksPerCycle*(1/FtxClk)); %MHz
                seq_length = 512;
                num_seq = 2*(2.5)-1;
                regFact = 0.001; %larger value allows more stable solution but reduces accuracy.
                symbolType = 31; %defines Viterbi model: see "viterbi_symbol2.m"
            case {'dpss.L22-14.wb.2'}%testsig-----------------------
                xdcr='L22-14vX'; %
                time_halfbandwidth = 90.0; %wider bandwidth
                %define Fc of probe so we can make the modulation band for our reference waveform.
                NtxClocksPerCycle = 14;
                FcProbe=1/(NtxClocksPerCycle*(1/FtxClk)); %MHz
                seq_length = 512;
                num_seq = 2*(2.5)-1;
                regFact = 0.001; %larger value allows more stable solution but reduces accuracy.
                symbolType = 33; %defines Viterbi model: see "viterbi_symbol2.m"
            case {'dpss.L22-14.wb.3'}%testsig-----------------------
                xdcr='L22-14vX'; %
                time_halfbandwidth = 50.0; %wider bandwidth
                %define Fc of probe so we can make the modulation band for our reference waveform.
                NtxClocksPerCycle = 16;
                FcProbe=1/(NtxClocksPerCycle*(1/FtxClk)); %MHz
                seq_length = 512;
                num_seq = 2*(2.5)-1;
                regFact = 0.001; %larger value allows more stable solution but reduces accuracy.
                symbolType = 33; %defines Viterbi model: see "viterbi_symbol2.m"
            case {'dpss.L22-14.wb.4'}%testsig-----------------------
                xdcr='L22-14vX'; %
                time_halfbandwidth =50.0; %wider bandwidth
                %define Fc of probe so we can make the modulation band for our reference waveform.
                NtxClocksPerCycle = 16;
                FcProbe=1/(NtxClocksPerCycle*(1/FtxClk)); %MHz
                seq_length = 512;
                num_seq = 2*(2.5)-1;
                regFact = 0.001; %larger value allows more stable solution but reduces accuracy.
                symbolType = 31; %defines Viterbi model: see "viterbi_symbol2.m"
            case {'dpss.L22-14.wb.4b'}%testsig-----------------------
                xdcr='L22-14vX'; %
                time_halfbandwidth =50.0; %wider bandwidth
                %define Fc of probe so we can make the modulation band for our reference waveform.
                NtxClocksPerCycle = 16;
                FcProbe=1/(NtxClocksPerCycle*(1/FtxClk)); %MHz
                seq_length = 512;
                num_seq = 2*(2.5)-1;
                regFact = 0.001; %larger value allows more stable solution but reduces accuracy.
                symbolType = 34; %defines Viterbi model: see "viterbi_symbol2.m"
            case {'dpss.L22-14.wb.4c'}%testsig-----------------------
                xdcr='L22-14vX'; %
                time_halfbandwidth =50.0; %wider bandwidth
                %define Fc of probe so we can make the modulation band for our reference waveform.
                NtxClocksPerCycle = 16;
                FcProbe=1/(NtxClocksPerCycle*(1/FtxClk)); %MHz
                seq_length = 512;
                num_seq = 2*(2.5)-1;
                regFact = 0.001; %larger value allows more stable solution but reduces accuracy.
                symbolType = 35; %defines Viterbi model: see "viterbi_symbol2.m"
            case {'dpss.L22-14.wb.4d'}%testsig-----------------------
                xdcr='L22-14vX'; %
                time_halfbandwidth =50.0; %wider bandwidth
                %define Fc of probe so we can make the modulation band for our reference waveform.
                NtxClocksPerCycle = 16;
                FcProbe=1/(NtxClocksPerCycle*(1/FtxClk)); %MHz
                seq_length = 512;
                num_seq = 2*(2.5)-1;
                regFact = 0.001; %larger value allows more stable solution but reduces accuracy.
                symbolType = 32; %defines Viterbi model: see "viterbi_symbol2.m"
            case {'dpss.L22-14.wb.5'}%testsig-----------------------
                xdcr='L22-14vX'; %
                time_halfbandwidth =35.0; %narrower bandwidth
                %define Fc of probe so we can make the modulation band for our reference waveform.
                NtxClocksPerCycle = 16;
                FcProbe=1/(NtxClocksPerCycle*(1/FtxClk)); %MHz
                seq_length = 512;
                num_seq = 2*(2.5)-1;
                regFact = 0.001; %larger value allows more stable solution but reduces accuracy.
                symbolType = 32; %defines Viterbi model: see "viterbi_symbol2.m"
            case {'lfm.L7-4.1' ,'lfm.L7-4.2' ,'lfm.L7-4.3' } %testsig %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                LFM_test_signal_setup;   %external script
            case {'file-ir.L22-14.c#18'} %testsig --------------------------------------
                xdcr='L22-14vX'; %
                testSigFileDefinition = struct;
                testSigFileDefinition.filename = 'ir_L2214vx_3005_c.mat';
                NtxClocksPerCycle = 16; %known from file
                FcProbe=1/(NtxClocksPerCycle*(1/FtxClk)); %MHz
                [~,tempstrXXX] = strtok(testsig,'#');
                testSigFileDefinition.channel = str2num(tempstrXXX(2:end));
                testSigFileDefinition.irVarName = 'Hall';
                testSigFileDefinition.Fs = 4*FcProbe;
                regFact = 0.001; %larger value allows more stable solution but reduces accuracy.
                symbolType = 32; %defines Viterbi model: see "viterbi_symbol2.m"

            case {'file-ir.L22-14.d#18'} %testsig --------------------------------------
                %      xdcr='L22-14vX'; %
                xdcr='L22-14vX:3005d'; %

                testSigFileDefinition = struct;
                testSigFileDefinition.filename = 'ir_L2214vx_3005_d.mat';
                NtxClocksPerCycle = 16; %known from file
                FcProbe=1/(NtxClocksPerCycle*(1/FtxClk)); %MHz
                [~,tempstrXXX] = strtok(testsig,'#');
                testSigFileDefinition.channel = str2num(tempstrXXX(2:end));
                testSigFileDefinition.irVarName = 'Hall';
                testSigFileDefinition.Fs = 4*FcProbe;
                regFact = 0.001; %larger value allows more stable solution but reduces accuracy.
                symbolType = 32; %defines Viterbi model: see "viterbi_symbol2.m"

            otherwise%testsig-----------------------
                error('bad ref. signal switch')
        end%testsig-----------------------

        [testSignalType,transducerSpecString]=strtok(testsig ,'.');
        refSigSpec = testSignalType; %legacy variable
        disp([mfilename,': xdcr: ',xdcr])
        switch testSignalType

            case {'dpss'}% %------------------
                %Obtain DPSS
                [dps_seq,lambda] = dpss(seq_length,time_halfbandwidth,num_seq);
                %modulate the Slepian sequences
                Fs = FtxClk;
                ttt=[1:size(dps_seq,1)];ttt=ttt-mean(ttt);
                sigEnvelope=dps_seq(1:end,1);
                sigCarrier=exp(2*pi*1i*ttt'*FcProbe/Fs);
                p2=sigEnvelope.*sigCarrier;
                %     freqz(p2,1,linspace(0,13,2222),Fs)
                x250_normalized = real(p2)/max(abs(p2));
                x250=signalLevel*x250_normalized;

            case {'cex-awe'} % ------------------------
                disp([mfilename,': CEX-AWE coded excitation+arb.wave.'])
                testSigFileContents=load(testSigFileDefinition.filename);
                testSigFileContents= testSigFileContents.sigCEXAWE;
                %test signal:
                p2 = testSigFileContents(...
                    testSigFileDefinition.TXchannel ...
                    ).modulatedSignal{testSigFileDefinition.PRI} ;

                NtxClocksPerCycle =  ...
                    testSigFileContents(1).symbolWaveform.parameters.LsymClocks ;
                FcProbe=1/(NtxClocksPerCycle*(1/FtxClk)); %MHz
                Fs = FtxClk;
                FsFile = round(1/diff(...
                    testSigFileContents(1).waveshape.parameters.support([1 2]...
                    )));
                USfactor = FtxClk/FsFile;  %need to upsample to get to 250
                x250_normalized = resample(p2,Fs,FsFile);

                x250 = signalLevel*x250_normalized ;

            case {'lfm'} % ------------------------
                disp([mfilename,': LFM setup'])
                Fs = FtxClk;

            case {'file-ir'} %------------------
                %use the transducer impulse response (from file) sas the signal to
                %generate (therefore ensuring at least one scaling will
                %produce perfect synthesis, using a single excitation pulse)
                testSigFileContents=load(testSigFileDefinition.filename);
                hall = testSigFileContents.(testSigFileDefinition.irVarName);
                channelXXXX = testSigFileDefinition.channel;
                p2 = hall(:,channelXXXX);
                p2  = double(real(p2)/max(abs(p2)));
                USfactor = FtxClk/testSigFileDefinition.Fs;  %need to upsample to get to 250
                x250_normalized = resample(p2,USfactor,1);
                x250 = signalLevel*x250_normalized ;

            otherwise
                error('unknown test signal type')
        end

        %setup reference signal zero padding:
        refSignal.x250 = x250 ;
        switch refSigPad
            case 1
                refSignal.offsetFactor = 0.25;
                refSignal.endpadFactor =   2;
            case 2
                refSignal.offsetFactor = 0.5;
                refSignal.endpadFactor =   4; %
            case {3,'short.pulse.1'}
                refSignal.offsetFactor = 1.5;
                refSignal.endpadFactor =   8;
            case {4,'short.pulse.2'}
                refSignal.offsetFactor = 3;
                refSignal.endpadFactor =   10;

            case {4.2,'short.pulse.3'}
                refSignal.offsetFactor = 8;
                refSignal.endpadFactor =   16;
            case {5,'centered.pulse.1'}
                refSignal.offsetFactor = 4;
                refSignal.endpadFactor =   6;
            case {6,'centered.pulse.2'}
                refSignal.offsetFactor = 1;
                refSignal.endpadFactor =   2;
            case {10.1,'imaging.pulse.1','no-zpad','none'}
                refSignal.offsetFactor = 0;
                refSignal.endpadFactor = 0;
            case {10.1,'imaging.pulse.1','LTBWP'}
                refSignal.offsetFactor = 0;
                refSignal.endpadFactor = 1;
            otherwise
                error('bad switch')
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %set factory parameter adjustments: These parameters are "not user-servicable"
        %zero padding method
        factMeth = 'backslash' ; %computation method
        viterbiMethod=   'approximate.resample' ; %legacy, symbol rate equalizer setup/usage

        %put parameters into structure
        vitParams=struct;
        vitParams.method =  viterbiMethod;
        vitParams.symbolRateSetup = symbolType; %defines clocks per pulse (see viterbi_symbol2.m)
        vitParams.refSignal = refSignal;  %waveform definition
        vitParams.amblesCode = prePostAmbMeth; %preamble and postamble
        vitParams.regFact = regFact;  %regularization factor (set to << 1.0)
        vitParams.gVec = gainSearchVec; %output scaling/gain (nuisance variable)
        vitParams.fMatchRatio = voltRat; % Vout/peakImpulseResponse
        vitParams.solverMeth = factMeth; % factorization method
        vitParams.cirFileCode = xdcr;  % transducer impulse response file name
        vitParams.reductionLengthFactor = RedLenFact;  % transducer impulse response file name
        vitParams.Nfir = []; %length of viterbi model/factorized pulse
        %insert errors in IR assumption:
        vitParams.vitEvalIRSource = vitEvalIRSource; %allows perf. eval. of different IR from that used for VA computation

        timeStart = datestr(now);
        executionTime.start = timeStart;

        %CALL encoder and preprocessing:
        vitResult= viterbi_factor2(vitParams); %CALL VA %%%%%%%%%%%%%%

        executionTime.stop = datestr(now);
        vitResult.executionTimeTotal = executionTime;

        %        resultsSaveName = 'encode_dpss_result';
        resultsSaveName = 'encode_va_result';  %changing name b.c. no longer only DPSS test waveform.
        if figureCleanup
            close all
            disp([mfilename,':figureCleanup: clearing functions: '])
            funcClearList = {'viterbi_factor2'}
            for kcf=1:length(funcClearList)
                clearingThisFunc = funcClearList{kcf};
                disp([mfilename,':clearing function: ',clearingThisFunc])
                eval(['clear ',clearingThisFunc])
            end
        end

    otherwise %testsig %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp([mfilename,': bad switch: testsig:  ',testsig])
        error('bad switch' )

end %testsig %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%post-processing and save results:

txTri250 = vitResult.txTriState250;  %transmitter input sequence


%The synthesized waveform:
ArbWav=filter(vitResult.paramSet.h250,1,double(txTri250)); %Achieved waveform
refZPadded = refSignal.x250;
if length(ArbWav)>length(refZPadded)
    refZPadded(length(ArbWav))=0;
end
[errS]=nrmseVec_VS(ArbWav,refZPadded);  %align and scale to compare with original.

if showFigures
    %show reconstructed signal (after transducer filtering):
    cfhand0Save = gcf;
    cfhand = figure(2000);

    hold off
    plot(errS.alignment.xa,'r')
    hold on
    plot(errS.alignment.xerr*errS.alignment.gain,'k')
    title('Waveform encoding')
    xlabel('TX clock sample')
    legend('Encoder Result','Error')

    %show performance NRMSE vs. gain:
    figure(cfhand.Number+1)
    plot(...
        vitResult.paramSet.gainVec,vitResult.paramSet.nrmseVec,'.-b', ...
        vitResult.paramSet.gainVec,vitResult.paramSet.nrmseVecAdj,'o-r')
    grid on
    title('Encoder Performance vs. Voltage Gain')
    legend('NRMSE','NRMSE(adjusted)','Location','NorthEast')
    xlabel('Relative Voltage Gain')
    ylabel('Normalized RMS Error (dB)')
    disp([mfilename,': restoring previous figure...'])
    figure(cfhand0Save)

end

%capture parameters and results
disp([mfilename,': save results ...'])
disp([mfilename,'Encoded Transmitter sequence: txTri250'])
disp([mfilename,'Synthesized waveform: ArbWav'])
disp([mfilename,': resultsSaveName: ',resultsSaveName])

save(resultsSaveName)

disp([mfilename,':',' done.'])


