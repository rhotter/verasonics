% MultiTWPulseCode.m is used to generate PulseCode for multiple waveforms.
% For instance, if the user would like to transmit different waveform
% between channels, the user can save the waveform into a .mat file with the
% format [waveform, number of waveforms]. This file will generate the pulse
% code for each waveform

% 1. Set important parameters, including
% Transducer: enter transducer name
% IR_Option: Three impulse response otpions. Only L7-4 has factory setting
% Symbol period, Encoder model, VS1 and VS2 search
% 2. Load waveforms
% 3. Generate Pulse Code

%% Parameters need to be modified

Trans.name = 'L7-4';
Trans.units = 'mm';
Trans = computeTrans(Trans);

bandwidth = 0.6;  % Bandwidth

IR_Option = 1;  % 1: nominal butterworth, 2: customized IR, 3. Factory Setting
EncoderModel = 0; % 0:Probe compensation, 1: DAC synthesis

%VS1 encoder level search --------------------
VS1Level = [0.2 2];
VS_Search = linspace(VS1Level(1),VS1Level(2),200);

postProcessing = 1; % 0: desable 1: enable the VS2 post processing
%--Parameters for VS2 post processing
VS2Level = [0.9 1.1];
maxIter = 5;
Nkv = 3;
%--End of post processing parameters

% Pre allocate the length for tri-state pulse and pulsecode
triStateLength = 5000;
pulseCodeLength = 500;

%% Load multiple waveforms

[FileName,PathName,FilterIndex] = uigetfile({'*.mat;*.txt'},'Select the customized waveform (*.mat) file');
if PathName == 0 %if the user pressed cancelled, then we exit this callback
    disp('Waveform file is not selected');
    return
else
    WaveformPath = [PathName,FileName];
    filetype = FileName(end-2:end);

    switch filetype

        case 'mat'
            S = load(WaveformPath);
            fn = fieldnames(S);
            if isnumeric(S.(fn{1}))
                multiWaveforms = S.(fn{1});
            else
                display('Incorrect waveform file!');
                return
            end

        otherwise
            display('Not supported format! Only .mat file is supported');
            return
    end
end

%% Parameters for pulsecode generation
setVerbose;
if postProcessing == 0
    enablePP = 0;
else
    enablePP = ['maxIter =',num2str(maxIter),';Nkv =',num2str(Nkv),...
        ';RVMod = linspace(',num2str(VS2Level(1)),',',num2str(VS2Level(2)),',Nkv);'];
end

Fc = Trans.frequency;

switch IR_Option

    case 1
        % -------------------------------------
        %make impulse response (nominal example) or use experimentally
        %determined impulse response for transducer element.

        BW = Fc * 0.6;

        Ncycles = 5;
        delaySpreadm30db_microSec = Ncycles/Fc;

        buttOrd = 2;
        BL_MHz = Fc+[-1 1]*BW/2;
        [btrans,atrans]=butter(buttOrd, 2*BL_MHz/250);  % coeff of butterworth filter

        Nsamp = delaySpreadm30db_microSec*250;
        htrans250 = impz(btrans,atrans,Nsamp);   % htrans250 is nominal impulse response

        %put into cir file format:
        cir_file1 =  [Trans.name,'_hfiletemp'];
        makeUserCIRFile('ieee2013',cir_file1,htrans250,250,[Trans.name,'.nominal']);  %writes the cir file used later

    case 2
        % -------------------------------------
        % Customized impulse response

        if isfield(Trans,'bandwidth')
            Trans = rmfield(Trans,'bandwidth');
        end
        BW = [];

        [FileName,PathName,~] = uigetfile({'*.mat'},'Select the customized impulse response (*.mat) file');
        if PathName == 0 %if the user pressed cancelled, then we exit this callback
            msgbox('no file is selected, PulseCode is NOT generated!');
            return
        else
            S = load([PathName,FileName]);
            fn = fieldnames(S);
            htrans250 = S.(fn{1});
        end

        %put into cir file format:
        cir_file1 =  [Trans.name,'_hfiletemp'];
        makeUserCIRFile('ieee2013',cir_file1,htrans250,250,[Trans.name,'.nominal']);  %writes the cir file used later

        cirS = load(cir_file1,'cirS');cirS = cirS.cirS;
        delaySpreadm30db_microSec =  cirS.dimensions.Nparam/ cirS.timing.TXclock;

    case 3
        % -------------------------------------
        % Factory impulse response, only L7-4 is supported
        if isfield(Trans,'bandwidth')
            Trans = rmfield(Trans,'bandwidth');
        end
        BW = [];

        ProbeName = Trans.name;
        switch ProbeName
            case 'L7-4'
                cir_file1 = 'cirest0701_code_acoustic-L74-23-tiny';
            otherwise
                msgbox('The impulse response of the probe you select does not exist, pleae select ButterWorth or use customized impulse response')
                return
        end

        cirS = load(cir_file1,'cirS');cirS = cirS.cirS;
        delaySpreadm30db_microSec =  cirS.dimensions.Nparam/cirS.timing.TXclock;
        htrans250 = cirS.CIR.gateSignalModel;

end

% Symbol period calculation
SymbolPeriod = round(250/Trans.frequency/4);
if SymbolPeriod < 10
    SymbolPeriod = 10;
end

%transducer specification
transSpec.Fc  = Fc;
transSpec.cirEstModelName  = '';
transSpec.BW  = BW;
transSpec.delaySpreadm20db_microSec  =[];
transSpec.delaySpreadm30db_microSec  = delaySpreadm30db_microSec; % may overwrite below
transSpec.windowIndexVec = [];
transSpec.modelName = 'external'; %defined in a file
transSpec.ACcouple = 3; %1=mean;2=lin;3 = quadratic, etc.
transSpec.tag = 'test';
transSpec.cirFileName = cir_file1;
transSpec.cirEstModelName  = 'gateSignalModel';

%transducer signal: used by wavesim.m
transSig.transducerSignalName = 'signal.sampled';
transSig.transducerSpec = transSpec;

%symbol set design
pulseSetLFM = 'nonaligned.1.augmented';

%encoding file contents (variable names):
saveVarList = {'pulseSeq','bestRSSI','optParm','wsInfo','pulseLFM','pc','transSpec'};

%% Generate the pulsecode for each waveform
sizeOfWave = size(multiWaveforms);
numOfWaveform = sizeOfWave(2);

triState = zeros(triStateLength,numOfWaveform);
PulseCode = zeros(pulseCodeLength,5,numOfWaveform);

h = waitbar(0,'Generating PulseCode...');
PClengthMax = 0;

for num = 1:numOfWaveform
    transSig.data = multiWaveforms(:,num);
    if nnz(transSig.data) > 0
        fprintf('\n')
        disp('---------------------------------')
        disp(['Waveform number: ',num2str(num)])
        setVerbose;
        try [pulseSeq,~,~,wsInfo]=wavesim('design.eq.1', transSig, 'gen3.b',...
                SymbolPeriod, EncoderModel, pulseSetLFM , enablePP , VS_Search);
        catch msg
            close(h)
            return
        end
        if postProcessing == 0
            triState(1:length(pulseSeq.MMSE),num) = pulseSeq.MMSE';
        else
            triState(1:length(pulseSeq.DFE),num) = pulseSeq.DFE';
        end

        [PC,~] = twgen('trinary2gen3',triState(:,num));

        if PClengthMax < length(PC)
            PClengthMax = length(PC);
        end

        PulseCode(1:length(PC),:,num) = PC;
        displayValue;
    end
    waitbar(num/numOfWaveform);
end

close(h);
clear TW; % clear existing TW structure
TW.type = 'pulseCode';
TW.PulseCode = PulseCode(1:PClengthMax,:,:);

% Check limit and other parameters calculation
Resource.Parameters.verbose = 0;
[~, ~, ~, ~, TW] = computeTWWaveform(TW);

fprintf('\n');
disp('Please save TW structure and load it in the SetUp script.');

