function filterTool(varargin)
% Copyright Verasonics, Inc.  All world-wide rights and remedies under all
% intellectual property laws and industrial property laws are reserved.
% Verasonics Registered U.S. Patent and Trademark Office.
%
% Notice:
%   This file is provided by Verasonics to end users as a programming
%   tool for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: filterTool.m
%
% Verasonics Provides a Filter Design Tool to
% 1) Visualize the default filter and design the customized filter in off-line mode
% 2) Program the decimRate, sampleMode, low-pass and band-pass filter
%    of the Vantage system, as well as visualize the effect on image in real-time
%
% 07-Jun-2019 --- Multiple fields should be used to determine the filter sets
% 07-Nov-2016 --- Support multiple filter sets based on the InputFilter

% If filterTool is a called by EventAnalysisTool, plot the LPF and BPF
% from a specfic rcvNum indicated in the Event Table.

% How many filter Sets? Default is 1
rcvInd = 1;
setInd = 1;
numOfSet = 1;
if evalin('base','exist(''Receive'',''var'')')
    Receive = evalin('base','Receive');
    setInd = ones(size(Receive,2),1);
    if isfield(Receive,'InputFilter')
        [~,rcvIndTemp,setIndTemp] = unique(cellfun(@(x)num2str(x),{Receive.InputFilter},'UniformOutput',false)');
        if length(rcvIndTemp)>length(rcvInd), rcvInd = rcvIndTemp; setInd = setIndTemp; end
    end
    if isfield(Receive,'LowPassCoef')
        [~,rcvIndTemp,setIndTemp] = unique(cellfun(@(x)num2str(x),{Receive.LowPassCoef},'UniformOutput',false)');
        if length(rcvIndTemp)>length(rcvInd), rcvInd = rcvIndTemp; setInd = setIndTemp; end
    end
    if isfield(Receive,'decimSampleRate')
        [~,rcvIndTemp,setIndTemp] = unique([Receive.decimSampleRate],'stable');
        if length(rcvIndTemp)>length(rcvInd), rcvInd = rcvIndTemp; setInd = setIndTemp; end
    end
    if isfield(Receive,'demodFrequency')
        [~,rcvIndTemp,setIndTemp] = unique([Receive.demodFrequency],'stable');
        if length(rcvIndTemp)>length(rcvInd), rcvInd = rcvIndTemp; setInd = setIndTemp; end
    end
    % sampleMode is required
    if isfield(Receive,'sampleMode')
        [~,rcvIndTemp,setIndTemp] = unique({Receive.sampleMode}','stable');
        if length(rcvIndTemp)>length(rcvInd), rcvInd = rcvIndTemp; setInd = setIndTemp; end
    else
        error('smapleMode is required for the Receive structure.');
    end
    numOfSet = length(rcvInd);
end

% Get rcvEvent from EventAnalysisTool?
% Check which set the rcvEvent belongs to.  Default is one
setValue = 1;
rcvEvent = 1;
if ~isempty(varargin) && (length(varargin) == 3)
    rcvEvent = varargin{3};
    setValue = setInd(rcvEvent); % determine which set will be displayed
end

% Design the filter on the fly? The decimSampleRate can't be adjusted with
% simulation mode on the fly
if ~isempty(findobj('tag','UI'))
    runningVSX = 1;
else
    runningVSX = 0;
end

% legend function gets changed since 2017a, V9.2
if verLessThan('Matlab','9.2')
    AutoUpdateLegend = 0;
else
    AutoUpdateLegend = 1;
end

%% Create supported frequency table for Vantage
VDASG3 = [];
VDASG3FourThird = [];
VDASInterleave = [];
DefaultLPF = [];
DefaultBPF = [];
setDefaultValue();

%% Retrive value from workspace for multiple sets of RcvParam
LPFParam = [];
BPFParam = [];
RcvParam = [];
retriveParam();

%% Create the GUI components
UI = [];
initializeUI();

%% Plot LPF and BPF
computeLPF;
computeBPF;

%% nested functions

    function setDefaultValue()
        SR = zeros(74,2); % Compute 74 entries from 250MHz to 2.5MHz.
        k = 1;
        for i = 1:100
            f = 250/i;
            for j=8:-1:1, if (i/j - floor(i/j)) < .0001, break, end, end
            if (i/j <= 25), SR(k,1) = f; SR(k,2) = j; k = k+1; end
        end
        % Correct the decimation factors in the first part of the table (easier to specify than compute).
        SR(1:25,2) = [1,1,1,1,1,1,1,2,1,2,1,3,1,2,3,4,1,3,1,5,3,2,1,6,5];

        Idx67BW = sort([find(SR(:,2)==1); find(SR(:,2)==2)]);
        SR67BW = SR(Idx67BW,:);
        ADCRate = SR(:,1).*SR(:,2);

        VDASG3.decimRate = SR(4:end,1);
        VDASG3FourThird.decimRate = SR67BW(4:end,1);
        VDASInterleave.decimRate = SR(2:4,1);

        VDASG3.decimFactor = SR(4:end,2);
        VDASG3FourThird.decimFactor =  SR67BW(4:end,2);
        VDASInterleave.decimFactor = SR(2:4,2);

        VDASG3.ADCRate = ADCRate(4:end);
        VDASG3FourThird.ADCRate = SR67BW(4:end,1).*SR67BW(4:end,2);
        VDASInterleave.ADCRate = ADCRate(2:4)/2;

        VDASG3.demodFreq = VDASG3.decimRate/4;
        VDASG3FourThird.demodFreq = VDASG3FourThird.decimRate/4*3;
        VDASInterleave.demodFreq = VDASInterleave.decimRate/4;

        Nfft = 1024;
        % ==== Default LPF parameters for Vantage ====
        DefaultLPF.cutoff = [2.00 1.83 1.81 1.72 1.79 1.76 1.80 1.89 0.69];
        DefaultLPF.beta   = [0.0 3.6 2.82 2.10 3.00 1.80 1.00 0.70 3.59];

        DefaultLPF.coef = [...
            [ +0.0000  +0.0000  +0.0000  +0.0000  +0.0000  +0.0000  +0.0000  +0.0000  +0.0000  +0.0000  +0.0000  +1.0000];...
            [ -0.0004  +0.0063  +0.0038  -0.0139  -0.0137  +0.0234  +0.0359  -0.0328  -0.0868  +0.0398  +0.3105  +0.4559];...
            [ -0.0057  -0.0005  +0.0115  +0.0197  +0.0095  -0.0208  -0.0497  -0.0411  +0.0285  +0.1444  +0.2545  +0.2997];...
            [ +0.0110  +0.0072  -0.0042  -0.0202  -0.0334  -0.0337  -0.0128  +0.0311  +0.0914  +0.1538  +0.2008  +0.2183];...
            [ -0.0006  -0.0058  -0.0130  -0.0191  -0.0192  -0.0084  +0.0161  +0.0531  +0.0972  +0.1394  +0.1699  +0.1810];...
            [ -0.0137  -0.0183  -0.0194  -0.0148  -0.0030  +0.0163  +0.0419  +0.0710  +0.1001  +0.1248  +0.1414  +0.1472];...
            [ -0.0212  -0.0198  -0.0139  -0.0031  +0.0124  +0.0315  +0.0528  +0.0744  +0.0942  +0.1101  +0.1205  +0.1240];...
            [ -0.0192  -0.0144  -0.0060  +0.0060  +0.0208  +0.0376  +0.0552  +0.0722  +0.0873  +0.0992  +0.1067  +0.1093];...
            [ +0.0030  +0.0034  -0.0093  -0.0068  +0.0214  +0.0106  -0.0441  -0.0141  +0.0932  +0.0166  -0.3135  +0.4820]];

        coefLPF = [DefaultLPF.coef fliplr(DefaultLPF.coef(:,1:end-1))];
        DefaultLPF.Fout = 20*log10(abs(fft(coefLPF',Nfft)))';

        % ==== BPF will be applied based on the samples per wavelength ====
        DefaultBPF.bndwdth = [1.0800   0.8300   0.3700   0.6000  4.0];
        DefaultBPF.xsnwdth = [0.4000   0.3400   0.3200   0.4000  0.4];

        DefaultBPF.coef = zeros(5,21);
        DefaultBPF.coef(1, :) = [ ...
            -0.00113 +0.00000 -0.00116 +0.00000 +0.00549 +0.00000 +0.00720 ...
            +0.00000 -0.01419 +0.00000 -0.02640 +0.00000 +0.02606 +0.00000 ...
            +0.07816 +0.00000 -0.03671 +0.00000 -0.30786 +0.00000 +0.54108 ];

        % BPF #2 (2 samples per wavelength, default coef row 2):
        DefaultBPF.coef(2, :) = [ ...
            +0.00034 +0.00000 +0.00244 +0.00000 -0.00629 +0.00000 -0.00333 ...
            +0.00000 +0.02188 +0.00000 -0.00897 +0.00000 -0.04745 +0.00000 ...
            +0.06076 +0.00000 +0.07294 +0.00000 -0.30048 +0.00000 +0.41632 ];

        % BPF #3 (1 sample per wavelength, default coef row 3):
        DefaultBPF.coef(3, :) = [ ...
            -0.00162 +0.00000 +0.00568 +0.00000 -0.01065 +0.00000 +0.01349 ...
            +0.00000 -0.00858 +0.00000 -0.00955 +0.00000 +0.04312 +0.00000 ...
            -0.08841 +0.00000 +0.13550 +0.00000 -0.17130 +0.00000 +0.18463 ];

        % BPF #4 (4/3 samples per wavelength aliased, default coef row 4):
        DefaultBPF.coef(4, :) = [ ...
            -0.00159 +0.00000 -0.00549 +0.00000 -0.01157 +0.00000 -0.02066 ...
            +0.00000 -0.03275 +0.00000 -0.04721 +0.00000 -0.06281 +0.00000 ...
            -0.07785 -0.00003 -0.09039 -0.00003 -0.09875 -0.00003 +0.89832 ];

        % BPF #5 (interleave, default coef row 5):
        DefaultBPF.coef(5, :) = [ ...
            +0.00000 +0.00214 +0.00000 -0.00409 +0.00000 +0.00693 +0.00000 ...
            -0.01093 +0.00000 +0.01654 +0.00000 -0.02457 +0.00000 +0.03665 ...
            +0.00000 -0.05713 +0.00000 +0.10217 +0.00000 -0.31735 +0.50067];

        coefBPF = [DefaultBPF.coef fliplr(DefaultBPF.coef(:,1:end-1))];
        DefaultBPF.Fout = 20*log10(abs(fft(coefBPF',Nfft)))';

    end

    function retriveParam()
        % Get RcvParam(setValue) name from base workspace. If not running with SetUp
        % script, the default Trans.frequency is 5.208 but can be modified in the GUI
        for n = 1:numOfSet
            if evalin('base','exist(''Trans'',''var'')')
                RcvParam(n).TransFreq = evalin('base','Trans.frequency');
            else
                RcvParam(n).TransFreq = 5.208;
                RcvParam(n).decimRate = 4*RcvParam(n).TransFreq;
            end

            % Get Receive from workspace, if LowPassCoef or InputFilter exists, the
            % filterTool will plot the frequency response of the filter. If not, the
            % filterTool will start from the default filter.

            RcvParam(n).interleave = 0;
            RcvParam(n).VDAS = VDASG3; % default

            if evalin('base','exist(''Receive'',''var'')')

                Receive = evalin('base','Receive');

                % Check for sampleMode or samplesPerWave attributes provided; if not found, set defaults
                if ~isfield(Receive(rcvInd(n)),'sampleMode')
                    if isfield(Receive(rcvInd(n)),'samplesPerWave') % for backwards compatibility
                        switch Receive(rcvInd(n)).samplesPerWave
                            case 4
                                Receive(rcvInd(n)).sampleMode = 'NS200BW';
                            case 2
                                Receive(rcvInd(n)).sampleMode = 'BS100BW';
                            case 4/3
                                Receive(rcvInd(n)).sampleMode = 'BS67BW';
                            case 1
                                Receive(rcvInd(n)).sampleMode = 'BS50BW';
                            otherwise
                                Receive(rcvInd(n)).sampleMode = 'custom';
                                Receive(rcvInd(n)).decimSampleRate = Receive(rcvInd(n)).samplesPerWave * RcvParam(n).TransFreq;
                        end
                    else
                        Receive(rcvInd(n)).sampleMode = 'NS200BW'; % default sampleMode if not provided.
                    end
                end

                RcvParam(n).sampleMode = Receive(rcvInd(n)).sampleMode;

                % Set target decimSampleRate - if provided, use it, otherwise use (4 or 4/3)*Trans.frequency.
                if isfield(Receive(rcvInd(n)),'decimSampleRate') && ~isempty(Receive(rcvInd(n)).decimSampleRate)
                    RcvParam(n).decimRate = Receive(rcvInd(n)).decimSampleRate;
                else
                    % if demodFrequency is provided, decimRate will be 4*demodFrequency
                    if isfield(Receive(rcvInd(n)),'demodFrequency') && ~isempty(Receive(rcvInd(n)).demodFrequency)
                        RcvParam(n).decimRate = 4*Receive(rcvInd(n)).demodFrequency;
                    elseif strcmp(Receive(rcvInd(n)).sampleMode,'BS67BW')
                        RcvParam(n).decimRate = (4/3)*RcvParam(n).TransFreq;
                    else
                        RcvParam(n).decimRate = 4*RcvParam(n).TransFreq;
                    end
                end

                switch RcvParam(n).sampleMode
                    case 'NS200BW'
                        RcvParam(n).BPFfilterNum = 1;
                        Receive(rcvInd(n)).samplesPerWave = RcvParam(n).decimRate/RcvParam(n).TransFreq;
                        RcvParam(n).sampleModeStr = {'NS200BW';'BS100BW';'BS50BW';'custom'};
                    case 'NS200BWI'
                        RcvParam(n).BPFfilterNum = 5;
                        RcvParam(n).interleave = 1;
                        RcvParam(n).VDAS = VDASInterleave;
                        Receive(rcvInd(n)).samplesPerWave = RcvParam(n).decimRate/(2*RcvParam(n).TransFreq);
                        RcvParam(n).sampleModeStr = {'NS200BWI';'custom'};
                    case 'BS100BW'
                        RcvParam(n).BPFfilterNum = 2;
                        Receive(rcvInd(n)).samplesPerWave = RcvParam(n).decimRate/(2*RcvParam(n).TransFreq);
                        RcvParam(n).sampleModeStr = {'NS200BW';'BS100BW';'BS50BW';'custom'};
                    case 'BS67BW'
                        RcvParam(n).VDAS = VDASG3FourThird;
                        RcvParam(n).BPFfilterNum = 4;
                        Receive(rcvInd(n)).samplesPerWave = RcvParam(n).decimRate/RcvParam(n).TransFreq;
                        RcvParam(n).sampleModeStr = {'BS67BW','custom'};
                    case 'BS50BW'
                        RcvParam(n).BPFfilterNum = 3;
                        Receive(rcvInd(n)).samplesPerWave = RcvParam(n).decimRate/(4*RcvParam(n).TransFreq);
                        RcvParam(n).sampleModeStr = {'NS200BW';'BS100BW';'BS50BW';'custom'};
                    case 'custom'
                        RcvParam(n).BPFfilterNum = 1;
                        Receive(rcvInd(n)).samplesPerWave = RcvParam(n).decimRate/RcvParam(n).TransFreq;
                        RcvParam(n).sampleModeStr = {'custom'};
                    otherwise
                        error('filterTool: Unrecognized Receive(%d).sampleMode.\n',rcvInd(n));
                end

                % Third, determind LPF and BPF
                if ~isfield(Receive(rcvInd(n)),'LowPassCoef')
                    LPFParam(n).LPFcoef = [];
                else
                    % user-defined filter exists; if of V1 size expand it to be
                    % compatible with Gen3

                    if size(Receive(rcvInd(n)).LowPassCoef, 2) == 6
                        Receive(rcvInd(n)).LowPassCoef = [zeros(1,6), Receive(rcvInd(n)).LowPassCoef];
                    end
                    LPFParam(n).LPFcoef = Receive(rcvInd(n)).LowPassCoef;
                end

                if ~isfield(Receive(rcvInd(n)),'InputFilter')
                    BPFParam(n).BPFcoef = [];
                else
                    if size(Receive(rcvInd(n)).InputFilter,2) == 6
                        % looks like a V1 filter provided by user so expand it to Gen3 size
                        Buffer = zeros(1, 21);
                        for ntap = 1:6
                            Buffer(1,9 + 2*ntap) = Receive(rcvInd(n)).InputFilter(:, ntap);
                        end
                        Receive(rcvInd(n)).InputFilter = Buffer;
                        clear Buffer
                    end

                    BPFParam(n).BPFcoef = Receive(rcvInd(n)).InputFilter;
                end

            else
                % default is 1 for NS200BW sampleMode
                RcvParam(n).BPFfilterNum = 1;
                RcvParam(n).sampleMode = 'NS200BW';
                RcvParam(n).sampleModeStr = {'NS200BW';'BS100BW';'BS50BW';'custom'};
                LPFParam(n).LPFcoef = [];
                BPFParam(n).BPFcoef = [];
            end

            RcvParam = freqCorrection(RcvParam,n);

            if ~runningVSX
                RcvParam(n).sampleModeStr = {'NS200BW';'NS200BWI';'BS100BW';'BS67BW';'BS50BW';'custom'};
            end

            % other parameters for LPF and BPF
            LPFParam(n).numTaps = 23;
            LPFParam(n).cutoffFreqMHz = RcvParam(n).demodFreq*DefaultLPF.cutoff(RcvParam(n).LPFfilterNum);
            LPFParam(n).kaiserBeta = DefaultLPF.beta(RcvParam(n).LPFfilterNum);
            LPFParam(n).DefaultCoef = DefaultLPF.coef(RcvParam(n).LPFfilterNum,:);
            LPFParam(n).DefaultFout = DefaultLPF.Fout(RcvParam(n).LPFfilterNum,:);
            LPFParam(n).modified = 0;

            BPFParam(n).numTaps = 41;
            BPFParam(n).allPass = 0;
            BPFParam(n).centerFreqMHz = RcvParam(n).demodFreq; if RcvParam(n).BPFfilterNum == 5, BPFParam(n).centerFreqMHz = RcvParam(n).demodFreq/2; end
            BPFParam(n).bndwdth    = DefaultBPF.bndwdth(RcvParam(n).BPFfilterNum); % initial value for each filter's bandwidth in units relative to Fc
            BPFParam(n).xsnwdth    = DefaultBPF.xsnwdth(RcvParam(n).BPFfilterNum); % initial value for width of each filter's transition bands, in units relative to Fc
            BPFParam(n).DCnulltap = [1 1 1 1]; % default tap position for DC zeroing
            BPFParam(n).Rangsel = [-90 10]; % Range to use for freq. response plots
            BPFParam(n).DefaultCoef = DefaultBPF.coef(RcvParam(n).BPFfilterNum,:);
            BPFParam(n).DefaultFout = DefaultBPF.Fout(RcvParam(n).BPFfilterNum,:);
            BPFParam(n).modified = 0;

        end
        assignin('base','RcvParam',RcvParam);
        assignin('base','LPFParam',LPFParam);
        assignin('base','BPFParam',BPFParam);
    end

    function initializeUI()
        % Define UIPos, which contains the default GUI positions. The origin
        % specified by UIPos is the lower left corner of a virtual box that encloses the control.

        ScrnSize = get(groot,'MonitorPositions');

        fWidth = 1300; fHighth = 600;

        if fWidth > ScrnSize(1,3) || fHighth > ScrnSize(1,4)
            fWidth = ScrnSize(1,3)*0.5;
            fHighth = ScrnSize(1,4)*0.45;
        end

        panelWidth = 0.15;
        panelHighth = 0.20;

        UIPos = zeros(8,2);
        UIPos(:,1) = 0.43;
        UIPos(:,2) = [0.82;0.74;0.56;0.49;0.40;0.27;0.18;0.10];

        % Each Pos means the left, button corner of
        % 1: Slider of TW frequency
        % 2: Slider of RcvParam(setValue).TransFreq
        % 3: Samples per wavelength
        % 4: Cutoff frequency
        % 5: Kaiser Beta
        % 6: Center frequency
        % 7: Relative bandwidth
        % 8: Filter Generation button

        % Define slider group offsets and sizes. All units are normalized.
        SG = struct('TO',[0.005,0.06],...   % title offset
            'TS',[0.14,0.03],...   % title size
            'TF',0.8,...           % title font size
            'SO',[0.01,0.014],...  % slider offset
            'SS',[0.08,0.035],...  % slider size
            'EO',[0.1,0.02],...   % edit box offset
            'ES',[0.04,0.032]);    % edit box size

        UI.figFont = 0.028;

        if ispc
            SG.TF = 0.75;
            SG.SO(2) = 0.018;
            SG.EO(2) = 0.018;
            SG.SS(2) = 0.028;
            SG.ES(2) = 0.030;
            UI.figFont = 0.025;
        end

        % Close any previously opened GUI windows.
        hf = findobj('tag','filterTool');
        if ishandle(hf)
            pos = get(hf,'Position');
            delete(hf);
        else
            pos = [(ScrnSize(1,3)-fWidth)/2,(ScrnSize(1,4)-fHighth)/2,fWidth,fHighth];
        end

        % Initialize and hide the GUI as it is being constructed.
        figFilter = figure('Visible','on',...  %'Units','normalized',...
            'Position',pos,... %'Position',[0.7,0.25,0.25,0.50],...
            'Name','Verasonics Filter Design Tool for Vantage system',...
            'NumberTitle','off',...
            'MenuBar','none', ...
            'Toolbar','figure',...
            'Resize','on', ...
            'tag','filterTool');

        % movegui(figFilter,'center');

        bkgrnd = get(figFilter,'Color');
        set(figFilter,'DefaultUicontrolBackgroundColor',bkgrnd)

        % Set UI control template of each txt, slider, and editbox set for each
        % position
        UItemp = repmat(struct('txt',1,'sldr',1,'edit',1),1,length(UIPos));
        for n = [1,4,5,6,7]%1:length(UIPos)-1
            UItemp(n).txt = uicontrol('Style','text',...
                'String','txt',...
                'Units','normalized',...
                'Position',[UIPos(n,:)+SG.TO,SG.TS],...
                'FontUnits','normalized',...
                'FontSize',SG.TF,...
                'FontWeight','bold');

            UItemp(n).sldr = uicontrol(figFilter,'Style','slider',...
                'Units','normalized',...
                'Position',[UIPos(n,:)+SG.SO,SG.SS],...
                'BackgroundColor',bkgrnd-0.05,...
                'Interruptible','off',...
                'BusyAction','cancel');

            UItemp(n).edit = uicontrol('Style','edit',...
                'Units','normalized',...
                'Position',[UIPos(n,:)+SG.EO,SG.ES],...
                'BackgroundColor',bkgrnd);
        end

        %% handle of two axes (LPF and BPF)
        if RcvParam(setValue).BPFfilterNum == 4
            txtLPF = 'HPF';
        else
            txtLPF = 'LPF';
        end

        if RcvParam(setValue).interleave
            txtBPF = 'HPF';
        else
            txtBPF = 'BPF';
        end

        UI.LPFtxt = uicontrol('Style','text',...
            'String',txtLPF,...
            'Units','normalized',...
            'Position',[0.05, 0.91, 0.05, 0.04],...
            'FontUnits','normalized',...
            'FontSize',0.8,...
            'Tag','UI.LPFtxt',...
            'FontWeight','bold');

        UI.LPFaxe = axes(...
            'Parent',figFilter,...
            'Units','normalized',...
            'Position',[0.05,0.1,0.35,0.8],...
            'NextPlot','replacechildren',...
            'FontUnits','normalized',...
            'FontSize',UI.figFont,...
            'Tag','figLPF');

        UI.BPFtxt = uicontrol('Style','text',...
            'String',txtBPF,...
            'Units','normalized',...
            'Position',[0.62, 0.91, 0.05, 0.04],...
            'FontUnits','normalized',...
            'FontSize',0.8,...
            'FontWeight','bold');

        UI.BPFaxe = axes(...
            'Parent',figFilter,...
            'Units','normalized',...
            'Position',[0.62,0.1,0.35,0.8],...
            'NextPlot','replacechildren',...
            'FontUnits','normalized',...
            'FontSize',UI.figFont,...
            'Tag','figBPF');

        %% dropdown menu for filter sets
        UI.filterSetTxt = uicontrol('Style','text',...
            'String','Filter Set',...
            'Units','normalized',...
            'Position',[UIPos(1,:)+[-0.013,0.12],0.1,SG.TS(2)],...
            'FontUnits','normalized',...
            'FontSize',SG.TF,...
            'FontWeight','bold');

        setIndStr = cell(size(rcvInd));
        for i = 1:size(rcvInd,1)
            setIndStr{i} = [blanks(10),int2str(i)];
        end

        UI.filterSetInd = uicontrol('Style','popupmenu',...
            'Units','normalized',...
            'Position',[UIPos(1,:)+[0.075,0.105],0.085,0.05],...
            'BackgroundColor',bkgrnd,...
            'String',setIndStr,...
            'Value',setValue,...
            'FontWeight','bold',...
            'Callback',@filterSetSelect);

        %% Pos(1): Slider for decimSampleRate
        UI.decimRatetxt = UItemp(1).txt;
        if strcmp(RcvParam(setValue).sampleMode,'NS200BWI')
            set(UI.decimRatetxt,'String','Interleave SampleRate');
        else
            set(UI.decimRatetxt,'String','decimSampleRate');
        end

        L = length(RcvParam(setValue).VDAS.decimRate);
        UI.decimRateSldr = UItemp(1).sldr;
        set(UI.decimRateSldr,'Max',-1,'Min',-L,'Value',-RcvParam(setValue).freqIdx,...
            'SliderStep',[1/(L-1), 5/(L-1)],...
            'Tag','UI.decimRateSldr',...
            'Callback',{@updateRate});

        UI.decimRateValue = UItemp(1).edit;
        set(UI.decimRateValue,'String',num2str(RcvParam(setValue).decimRate,'%2.3f'),...
            'Tag','UI.decimRateValue',...
            'Callback',{@updateRate});

        %% Pos(2) sampleMode
        UI.sampleModetxt = uicontrol('Style','text',...
            'String','sampleMode',...
            'Units','normalized',...
            'Position',[UIPos(2,:)+[-0.013,0.045],0.1,SG.TS(2)],...
            'FontUnits','normalized',...
            'FontSize',SG.TF,...
            'FontWeight','bold');

        UI.sampleMode = uicontrol('Style','popupmenu',...
            'Units','normalized',...
            'Position',[UIPos(2,:)+[0.075,0.03],0.085,0.05],...
            'BackgroundColor',bkgrnd,...
            'String',RcvParam(setValue).sampleModeStr,...
            'FontWeight','bold',...
            'Callback',@sampleModeSelect);

        % set default value based on Reeive structure
        set(UI.sampleMode,'Value',find(strcmp(RcvParam(setValue).sampleModeStr,RcvParam(setValue).sampleMode)));

        %% ADCRate, decimFactor and Trans.frequency
        Pos = UIPos(2,:); offset = -0.05;

        Pos = Pos + [0.005, offset+0.02];
        UI.ADCRatetxt = uicontrol('Style','text',...
            'String','ADCRate (Fs)',...
            'Units','normalized',...
            'Position',[Pos+[0,0.02],0.1,SG.TS(2)],...
            'FontUnits','normalized',...
            'FontSize',SG.TF,...
            'FontWeight','bold');

        UI.ADCRateValue = uicontrol('Style','edit',...
            'String',num2str(RcvParam(setValue).ADCRate,'%2.3f'),...
            'Units','normalized',...
            'Position',[Pos+SG.EO,SG.ES],...
            'BackgroundColor',bkgrnd,...
            'FontUnits','normalized',...
            'FontSize',0.6,...
            'Tag','ADCRate',...
            'Enable','off');

        Pos = Pos + [0, offset];
        UI.decimFactortxt = uicontrol('Style','text',...
            'String','decimFactor',...
            'Units','normalized',...
            'Position',[Pos+[0,0.02],0.1,SG.TS(2)],...
            'FontUnits','normalized',...
            'FontSize',SG.TF,...
            'FontWeight','bold');

        UI.decimFactorValue = uicontrol('Style','edit',...
            'String',num2str(RcvParam(setValue).decimFactor,'%2.0f'),...
            'Units','normalized',...
            'Position',[Pos+SG.EO,SG.ES],...
            'BackgroundColor',bkgrnd,...
            'FontUnits','normalized',...
            'FontSize',0.6,...
            'Tag','decimFactor',...
            'Enable','off');

        Pos = Pos + [0,offset];
        UI.TransFreqtxt = uicontrol('Style','text',...
            'String','Trans.frequency',...
            'Units','normalized',...
            'Position',[Pos+[0,0.02],0.1,SG.TS(2)],...
            'FontUnits','normalized',...
            'FontSize',SG.TF,...
            'FontWeight','bold');

        UI.TransFreqValue = uicontrol('Style','edit',...
            'String',num2str(RcvParam(setValue).TransFreq,'%2.3f'),...
            'Units','normalized',...
            'Position',[Pos+SG.EO,SG.ES],...
            'BackgroundColor',bkgrnd,...
            'FontUnits','normalized',...
            'FontSize',0.6,...
            'Tag','UI.TransFreqValue',...
            'Enable','off',...
            'Callback',{@changeTrans});
        %
        % if ~runningVSX
        %     set(UI.TransFreqValue,'Enable','on');
        % end
        %
        %     function changeTrans(hObject,~)
        %         transFreq = str2double(get(hObject,'String'));
        %         [~, freqIdx] = min(abs(VDASG3.demodFreq - transFreq));
        %         if evalin('base','exist(''Trans'',''var'')')
        %             Trans = evalin('base','Trans');
        %         end
        %         Trans.frequency = VDASG3.demodFreq(freqIdx);
        %         RcvParam(setValue).decimRate = 4*Trans.frequency;
        %
        %         if strcmp(RcvParam(setValue).sampleMode,'BS67BW')
        %             RcvParam(setValue).decimRate = (4/3)*Trans.frequency;
        %         end
        %
        %         set(UI.TransFreqValue,'String',num2str(Trans.frequency));
        %         assignin('base','Trans',Trans);
        %         RcvParam = freqCorrection(RcvParam);
        %         backToDefault
        %     end

        %% LPF panel
        UI.LPFpanel = uipanel('Title',[txtLPF,' panel'],...
            'FontUnits','normalized',...
            'FontSize',0.08,...
            'Position',[UIPos(5,:),panelWidth,panelHighth],...
            'BackgroundColor',bkgrnd + 0.02);

        % Pos(4), Slider for cutoff frequency
        UI.cutofftxt = UItemp(4).txt; uistack(UI.cutofftxt,'top');
        set(UI.cutofftxt,'String','Cutoff frequency (MHz)',...
            'BackgroundColor',bkgrnd + 0.02);

        UI.cutoffSldr = UItemp(4).sldr; uistack(UI.cutoffSldr,'top');
        set(UI.cutoffSldr,...
            'Max',62.5,'Min',0.5,'Value',LPFParam(setValue).cutoffFreqMHz,...
            'SliderStep',[0.01, 0.05],...
            'Tag','UI.cutoffSldr',...
            'Callback',{@cutoffSldr_Callback});

        UI.cutoffValue = UItemp(4).edit; uistack(UI.cutoffValue,'top');
        set(UI.cutoffValue,...
            'String',num2str(LPFParam(setValue).cutoffFreqMHz,'%2.3f'),...
            'Tag','UI.cutoffValue',...
            'Callback',{@cutoffValue_Callback});

        % Pos(5), Slider for cutoff frequency
        UI.betatxt = UItemp(5).txt; uistack(UI.betatxt,'top');
        set(UI.betatxt,'String','Kaiser Beta',...
            'BackgroundColor',bkgrnd + 0.02);

        UI.betaSldr = UItemp(5).sldr; uistack(UI.betaSldr,'top');
        set(UI.betaSldr,...
            'Max',5,'Min',0,'Value',LPFParam(setValue).kaiserBeta ,...
            'SliderStep',[0.01, 0.05],...
            'Callback',{@betaSldr_Callback});

        UI.betaValue = UItemp(5).edit; uistack(UI.betaValue,'top');
        set(UI.betaValue,...
            'String',num2str(LPFParam(setValue).kaiserBeta ,'%1.2f'),...
            'Callback',{@betaValue_Callback});

        %% BPF panel
        UI.BPFpanel = uipanel('Title',[txtBPF,' panel'],...
            'FontUnits','normalized',...
            'FontSize',0.08,...
            'Position',[UIPos(7,:),panelWidth,panelHighth],...
            'BackgroundColor',bkgrnd + 0.02);

        % Pos(6), Slider of Center frequency
        UI.CenterFreqtxt = UItemp(6).txt; uistack(UI.CenterFreqtxt,'top');
        set(UI.CenterFreqtxt,'String','Center frequency (MHz)',...
            'BackgroundColor',bkgrnd + 0.02);


        UI.CenterFreqSldr = UItemp(6).sldr; uistack(UI.CenterFreqSldr,'top');
        set(UI.CenterFreqSldr,...
            'Max',62.5,'Min',0.5,'Value',5,...
            'SliderStep',[0.1/(62.5-0.5), 0.5/(62.5-0.5)],...
            'Tag','CFsldr',...
            'Callback',{@CenterFreqSldr_Callback});

        UI.CenterFreqValue = UItemp(6).edit; uistack(UI.CenterFreqValue,'top');
        set(UI.CenterFreqValue,...
            'String',num2str(RcvParam(setValue).TransFreq,'%2.3f'),...
            'Tag','CFValue',...
            'Callback',{@CenterFreqValue_Callback});

        % Pos(7) Relative bandwidth
        bndwdth = 1.08;
        UI.bwtxt = UItemp(7).txt; uistack(UI.bwtxt,'top');
        set(UI.bwtxt,'String','Relative bandwidth',...
            'BackgroundColor',bkgrnd+0.02);

        UI.bwSldr = UItemp(7).sldr; uistack(UI.bwSldr,'top');
        set(UI.bwSldr,...
            'Max',2,'Min',0.1,'Value',bndwdth,...
            'SliderStep',[0.01/1.9, 0.05/1.9],...
            'Tag','UI.bwSldr',...
            'Callback',{@bwSldr_Callback});

        UI.bwValue = UItemp(7).edit; uistack(UI.bwValue,'top');
        set(UI.bwValue,...
            'String',num2str(bndwdth,'%1.2f'),...
            'Tag','UI.bwValue',...
            'Callback',{@bwValue_Callback});

        % interleave will change BPF to HPF
        if RcvParam(setValue).interleave
            set(UI.CenterFreqtxt,'String','Cutoff frequency (MHz)');
            set(UI.bwtxt,'String','Kaiser Beta');
            set(UI.bwSldr,'Max',5,'Min',0,'SliderStep',[0.01/5, 0.2/5]);
        end

        %% Filter Generation pushbutton
        Pos = UIPos(8,:);
        UI.FiltGeneration = uicontrol('Style','pushbutton',...
            'String','Filter Generation',...
            'Units','normalized',...
            'Position',[Pos,panelWidth,0.06],...
            'FontUnits','normalized',...
            'FontSize',0.5,...
            'FontWeight','bold',...
            'Callback',{@filterGen_Callback});

    end

    function filterSetSelect(hObject,~)
        setValue = get(hObject,'Value');
        set(UI.sampleMode,'Value',find(strcmp(RcvParam(setValue).sampleModeStr,RcvParam(setValue).sampleMode)));
        setUI();
        computeLPF;
        computeBPF;
    end

    function setUI
        set(UI.LPFtxt,'String','LPF');
        set(UI.LPFpanel,'Title','LPF panel');
        set(UI.BPFtxt,'String','BPF');
        set(UI.BPFpanel,'Title','BPF panel');
        set(UI.decimRatetxt,'String','decimSampleRate');
        set(UI.CenterFreqtxt,'String','Center frequency (MHz)');
        set(UI.bwtxt,'String','Relative bandwidth');

        switch RcvParam(setValue).sampleMode
            case 'NS200BWI'
                set(UI.BPFtxt,'String','HPF');
                set(UI.BPFpanel,'Title','HPF panel');
                set(UI.decimRatetxt,'String','Interleave SampleRate');
                set(UI.CenterFreqtxt,'String','Cutoff frequency (MHz)');
                set(UI.bwtxt,'String','Kaiser Beta');
                set(UI.bwSldr,'Max',5,'Min',0,'SliderStep',[0.01/5, 0.2/5]);
                if ~evalin('base','exist(''Trans'',''var'')')
                    set(UI.TransFreqValue,'String',num2str(RcvParam(setValue)));
                end
            case 'BS67BW'
                set(UI.LPFtxt,'String','HPF');
                set(UI.LPFpanel,'Title','HPF panel');
        end

        L = length(RcvParam(setValue).VDAS.decimRate);
        freqIndex = -round(get(UI.decimRateSldr,'Value'));
        set(UI.decimRateSldr,'Min',-L);
        set(UI.decimRateSldr,'SliderStep',[1/(L-1), 5/(L-1)]);
        if freqIndex > L
            set(UI.decimRateSldr,'Value',-L);
        end

        RcvParam = freqCorrection(RcvParam,setValue);
    end

    function updateRate(hObject,~)
        decimRateOld = RcvParam(setValue).decimRate;

        if evalin('base','exist(''Resource'',''var'')')
            Resource = evalin('base', 'Resource');
            if (Resource.Parameters.simulateMode == 1) && runningVSX
                RcvParam = freqCorrection(RcvParam,setValue);
                return;
            end
        end

        Cntrl = get(hObject,'Style');

        if strcmp(Cntrl,'slider')
            freqIndex = -round(get(hObject,'Value'));
            RcvParam(setValue).decimRate = RcvParam(setValue).VDAS.decimRate(freqIndex);
        else
            RcvParam(setValue).decimRate = str2double(get(hObject,'String'));
        end

        if ~strcmp(RcvParam(setValue).sampleMode,'BS67BW') && RcvParam(setValue).decimRate < 2 * RcvParam(setValue).TransFreq
            RcvParam(setValue).decimRate = decimRateOld;
        end

        RcvParam = freqCorrection(RcvParam,setValue);
        backToDefault;

    end

    function sampleModeSelect(hObject,~)
        sampleModeValue = get(hObject,'Value');

        RcvParam(setValue).interleave = 0;
        RcvParam(setValue).decimRate = 4*RcvParam(setValue).TransFreq;
        RcvParam(setValue).sampleMode = RcvParam(setValue).sampleModeStr{sampleModeValue};

        RcvParam(setValue).VDAS = VDASG3; % default

        switch RcvParam(setValue).sampleMode
            case 'NS200BW'
                RcvParam(setValue).BPFfilterNum = 1;
            case 'NS200BWI'
                % interleave will change BPF to HPF
                RcvParam(setValue).VDAS = VDASInterleave;
                RcvParam(setValue).BPFfilterNum = 5;
                RcvParam(setValue).interleave = 1;

                if ~evalin('base','exist(''Trans'',''var'')')
                    RcvParam(setValue).decimRate = RcvParam(setValue).TransFreq*4;
                end
            case 'BS100BW'
                RcvParam(setValue).BPFfilterNum = 2;
            case 'BS67BW'
                RcvParam(setValue).VDAS = VDASG3FourThird;
                RcvParam(setValue).decimRate = (4/3)*RcvParam(setValue).TransFreq;
                RcvParam(setValue).BPFfilterNum = 4;
            case 'BS50BW'
                RcvParam(setValue).BPFfilterNum = 3;
            case 'custom'
                RcvParam(setValue).BPFfilterNum = 1;
            otherwise
                error('filterTool: Unrecognized Receive(%d).sampleMode.\n',rcvEvent);
        end

        setUI;
        backToDefault;
    end

    function cutoffSldr_Callback(varargin)
        LPFParam(setValue).modified = 1;
        LPFParam(setValue).LPFcoef = [];
        LPFParam(setValue).cutoffFreqMHz = get(UI.cutoffSldr,'Value');
        computeLPF
        updateVDASfilter('LPF');
    end

    function cutoffValue_Callback(varargin)
        LPFParam(setValue).modified = 1;
        LPFParam(setValue).LPFcoef = [];
        LPFParam(setValue).cutoffFreqMHz = str2double(get(UI.cutoffValue,'String'));
        computeLPF;
        updateVDASfilter('LPF');
    end

    function betaSldr_Callback(varargin)
        LPFParam(setValue).modified = 1;
        LPFParam(setValue).LPFcoef = [];
        LPFParam(setValue).kaiserBeta = get(UI.betaSldr,'Value');
        set(UI.betaValue,'String',num2str(LPFParam(setValue).kaiserBeta,'%2.3f'));
        computeLPF;
        updateVDASfilter('LPF');
    end

    function betaValue_Callback(varargin)
        LPFParam(setValue).modified = 1;
        LPFParam(setValue).LPFcoef = [];
        LPFParam(setValue).kaiserBeta = str2double(get(UI.betaValue,'String'));
        set(UI.betaSldr,'Value',LPFParam(setValue).kaiserBeta);
        computeLPF;
        updateVDASfilter('LPF');
    end

    function CenterFreqSldr_Callback(varargin)
        BPFParam(setValue).modified = 1;
        BPFParam(setValue).BPFcoef = [];
        BPFParam(setValue).centerFreqMHz = get(UI.CenterFreqSldr,'Value');
        set(UI.CenterFreqValue,'String',num2str(BPFParam(setValue).centerFreqMHz,'%2.3f'));
        computeBPF;
        updateVDASfilter('BPF');
    end

    function CenterFreqValue_Callback(varargin)
        BPFParam(setValue).modified = 1;
        BPFParam(setValue).BPFcoef = [];
        BPFParam(setValue).centerFreqMHz = str2double(get(UI.CenterFreqValue,'String'));
        set(UI.CenterFreqSldr,'Value',BPFParam(setValue).centerFreqMHz);
        computeBPF;
        updateVDASfilter('BPF');
    end

    function bwSldr_Callback(varargin)
        BPFParam(setValue).modified = 1;
        BPFParam(setValue).BPFcoef = [];
        BPFParam(setValue).bndwdth = get(UI.bwSldr,'Value');
        set(UI.bwValue,'String',num2str(get(UI.bwSldr,'Value'),'%1.2f'));
        computeBPF;
        updateVDASfilter('BPF');
    end

    function bwValue_Callback(varargin)
        BPFParam(setValue).modified = 1;
        BPFParam(setValue).BPFcoef = [];
        bandwidth = str2double(get(UI.bwValue,'String'));
        BPFParam(setValue).bndwdth = bandwidth;
        set(UI.bwSldr,'Value',BPFParam(setValue).bndwdth);
        computeBPF;
        updateVDASfilter('BPF');
    end

    function filterGen_Callback(varargin)

        computeLPF;
        computeBPF;

        LPF = evalin('base','LPFParam');
        BPF = evalin('base','BPFParam');

        if ~isdeployed()
            % Give focus to the command window.
            % This function is not supported by the MATLAB Compiler.
            commandwindow;
        end

        disp('Receive.LowPassCoef =');
        disp(['[', num2str(LPF(setValue).LPFcoef(1:6),'%+9.5f'), '...']);
        disp([' ', num2str(LPF(setValue).LPFcoef(7:12),'%+9.5f'), '];']);
        disp(' ')
        disp('Receive.InputFilter =');
        disp(['[',  num2str(BPF(setValue).BPFcoef(1:7),'%+9.5f'), ' ...']);
        disp([' ',  num2str(BPF(setValue).BPFcoef(8:14),'%+9.5f'), ' ...']);
        disp([' ',  num2str(BPF(setValue).BPFcoef(15:21),'%+9.5f'), '];']);
        disp(' ')

    end

% backtoDefault is used when sampleMode or decimSampleRate is changed
    function backToDefault(varargin)

        % Start from default: no modification and no LPFcoef
        LPFParam(setValue).modified = 0;
        LPFParam(setValue).LPFcoef = [];
        BPFParam(setValue).modified = 0;
        BPFParam(setValue).BPFcoef = [];

        % LPF default
        LPFParam(setValue).DefaultCoef = DefaultLPF.coef(RcvParam(setValue).LPFfilterNum,:);
        LPFParam(setValue).DefaultFout = DefaultLPF.Fout(RcvParam(setValue).LPFfilterNum,:);
        LPFParam(setValue).cutoffFreqMHz = RcvParam(setValue).demodFreq*DefaultLPF.cutoff(RcvParam(setValue).LPFfilterNum);
        LPFParam(setValue).kaiserBeta = DefaultLPF.beta(RcvParam(setValue).LPFfilterNum);

        if RcvParam(setValue).interleave
            BPFParam(setValue).centerFreqMHz = RcvParam(setValue).demodFreq*0.5;
        else
            BPFParam(setValue).centerFreqMHz = RcvParam(setValue).demodFreq;
        end

        BPFParam(setValue).DefaultCoef = DefaultBPF.coef(RcvParam(setValue).BPFfilterNum,:);
        BPFParam(setValue).DefaultFout = DefaultBPF.Fout(RcvParam(setValue).BPFfilterNum,:);
        BPFParam(setValue).bndwdth = DefaultBPF.bndwdth(RcvParam(setValue).BPFfilterNum);
        BPFParam(setValue).xsnwdth = DefaultBPF.xsnwdth(RcvParam(setValue).BPFfilterNum);

        computeLPF;
        computeBPF;
        updateVDASfilter('ALL');

    end

% updateVDASfilter is used to update VDAS when the filter or decimSampleRate is
% modified
    function updateVDASfilter(updateCase)

        if evalin('base','exist(''Receive'',''var'')')

            % Update Receive structures.
            Receive = evalin('base', 'Receive');

            LPF = evalin('base','LPFParam');
            BPF = evalin('base','BPFParam');

            for i = 1:size(Receive,2)
                if isequal(setInd(i),setValue)
                    Receive(i).sampleMode = RcvParam(setValue).sampleMode;
                    Receive(i).decimSampleRate = RcvParam(setValue).decimRate;

                    switch updateCase
                        case 'LPF'
                            Receive(i).LowPassCoef = LPF(setInd(i)).LPFcoef;
                        case 'BPF'
                            Receive(i).InputFilter = BPF(setInd(i)).BPFcoef;
                        case 'ALL'
                            Receive(i).LowPassCoef = LPF(setInd(i)).LPFcoef;
                            Receive(i).InputFilter = BPF(setInd(i)).BPFcoef;
                    end
                end
            end
            assignin('base','Receive',Receive);

            %             if runningVSX
            %                 evalin('base','VsUpdate(''Receive'')');
            %             end

            if evalin('base','exist(''Control'',''var'')')
                Control = evalin('base','Control');
                if isempty(Control(1).Command)
                    n=1;
                else
                    n=length(Control)+1;
                end
                Control(n).Command = 'update&Run';
                Control(n).Parameters = {'Receive','Recon'};
                assignin('base','Control', Control);
            end
        end

    end

% RcvParam = freqCorrection(RcvParam) is used to obtain correct
% decimSampleRate supported by VDAS
    function RcvParam = freqCorrection(RcvParam,setNum)

        % - Determine the Receive.ADCRate and Receive.decimFactor based on RcvParam(setNum).TransFreq
        % used and return to base workspace
        [~, RcvParam(setNum).freqIdx] = min(abs(RcvParam(setNum).VDAS.decimRate - RcvParam(setNum).decimRate));
        RcvParam(setNum).ADCRate = RcvParam(setNum).VDAS.ADCRate(RcvParam(setNum).freqIdx);
        RcvParam(setNum).decimRate = RcvParam(setNum).VDAS.decimRate(RcvParam(setNum).freqIdx);
        RcvParam(setNum).decimFactor = RcvParam(setNum).VDAS.decimFactor(RcvParam(setNum).freqIdx);
        RcvParam(setNum).demodFreq = RcvParam(setNum).VDAS.demodFreq(RcvParam(setNum).freqIdx);

        if ~isempty(findobj('Tag','filterTool'))
            set(findobj('Tag','UI.decimRateSldr'),'Value',-RcvParam(setNum).freqIdx);
            set(findobj('Tag','UI.decimRateValue'),'String',num2str(RcvParam(setNum).decimRate,'%2.3f'));
            set(findobj('Tag','ADCRate'),'String',num2str(RcvParam(setNum).ADCRate,'%2.3f'));
            set(findobj('Tag','decimFactor'),'String',num2str(RcvParam(setNum).decimFactor,'%2.0f'));
        end

        if strcmp(RcvParam(setNum).sampleMode,'BS67BW') && isequal(RcvParam(setNum).decimFactor,2)
            RcvParam(setNum).LPFfilterNum = 9;
        else
            RcvParam(setNum).LPFfilterNum = RcvParam(setNum).decimFactor;
        end

    end

    function computeLPF

        % Define coefficients for CGD lowpass filter and plot the response,
        % using parameter values in LPFParam(setValue) structure.

        if evalin('base','~exist(''LPFParam'',''var'')')
            assignin('base','LPFParam',LPFParam);
            assignin('base','BPFParam',BPFParam);
        end

        nQ = 15; % number of magnitude bits for coefficient quantization
        Nfilt = RcvParam(setValue).decimFactor; % subsampling used in VDAS
        LPFxlim = RcvParam(setValue).ADCRate/2;

        % If Fs is 4/3 samples per wavelength, limit will be changed to HPF
        if RcvParam(setValue).BPFfilterNum == 4
            win = 'high';
        else % cutoff frequency must less than Nyquist Frequency
            win = 'low';
        end

        ctffFreq = LPFParam(setValue).cutoffFreqMHz/LPFxlim; % Needs to be normalized to Nyquist frequency

        % Set up parameters for calculating the frequency response:
        Nfft = 1024; % size of the fft
        Fsteps = LPFxlim*(0:2/Nfft:1);
        step = 4*Nfilt/Nfft; % size of one step in our normalized-to-Fc frequency units

        showTitle = 1;

        if ~isempty(LPFParam(setValue).LPFcoef)

            fcoef = [LPFParam(setValue).LPFcoef, fliplr(LPFParam(setValue).LPFcoef(1:end-1))];

            if isequal(LPFParam(setValue).modified,0) && max(abs(LPFParam(setValue).LPFcoef-DefaultLPF.coef(RcvParam(setValue).decimFactor,:))) > 1e-4
                showTitle = 0;
                set(UI.cutoffValue,'String','Custom');
                set(UI.betaValue,'String','Custom');
            end

        else

            if ctffFreq >= 0.99

                LPFParam(setValue).cutoffFreqMHz = LPFxlim;
                Fout = zeros(1,Nfft);
                LPFParam(setValue).LPFcoef = [zeros(1,11), 1];
                fcoef = [LPFParam(setValue).LPFcoef, fliplr(LPFParam(setValue).LPFcoef(1:end-1))];

                showTitle = 0;
                set(findobj('Tag','UI.cutoffSldr'),'Value',LPFParam(setValue).cutoffFreqMHz);
                set(UI.cutoffValue,'String','All Pass');
                set(UI.betaValue,'String','All Pass');

            else

                %% Coef calculation
                % if enabled this section calculates coefficients using Kaiser window function:
                % Compute coefficients using Kaiser window:

                kbeta = LPFParam(setValue).kaiserBeta;
                wndw  = kaiser(LPFParam(setValue).numTaps,kbeta);
                fcoef = fir1(LPFParam(setValue).numTaps-1, ctffFreq, win, wndw );

                % Make sure coefficient values are symmetric
                fcoef = (fcoef + fliplr(fcoef))/2;

                % Now quantize to coefNbits resolution
                fcoef = (round(fcoef*2^nQ))/(2^nQ); % 16 bits for Gen3 CGD design
                LPFParam(setValue).LPFcoef = fcoef(1:12);

                set(UI.cutoffSldr,'Value',LPFParam(setValue).cutoffFreqMHz);
                set(UI.cutoffValue,'String',num2str(LPFParam(setValue).cutoffFreqMHz,'%2.3f'));
                set(UI.betaValue,'String',num2str(get(UI.betaSldr,'Value'),'%2.3f'));
                set(UI.betaSldr,'Value',LPFParam(setValue).kaiserBeta);

            end

        end



        %% now find passband and stopband for selected filter
        if ~isequal(fcoef,LPFParam(setValue).DefaultCoef(1,:))
            Fout = 20*log10(abs(fft(fcoef,Nfft)));

            acceptRange = 0.01;
            ind1 = find(abs((Fout+3)./Fout)<acceptRange);
            ind2 = find(abs((Fout+20)./Fout)<acceptRange);

            status = 0;
            while status ~= 1
                if isempty(ind1) || isempty(ind2)
                    acceptRange = acceptRange + 0.01;
                    ind1 = find(abs((Fout+3)./Fout)<acceptRange);
                    ind2 = find(abs((Fout+20)./Fout)<acceptRange);
                    if acceptRange > 0.05, status = 1; showTitle = 0; end
                else
                    status = 1;
                end
            end
        end
        LPFParam(setValue).Fout = Fout;

        %% plot LPF or HPF

        % If Fs is 4/3 samples per wavelength with decimFactor = 2, Nfilt will be 9
        % for HPF

        FoutDefault = LPFParam(setValue).DefaultFout;

        axeLPF = findobj('Tag','figLPF');
        axes(axeLPF),
        plot(Fsteps,FoutDefault((1:Nfft/2+1)),Fsteps,Fout(1:Nfft/2+1),'r');
        if AutoUpdateLegend
            legend('Default','Design','AutoUpdate','off');
        else
            legend('Default','Design');
        end

        ylim([-60, 10]);
        xlim([0 LPFxlim]);
        xlabel('Frequency (MHz)','FontUnits','normalized','FontSize',UI.figFont);
        ylabel('Amplitude Response in dB','FontUnits','normalized','FontSize',UI.figFont);

        if showTitle

            dBfreq1 = Fsteps(ind1(1));
            dBfreq2 = Fsteps(ind2(1));
            fbw = dBfreq1/RcvParam(setValue).TransFreq*100;

            nbw = 0.8; % fractional bandwidth over which we will count the noise
            nsum = 1;
            ncount = 1;
            if Nfilt>1  % no subsampling and thus no aliasing if nfilt = 1
                nsum = 0;
                ncount = 0;
                for i=3:2:2*Nfilt-1 % only the odd harmonics will alias to Fc
                    for j=round((i-nbw/2)/step+1):round((i+nbw/2)/step+1)
                        nsum = nsum + 10^(Fout(j)/10); % sum noise power from each fft spectral point
                        ncount = ncount + 1; % cumulative count of points summed, for normalization
                    end
                end
            end
            aliasnoise = 10*log(nsum/ncount); % normalized total aliased noise power in dB

            title({['-3 dB bandwidth ' num2str(fbw,'%3.0f') ' %  with cutoff frequency ',num2str(dBfreq1,'%1.2f'),' MHz'];...
                ['Stopband -20 dB point ' num2str(dBfreq2,'%1.2f'),' MHz'];...
                ['aliased noise power ' num2str(aliasnoise, '%.1f') ' dB using ' num2str(nbw, '%.2f')]},'FontUnits','normalized','FontSize',UI.figFont);
        else
            title('')
        end
        line([0 LPFxlim], [3 3],'color','black','LineStyle','--');
        line([0 LPFxlim], [-3 -3],'color','black','LineStyle','--');
        line([0 LPFxlim], [-20 -20],'color','black','LineStyle','--');

        % freqFilt is the frequency limit after decimation, so customer will be
        % able to know the actual frequency range after LPF and subsampling

        if 1 < Nfilt % && (Nfilt<9)
            line([RcvParam(setValue).decimRate/2 RcvParam(setValue).decimRate/2], [-60 10],'color','blue','LineStyle','--');
            text(RcvParam(setValue).decimRate/2,-15,'\leftarrow','FontSize',10);
            text(RcvParam(setValue).decimRate/2,-15.4,'      Frequency limit after subsampling','FontUnits','normalized','FontSize',UI.figFont);
            if showTitle, line([LPFParam(setValue).cutoffFreqMHz LPFParam(setValue).cutoffFreqMHz], [-6 10],'color','red','LineStyle','--');end

        end

        assignin('base','LPFParam',LPFParam);
        if isequal(LPFParam(setValue).modified,1)
            LPFParam(setValue).LPFcoef = [];
        end

    end

    function computeBPF

        % Function to compute the CGD bandpass filter coefficients and plot the
        % filter's frequency response.  Coefficient values are saved in the
        % BPFParam(setValue) structure, for importing into VSX after exiting the coef
        % development script.

        % Because the sample rate at the point right before bandpass filter
        % (Receive.InputFilterCoefs) is 4*RcvParam(setValue).TransFreq, the Nyquist for the
        % bandpass filter will be always 2*RcvParam(setValue).TransFreq. Therefore, if fmax is
        % awlays 2, Fc will be center freq/Transquency

        % In the case of NS200BWI smapleMode, BPF will become a HPF.

        if evalin('base','~exist(''BPFParam'',''var'')')
            assignin('base','BPFParam',BPFParam);
            assignin('base','LPFParam',LPFParam);
        end

        showTitle = 1;
        set(findobj('tag','UI.bwValue'),'String',num2str(BPFParam(setValue).bndwdth,'%1.2f'));
        set(findobj('tag','UI.bwSldr'),'Value',BPFParam(setValue).bndwdth);
        set(findobj('tag','CFValue'),'String',num2str(BPFParam(setValue).centerFreqMHz,'%2.3f'));
        set(findobj('tag','CFsldr'),'Value',BPFParam(setValue).centerFreqMHz);

        if RcvParam(setValue).interleave % HPF calculation

            nQ = 15; % number of magnitude bits for coefficient quantization
            Nfft = 1024;

            if isempty(BPFParam(setValue).BPFcoef)

                ctffFreq = BPFParam(setValue).centerFreqMHz/(RcvParam(setValue).decimRate/4);

                if ctffFreq > 0.99
                    ctffFreq = 0.99;
                    BPFParam(setValue).centerFreqMHz = ctffFreq * RcvParam(setValue).decimRate/4;
                end


                kbeta = BPFParam(setValue).bndwdth;
                wndw  = kaiser(BPFParam(setValue).numTaps,kbeta);
                fcoef = fir1(BPFParam(setValue).numTaps-1, ctffFreq, 'high', wndw);

                % Make sure coefficient values are symmetric
                fcoef = (fcoef + fliplr(fcoef))/2;

            else
                fcoef = [BPFParam(setValue).BPFcoef, fliplr(BPFParam(setValue).BPFcoef(1:end-1))];
                if isequal(BPFParam(setValue).modified,0) && max(abs(BPFParam(setValue).BPFcoef-BPFParam(setValue).DefaultCoef)) > 1e-4
                    showTitle = 0;
                    set(UI.CenterFreqValue,'String','Custom');
                    set(UI.bwValue,'String','Custom');
                end
            end

            % Now quantize to coefNbits resolution
            fcoef = (round(fcoef*2^nQ))/(2^nQ); % 16 bits for Gen3 CGD design
            BPFParam(setValue).BPFcoef = fcoef(1:(BPFParam(setValue).numTaps+1)/2); %
            Fout = 20*log10(abs(fft(fcoef,Nfft)));
            BPFParam(setValue).Fout = Fout;
            BPFxlim = RcvParam(setValue).decimRate/2;
            Fsteps = BPFxlim*(1/Nfft:1/Nfft:1);

            if ~isequal(fcoef,LPFParam(setValue).DefaultCoef(1,:))
                Fout = 20*log10(abs(fft(fcoef,Nfft)));

                acceptRange = 0.01;
                ind1 = find(abs((Fout+3)./Fout)<acceptRange);
                ind2 = find(abs((Fout+20)./Fout)<acceptRange);

                status = 0;
                while status ~= 1
                    if isempty(ind1) || isempty(ind2)
                        acceptRange = acceptRange + 0.01;
                        ind1 = find(abs((Fout+3)./Fout)<acceptRange);
                        ind2 = find(abs((Fout+20)./Fout)<acceptRange);
                        if acceptRange > 0.05, status = 1; showTitle = 0; end
                    else
                        status = 1;
                    end
                end
            end

            arrow1 = findall(0,'Tag','arrow1'); if ishandle(arrow1),delete(arrow1); end
            arrow2 = findall(0,'Tag','arrow2'); if ishandle(arrow2),delete(arrow2); end

            axeBPF = findobj('Tag','figBPF');
            axes(axeBPF); % select the display figure for updating

            % 0.62 is the bottom corner of the BPF
            plot(Fsteps,BPFParam(setValue).DefaultFout,Fsteps,Fout,'r');
            xlim([0 BPFxlim]),ylim(BPFParam(setValue).Rangsel); title('');
            if AutoUpdateLegend
                legend('Default','Design','AutoUpdate','off');
            else
                legend('Default','Design');
            end
            x1 = 0.62 + 0.35 * (RcvParam(setValue).TransFreq/BPFxlim);
            y1 = 0.15;
            y2 = 0.1;
            annotation('Textarrow',[x1,x1],[y1,y2],'Color','b','LineWidth',1,...
                'String','Trans.frequency','TextColor','k','FontUnits','normalized','FontSize',UI.figFont,...
                'TextMargin',65,'Tag','arrow1');
            line([BPFxlim/2 BPFxlim/2], [-70 5],'color','b','LineStyle','--');

            if showTitle
                dBfreq1 = Fsteps(ind1(1));
                dBfreq2 = Fsteps(ind2(1));
                fbw = dBfreq1/RcvParam(setValue).TransFreq*100;
                nbw = 0.8; % fractional bandwidth over which we will count the noise

                nsum = 0;
                ncount = 0;
                for i=3:2:1 % only the odd harmonics will alias to Fc
                    for j=round((i-nbw/2)/step+1):round((i+nbw/2)/step+1)
                        nsum = nsum + 10^(Fout(j)/10); % sum noise power from each fft spectral point
                        ncount = ncount + 1; % cumulative count of points summed, for normalization
                    end
                end
                aliasnoise = 10*log(nsum/ncount); % normalized total aliased noise power in dB
                title({['-3 dB bandwidth ' num2str(fbw,'%3.0f') ' %  with cutoff frequency ',num2str(dBfreq1,'%1.2f'),' MHz'];...
                    ['Stopband -20 dB point ' num2str(dBfreq2,'%1.2f'),' MHz'];...
                    ['aliased noise power ' num2str(aliasnoise, '%.1f') ' dB using ' num2str(nbw, '%.2f')]},'FontUnits','normalized','FontSize',UI.figFont);
            else
                title('')
            end

            line([0 BPFxlim], [3 3],'color','black','LineStyle','--');
            line([0 BPFxlim], [-3 -3],'color','black','LineStyle','--');
            line([0 BPFxlim], [-20 -20],'color','black','LineStyle','--');


        else % Regular BPF calculation and plot

            %% step 1: Initialize parameters

            nQ = 15; % number of magnitude bits for coefficient quantization
            Nfilt = RcvParam(setValue).BPFfilterNum; % find which filter we're processing

            % Set up parameters for calculating the frequency response:
            Nfft = 1024; % size of the fft
            Fout= -100*ones(1,Nfft); % array for fft results
            FoutDefault = -100*ones(1,Nfft);

            bndwdth = BPFParam(setValue).bndwdth;  % The value has been modified to be within [0.1 2] range
            BWrange = [2;1;0.5;0.6];
            if bndwdth > BWrange(Nfilt), bndwdth = BWrange(Nfilt);
            elseif bndwdth < 0.1, bndwdth = 0.1; end

            Freqrange = [0.1 1.9; 0.5 1.5; 0.75 1.25; 2.1 3.9]; % allowed range for center frequency

            if Nfilt == 4
                % in this case, 2<cetFreq<4, and 0.1<bandwidth<0.6
                BPFxlim = RcvParam(setValue).decimRate;
                ctrFreq = BPFParam(setValue).centerFreqMHz/(BPFxlim/4);
                Fsteps = BPFxlim*(1/Nfft:1/Nfft:1); % scale the frequency steps so Fc is always 1 (and Nyquist limit is at BPFxlim)
            else
                BPFxlim = RcvParam(setValue).decimRate/2;
                ctrFreq = BPFParam(setValue).centerFreqMHz/(BPFxlim/2);
                Fsteps = BPFxlim*(0:2/Nfft:1); % scale the frequency steps so Fc is always 1 (and Nyquist limit is at BPFxlim)
            end

            % Check the Freq is within correct range
            if ctrFreq*(1 + bndwdth/2) > Freqrange(Nfilt,2)
                ctrFreq = Freqrange(Nfilt,2)/(1 + bndwdth/2);
                if ctrFreq*(1 - bndwdth/2) < Freqrange(Nfilt,1)
                    ctrFreq = sum(Freqrange(Nfilt,:))/2;
                end

            elseif ctrFreq*(1 - bndwdth/2) < Freqrange(Nfilt,1)
                ctrFreq = Freqrange(Nfilt,1)/(1 - bndwdth/2);
                if ctrFreq*(1 + bndwdth/2) > Freqrange(Nfilt,2)
                    ctrFreq = sum(Freqrange(Nfilt,:))/2;
                end
            end

            % Have correct value shown in the GUI
            if Nfilt == 4
                BPFParam(setValue).centerFreqMHz = BPFxlim*ctrFreq/4;
                BPFParam(setValue).bndwdth  = bndwdth;
                bndwdth = (ctrFreq/(4-ctrFreq))*bndwdth;

                ctrFreq = 4-ctrFreq; % Aliasing to 0-2 range for filter coef determination
                Freqrange(Nfilt,:) = Freqrange(1,:);
                BPFParam(setValue).centerFreq  = ctrFreq;
            else
                BPFParam(setValue).centerFreq  = ctrFreq;
                BPFParam(setValue).bndwdth  = bndwdth;
                BPFParam(setValue).centerFreqMHz = ctrFreq*BPFxlim/2;
            end

            %% step 2: Create frequency point array for firpm function

            hxsn=0.5*BPFParam(setValue).xsnwdth;
            apm = [0 0 1 1 0 0 ]; % amplitude points for the firpm function
            fpm = [0 .1 .4 .5 .9 1]; % frequency points for the firpm function; actual values will be set below

            FL = (1 - bndwdth/2)*ctrFreq; % normalized band edges
            FH = (1 + bndwdth/2)*ctrFreq;
            if FH > Freqrange(Nfilt,2)
                FH = Freqrange(Nfilt,2);
            end
            if FL < Freqrange(Nfilt,1)
                FL = Freqrange(Nfilt,1);
            end

            FL=FL/2; % scale to units used by fft for 4X sampling
            FH=FH/2;
            hxsn=hxsn/2;

            % set nominal values  for fpm vector [0 FL-hxsn FL+hxsn FH-hxsn FH+hxsn 1] but check
            % for overlaps and restrict transition width if necessary:
            fpm(2) = max((FL-hxsn),.001); % don't let low edge of transition band go to zero or below
            fpm(5) = min((FH+hxsn),.999); % don't let high edge of transition band go to one or above
            if (FL+hxsn)<(FH-hxsn) % make sure passband edges are in ascending order
                fpm(3) = max(fpm(2)+.001,FL+hxsn); % keep transition band edges in correct order too
                fpm(4) = min(fpm(5)-.001,FH-hxsn);
            else
                fpm(3) = (FL+FH)/2-.001; % stay in the middle if hxsn is too big
                fpm(4) = (FL+FH)/2+.001;
            end


            %% step 3: get firpm coeff's
            if isempty(BPFParam(setValue).BPFcoef)

                fcoef = firpm((BPFParam(setValue).numTaps-1),fpm,apm); % compute FIR coef's using firpm function

                % Make sure coefficient values are symmetric
                fcoef = (fcoef + fliplr(fcoef))/2;

                % Now null out the DC term
                DC=sum(fcoef);
                fcoef(BPFParam(setValue).DCnulltap ) = fcoef(BPFParam(setValue).DCnulltap ) - DC/2; % adjust coef to null DC
                fcoef = fliplr(fcoef); % flip to do it again on the other end
                fcoef(BPFParam(setValue).DCnulltap ) = fcoef(BPFParam(setValue).DCnulltap ) - DC/2;

                % Now normalize gain at Fc, quantize, etc.
                testresp = abs(fft(fcoef,Nfft));

                FLindx = round(fpm(3)*Nfft/2 + 1); % bandwidth limits for gain normalization
                FHindx = round(fpm(4)*Nfft/2 + 1);

                % compute average gain over passband, but don't let it go below .01
                Fcgain =max(sum(testresp(FLindx:FHindx))/(FHindx-FLindx+1),.01);
                fcoef = fcoef/Fcgain; % note scaling preserves zero at DC

                % Now quantize to coefNbits resolution
                fcoef = round(fcoef*2^nQ);
                DC=sum(fcoef); % remove any DC offset from the rounding
                ctrtap = round(BPFParam(setValue).numTaps(1)/2);
                fcoef(ctrtap) = fcoef(ctrtap) - DC; % offset center coef, to preserve integer steps
                fcoef = fcoef/(2^nQ);

                BPFParam(setValue).BPFcoef = fcoef(1:(BPFParam(setValue).numTaps+1)/2); %
            else
                fcoef = [BPFParam(setValue).BPFcoef, fliplr(BPFParam(setValue).BPFcoef(1:end-1))];
                if isequal(BPFParam(setValue).modified,0) && max(abs(BPFParam(setValue).BPFcoef-BPFParam(setValue).DefaultCoef)) > 1e-4
                    %                     ~isequal(BPFParam(setValue).BPFcoef,BPFParam(setValue).DefaultCoef)
                    set(UI.CenterFreqValue,'String','Custom');
                    set(UI.bwValue,'String','Custom');
                    showTitle = 0;
                end
            end

            %% step 4: compute frequency response and find passband and stopband
            Fout(:) = 20*log10(abs(fft(fcoef,Nfft)));

            if Nfilt == 4
                FoutLeft = Fout;
                Fout(1:Nfft/2+1) = -100;
                FoutLeft(Nfft/2+1:end) = -100;
                FoutDefault(1:Nfft/2+1) = BPFParam(setValue).DefaultFout(1:Nfft/2+1);
            else
                Fout = Fout(1:length(Fsteps));
                FoutDefault = BPFParam(setValue).DefaultFout;
            end

            ind1 = find(Fout+3>0);
            ind2 = find(Fout+20>0);
            ind3 = find(Fout+6>0);

            pbL = Fsteps(ind1(1)); % -3 dB point at low end in relative frequency units
            pbH = Fsteps(ind1(end));
            fbw = 100*(pbH - pbL)/BPFParam(setValue).centerFreqMHz;

            stpL = Fsteps(ind2(1)); % stpbd is the frequency at which response first goes above -20 dB; note i=1 means DC or F=0
            stpH = Fsteps(ind2(end)); % stpbd is the frequency at which response first goes above -20 dB; note i=1 means DC or F=0


            %% Plot BPF

            arrow1 = findall(0,'Tag','arrow1'); if ishandle(arrow1),delete(arrow1); end
            arrow2 = findall(0,'Tag','arrow2'); if ishandle(arrow2),delete(arrow2); end

            axeBPF = findobj('Tag','figBPF');
            axes(axeBPF); % select the display figure for updating

            % 0.62 is the bottom corner of the BPF
            if Nfilt == 4
                plot(Fsteps,FoutDefault,Fsteps,Fout,'r',Fsteps,FoutLeft,'r--');
                x1 = 0.62 + 0.35 * (BPFParam(setValue).centerFreqMHz/BPFxlim);
                x2 = 0.62 + 0.35 * (1-BPFParam(setValue).centerFreqMHz/BPFxlim);
                y1 = 0.15;
                y2 = 0.1;
                annotation('Textarrow',[x1,x1],[y1,y2],'Color','b','LineWidth',1,...
                    'String','Center Frequency','TextColor','k','FontUnits','normalized','FontSize',UI.figFont,...
                    'TextMargin',65,'Tag','arrow1');
                annotation('Textarrow',[x2,x2],[y1,y2],'Color','b','LineWidth',1,...
                    'String','Aliasing Center Frequency','TextColor','k','FontUnits','normalized','FontSize',UI.figFont,...
                    'TextMargin',85,'Tag','arrow2');
            else
                plot(Fsteps,FoutDefault(1:length(Fsteps)),Fsteps,Fout,'r');
                x1 = 0.62 + 0.35 * (RcvParam(setValue).TransFreq/BPFxlim);
                y1 = 0.15;
                y2 = 0.1;
                annotation('Textarrow',[x1,x1],[y1,y2],'Color','b','LineWidth',1,...
                    'String','Trans.frequency','TextColor','k','FontUnits','normalized','FontSize',UI.figFont,...
                    'TextMargin',65,'Tag','arrow1');
            end

            if AutoUpdateLegend
                legend('Default','Design','AutoUpdate','off');
            else
                legend('Default','Design');
            end

            if Nfilt == 2 || Nfilt == 3
                titleStr = 'Green lines show Nyquist limit.';
            else
                titleStr = '';
            end

            ylim(BPFParam(setValue).Rangsel);
            xlim([0 BPFxlim]);
            xlabel('Frequency (MHz)','FontUnits','Normalized','FontSize',UI.figFont);
            ylabel('Amplitude Response in dB','FontUnits','Normalized','FontSize',UI.figFont);
            if showTitle
                title({['-3 dB bandwidth ' num2str(fbw,'%3.0f') ' % (from ' num2str(pbL,'%1.2f') '  to '  num2str(pbH,'%1.2f'), ' MHz'];...
                    ['Stopband -20 dB points ' num2str(stpL,'%1.2f') '  to  ' num2str(stpH,'%1.2f'),' MHz' ];...
                    titleStr},'FontUnits','Normalized','FontSize',UI.figFont);
            end
            line([0 BPFxlim], [3 3],'color','black','LineStyle','--');
            line([0 BPFxlim], [-3 -3],'color','black','LineStyle','--');
            line([0 BPFxlim], [-20 -20],'color','black','LineStyle','--');

            freq = RcvParam(setValue).decimRate/4;

            if Nfilt == 2
                line([0.5*freq 0.5*freq], [-60 -3],'color','g','LineStyle','--','LineWidth',1.5);
                line([1.5*freq 1.5*freq], [-60 -3],'color','g','LineStyle','--','LineWidth',1.5);
            elseif Nfilt == 3
                line([0.75*freq 0.75*freq], [-60 -3],'color','g','LineStyle','--','LineWidth',1.5);
                line([1.25*freq 1.25*freq], [-60 -3],'color','g','LineStyle','--','LineWidth',1.5);
            end

            if ~isequal(length(ind3),length(Fsteps))
                line([Fsteps(ind3(1)) Fsteps(ind3(1))], [-6 10],'color','red','LineStyle','--');
                line([Fsteps(ind3(end)+1) Fsteps(ind3(end)+1)], [-6 10],'color','red','LineStyle','--');
            end

        end

        assignin('base','BPFParam', BPFParam);
        if isequal(BPFParam(setValue).modified,1)
            BPFParam(setValue).BPFcoef = [];
        end
    end

end
