function showTXPD(varargin)
% showTXPD visualize transmit TXPD data with cutoffs
%   Call with no arguments to process all TX structures with first PData struct.
%   If a range of TX values is desired, set the first varargin to the range, e.g. [1:10].
%   If a PData struct other than one is to be specified, set the range of TX values and
%      then specify the PData number in the 2nd argument.
%
% Last Update: 1-10-2017

if evalin('base','exist(''TX'',''var'')')
    TX = evalin('base','TX');
else
    error('ShowTXPD: No TX structure found in base workspace.');
end

if evalin('base','exist(''PData'',''var'')')
    PDataBase = evalin('base','PData');
    PDsldrMax = length(PDataBase);
    PDsldrMin = 0.99;
    if length(PDataBase) > 1,
        PDsldrSS  = 1/max(1,length(PDataBase)-1);
    else PDsldrSS = 1;
    end
else
    error('ShowTXPD: No PData structure found in base workspace.');
end

if evalin('base','exist(''Trans'',''var'')')
    Trans = evalin('base','Trans');
else
    error('ShowTXPD: No Trans structure found in base workspace.');
end

% default value
switch nargin
    case 0
        TXIn = 1:size(TX,2);
        PDIn = 1;
    case 1
        TXIn = varargin{1};
        PDIn = 1;
    case 2
        TXIn = varargin{1};
        PDIn = varargin{2};
        if length(PDIn) > 1, error('Multiple PData numbers is not allowed in the 2nd argument'); end
        PDsldrMax = PDIn;
        PDsldrMin = PDIn-0.01;
    otherwise
        error('Too many input arguments.');
end

pdelta=[]; pdeltaX=[]; pdeltaY=[]; pdeltaZ=[]; pdeltaR=[]; pdeltaT=[];
%% Get the default value, PDIn is pdata num, txn is tx num
txn = TXIn(1);
PData = PDataBase(PDIn);

% Check for presence of the TXPD data, and compute if not already computed.
if ~isfield(TX,'TXPD')
    h = waitbar(0,'Calculate TPXD, please wait!');
    for i=1:size(TXIn,2)
        TX(TXIn(i)).TXPD = computeTXPD(TX(TXIn(i)),PData);
        waitbar(i/size(TXIn,2))
    end
    [TX.peakCutOff] = deal(1.0);
    [TX.peakBLMax] = deal(4.0);
    close(h)
else
    for i=1:size(TXIn,2)
        if isempty(TX(TXIn(i)).TXPD)
            TX(TXIn(i)).TXPD = computeTXPD(TX(TXIn(i)),PData);
            TX(TXIn(i)).peakCutOff = 1.0;
            TX(TXIn(i)).peakBLMax = 4.0;
        end
    end
end

% update the pdelta setting for showTXPD plot
Delta = updatePData(PData);
deltaFieldsName = fieldnames(Delta);
for i = 1:length(deltaFieldsName)
    eval(['pdelta',deltaFieldsName{i},'= Delta.',deltaFieldsName{i},';']);
end
Nrow = PData.Size(1); Ncol = PData.Size(2); Nsec = PData.Size(3);

    function Delta = updatePData(PData)

        if isfield(PData,'PDelta') && ~isempty(PData.PDelta)
            if isfield(PData,'Coord')
                switch PData.Coord
                    case 'rectangular'
                        if size(PData.PDelta,2) ~= 3
                            error('computeRegions: PData(x).PDelta number of columns must be 3.');
                        end
                        Delta.X = PData.PDelta(1);
                        Delta.Y = PData.PDelta(2);
                        Delta.Z = PData.PDelta(3);
                    case 'polar'
                        if size(PData.PDelta,2) ~= 3
                            error('computeRegions: PData(x).PDelta number of columns must be 3.');
                        end
                        Delta.T = PData.PDelta(1);
                        Delta.R = PData.PDelta(2);
                        Delta.Z = PData.PDelta(3);
                    case 'cylindrical'
                        error('computeRegions: PData(x).Coor = ''cylindrical'' is not yet supported.');
                    otherwise
                        error('computeRegions: Unrecognized PData(x).Coor string.');
                end
            else
                % Default to rectangular coordinates if no Coord attribute
                if size(PData.PDelta,2) ~= 3
                    error('computeRegions: PData(x).PDelta number of columns must be 3.');
                end
                Delta.X = PData.PDelta(1);
                Delta.Y = PData.PDelta(2);
                Delta.Z = PData.PDelta(3);
            end
        elseif isfield(PData,'pdelta') && ~isempty(PData.pdelta)
            % if no PDelta and PData.pdelta is set, it implies all non-specified pdelta? are equal to its value.
            Delta.pdelta = PData.pdelta;
            Delta.X = pdelta;
            Delta.Y = pdelta;
            Delta.Z = pdelta;
            Delta.R = pdelta;
            Delta.T = pdelta/64;  % pdeltaT is in radians, for example if pdelta=0.5, pdeltaT = .0078 or
            % about 200 radial lines in 90 degrees
        else
            % if PDelta and pdelta are missing, at the least, pdeltaX and pdeltaZ must be specified.
            if ~isfield(PData,'pdeltaX') || isempty(PData.pdeltaX)
                error('computeRegions: PData(x).pdeltaX must be specified if PData(x).PDelta or pdelta missing.');
            elseif ~isfield(PData,'pdeltaZ') || isempty(PData.pdeltaZ)
                error('computeRegions: PData(x).pdeltaZ must be specified if PData(x).PDelta or pdelta missing.');
            end
            % Read specified pdeltaX,Y,Z,R or T.
            Delta.X = PData.pdeltaX;
            if isfield(PData,'pdeltaY')&&(~isempty(PData.pdeltaY)), Delta.Y = PData.pdeltaY; end
            Delta.Z = PData.pdeltaZ;
            if isfield(PData,'pdeltaR')&&(~isempty(PData.pdeltaR)), Delta.R = PData.pdeltaR; end
            if isfield(PData,'pdeltaT')&&(~isempty(PData.pdeltaT)), Delta.T = PData.pdeltaT; end
        end
    end

%%
% Create the figure window
ScrnSize = get(0,'ScreenSize');
if isfield(PData,'Coord')&&(strcmp(PData.Coord,'polar'))
    aspect = 150*PData.Size(2)*pdeltaT/(PData.Size(1)*pdeltaR);
else
    aspect = PData.Size(2)*pdeltaX/(PData.Size(1)*pdeltaZ);
end
imgHeight = 500; %ScrnSize(4)/3;  % set image height (in pixels) to 1/3 of screen height.
imgWidth = imgHeight * aspect;
figHeight = imgHeight + 275;
figWidth = imgWidth + 100;
bkgrnd = [0.8 0.8 0.8];
hf = figure('Visible','off',...
    'Position',[(ScrnSize(3)-figWidth)/2,(ScrnSize(4)-figHeight)/2,figWidth,figHeight],...
    'Name','showTXPD',...
    'Color',bkgrnd,...
    'NumberTitle','off',...
    'Toolbar','figure', ...
    'Resize','on',...
    'tag','TXPD');
set(hf,'CloseRequestFcn',{@closefunc});
set(hf,'DefaultUicontrolBackgroundColor',bkgrnd)

% Determine limits for axes. Default is for 2D display, will be changed by
% 3D display
if isfield(PData,'Coord')&&(strcmp(PData.Coord,'polar'))
    XLim = [pdeltaT*(1-PData.Size(2))/2, pdeltaT*(PData.Size(2)-1)/2];
    YLim = [0, pdeltaR * PData.Size(1)];
else
    XLim = [PData.Origin(1), PData.Origin(1) + pdeltaX * PData.Size(2)];
    YLim = [PData.Origin(3), PData.Origin(3) + pdeltaZ * PData.Size(1)];
end

% Draw the axes for the image plot
ha = axes('Parent',hf,...
    'Units','normalized',...
    'OuterPosition',[0,0.24,1,0.76],...
    'YDir','reverse',...
    'XLim',XLim,...
    'YLim',YLim);

% Define offsets for creating UI controls with normalized units.
% Slider group offsets and sizes.
SG = struct('TO',[0.0,0.07],...    % title offset
    'TS',[0.25,0.03],...    % title size
    'TF',0.5,...            % title font size
    'SO',[0.0,0.05],...     % slider offset
    'SS',[0.25,0.025],...   % slider size
    'EO',[0.075,0.025],...  % edit box offset
    'ES',[0.11,0.025]);     % edit box size
% buttonGroup offsets and sizes.
BG = struct('TFS',0.1,...           % button group title font size
    'BGO',[0.0,0.0],...     % button group offset
    'BGS',[0.28,0.14],...   % button group size
    'BI',0.28,...           % button increment
    'BO',[0.1,0.05],...     % button offset (units relative to BGS box)
    'BS',[0.9,0.36],...     % button size      "    "    "    "     "
    'BFS',0.27);            % button font size
% Push button offsets and sizes.
PB = struct('BS',[0.14,0.045],...   % button size
    'FS',0.30);

%% UIbutton group for plot selection.
bgx = (1 - 3*0.28)/4;
bgy = 0.10;
hbg = uibuttongroup('Units','normalized',...
    'Position',[[bgx bgy]+BG.BGO BG.BGS],...
    'FontUnits','normalized',...
    'FontSize',BG.TFS,...
    'FontWeight','bold',...
    'Title','TXPD Param.',...
    'BackgroundColor',[0.8,0.8,0.8],...
    'Interruptible','off',...
    'SelectionChangeFcn',{@parameter_Callback});
bIntens = uicontrol(hbg,'Style','radiobutton',...
    'Units','normalized',...
    'Position', [BG.BO(1) BG.BO(2)+2*BG.BI BG.BS],...
    'FontUnits','normalized',...
    'FontSize',BG.BFS,...
    'String','Intensity',...
    'Tag','rbIntens');
bTime = uicontrol(hbg,'Style','radiobutton','String','Peak Time',...
    'Units','normalized',...
    'Position', [BG.BO(1) BG.BO(2)+BG.BI BG.BS],...
    'FontUnits','normalized',...
    'FontSize',BG.BFS,...
    'Tag','rbTime');
bBL = uicontrol(hbg,'Style','radiobutton',...
    'Units','normalized',...
    'Position', [BG.BO BG.BS],...
    'FontUnits','normalized',...
    'FontSize',BG.BFS,...
    'String','Burst Length',...
    'Tag','rbBL');

% Read TX control - rereads the TX structure from the base workspace.
readTXx = (1 - 3*0.28)/4;
readTXy = 0.044;
hreadTX = uicontrol('Style','pushbutton',...
    'String','Read TX',...
    'Units','normalized',...
    'Position',[readTXx,readTXy,PB.BS],...
    'FontUnits','normalized',...
    'FontSize',PB.FS,...
    'Interruptible','off',...
    'BusyAction','cancel',...
    'tag','readTX',...
    'Callback',{@readTX_Callback});

% Write TX control - puts the modified TX structure into the base workspace.
writeTXx = (1 - 3*0.28)/4 + 0.14;
writeTXy = 0.044;
hwriteTX = uicontrol('Style','pushbutton',...
    'String','Write TX',...
    'Units','normalized',...
    'Position',[writeTXx,writeTXy,PB.BS],...
    'FontUnits','normalized',...
    'FontSize',PB.FS,...
    'Interruptible','off',...
    'BusyAction','cancel',...
    'tag','writeTX',...
    'Callback',{@writeTX_Callback});

% TX index number UI control.
txNumx = 2*(1-3*0.28)/4 + 0.28;  % 2*(space between 3 controls of width 0.28) + 1 control
txNumy = 0.16;
Pos = [txNumx,txNumy];
stps = max(1,size(TXIn,2)-1);
txNumtxt = uicontrol(hf,'Style','text',...
    'String','TX Index',...
    'Units','normalized',...
    'Position',[Pos+SG.TO,SG.TS],...
    'FontUnits','normalized',...
    'FontSize',SG.TF,...
    'FontWeight','bold');
txNumSldr = uicontrol(hf,'Style','slider',...
    'Max',TXIn(end),'Min',TXIn(1)-.01,'Value',txn,...
    'SliderStep',[1/stps,1/stps],...
    'Units','normalized',...
    'Position',[Pos+SG.SO,SG.SS],...
    'BackgroundColor',bkgrnd-0.05,...
    'Callback',{@txNum_Callback});
txNumValue = uicontrol('Style','edit','String',num2str(txn,'%3.0f'),...
    'Units','normalized',...
    'Position',[Pos+SG.EO,SG.ES],...
    'BackgroundColor',bkgrnd+0.1,...
    'Callback',{@txNum_Callback});

% TX peak cutoff slider
txCutoffx = 2*(1-3*0.28)/4 + 0.28;  % 2*(space between 3 controls of width 0.28) + 1 control
txCutoffy = 0.08;
Pos = [txCutoffx,txCutoffy];
peakMax = 40;
if ~isfield(TX,'peakCutOff')||isempty(TX(txn).peakCutOff), TX(txn).peakCutOff = 0; end
txCutofftxt = uicontrol(hf,'Style','text',...
    'String','TX Intens. Cutoff',...
    'Units','normalized',...
    'Position',[Pos+SG.TO,SG.TS],...
    'FontUnits','normalized',...
    'FontSize',SG.TF,...
    'FontWeight','bold');
txCutoffSldr = uicontrol(hf,'Style','slider',...
    'Max',peakMax,'Min',0,'Value',TX(txn).peakCutOff,...
    'SliderStep',[1/400 1/80],...
    'Units','normalized',...
    'Position',[Pos+SG.SO,SG.SS],...
    'BackgroundColor',bkgrnd-0.05,...
    'Interruptible','off',...
    'BusyAction','cancel',...
    'Callback',{@txCutoffSldr_Callback});
txCutoffValue = uicontrol(hf,'Style','edit',...
    'String',num2str(TX(txn).peakCutOff,'%2.2f'),...
    'Units','normalized',...
    'Position',[Pos+SG.EO,SG.ES],...
    'BackgroundColor',bkgrnd+0.1,...
    'Callback',{@txCutoffValue_Callback});

% TX burst length max. slider
txBLMaxx = 2*(1-3*0.28)/4 + 0.28;  % 2*(space between 3 controls of width 0.28) + 1 control
txBLMaxy = 0.0;
Pos = [txBLMaxx,txBLMaxy];
BLMax = 40.5;
BLMin = 0.5;
if ~isfield(TX,'peakBLMax')||isempty(TX(txn).peakBLMax), TX(txn).peakBLMax = BLMax; end
txBLMaxtxt = uicontrol(hf,'Style','text',...
    'String','Max. Burst Length',...
    'Units','normalized',...
    'Position',[Pos+SG.TO,SG.TS],...
    'FontUnits','normalized',...
    'FontSize',SG.TF,...
    'FontWeight','bold');
txBLMaxSldr = uicontrol(hf,'Style','slider',...
    'Max',BLMax,'Min',0.5,'Value',TX(txn).peakBLMax,...
    'SliderStep',[0.5/(BLMax-BLMin) 2/(BLMax-BLMin)],...
    'Units','normalized',...
    'Position',[Pos+SG.SO,SG.SS],...
    'BackgroundColor',bkgrnd-0.05,...
    'Interruptible','off',...
    'BusyAction','cancel',...
    'Callback',{@txBLMaxSldr_Callback});
txBLMaxValue = uicontrol(hf,'Style','edit',...
    'String',num2str(TX(txn).peakBLMax,'%2.2f'),...
    'Units','normalized',...
    'Position',[Pos+SG.EO,SG.ES],...
    'BackgroundColor',bkgrnd+0.1,...
    'Callback',{@txBLMaxValue_Callback});

% Section number slider appears only with 3D data
if PData.Size(3) > 1
    nsecx = 3*(1-3*0.28)/4 + 0.56;  % 3*(space between 3 controls of width 0.28) + 2 control
    nsecy = 0.16;
    Pos = [nsecx,nsecy];

    sectNumX = round(PData.Size(2)/2);
    sectNumY = round(PData.Size(1)/2);
    sectNumZ = round(PData.Size(3)/2);

    % default orientation is xy
    shapeOrientation = 'xy';
    numSections = PData.Size(3);
    sectNum = sectNumZ;
    sliceLA = sub2ind(PData.Size, (ones(Ncol,1)*(1:1:Nrow))',ones(Nrow,1)*(1:1:Ncol),sectNum*ones(Nrow,Ncol));

    if numSections>1, stps = numSections-1; else stps = 2; end
    nsectxt = uicontrol(hf,'Style','text',...
        'String','Section No.',...
        'Units','normalized',...
        'Position',[Pos+SG.TO,SG.TS],...
        'FontUnits','normalized',...
        'FontSize',SG.TF,...
        'FontWeight','bold');
    nsecSldr = uicontrol(hf,'Style','slider',...
        'Max',numSections,'Min',0.99,'Value',sectNum,...
        'SliderStep',[1/stps 2/stps],...
        'Units','normalized',...
        'Position',[Pos+SG.SO,SG.SS],...
        'BackgroundColor',bkgrnd-0.05,...
        'Interruptible','off',...
        'BusyAction','cancel',...
        'Callback',{@updateSect_Callback});
    nsecValue = uicontrol(hf,'Style','edit',...
        'String',num2str(sectNum,'%3.0f'),...
        'Units','normalized',...
        'Position',[Pos+SG.EO,SG.ES],...
        'BackgroundColor',bkgrnd+0.1,...
        'Callback',{@updateSect_Callback});
    changeView = uicontrol('Style','popupmenu',...
        'String',{'XZ View','YZ View', 'XY View'},...
        'Units','normalized',...
        'Position',[0.75,0.32,0.15,0.04],...
        'FontUnits','normalized',...
        'FontSize',PB.FS+0.08,...
        'Interruptible','off',...
        'BusyAction','cancel',...
        'tag','viewButton',...
        'Callback',{@changeView_Callback});
    set(changeView,'Value',3);
end

% PData number slider appers only with multiple PData settings
if length(PDataBase) > 1
    nsecx = 3*(1-3*0.28)/4 + 0.56;  % 3*(space between 3 controls of width 0.28) + 2 control
    nsecy = 0.16;
    Pos = [nsecx,nsecy];

    pdNumtxt = uicontrol(hf,'Style','text',...
        'String','PData No.',...
        'Units','normalized',...
        'Position',[Pos+SG.TO,SG.TS],...
        'FontUnits','normalized',...
        'FontSize',SG.TF,...
        'FontWeight','bold');
    pdNumSldr = uicontrol(hf,'Style','slider',...
        'Max',PDsldrMax,'Min',PDsldrMin,'Value',PDIn,...
        'SliderStep',[PDsldrSS PDsldrSS],...
        'Units','normalized',...
        'Position',[Pos+SG.SO,SG.SS],...
        'BackgroundColor',bkgrnd-0.05,...
        'Interruptible','off',...
        'BusyAction','cancel',...
        'Callback',{@pdNum_Callback});
    pdNumValue = uicontrol(hf,'Style','edit',...
        'String',num2str(PDIn,'%2.0f'),...
        'Units','normalized',...
        'Position',[Pos+SG.EO,SG.ES],...
        'BackgroundColor',bkgrnd+0.1,...
        'Callback',{@pdNum_Callback});
end


% Region number slider
nregx = 3*(1-3*0.28)/4 + 0.56;  % 3*(space between 3 controls of width 0.28) + 2 control
nregy = 0.08;
Pos = [nregx,nregy];
if ~isfield(PData,'Region')||isempty(PData.Region)
    evalin('base',['[PData(',int2str(PDIn),').Region] = computeRegions(PData(',int2str(PDIn),'));']);
    PData = evalin('base',['PData(',int2str(PDIn),')']);
end
regNum = 1;
numRegions = size(PData.Region,2);
if numRegions>1, stps = numRegions-1; else stps = 2; end
nregtxt = uicontrol(hf,'Style','text',...
    'String','Region No.',...
    'Units','normalized',...
    'Position',[Pos+SG.TO,SG.TS],...
    'FontUnits','normalized',...
    'FontSize',SG.TF,...
    'FontWeight','bold');
nregSldr = uicontrol(hf,'Style','slider',...
    'Max',numRegions,'Min',0.99,'Value',regNum,...
    'SliderStep',[1/stps 2/stps],...
    'Units','normalized',...
    'Position',[Pos+SG.SO,SG.SS],...
    'BackgroundColor',bkgrnd-0.05,...
    'Interruptible','off',...
    'BusyAction','cancel',...
    'Callback',{@nregSldr_Callback});
nregValue = uicontrol(hf,'Style','edit',...
    'String',num2str(regNum,'%3.0f'),...
    'Units','normalized',...
    'Position',[Pos+SG.EO,SG.ES],...
    'BackgroundColor',bkgrnd+0.1,...
    'Callback',{@nregValue_Callback});

% Show Regions push button.
showRegx = 3*(1-3*0.28)/4 + 0.55 + (0.28 - PB.BS(1))/2;
showRegy = 0.044;
showReg = uicontrol('Style','togglebutton',...
    'String','Show Reg.',...
    'Max',1,'Min',0,'Value',0,...
    'Units','normalized',...
    'Position',[showRegx,showRegy,PB.BS],...
    'FontUnits','normalized',...
    'FontSize',PB.FS,...
    'Interruptible','off',...
    'BusyAction','cancel',...
    'tag','showReg',...
    'Callback',{@showReg_Callback});

% Lock Regions and TX, only happens while numRegions is the same as num of
% TX
lockSx = 125/figWidth;
lockSy = SG.SS(2);
lockReg = 0;
if isequal(numRegions,size(TXIn,2)) && numRegions > 1
    lockTxReg = uicontrol('Style','checkbox',...
        'String','Lock TX & Region',...
        'Units','normalized',...
        'Position',[showRegx+(PB.BS(1)-lockSx)/2,0.02,lockSx,lockSy],...
        'FontUnits','normalized',...
        'FontSize',SG.TF,...
        'HorizontalAlignment','right',...
        'Interruptible','off',...
        'BusyAction','cancel',...
        'tag','lockTxReg',...
        'Callback',{@lockTxReg_Callback});
    set(showReg,'Position',[showRegx,showRegy+0.01,PB.BS]);
end

Uwidth = 120/figWidth;
AxesPos = get(ha,'Position');

changeUnit= uicontrol('Style','popupmenu',...
    'String',{'   wavelength','         mm   '},...
    'Units','normalized',...
    'Position',[0.525-Uwidth/2,0.275,Uwidth,0.02],...
    'FontUnits','normalized',...
    'FontSize',0.7,...
    'Interruptible','off',...
    'BusyAction','cancel',...
    'tag','changeUnit',...
    'Callback',{@changeUnit_Callback});
AxesUnit = 'wavelength';

speedOfSound = 1540;
if evalin('base','exist(''Resource'',''var'')')
    Resource = evalin('base','Resource');
    if isfield(Resource,'Parameters')
        if isfield(Resource.Parameters,'speedOfSound')
            speedOfSound = Resource.Parameters.speedOfSound/1000; % speed of sound in mm/usec
        end
    end

    if isfield(Resource,'DisplayWindow')
        if isfield(Resource.DisplayWindow,'Title')
            hf.Name = ['showTXPD - ',Resource.DisplayWindow.Title];
        end
    end

else
    disp('Speed of sound is not defined, use 1540 m/s as the default value');
end

% Make the GUI visible.
set(hf,'Visible','on');

% Calculate the maximum intensity of the TXPD data.
maxIntensity = 0;
for i = 1:size(TXIn,2)
    for j = 1:Nsec
        TXPD = double(TX(TXIn(i)).TXPD(:,:,j,1));
        peak = max(max(TXPD));
        if peak>maxIntensity, maxIntensity = peak; end
        clear TXPD
    end
end
maxIntensity = maxIntensity/256;
plotCode = 'intensity';  % default plotCode
[Data,clims] = updateImage;

%% Nested functions
    function [Data,clims] = updateImage()
        Nrow = PData.Size(1); Ncol = PData.Size(2); Nsec = PData.Size(3);
        setYDir = 0;
        if Nsec > 1
            switch shapeOrientation
                case 'xz'
                    TXPD = zeros(Nsec,Ncol,3);
                    TXPD(:,:,1) = permute(double(TX(txn).TXPD(sectNum,:,:,1)),[3 2 1])/256;
                    TXPD(:,:,2) = permute(double(TX(txn).TXPD(sectNum,:,:,2)),[3 2 1])/16;
                    TXPD(:,:,3) = permute(double(TX(txn).TXPD(sectNum,:,:,3)),[3 2 1])/16;
                    XLim = [-PData.Origin(2), -PData.Origin(2) + pdeltaY * PData.Size(1)];
                    YLim = [PData.Origin(3), PData.Origin(3) + pdeltaZ * PData.Size(3)];
                    sliceLA = sub2ind(PData.Size, sectNum*(ones(Nsec,Ncol)),ones(Nsec,1)*(1:1:Ncol),(ones(Ncol,1)*(1:1:Nsec))');
                case 'yz'
                    TXPD = zeros(Nsec,Nrow,3);
                    TXPD(:,:,1) = permute(double(TX(txn).TXPD(:,sectNum,:,1)),[3 1 2])/256;
                    TXPD(:,:,2) = permute(double(TX(txn).TXPD(:,sectNum,:,2)),[3 1 2])/16;
                    TXPD(:,:,3) = permute(double(TX(txn).TXPD(:,sectNum,:,3)),[3 1 2])/16;
                    XLim = [PData.Origin(1), PData.Origin(1) + pdeltaX * PData.Size(2)];
                    YLim = [PData.Origin(3), PData.Origin(3) + pdeltaZ * PData.Size(3)];
                    sliceLA = sub2ind(PData.Size, ones(Nsec,1)*(1:1:Nrow),sectNum*ones(Nsec,Nrow),(ones(Nrow,1)*(1:1:Nsec))');
                case 'xy'
                    TXPD = zeros(Nrow,Ncol,3);
                    TXPD(:,:,1) = (double(TX(txn).TXPD(:,:,sectNum,1)))/256;
                    TXPD(:,:,2) = (double(TX(txn).TXPD(:,:,sectNum,2)))/16;
                    TXPD(:,:,3) = (double(TX(txn).TXPD(:,:,sectNum,3)))/16;
                    XLim = [PData.Origin(1), PData.Origin(1) + pdeltaX * PData.Size(2)];
                    YLim = [PData.Origin(2), PData.Origin(2) - pdeltaY * PData.Size(1)];
                    sliceLA = sub2ind(PData.Size, (ones(Ncol,1)*(1:1:Nrow))',ones(Nrow,1)*(1:1:Ncol),sectNum*ones(Nrow,Ncol));
                    setYDir = 1;
            end
        else
            TXPD = double(TX(txn).TXPD);
            TXPD(:,:,1) = TXPD(:,:,1)/256;
            TXPD(:,:,2:3) = TXPD(:,:,2:3)/16;
            if isfield(PData,'Coord')&&(strcmp(PData.Coord,'polar'))
                XLim = [pdeltaT*(1-PData.Size(2))/2, pdeltaT*(PData.Size(2)-1)/2];
                YLim = [0, pdeltaR * PData.Size(1)];
            else
                XLim = [PData.Origin(1), PData.Origin(1) + pdeltaX * PData.Size(2)];
                YLim = [PData.Origin(3), PData.Origin(3) + pdeltaZ * PData.Size(1)];
            end
            sliceLA = reshape(find(ones(Nrow,Ncol)),Nrow,Ncol);
        end

        P = (TXPD(:,:,1)>TX(txn).peakCutOff)&(TXPD(:,:,3)<TX(txn).peakBLMax);
        P = +P;  % convert logical to double.
        switch plotCode
            case 'intensity'
                Data = TXPD(:,:,1) .* P;
                clims = [0 maxIntensity];
            case 'peakTime'
                Data = TXPD(:,:,2) .* P;
                clims = [0 max(max(TXPD(:,:,2)))];
            case 'burstLength'
                Data = TXPD(:,:,3) .* P;
                clims = [0 max(0.5,max(max(Data)))]; % don't allow max to go to 0
        end

        if isequal(clims(1),clims(2)), clims(2) = clims(1)+0.01; end % prevent error message from caxis

        if (get(showReg,'Value') == 1)&&(regNum~=0)
            regNum = round(get(nregSldr,'Value'));
            if ~isfield(PData,'Region')||isempty(PData.Region), return, end
            if ~isfield(PData.Region(regNum),'numPixels')||isempty(PData.Region(regNum).numPixels), return, end
            ci = clims(2)/4; % intensity of Region pixels
            [row,col] = find(ismember(sliceLA,PData.Region(regNum).PixelsLA+1));
            Data(sub2ind(size(Data),row(iseven(row+col)),col(iseven(row+col)))) = ci;
        end

        axes(ha),
        if strcmp(AxesUnit,'mm')
            XLim = XLim * speedOfSound/Trans.frequency;
            YLim = YLim * speedOfSound/Trans.frequency;
        end

        imagesc(XLim,YLim,Data,clims); if setYDir, set(ha,'YDir','normal'); end
        if isfield(PData,'Coord')&&(strcmp(PData.Coord,'polar'))
            xlabel('Angle in Radians');
            ylabel('Depth in Wavelengts');
            set(changeUnit,'Visible','off');
        end

        % plot Transducer elements
        if Trans.type ~= 2
            scaleToWvl = Trans.frequency/speedOfSound; % conversion factor from mm to wavelengths

            if strcmp(Trans.units, 'mm') && strcmp(AxesUnit,'wavelength')
                ElementPos = Trans.ElementPos * scaleToWvl;
            elseif strcmp(Trans.units, 'wavelengths') && strcmp(AxesUnit,'mm')
                ElementPos = Trans.ElementPos / scaleToWvl;
            else
                ElementPos = Trans.ElementPos;
            end
            if isfield(PData,'Coord')&&(strcmp(PData.Coord,'polar'))
                ElementPos(:,1) = linspace(XLim(1),XLim(2),Trans.numelements);
                switch Trans.type
                    case 0 % phase array
                        if isfield(PData.Region(regNum).Shape,'Position')
                            z = PData.Region(regNum).Shape.Position(3);
                        else
                            z = PData.Origin(3); % entire PData
                        end
                        ElementPos(:,3) = abs(z./cos(ElementPos(:,1)));
                    case 1 % curved array
                        ElementPos(:,3) = Trans.radius;
                end
            end
            hold on;
            if isfield(Trans,'HVMux')
                if isfield(TX(1,txn),'aperture') == 1
                    aperture = TX(1,txn).aperture;  %Import aperture Data
                else
                    aperture = 1;
                end
                if isfield(Trans.HVMux,'ApertureES')
                    ActEle = Trans.HVMux.ApertureES(:,aperture);
                else
                    ActEle = Trans.HVMux.Aperture(:,aperture);
                end
                ActEle(ActEle ~= 0) = 1;
                ActEle(ActEle == 0) = NaN;
            else
                ActEle = ones(1,Trans.numelements);
            end
            plot(ElementPos(ActEle==1,1),ElementPos(ActEle==1,3),'rs');
            hold off
        end
    end

    function parameter_Callback(~,eventdata)
        switch get(eventdata.NewValue,'Tag')
            case 'rbIntens'
                plotCode = 'intensity';
                [Data,clims] = updateImage;
            case 'rbTime'
                plotCode = 'peakTime';
                [Data,clims] = updateImage;
            case 'rbBL'
                plotCode = 'burstLength';
                [Data,clims] = updateImage;
        end
    end

    function txNum_Callback(hObject,~)
        Cntrl = get(hObject,'Style');
        if strcmp(Cntrl,'slider')
            txn = round(get(hObject,'Value'));
        else
            txn = round(str2double(get(hObject,'String')));
            txn = max(min(length(TX),txn),1);
        end

        set(txNumSldr,'Value',txn);
        set(txNumValue,'String',num2str(txn,'%3.0f'));
        if lockReg
            regNum = find(TXIn==txn);
            set(nregSldr,'Value',regNum);
            set(nregValue,'String',num2str(regNum,'%3.0f'));
        end

        if isequal(PData.Size(1:2), size(TX(txn).TXPD(:,:,1)))
            [Data,clims] = updateImage;
        else
            dismatchMsg(txn,PDIn);
        end

    end

    function txCutoffSldr_Callback(~,~)
        peakCutOff = get(txCutoffSldr,'Value');
        for nTX = 1:size(TX,2)
            TX(nTX).peakCutOff = peakCutOff;
        end
        [Data,clims] = updateImage;
        set(txCutoffValue,'String',num2str(peakCutOff,'%2.2f'));
    end

    function txCutoffValue_Callback(~,~)
        peakCutOff = str2double(get(txCutoffValue,'String'));
        peakCutOff = min(max(0,peakCutOff),peakMax);
        set(txCutoffSldr,'Value',peakCutOff);
        set(txCutoffValue,'String',num2str(peakCutOff,'%2.2f'));
        for nTX = 1:size(TX,2)
            TX(nTX).peakCutOff = peakCutOff;
        end
        [Data,clims] = updateImage;
    end

    function txBLMaxSldr_Callback(~,~)
        peakBLMax = get(txBLMaxSldr,'Value');
        set(txBLMaxValue,'String',num2str(peakBLMax,'%2.2f'));
        for nTX = 1:size(TX,2)
            TX(nTX).peakBLMax = peakBLMax;
        end
        [Data,clims] = updateImage;
    end

    function txBLMaxValue_Callback(~,~)
        peakBLMax = str2double(get(txBLMaxValue,'String'));
        peakBLMax = min(max(0,peakBLMax),BLMax);
        set(txBLMaxSldr,'Value',peakBLMax);
        set(txBLMaxValue,'String',num2str(peakBLMax,'%2.2f'));
        for nTX = 1:size(TX,2)
            TX(nTX).peakBLMax = peakBLMax;
        end
        [Data,clims] = updateImage;
    end

    function readTX_Callback(~,~)
        TX = evalin('base','TX');
        % Check for presence of the TXPD data, and compute if not already computed.
        if ~isfield(TX,'TXPD')
            for i=1:size(TX,2)
                TX(TXIn(i)).TXPD = computeTXPD(TX(TXIn(i)),PData);
            end
            [TX.peakCutOff] = deal(1.0);
            [TX.peakBLMax] = deal(4.0);
        else
            for i=1:size(TXIn,2)
                if isempty(TX(TXIn(i)).TXPD)
                    TX(TXIn(i)).TXPD = computeTXPD(TX(TXIn(i)),PData);
                end
            end
            if ~isfield(TX,'peakCutOff'), [TX.peakCutOff] = deal(1.0); end
            if ~isfield(TX,'peakBLMax'), [TX.peakBLMax] = deal(4.0); end
        end

        maxIntensity = 0;
        TXPD = zeros(Nrow,Ncol);
        for i = 1:size(TXIn,2)
            for j = 1:Nsec
                TXPD(:,:) = double(TX(TXIn(i)).TXPD(:,:,j,1));
                peak = max(max(TXPD));
                if peak>maxIntensity, maxIntensity = peak; end
            end
        end
        maxIntensity = maxIntensity/256;
        if isempty(TX(txn).peakCutOff), TX(txn).peakCutOff = 1.0; end
        set(txCutoffSldr,'Value',TX(txn).peakCutOff);
        set(txCutoffValue,'String',num2str(TX(txn).peakCutOff,'%2.2f'));
        if isempty(TX(txn).peakBLMax), TX(txn).peakBLMax = 4.0; end
        set(txBLMaxSldr,'Value',TX(txn).peakBLMax);
        set(txBLMaxValue,'String',num2str(TX(txn).peakBLMax,'%2.2f'));
        [Data,clims] = updateImage;
    end

    function writeTX_Callback(~,~)
        assignin('base', 'TX', TX);
        if evalin('base','exist(''Control'',''var'')')
            Control = evalin('base','Control');
            Control.Command = 'update&Run';
            Control.Parameters = {'TX','Recon'};
            assignin('base','Control', Control);
        end
    end

    function updateSect_Callback(hObject,~)

        Cntrl = get(hObject,'Style');
        if strcmp(Cntrl,'slider')
            sectNum = round(get(hObject,'Value'));
        else
            sectNum = round(str2double(get(hObject,'String')));
            sectNum = max(min(numSections,sectNum),0.99);
        end
        set(nsecSldr,'Value',sectNum);
        set(nsecValue,'String',num2str(sectNum,'%3.0f'));

        switch shapeOrientation
            case 'xz'
                sectNumY = sectNum;
            case 'yz'
                sectNumX = sectNum;
            case 'xy'
                sectNumZ = sectNum;
        end

        [Data,clims] = updateImage;
    end

    function pdNum_Callback(hObject,~)

        Cntrl = get(hObject,'Style');
        if strcmp(Cntrl,'slider')
            PDIn = round(get(hObject,'Value'));
        else
            PDIn = round(str2double(get(hObject,'String')));
            PDIn = max(min(length(PDataBase),PDIn),1);
        end

        set(pdNumSldr,'Value',PDIn);
        set(pdNumValue,'String',num2str(PDIn,'%2.0f'));

        % if PDdata number is changed, new PDIn and Region slider should be
        % assigned

        PData = PDataBase(PDIn);
        Delta = updatePData(PData);
        deltaFieldsName = fieldnames(Delta);
        for i = 1:length(deltaFieldsName)
            eval(['pdelta',deltaFieldsName{i},'= Delta.',deltaFieldsName{i},';']);
        end

        if ~isfield(PData,'Region')||isempty(PData.Region)
            evalin('base',['[PData(',int2str(PDIn),').Region] = computeRegions(PData(',int2str(PDIn),'));']);
            PData = evalin('base',['PData(',int2str(PDIn),')']);
        end
        regNum = 1;
        numRegions = size(PData.Region,2);
        if numRegions>1, stps = numRegions-1; else stps = 2; end
        set(nregSldr,'Max',numRegions);
        set(nregSldr,'Value',regNum);
        set(nregSldr,'SliderStep',[1/stps 2/stps]);
        set(nregValue,'String',num2str(regNum,'%3.0f'));

        if isequal(PData.Size(1:2), size(TX(txn).TXPD(:,:,1)))
            [Data,clims] = updateImage;
        else
            dismatchMsg(txn,PDIn);
        end
    end

    function nregSldr_Callback(~,~)
        regNum = round(get(nregSldr,'Value'));
        set(nregValue,'String',num2str(regNum,'%3.0f'));
        if lockReg
            txn = TXIn(regNum);
            set(txNumSldr,'Value',txn);
            set(txNumValue,'String',num2str(txn,'%3.0f'));
        end
        if isequal(PData.Size(1:2), size(TX(txn).TXPD(:,:,1)))
            [Data,clims] = updateImage;
        else
            dismatchMsg(txn,PDIn);
        end
    end

    function nregValue_Callback(~,~)
        nrv = str2double(get(nregValue,'String'));
        if (nrv>numRegions)||(nrv<1)
            set(nregValue,'String',num2str(round(get(nregSldr,'Value')),'%3.0f'));
            return
        end
        regNum = round(nrv);
        set(nregSldr,'Value',regNum);
        if lockReg
            txn = TXIn(regNum);
            set(txNumSldr,'Value',txn);
            set(txNumValue,'String',num2str(txn,'%3.0f'));
        end
        if isequal(PData.Size(1:2), size(TX(txn).TXPD(:,:,1)))
            [Data,clims] = updateImage;
        else
            dismatchMsg(txn,PDIn);
        end
    end

    function showReg_Callback(~,~)
        if (get(showReg,'Value') == 1)&&(regNum~=0)
            if isequal(PData.Size(1:2), size(TX(txn).TXPD(:,:,1)))
                [Data,clims] = updateImage;
            else
                dismatchMsg(txn,PDIn);
            end
        else
            [Data,clims] = updateImage;
        end
    end

    function lockTxReg_Callback(~,~)
        lockReg = get(lockTxReg,'Value');
        regNum = find(TXIn==txn);
        set(nregSldr,'Value',regNum);
        set(nregValue,'String',num2str(regNum,'%3.0f'));
        [Data,clims] = updateImage;
    end

    function changeView_Callback(hObject,~)
        planeValue = get(hObject,'Value');
        switch planeValue
            case 1 % xz
                shapeOrientation = 'xz';
                numSections = PData.Size(1);
                sectNum = sectNumY;
                numSections = Nrow;
            case 2 % yz
                shapeOrientation = 'yz';
                numSections = PData.Size(2);
                sectNum = sectNumX;
                numSections = Ncol;
            case 3 % xy
                shapeOrientation = 'xy';
                numSections = PData.Size(3);
                sectNum = sectNumZ;
                numSections = Nsec;
        end
        stp = numSections - 1;
        set(nsecSldr,'Max',numSections,'Value',sectNum,'SliderStep',[1/stp 2/stp]);
        set(nsecValue,'String',num2str(sectNum,'%3.0f'));
        [Data,clims] = updateImage;
    end

    function changeUnit_Callback(varargin)
        if isequal(get(changeUnit,'Value'),1)
            AxesUnit = 'wavelength';
        else
            AxesUnit = 'mm';
        end
        [Data,clims] = updateImage;
    end

    function dismatchMsg(TXnum, PDnum)

        if (get(showReg,'Value') == 1)&&(regNum~=0) % remove region display

            P = (TXPD(:,:,1)>TX(txn).peakCutOff)&(TXPD(:,:,3)<TX(txn).peakBLMax);
            P = +P;  % convert logical to double.

            switch plotCode
                case 'intensity'
                    Data = TXPD(:,:,1) .* P;
                    clims = [0 maxIntensity];
                case 'peakTime'
                    Data = TXPD(:,:,2) .* P;
                    clims = [0 max(max(TXPD(:,:,2)))];
                case 'burstLength'
                    Data = TXPD(:,:,3) .* P;
                    clims = [0 max(0.5,max(max(Data)))]; % don't allow max to go to 0
            end

            axes(ha),
            if strcmp(AxesUnit,'mm')
                XLim = XLim * speedOfSound/Trans.frequency;
                YLim = YLim * speedOfSound/Trans.frequency;
            end

            imagesc(XLim,YLim,Data,clims);
            if isfield(PData,'Coord')&&(strcmp(PData.Coord,'polar'))
                xlabel('Angle in Radians');
                set(changeUnit,'Visible','off');
            end

            % plot Transducer elements in wavelength
            if Trans.type ~= 2
                scaleToWvl = Trans.frequency/speedOfSound; % conversion factor from mm to wavelengths
                if strcmp(Trans.units, 'mm') && strcmp(AxesUnit,'wavelength')
                    Trans.ElementPos = Trans.ElementPos * scaleToWvl;
                elseif strcmp(Trans.units, 'wavelengths') && strcmp(AxesUnit,'mm')
                    Trans.ElementPos = Trans.ElementPos / scaleToWvl;
                end
                if isfield(PData,'Coord')&&(strcmp(PData.Coord,'polar'))
                    Trans.ElementPos(:,1) = linspace(XLim(1),XLim(2),Trans.numelements);
                    switch Trans.type
                        case 0 % phase array
                            z = PData.Region(1).Shape.Position(3);
                            Trans.ElementPos(:,3) = abs(z./cos(Trans.ElementPos(:,1)));
                        case 1 % curved array
                            Trans.ElementPos(:,3) = Trans.radius;
                    end
                end

                hold on;plot(Trans.ElementPos(:,1), Trans.ElementPos(:,3),'rs'); hold off
            end
        end
        warndlg(['TX(',num2str(TXnum),').TXPD does not match PData(',num2str(PDnum),')!'],'Warning');
        set(showReg,'Value',0);
    end

    function closefunc(~,~)
        hFigure = findobj('Type','figure');
        if ~isempty(hFigure)
            close(hFigure(contains(get(hFigure,'Tag'),'Msgbox')));
        end
        delete(hf);
    end

end

function ise = iseven( val )
    ise = vsv.math.EquMath.iseven(val);
end
